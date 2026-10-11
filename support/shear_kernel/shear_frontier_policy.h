#ifndef CBALLS_SHEAR_FRONTIER_POLICY_H
#define CBALLS_SHEAR_FRONTIER_POLICY_H

#include <errno.h>
#include <limits.h>

/* A rank owns only its own task histograms. Catalogs and schedules remain
 * replicated. Keep the previous OpenMP target unless explicitly overridden. */
static int shear_frontier_integer(struct cmdline_data *cmd, const char *name,
                                  size_t fallback, size_t minimum,
                                  size_t maximum, size_t *result)
{
    const char *text = getenv(name);
    char *end;
    unsigned long long value;
    *result = fallback;
    if (text == NULL) return SUCCESS;
    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno || end == text || *end || text[0] == '-'
        || value < minimum || value > maximum) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s must be an integer in [%zu,%zu]", name, minimum, maximum);
        return FAILURE;
    }
    *result = (size_t)value;
    return SUCCESS;
}

static int shear_frontier_policy(struct cmdline_data *cmd, int threads,
        int ranks, size_t stride, size_t serial_target,
        INTEGER *target, size_t *rank_limit)
{
    const size_t maximum = (size_t)INT_MAX - 8;
    size_t requested, memory_mb, desired = serial_target;
    size_t bytes, budget, global_limit;
    ranks = MAX(1, ranks);
    if (shear_frontier_integer(cmd, "CBALLS_SHEAR_TASKS", 0, 0,
                              maximum, &requested) == FAILURE
        || shear_frontier_integer(cmd, "CBALLS_SHEAR_TASK_MEMORY_MB", 256,
                     1, MIN((size_t)1048576, SIZE_MAX >> 20), &memory_mb)
           == FAILURE)
        return FAILURE;
    if (!stride || stride > SIZE_MAX/sizeof(real)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid task histogram size");
        return FAILURE;
    }
    bytes = stride*sizeof(real);
    budget = memory_mb << 20;
    /* One histogram is the unavoidable minimum, even if it exceeds budget. */
    *rank_limit = MAX((size_t)1, budget/bytes);
    *rank_limit = MIN(*rank_limit, maximum);
    global_limit = *rank_limit > maximum/(size_t)ranks
        ? maximum : *rank_limit*(size_t)ranks;
#ifdef SHEAR_MPI_ENABLED
    size_t workers = (size_t)MAX(1, threads);
    desired = workers > maximum/8/(size_t)ranks
        ? maximum : MAX((size_t)64, 8*workers*(size_t)ranks);
    /* Retain the native octree's existing minimum frontier size. */
    desired = MAX(desired, serial_target);
#else
    (void)threads;
#endif
    if (requested) desired = requested;
    *target = (INTEGER)MIN(desired, global_limit);
    return SUCCESS;
}

typedef struct {
    INTEGER task;
    long double cost;
    int owner;
    INTEGER local_index;
} shear_frontier_priority;

static int shear_frontier_compare(const void *a, const void *b)
{
    const shear_frontier_priority *x = a, *y = b;
    if (x->cost != y->cost) return x->cost > y->cost ? -1 : 1;
    return (x->task > y->task) - (x->task < y->task);
}

/* Longest estimated task first, deterministic ties, bounded rank scratch.
 * Sorting avoids the previous quadratic scan for every assignment. */
static int shear_frontier_assign(struct cmdline_data *cmd,
        shear_frontier_priority *order, INTEGER count, int ranks,
        size_t rank_limit)
{
    long double *work = calloc((size_t)ranks, sizeof(*work));
    INTEGER *tasks = calloc((size_t)ranks, sizeof(*tasks));
    int status = FAILURE;
    if (work == NULL || tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": rank-work allocation failed");
        goto cleanup;
    }
    qsort(order, (size_t)count, sizeof(*order), shear_frontier_compare);
    for (INTEGER i = 0; i < count; i++) {
        int owner = -1;
        for (int rank = 0; rank < ranks; rank++) {
            if ((size_t)tasks[rank] >= rank_limit) continue;
            if (owner < 0 || work[rank] < work[owner]
                || (work[rank] == work[owner] && tasks[rank] < tasks[owner]))
                owner = rank;
        }
        if (owner < 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": pivot tasks exceed rank scratch budget");
            goto cleanup;
        }
        order[i].owner = owner;
        order[i].local_index = tasks[owner]++;
        work[owner] += order[i].cost;
    }
    status = SUCCESS;
cleanup:
    free(tasks);
    free(work);
    return status;
}

#endif
