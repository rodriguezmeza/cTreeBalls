/* FCFC-style dynamic MPI scheduling for the octree-GGG estimator. */

#if defined(OCTREE2BALLS_GGG_MPI_COMPAT) && !defined(OCTREEGGGMPI)
#define OCTREEGGGMPI
#endif

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_octree_ggg_mpi.h"

#ifdef OCTREEGGGMPI

#define FCFC_MPI_ROOT 0

#define FCFC_MPI_REAL MPI_DOUBLE

#ifdef LONGINT
#define FCFC_MPI_INTEGER MPI_LONG
#else
#define FCFC_MPI_INTEGER MPI_INT
#endif

/* Per-context rank view; process MPI ownership lives in mpi_runtime.c. */
#define CBALLS_MPI_STATE (cballs_runtime_current()->mpi[CBALLS_MPI_FCFC_OCTREE_GGG_MPI])
#define mpi_active (CBALLS_MPI_STATE.active)
#define mpi_rank (CBALLS_MPI_STATE.rank)
#define mpi_size (CBALLS_MPI_STATE.size ? CBALLS_MPI_STATE.size : 1)

static bool fcfc_octree_ggg_mpi_selected(const struct cmdline_data *cmd)
{
    if (cmd == NULL || cmd->searchMethod == NULL)
        return false;
    if (strcmp(cmd->searchMethod, "octree-ggg-mpi") == 0)
        return true;
#ifdef OCTREE2BALLSMPI
    return strcmp(cmd->searchMethod, "octree-2balls-mpi") == 0
        && cballs_opt_legacy_one_ball(cmd);
#else
    return false;
#endif
}



static int mpi_error(struct cmdline_data *cmd, const char *operation,
                     int status)
{
    return cballs_mpi_shared_error(cmd, operation, status);
}

int fcfc_octree_ggg_mpi_prepare(struct cmdline_data *cmd,
                                struct global_data *gd)
{
    return cballs_mpi_shared_prepare(&CBALLS_MPI_STATE, cmd, gd, !(!fcfc_octree_ggg_mpi_selected(cmd)));
}

int fcfc_octree_ggg_mpi_finalize(struct cmdline_data *cmd)
{
    return cballs_mpi_shared_finalize(&CBALLS_MPI_STATE, cmd);
}

int fcfc_octree_ggg_mpi_active(void) { return mpi_active; }
int fcfc_octree_ggg_mpi_is_root(void) {
    return !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}
int fcfc_octree_ggg_mpi_rank(void) { return mpi_rank; }
int fcfc_octree_ggg_mpi_size(void) { return mpi_size; }

int fcfc_octree_ggg_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !fcfc_octree_ggg_mpi_selected(cmd)
        || !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}

int fcfc_octree_ggg_mpi_consensus(struct cmdline_data *cmd,
                                  int local_status,
                                  const char *operation)
{
    return cballs_mpi_shared_consensus(&CBALLS_MPI_STATE, cmd, !(!fcfc_octree_ggg_mpi_selected(cmd)), local_status, operation);
}

int fcfc_octree_ggg_mpi_scheduler_init(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler,
    INTEGER task_count, INTEGER step)
{
    MPI_Aint bytes = 0;
    int root_status = SUCCESS;
    int status;

    memset(scheduler, 0, sizeof(*scheduler));
    scheduler->window = MPI_WIN_NULL;
    scheduler->task_count = (uint64_t)task_count;
    scheduler->step = step > 0 ? (uint64_t)step : 1;
    if (mpi_size == 1) {
        scheduler->ready = TRUE;
        return SUCCESS;
    }
    if (mpi_rank == FCFC_MPI_ROOT) {
        bytes = sizeof(*scheduler->counter);
        if ((status = MPI_Alloc_mem(bytes, MPI_INFO_NULL,
                                    &scheduler->counter)) != MPI_SUCCESS)
            root_status = mpi_error(cmd, "MPI_Alloc_mem", status);
        else
            *scheduler->counter = 0;
    }
    if ((status = MPI_Bcast(&root_status, 1, MPI_INT, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI scheduler-status broadcast", status);
    if (root_status == FAILURE) return FAILURE;
    status = MPI_Win_create(scheduler->counter, bytes, sizeof(uint64_t),
                            MPI_INFO_NULL, MPI_COMM_WORLD, &scheduler->window);
    if (status != MPI_SUCCESS) {
        if (mpi_rank == FCFC_MPI_ROOT) MPI_Free_mem(scheduler->counter);
        scheduler->counter = NULL;
        return mpi_error(cmd, "MPI_Win_create", status);
    }
    scheduler->ready = TRUE;
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_scheduler_claim(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler,
    INTEGER *first, INTEGER *last)
{
    uint64_t start;
    int status;

    if (!scheduler->ready) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-mpi scheduler is not initialized");
        return FAILURE;
    }
    if (mpi_size == 1) {
        start = scheduler->local_next;
        scheduler->local_next += scheduler->step;
    } else {
        int locked = FALSE;
        status = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, FCFC_MPI_ROOT, 0,
                              scheduler->window);
        if (status == MPI_SUCCESS) {
            locked = TRUE;
            status = MPI_Fetch_and_op(&scheduler->step, &start, MPI_UINT64_T,
                                      FCFC_MPI_ROOT, 0, MPI_SUM,
                                      scheduler->window);
        }
        if (locked) {
            const int unlock_status = MPI_Win_unlock(FCFC_MPI_ROOT,
                                                      scheduler->window);
            if (status == MPI_SUCCESS) status = unlock_status;
        }
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI RMA task claim", status);
    }
    if (start >= scheduler->task_count) {
        *first = *last = (INTEGER)scheduler->task_count;
        return SUCCESS;
    }
    uint64_t end = start + scheduler->step;
    if (end > scheduler->task_count) end = scheduler->task_count;
    *first = (INTEGER)start;
    *last = (INTEGER)end;
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_scheduler_destroy(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler)
{
    int status = MPI_SUCCESS;
    int free_status = MPI_SUCCESS;

    if (!scheduler->ready) return SUCCESS;
    if (mpi_size > 1) {
        status = MPI_Win_free(&scheduler->window);
        if (mpi_rank == FCFC_MPI_ROOT && scheduler->counter != NULL)
            free_status = MPI_Free_mem(scheduler->counter);
    }
    memset(scheduler, 0, sizeof(*scheduler));
    if (status != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Win_free", status);
    if (free_status != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Free_mem", free_status);
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_reduce_reals(struct cmdline_data *cmd,
                                     real *values, size_t count)
{
    size_t offset = 0;
    while (offset < count) {
        int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;
        if (mpi_rank == FCFC_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                FCFC_MPI_REAL, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk, FCFC_MPI_REAL,
                                MPI_SUM, FCFC_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI histogram reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_reduce_integers(struct cmdline_data *cmd,
                                        INTEGER *values, size_t count)
{
    size_t offset = 0;
    while (offset < count) {
        int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;
        if (mpi_rank == FCFC_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                FCFC_MPI_INTEGER, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                FCFC_MPI_INTEGER, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI counter reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* OCTREEGGGMPI */
