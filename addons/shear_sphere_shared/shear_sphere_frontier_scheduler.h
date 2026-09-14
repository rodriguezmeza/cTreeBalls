#ifndef CTREEBALLS_SHEAR_SPHERE_FRONTIER_SCHEDULER_H
#define CTREEBALLS_SHEAR_SPHERE_FRONTIER_SCHEDULER_H

/* Adaptive exact-pivot scheduling for the native spherical octree. */
#define SHEAR_NATIVE_PIVOT_TASK_TARGET ((INTEGER)256)
#define SHEAR_NATIVE_NEIGHBOR_FRONTIER_TARGET ((INTEGER)256)
#define SHEAR_NATIVE_TASK_HISTOGRAM_MEMORY ((size_t)256 << 20)

typedef struct {
    nodeptr pivot_node;
    size_t pivot_offset;
    INTEGER pivot_count;
    int owner;
    INTEGER local_index;
    size_t first_offset;
    INTEGER first_count;
    size_t second_offset;
    INTEGER second_count;
    long double estimated_work;
} shear_native_pivot_task;

typedef struct {
    shear_native_pivot_task *tasks;
    INTEGER task_count;
    bodyptr *pivots;
    size_t pivot_count;
    nodeptr *first_nodes;
    nodeptr *second_nodes;
    size_t first_node_count;
    size_t second_node_count;
} shear_native_frontier_schedule;

static INTEGER shear_native_node_point_count(nodeptr node)
{
    return Type(node) == CELL ? Nb(node) : 1;
}

static bool shear_native_pivot_is_active(
        struct cmdline_data *cmd, bodyptr pivot, bodyptr first,
        bodyptr finish, bool read_mask)
{
    if (pivot < first || pivot >= finish
        || (read_mask && Mask(pivot) != MASK_NODE_VALID))
        return FALSE;
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd) && !Update(pivot))
        return FALSE;
#else
    (void)cmd;
#endif
    return TRUE;
}

static int shear_native_collect_pivots(
        struct cmdline_data *cmd, nodeptr node, bodyptr first,
        bodyptr finish, bool read_mask, bodyptr *result,
        size_t capacity, size_t *count)
{
    nodeptr child;

    if (node == NULL || (read_mask && Mask(node) == MASK_NODE_MASKED))
        return SUCCESS;
    if (Type(node) == CELL) {
        for (child = More(node); child != Next(node); child = Next(child))
            if (shear_native_collect_pivots(
                    cmd, child, first, finish, read_mask,
                    result, capacity, count) == FAILURE)
                return FAILURE;
        return SUCCESS;
    }
    if ((Type(node) != BODY && Type(node) != BODY3)
        || !shear_native_pivot_is_active(
               cmd, (bodyptr)node, first, finish, read_mask))
        return SUCCESS;
    if (*count == SIZE_MAX || (result != NULL && *count >= capacity)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": pivot-frontier body overflow");
        return FAILURE;
    }
    if (result != NULL) result[*count] = (bodyptr)node;
    (*count)++;
    return SUCCESS;
}

static int shear_native_build_node_frontier(
        const shear_sphere_pair_context *context, nodeptr root,
        INTEGER target, nodeptr **result, INTEGER *result_count)
{
    const size_t capacity = (size_t)target + NSUB;
    nodeptr *frontier;
    INTEGER count = 1;

    *result = NULL;
    *result_count = 0;
    if (root == NULL || target < 1
        || capacity < (size_t)target
        || capacity > SIZE_MAX/sizeof(*frontier)
        || (frontier = malloc(capacity*sizeof(*frontier))) == NULL) {
        snprintf(context->cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": neighbor-frontier allocation failed");
        return FAILURE;
    }
    frontier[0] = shear_sphere_collapse_single_child(context, root);
    while (count < target) {
        nodeptr selected_children[NSUB];
        INTEGER selected = -1;
        INTEGER largest = -1;
        int selected_count = 0;

        for (INTEGER index = 0; index < count; index++) {
            nodeptr children[NSUB];
            const int child_count = shear_sphere_children(
                context, frontier[index], children);

            if (child_count >= 2
                && shear_native_node_point_count(frontier[index]) > largest) {
                selected = index;
                largest = shear_native_node_point_count(frontier[index]);
                selected_count = child_count;
                for (int child = 0; child < child_count; child++)
                    selected_children[child] =
                        shear_sphere_collapse_single_child(
                            context, children[child]);
            }
        }
        if (selected < 0) break;
        memmove(frontier + selected + selected_count,
                frontier + selected + 1,
                (size_t)(count - selected - 1)*sizeof(*frontier));
        memcpy(frontier + selected, selected_children,
               (size_t)selected_count*sizeof(*frontier));
        count += selected_count - 1;
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}

static long double shear_native_estimate_task(
        const shear_sphere_pair_context *context, nodeptr pivot,
        size_t pivot_points, nodeptr *first_frontier,
        INTEGER first_count, nodeptr *second_frontier,
        INTEGER second_count)
{
    long double neighbor_points = 0.0L;

    for (INTEGER index = 0; index < first_count; index++)
        if (shear_sphere_pair_frontier_overlap(
                context, pivot, first_frontier[index]))
            neighbor_points += (long double)shear_native_node_point_count(
                first_frontier[index]);
    for (INTEGER index = 0; index < second_count; index++)
        if (shear_sphere_pair_frontier_overlap(
                context, pivot, second_frontier[index]))
            neighbor_points += (long double)shear_native_node_point_count(
                second_frontier[index]);
    return (long double)pivot_points*neighbor_points;
}

static shear_native_pivot_task shear_native_make_task(
        const shear_sphere_pair_context *context, nodeptr pivot,
        nodeptr *first_frontier, INTEGER first_count,
        nodeptr *second_frontier, INTEGER second_count)
{
    shear_native_pivot_task task;

    memset(&task, 0, sizeof(task));
    task.pivot_node = pivot;
    task.owner = -1;
    task.estimated_work = shear_native_estimate_task(
        context, pivot, (size_t)shear_native_node_point_count(pivot),
        first_frontier, first_count, second_frontier, second_count);
    return task;
}

static int shear_native_assign_tasks(
        struct cmdline_data *cmd, shear_native_pivot_task *tasks,
        INTEGER count, int rank_count)
{
    long double *rank_work;
    INTEGER *rank_tasks;

    rank_count = MAX(1, rank_count);
    rank_work = calloc((size_t)rank_count, sizeof(*rank_work));
    rank_tasks = calloc((size_t)rank_count, sizeof(*rank_tasks));
    if (rank_work == NULL || rank_tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": rank-work allocation failed");
        free(rank_tasks);
        free(rank_work);
        return FAILURE;
    }
    for (INTEGER task = 0; task < count; task++) tasks[task].owner = -1;
    for (INTEGER assigned = 0; assigned < count; assigned++) {
        INTEGER largest = -1;
        int owner = 0;

        for (INTEGER task = 0; task < count; task++)
            if (tasks[task].owner < 0
                && (largest < 0
                    || tasks[task].estimated_work
                       > tasks[largest].estimated_work))
                largest = task;
        for (int rank = 1; rank < rank_count; rank++)
            if (rank_work[rank] < rank_work[owner]) owner = rank;
        tasks[largest].owner = owner;
        tasks[largest].local_index = rank_tasks[owner]++;
        rank_work[owner] += tasks[largest].estimated_work;
    }
    free(rank_tasks);
    free(rank_work);
    return SUCCESS;
}

static size_t shear_native_count_neighbors(
        const shear_sphere_pair_context *context, nodeptr pivot,
        nodeptr *frontier, INTEGER count)
{
    size_t result = 0;

    for (INTEGER index = 0; index < count; index++)
        if (shear_sphere_pair_frontier_overlap(
                context, pivot, frontier[index]))
            result++;
    return result;
}

static void shear_native_fill_neighbors(
        const shear_sphere_pair_context *context, nodeptr pivot,
        nodeptr *frontier, INTEGER count, nodeptr *result)
{
    size_t output = 0;

    for (INTEGER index = 0; index < count; index++)
        if (shear_sphere_pair_frontier_overlap(
                context, pivot, frontier[index]))
            result[output++] = frontier[index];
}

static void shear_native_release_frontier_schedule(
        shear_native_frontier_schedule *schedule)
{
    if (schedule == NULL) return;
    free(schedule->second_nodes);
    free(schedule->first_nodes);
    free(schedule->pivots);
    free(schedule->tasks);
    memset(schedule, 0, sizeof(*schedule));
}

static int shear_native_build_frontier_schedule(
        struct cmdline_data *cmd, struct global_data *gd,
        nodeptr pivot_root, bodyptr pivot_first, bodyptr pivot_finish,
        nodeptr first_root, nodeptr second_root, INTEGER target,
        int rank_count, shear_native_frontier_schedule *schedule)
{
    shear_sphere_pair_context context;
    nodeptr *first_frontier = NULL;
    nodeptr *second_frontier = NULL;
    INTEGER first_frontier_count = 0;
    INTEGER second_frontier_count = 0;
    const size_t task_capacity = (size_t)target + NSUB;
    INTEGER count = 1;
    size_t pivot_total = 0;
    size_t first_total = 0;
    size_t second_total = 0;
    size_t expected = 0;
    const bool read_mask = cballs_opt_read_mask(cmd);
    int status = FAILURE;

    memset(schedule, 0, sizeof(*schedule));
    context.cmd = cmd;
    context.gd = gd;
    context.bins = cmd->sizeHistN;
    if (pivot_root == NULL || first_root == NULL || target < 1
        || task_capacity < (size_t)target
        || task_capacity > SIZE_MAX/sizeof(*schedule->tasks)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid native pivot frontier");
        return FAILURE;
    }
    if (shear_native_build_node_frontier(
            &context, first_root, SHEAR_NATIVE_NEIGHBOR_FRONTIER_TARGET,
            &first_frontier, &first_frontier_count) == FAILURE)
        goto cleanup;
    if (second_root != NULL
        && shear_native_build_node_frontier(
               &context, second_root,
               SHEAR_NATIVE_NEIGHBOR_FRONTIER_TARGET,
               &second_frontier, &second_frontier_count) == FAILURE)
        goto cleanup;
    schedule->tasks = calloc(task_capacity, sizeof(*schedule->tasks));
    if (schedule->tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": native pivot-frontier allocation failed");
        goto cleanup;
    }
    schedule->tasks[0] = shear_native_make_task(
        &context, shear_sphere_collapse_single_child(&context, pivot_root),
        first_frontier, first_frontier_count,
        second_frontier, second_frontier_count);
    while (count < target) {
        INTEGER best = -1;
        long double largest_work = -1.0L;
        nodeptr best_children[NSUB];
        int best_child_count = 0;

        for (INTEGER task = 0; task < count; task++) {
            nodeptr children[NSUB];
            const int child_count = shear_sphere_children(
                &context, schedule->tasks[task].pivot_node, children);

            if (child_count >= 2
                && schedule->tasks[task].estimated_work > largest_work) {
                best = task;
                largest_work = schedule->tasks[task].estimated_work;
                best_child_count = child_count;
                memcpy(best_children, children,
                       (size_t)child_count*sizeof(*children));
            }
        }
        if (best < 0) break;
        memmove(schedule->tasks + best + best_child_count,
                schedule->tasks + best + 1,
                (size_t)(count - best - 1)*sizeof(*schedule->tasks));
        for (int child = 0; child < best_child_count; child++) {
            nodeptr child_root = shear_sphere_collapse_single_child(
                &context, best_children[child]);

            schedule->tasks[best + child] = shear_native_make_task(
                &context, child_root,
                first_frontier, first_frontier_count,
                second_frontier, second_frontier_count);
        }
        count += best_child_count - 1;
    }

    {
        INTEGER active_tasks = 0;

        for (INTEGER task = 0; task < count; task++) {
            size_t task_pivots = 0;

            if (shear_native_collect_pivots(
                    cmd, schedule->tasks[task].pivot_node,
                    pivot_first, pivot_finish, read_mask,
                    NULL, 0, &task_pivots) == FAILURE)
                goto cleanup;
            if (task_pivots == 0) continue;
            if (task_pivots > (size_t)
#ifdef LONGINT
                    LONG_MAX
#else
                    INT_MAX
#endif
                ) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME ": task pivot count overflow");
                goto cleanup;
            }
            schedule->tasks[active_tasks] = schedule->tasks[task];
            schedule->tasks[active_tasks].pivot_offset = pivot_total;
            schedule->tasks[active_tasks].pivot_count = (INTEGER)task_pivots;
            schedule->tasks[active_tasks].estimated_work =
                shear_native_estimate_task(
                    &context, schedule->tasks[active_tasks].pivot_node,
                    task_pivots, first_frontier, first_frontier_count,
                    second_frontier, second_frontier_count);
            if (task_pivots > SIZE_MAX - pivot_total) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME ": pivot count overflow");
                goto cleanup;
            }
            pivot_total += task_pivots;
            active_tasks++;
        }
        count = active_tasks;
    }
    for (bodyptr pivot = pivot_first; pivot < pivot_finish; pivot++)
        if (shear_native_pivot_is_active(
                cmd, pivot, pivot_first, pivot_finish, read_mask))
            expected++;
    if (count <= 0 || pivot_total != expected
#ifdef LONGINT
        || pivot_total > (size_t)LONG_MAX
#else
        || pivot_total > (size_t)INT_MAX
#endif
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME
                 ": native pivot frontier covers %zu pivots, expected %zu",
                 pivot_total, expected);
        goto cleanup;
    }
    schedule->pivots = malloc(pivot_total*sizeof(*schedule->pivots));
    if (schedule->pivots == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": native pivot-list allocation failed");
        goto cleanup;
    }
    for (INTEGER task = 0; task < count; task++) {
        size_t filled = schedule->tasks[task].pivot_offset;
        const size_t finish = filled + (size_t)schedule->tasks[task].pivot_count;

        if (shear_native_collect_pivots(
                cmd, schedule->tasks[task].pivot_node,
                pivot_first, pivot_finish, read_mask,
                schedule->pivots, pivot_total, &filled) == FAILURE
            || filled != finish) {
            if (cmd->error_message[0] == '\0')
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME ": native pivot fill mismatch");
            goto cleanup;
        }
    }
    if (shear_native_assign_tasks(
            cmd, schedule->tasks, count, rank_count) == FAILURE)
        goto cleanup;

    for (INTEGER task = 0; task < count; task++) {
        const size_t first_count = shear_native_count_neighbors(
            &context, schedule->tasks[task].pivot_node,
            first_frontier, first_frontier_count);
        const size_t second_count = second_root != NULL
            ? shear_native_count_neighbors(
                  &context, schedule->tasks[task].pivot_node,
                  second_frontier, second_frontier_count)
            : 0;

        if (first_count > SIZE_MAX - first_total
            || second_count > SIZE_MAX - second_total
#ifdef LONGINT
            || first_count > (size_t)LONG_MAX
            || second_count > (size_t)LONG_MAX
#else
            || first_count > (size_t)INT_MAX
            || second_count > (size_t)INT_MAX
#endif
            ) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": native neighbor count overflow");
            goto cleanup;
        }
        schedule->tasks[task].first_offset = first_total;
        schedule->tasks[task].first_count = (INTEGER)first_count;
        schedule->tasks[task].second_offset = second_total;
        schedule->tasks[task].second_count = (INTEGER)second_count;
        first_total += first_count;
        second_total += second_count;
    }
    if ((first_total > 0
         && first_total > SIZE_MAX/sizeof(*schedule->first_nodes))
        || (second_total > 0
            && second_total > SIZE_MAX/sizeof(*schedule->second_nodes))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": native neighbor allocation overflow");
        goto cleanup;
    }
    if (first_total > 0)
        schedule->first_nodes = malloc(
            first_total*sizeof(*schedule->first_nodes));
    if (second_total > 0)
        schedule->second_nodes = malloc(
            second_total*sizeof(*schedule->second_nodes));
    if ((first_total > 0 && schedule->first_nodes == NULL)
        || (second_total > 0 && schedule->second_nodes == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": native neighbor allocation failed");
        goto cleanup;
    }
    for (INTEGER task = 0; task < count; task++) {
        if (schedule->tasks[task].first_count > 0)
            shear_native_fill_neighbors(
                &context, schedule->tasks[task].pivot_node,
                first_frontier, first_frontier_count,
                schedule->first_nodes + schedule->tasks[task].first_offset);
        if (schedule->tasks[task].second_count > 0)
            shear_native_fill_neighbors(
                &context, schedule->tasks[task].pivot_node,
                second_frontier, second_frontier_count,
                schedule->second_nodes + schedule->tasks[task].second_offset);
    }
    schedule->task_count = count;
    schedule->pivot_count = pivot_total;
    schedule->first_node_count = first_total;
    schedule->second_node_count = second_total;
    status = SUCCESS;

cleanup:
    free(second_frontier);
    free(first_frontier);
    if (status == FAILURE) shear_native_release_frontier_schedule(schedule);
    return status;
}

static INTEGER shear_native_pivot_task_target(
        int threads, int ranks, size_t accumulator_stride)
{
    size_t target = (size_t)SHEAR_NATIVE_PIVOT_TASK_TARGET;
    size_t bytes_per_task;
    size_t memory_target;

    (void)threads;
    (void)ranks;
    if (accumulator_stride == 0
        || accumulator_stride > SIZE_MAX/sizeof(real))
        return 1;
    bytes_per_task = accumulator_stride*sizeof(real);
    memory_target = SHEAR_NATIVE_TASK_HISTOGRAM_MEMORY/bytes_per_task;
    target = MIN(target, MAX((size_t)1, memory_target));
    return (INTEGER)MAX((size_t)1, target);
}

#undef SHEAR_NATIVE_PIVOT_TASK_TARGET
#undef SHEAR_NATIVE_NEIGHBOR_FRONTIER_TARGET
#undef SHEAR_NATIVE_TASK_HISTOGRAM_MEMORY
#endif
