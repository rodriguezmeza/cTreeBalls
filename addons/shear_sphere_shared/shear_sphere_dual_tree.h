#ifndef CTREEBALLS_SHEAR_SPHERE_DUAL_TREE_H
#define CTREEBALLS_SHEAR_SPHERE_DUAL_TREE_H

/* This file is included by search_octree_shear_omp.c after its spin helpers. */

#define SHEAR_SPHERE_PAIR_MAX_FRONTIER ((INTEGER)256)
#define SHEAR_SPHERE_PAIR_FRONTIER ((INTEGER)256)
#define SHEAR_SPHERE_DUAL_NODE_SPLIT_FACTOR ((real)0.585)

typedef struct {
    real *xi_plus_re;
    real *xi_plus_im;
    real *xi_minus_re;
    real *xi_minus_im;
    real *weight;
    INTEGER body_pairs;
    INTEGER cell_pairs;
    shear_profile_counters *profile;
} shear_sphere_pair_histogram;

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    int bins;
} shear_sphere_pair_context;

static bool shear_sphere_live_node(const shear_sphere_pair_context *context,
                                   nodeptr q)
{
    if (q == NULL)
        return FALSE;
    if (cballs_opt_read_mask(context->cmd)
        && Mask(q) == MASK_NODE_MASKED)
        return FALSE;
    if (Type(q) == CELL)
        return Nb(q) > 0;
    return Type(q) == BODY || Type(q) == BODY3;
}

static int shear_sphere_children(const shear_sphere_pair_context *context,
                                 nodeptr q, nodeptr children[NSUB])
{
    int count = 0;
    int index;

    if (q == NULL || Type(q) != CELL)
        return 0;
    for (index = 0; index < NSUB; index++) {
        nodeptr child = Subp(q)[index];

        if (shear_sphere_live_node(context, child))
            children[count++] = child;
    }
    return count;
}

static nodeptr shear_sphere_collapse_single_child(
        const shear_sphere_pair_context *context, nodeptr q)
{
    while (q != NULL && Type(q) == CELL) {
        nodeptr children[NSUB];
        const int child_count = shear_sphere_children(context, q, children);

        if (child_count != 1)
            break;
        q = children[0];
    }
    return q;
}

static bool shear_sphere_pair_node_frame(
        const shear_sphere_pair_context *context, nodeptr q,
        compute_vector unit, compute_vector east, compute_vector north,
        real *radius, real *angle)
{
    const real qsize = Type(q) == CELL ? (real)Size(q) : 0.0;

    return shear_spherical_node_frame(context->cmd, q, qsize,
                                      unit, east, north, radius, angle);
}

static real shear_sphere_node_weight(nodeptr q)
{
    return Type(q) == CELL ? ShearWeightSum(q) : Weight(q);
}

static shear_complex shear_sphere_node_weighted_gamma(nodeptr q)
{
    const real weight = shear_sphere_node_weight(q);

    return shear_scale(shear_make(Gamma1(q), Gamma2(q)), weight);
}

static bool shear_sphere_pair_phase(
        const real *pivot_unit, const real *pivot_east,
        const real *pivot_north, const real *neighbor_unit,
        shear_complex *phase)
{
    compute_vector tangent;
    real projection = shear_dot3(pivot_unit, neighbor_unit);
    real norm;

    tangent[0] = neighbor_unit[0] - projection*pivot_unit[0];
    tangent[1] = neighbor_unit[1] - projection*pivot_unit[1];
    tangent[2] = neighbor_unit[2] - projection*pivot_unit[2];
    norm = rsqrt(shear_dot3(tangent, tangent));
    if (!(norm > 64.0*DBL_EPSILON) || !isfinite(norm))
        return FALSE;
    phase->re = shear_dot3(tangent, pivot_east)/norm;
    phase->im = shear_dot3(tangent, pivot_north)/norm;
    return isfinite(phase->re) && isfinite(phase->im);
}

/* Return -2 outside the histogram, -1 when the pair must split, or a bin. */
static int shear_sphere_pair_bin(
        const shear_sphere_pair_context *context, nodeptr first,
        nodeptr second, bool permit_cells,
        shear_sphere_pair_histogram *hist)
{
    compute_vector first_unit;
    compute_vector first_east;
    compute_vector first_north;
    compute_vector second_unit;
    compute_vector second_east;
    compute_vector second_north;
    compute_vector difference;
    real first_radius;
    real second_radius;
    real first_angle;
    real second_angle;
    real distance2;
    real distance;
    real size;
#ifndef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    real lower;
    real upper;
#endif
    real tolerance;
    real separation;
    real error;
    int bin;

    if (!shear_sphere_pair_node_frame(
            context, first, first_unit, first_east, first_north,
            &first_radius, &first_angle)
        || !shear_sphere_pair_node_frame(
            context, second, second_unit, second_east, second_north,
            &second_radius, &second_angle))
        return -1;
    (void)first_east;
    (void)first_north;
    (void)second_east;
    (void)second_north;
    DOTPSUBV(distance2, difference, first_unit, second_unit);
    if (!(distance2 > 0.0) || !isfinite(distance2))
        return -2;
    distance = rsqrt(distance2);
    size = first_radius + second_radius;
    if (distance + size <= context->cmd->rminHist
        || distance - size >= context->cmd->rangeN)
        return -2;
    if (hist->profile == NULL) {
        bin = shear_radial_bin(context->cmd, context->gd, distance);
    } else {
        double radial_started = 0.0;
        bool radial_active = shear_profile_sample_begin(
            &hist->profile->radial, &radial_started);
        bin = shear_radial_bin(context->cmd, context->gd, distance);
        shear_profile_sample_end(
            &hist->profile->radial, radial_started, radial_active);
    }
    if (Type(first) != CELL && Type(second) != CELL)
        return bin < 0 ? -2 : bin;
    if (!permit_cells || !(context->cmd->theta > 0.0)
#ifdef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
        || cballs_opt_no_two_balls(context->cmd)
        || cballs_opt_no_one_ball(context->cmd)
#endif
        || (cballs_opt_read_mask(context->cmd)
            && ((Type(first) == CELL && Mask(first) != MASK_NODE_VALID)
                || (Type(second) == CELL
                    && Mask(second) != MASK_NODE_VALID)))
        || (Type(first) == CELL
            && !isfinite(ShearTransportError(first)))
        || (Type(second) == CELL
            && !isfinite(ShearTransportError(second))))
        return -1;
#ifdef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    if (bin < 0)
        /* The center can lie outside while the node-radius interval still
         * overlaps the histogram. Split until the descendants can be
         * classified; only the interval test above is allowed to prune. */
        return -1;
    if (context->cmd->useLogHist) {
        const real logarithmic_bin_width = context->cmd->rminHist > 0.0
            ? rlog(10.0)*context->gd->deltaR
            : rlog(10.0)/(real)context->cmd->logHistBinsPD;

        if (!(logarithmic_bin_width > 0.0)
            || size/distance
                 > context->cmd->theta*logarithmic_bin_width)
            return -1;
    } else if (!(context->gd->deltaR > 0.0)
               || size > context->cmd->theta*context->gd->deltaR) {
        return -1;
    }
#else
    lower = distance - size;
    upper = distance + size;
    if (!(lower > context->cmd->rminHist
          && upper < context->cmd->rangeN)
        || bin < 0
        || shear_radial_bin(context->cmd, context->gd, lower) != bin
        || shear_radial_bin(context->cmd, context->gd, upper) != bin)
        return -1;
#endif

    /* xi-minus carries phase order four. The same tolerance controls the
     * center-bearing approximation and spherical transport holonomy. */
    tolerance = MIN(0.5*PI, context->cmd->theta*PI/9.0);
    if (size/distance > rsin(MAX(0.0, tolerance)))
        return -1;
    separation = 2.0*rasin(MIN(1.0, 0.5*distance));
    error = (Type(first) == CELL ? ShearTransportError(first) : 0.0)
          + (Type(second) == CELL ? ShearTransportError(second) : 0.0)
          + 2.0*(first_angle + second_angle)*separation;
    return isfinite(error) && error <= tolerance ? bin : -1;
}

typedef struct {
    compute_vector unit;
    compute_vector east;
    compute_vector north;
    shear_complex weighted_gamma;
    real weight;
} shear_sphere_pair_frame;

static bool shear_sphere_prepare_pair_frame(
        const shear_sphere_pair_context *context, nodeptr node,
        shear_sphere_pair_frame *frame)
{
    real radius;
    real angle;

    if (!shear_sphere_pair_node_frame(
            context, node, frame->unit, frame->east, frame->north,
            &radius, &angle))
        return FALSE;
    frame->weighted_gamma = shear_sphere_node_weighted_gamma(node);
    frame->weight = shear_sphere_node_weight(node);
    return TRUE;
}

static void shear_sphere_accumulate_oriented_frames(
        shear_sphere_pair_histogram *hist,
        const shear_sphere_pair_frame *pivot,
        const shear_sphere_pair_frame *neighbor,
        shear_complex rotation, int bin)
{
    shear_complex phase;
    shear_complex neighbor_gamma;
    shear_complex z2;
    shear_complex z4;
    shear_complex value;

    if (!shear_sphere_pair_phase(
            pivot->unit, pivot->east, pivot->north, neighbor->unit, &phase))
        return;
    neighbor_gamma = shear_mul(neighbor->weighted_gamma, rotation);
    value = shear_mul(pivot->weighted_gamma, shear_conj(neighbor_gamma));
    hist->xi_plus_re[bin] += value.re;
    hist->xi_plus_im[bin] += value.im;
    z2 = shear_mul(phase, phase);
    z4 = shear_mul(z2, z2);
    value = shear_mul(shear_mul(pivot->weighted_gamma, neighbor_gamma),
                      shear_conj(z4));
    hist->xi_minus_re[bin] += value.re;
    hist->xi_minus_im[bin] += value.im;
    hist->weight[bin] += pivot->weight*neighbor->weight;
}

static void shear_sphere_accumulate_pair(
        const shear_sphere_pair_context *context,
        shear_sphere_pair_histogram *hist, nodeptr first, nodeptr second,
        int bin, bool bidirectional)
{
    shear_sphere_pair_frame first_frame;
    shear_sphere_pair_frame second_frame;
    shear_complex rotation;
    bool rotation_valid;

    if (bin < 0 || bin >= context->bins
        || !shear_sphere_prepare_pair_frame(context, first, &first_frame)
        || !shear_sphere_prepare_pair_frame(context, second, &second_frame))
        return;
    if (hist->profile == NULL) {
        rotation_valid = shear_transport_rotation_between_frames(
            first_frame.unit, first_frame.east, first_frame.north,
            second_frame.unit, second_frame.east, &rotation);
    } else {
        double transport_started = 0.0;
        bool transport_active = shear_profile_sample_begin(
            &hist->profile->transport, &transport_started);
        rotation_valid = shear_transport_rotation_between_frames(
            first_frame.unit, first_frame.east, first_frame.north,
            second_frame.unit, second_frame.east, &rotation);
        shear_profile_sample_end(
            &hist->profile->transport, transport_started, transport_active);
    }
    if (!rotation_valid) return;
    shear_sphere_accumulate_oriented_frames(
        hist, &first_frame, &second_frame, rotation, bin);
    if (bidirectional)
        shear_sphere_accumulate_oriented_frames(
            hist, &second_frame, &first_frame, shear_conj(rotation), bin);
}

static void shear_sphere_process_pair(
        const shear_sphere_pair_context *context,
        shear_sphere_pair_histogram *hist, nodeptr first, nodeptr second,
        bool bidirectional)
{
#ifndef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    nodeptr children[NSUB];
#endif
    int bin;
#ifndef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    int child_count;
    int child;
#endif

    if (!shear_sphere_live_node(context, first)
        || !shear_sphere_live_node(context, second))
        return;
    bin = shear_sphere_pair_bin(context, first, second, TRUE, hist);
    if (bin == -2)
        return;
    if (bin >= 0) {
        shear_sphere_accumulate_pair(
            context, hist, first, second, bin, bidirectional);
        if (Type(first) == CELL || Type(second) == CELL)
            hist->cell_pairs++;
        else
            hist->body_pairs++;
        return;
    }

    if (Type(first) != CELL && Type(second) != CELL)
        return;
#ifdef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    {
        bool split_first = Type(first) == CELL;
        bool split_second = Type(second) == CELL;
        nodeptr first_children[NSUB];
        nodeptr second_children[NSUB];
        int first_count = 0;
        int second_count = 0;

        if (split_first && split_second) {
            compute_vector first_unit;
            compute_vector second_unit;
            compute_vector east;
            compute_vector north;
            compute_vector delta;
            real first_radius = 0.0;
            real second_radius = 0.0;
            real ignored_angle;
            real distance2;
            real distance;
            real effective_width;

            if (!shear_sphere_pair_node_frame(
                    context, first, first_unit, east, north,
                    &first_radius, &ignored_angle)
                || !shear_sphere_pair_node_frame(
                    context, second, second_unit, east, north,
                    &second_radius, &ignored_angle))
                return;
            DOTPSUBV(distance2, delta, first_unit, second_unit);
            distance = distance2 > 0.0 ? rsqrt(distance2) : 0.0;
            effective_width = context->cmd->useLogHist
                ? context->cmd->theta
                    * (context->cmd->rminHist > 0.0
                       ? rlog(10.0)*context->gd->deltaR
                       : rlog(10.0)/(real)context->cmd->logHistBinsPD)
                    * distance
                : context->cmd->theta*context->gd->deltaR;

            if (second_radius > first_radius) {
                split_first = !(second_radius > 2.0*first_radius)
                    && first_radius
                         > SHEAR_SPHERE_DUAL_NODE_SPLIT_FACTOR*effective_width;
                split_second = TRUE;
            } else {
                split_first = TRUE;
                split_second = !(first_radius > 2.0*second_radius)
                    && second_radius
                         > SHEAR_SPHERE_DUAL_NODE_SPLIT_FACTOR*effective_width;
            }
        }
        if (split_first)
            first_count = shear_sphere_children(
                context, first, first_children);
        if (split_second)
            second_count = shear_sphere_children(
                context, second, second_children);
        if (split_first && split_second) {
            for (int first_child = 0; first_child < first_count; first_child++)
                for (int second_child = 0;
                     second_child < second_count; second_child++)
                    shear_sphere_process_pair(
                        context, hist, first_children[first_child],
                        second_children[second_child], bidirectional);
        } else if (split_first) {
            for (int child_index = 0; child_index < first_count; child_index++)
                shear_sphere_process_pair(
                    context, hist, first_children[child_index], second,
                    bidirectional);
        } else {
            for (int child_index = 0; child_index < second_count; child_index++)
                shear_sphere_process_pair(
                    context, hist, first, second_children[child_index],
                    bidirectional);
        }
        return;
    }
#else
    if (Type(second) != CELL
        || (Type(first) == CELL
            && ((real)Size(first) > (real)Size(second)
                || ((real)Size(first) == (real)Size(second)
                    && Nb(first) >= Nb(second))))) {
        child_count = shear_sphere_children(context, first, children);
        for (child = 0; child < child_count; child++)
            shear_sphere_process_pair(
                context, hist, children[child], second, bidirectional);
    } else {
        child_count = shear_sphere_children(context, second, children);
        for (child = 0; child < child_count; child++)
            shear_sphere_process_pair(
                context, hist, first, children[child], bidirectional);
    }
#endif
}

static void shear_sphere_process_auto(
        const shear_sphere_pair_context *context,
        shear_sphere_pair_histogram *hist, nodeptr q)
{
    nodeptr children[NSUB];
    int child_count;
    int first;
    int second;

    if (q == NULL || Type(q) != CELL)
        return;
    child_count = shear_sphere_children(context, q, children);
    for (first = 0; first < child_count; first++) {
        shear_sphere_process_auto(context, hist, children[first]);
        for (second = first + 1; second < child_count; second++)
            shear_sphere_process_pair(
                context, hist, children[first], children[second], TRUE);
    }
}

static INTEGER shear_sphere_frontier_target(struct cmdline_data *cmd)
{
    (void)cmd;
    /* A thread-count-independent decomposition is part of the numerical
     * contract: changing OpenMP workers must not change accepted node pairs. */
    return MIN(SHEAR_SPHERE_PAIR_MAX_FRONTIER,
               SHEAR_SPHERE_PAIR_FRONTIER);
}

#define SHEAR_SPHERE_PAIR_TASK_CHUNK ((INTEGER)64)

static int shear_sphere_build_frontier(
        const shear_sphere_pair_context *context, nodeptr root,
        nodeptr **result, INTEGER *result_count)
{
    const INTEGER target = shear_sphere_frontier_target(context->cmd);
    const INTEGER capacity = target + NSUB;
    nodeptr *frontier = NULL;
    INTEGER count = 1;

    *result = NULL;
    *result_count = 0;
    if ((size_t)capacity > SIZE_MAX/sizeof(*frontier)
        || (frontier = malloc((size_t)capacity*sizeof(*frontier))) == NULL) {
        snprintf(context->cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree frontier allocation failed");
        return FAILURE;
    }
    frontier[0] = shear_sphere_collapse_single_child(context, root);
    while (count < target) {
        nodeptr children[NSUB];
        INTEGER selected = -1;
        INTEGER largest = -1;
        int child_count = 0;
        INTEGER index;

        for (index = 0; index < count; index++) {
            int available = shear_sphere_children(
                context, frontier[index], children);
            if (available >= 2 && Nb(frontier[index]) > largest) {
                selected = index;
                largest = Nb(frontier[index]);
                child_count = available;
            }
        }
        if (selected < 0)
            break;
        child_count = shear_sphere_children(
            context, frontier[selected], children);
        frontier[selected] = shear_sphere_collapse_single_child(
            context, children[0]);
        for (int child = 1; child < child_count; child++)
            frontier[count++] = shear_sphere_collapse_single_child(
                context, children[child]);
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}

static bool shear_sphere_dual_tree_eligible(
        const shear_sphere_pair_context *context, bodyptr pivots,
        INTEGER nbody, INTEGER ipmin, INTEGER ipmax)
{
    INTEGER index;

    if (!(context->cmd->theta > 0.0)
#ifndef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
        || cballs_opt_no_one_ball(context->cmd)
#endif
        || ipmin != 1 || ipmax != nbody)
        return FALSE;
    for (index = 0; index < nbody; index++) {
        bodyptr pivot = pivots + index;

        if ((!cballs_opt_read_mask(context->cmd)
             || Mask(pivot) == MASK_NODE_VALID)
            && !Update(pivot))
            return FALSE;
    }
    return TRUE;
}

static INTEGER shear_sphere_pair_point_count(nodeptr node)
{
    return Type(node) == CELL ? Nb(node) : 1;
}

static bool shear_sphere_pair_frontier_overlap(
        const shear_sphere_pair_context *context,
        nodeptr first, nodeptr second)
{
    compute_vector first_unit;
    compute_vector first_east;
    compute_vector first_north;
    compute_vector second_unit;
    compute_vector second_east;
    compute_vector second_north;
    compute_vector delta;
    real first_radius;
    real second_radius;
    real first_angle;
    real second_angle;
    real distance2;
    real distance;

    if (!shear_sphere_pair_node_frame(
            context, first, first_unit, first_east, first_north,
            &first_radius, &first_angle)
        || !shear_sphere_pair_node_frame(
            context, second, second_unit, second_east, second_north,
            &second_radius, &second_angle))
        return TRUE;
    DOTPSUBV(distance2, delta, first_unit, second_unit);
    if (!(distance2 >= 0.0) || !isfinite(distance2)) return TRUE;
    distance = rsqrt(distance2);
    return distance + first_radius + second_radius > context->cmd->rminHist
        && distance - first_radius - second_radius < context->cmd->rangeN;
}

static long double shear_sphere_pair_task_work(
        const shear_sphere_pair_context *context,
        nodeptr *frontier1, nodeptr *frontier2,
        INTEGER first_index, INTEGER second_index, INTEGER last_index,
        bool auto_correlation)
{
    const long double first_points =
        (long double)shear_sphere_pair_point_count(frontier1[first_index]);
    long double work = 0.0L;

    if (auto_correlation && first_index == second_index)
        return 0.5L*first_points*MAX(0.0L, first_points - 1.0L);
    for (INTEGER second = second_index; second <= last_index; second++)
        if (shear_sphere_pair_frontier_overlap(
                context, frontier1[first_index], frontier2[second]))
            work += first_points
                * (long double)shear_sphere_pair_point_count(
                    frontier2[second]);
    return work;
}

static int shear_sphere_dual_tree_2pcf(
        struct cmdline_data *cmd, struct global_data *gd,
        nodeptr root1, nodeptr root2, bool auto_correlation, int bins)
{
    shear_sphere_pair_context context;
    nodeptr *frontier1 = NULL;
    nodeptr *frontier2 = NULL;
    INTEGER count1 = 0;
    INTEGER count2 = 0;
    INTEGER task_count = 0;
    INTEGER local_task_count = 0;
    INTEGER *task_first = NULL;
    INTEGER *task_second = NULL;
    INTEGER *task_last = NULL;
    INTEGER *task_local_index = NULL;
    INTEGER *local_tasks = NULL;
    int *task_owner = NULL;
    long double *task_work = NULL;
    long double *rank_work = NULL;
    real *task_storage = NULL;
    INTEGER *task_body_pairs = NULL;
    INTEGER *task_cell_pairs = NULL;
    shear_profile_counters *profiles = NULL;
    size_t stride;
    size_t values;
    int threads = 1;
    int status = FAILURE;

#ifdef OPENMPCODE
    threads = omp_get_max_threads();
#endif
    if (threads < 1) threads = 1;
    context.cmd = cmd;
    context.gd = gd;
    context.bins = bins;
    if (shear_sphere_build_frontier(
            &context, root1, &frontier1, &count1) == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        frontier2 = frontier1;
        count2 = count1;
    } else if (shear_sphere_build_frontier(
                   &context, root2, &frontier2, &count2) == FAILURE) {
        goto cleanup;
    }
    if (count1 <= 0 || count2 <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": empty dual-tree frontier");
        goto cleanup;
    }
    for (INTEGER first = 0; first < count1; first++) {
        const INTEGER pair_count = auto_correlation
            ? count2 - first - 1 : count2;
        const INTEGER chunks = pair_count > 0
            ? 1 + (pair_count - 1)/SHEAR_SPHERE_PAIR_TASK_CHUNK : 0;
        task_count += chunks + (auto_correlation ? 1 : 0);
    }
    stride = (size_t)bins;
    if (task_count <= 0 || stride == 0 || stride > SIZE_MAX/5
        || (size_t)task_count > SIZE_MAX/sizeof(*task_first)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result size overflow");
        goto cleanup;
    }
    task_first = malloc((size_t)task_count*sizeof(*task_first));
    task_second = malloc((size_t)task_count*sizeof(*task_second));
    task_last = malloc((size_t)task_count*sizeof(*task_last));
    task_local_index = malloc(
        (size_t)task_count*sizeof(*task_local_index));
    if (task_first == NULL || task_second == NULL || task_last == NULL
        || task_local_index == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree task allocation failed");
        goto cleanup;
    }
    {
        INTEGER task = 0;

        for (INTEGER first = 0; first < count1; first++) {
            INTEGER second = auto_correlation ? first + 1 : 0;

            if (auto_correlation) {
                task_first[task] = first;
                task_second[task] = first;
                task_last[task] = first;
                task++;
            }
            while (second < count2) {
                task_first[task] = first;
                task_second[task] = second;
                task_last[task] = MIN(
                    count2 - 1,
                    second + SHEAR_SPHERE_PAIR_TASK_CHUNK - 1);
                second = task_last[task] + 1;
                task++;
            }
        }
#ifdef SHEAR_MPI_ENABLED
        const int rank_count = MAX(1, SHEAR_MPI_SIZE());
        const int rank = SHEAR_MPI_RANK();

        task_owner = malloc((size_t)task_count*sizeof(*task_owner));
        task_work = malloc((size_t)task_count*sizeof(*task_work));
        rank_work = calloc((size_t)rank_count, sizeof(*rank_work));
        if (task_owner == NULL || task_work == NULL || rank_work == NULL) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME
                     ": dual-tree rank scheduling allocation failed");
            goto cleanup;
        }
        for (task = 0; task < task_count; task++) {
            task_owner[task] = -1;
            task_work[task] = shear_sphere_pair_task_work(
                &context, frontier1, frontier2,
                task_first[task], task_second[task], task_last[task],
                auto_correlation);
        }
        for (INTEGER assigned = 0; assigned < task_count; assigned++) {
            INTEGER largest = -1;
            int owner = 0;

            for (task = 0; task < task_count; task++)
                if (task_owner[task] < 0
                    && (largest < 0 || task_work[task] > task_work[largest]))
                    largest = task;
            for (int candidate = 1; candidate < rank_count; candidate++)
                if (rank_work[candidate] < rank_work[owner]) owner = candidate;
            task_owner[largest] = owner;
            rank_work[owner] += task_work[largest];
        }
        for (task = 0; task < task_count; task++) {
            if (task_owner[task] == rank)
                task_local_index[task] = local_task_count++;
            else
                task_local_index[task] = -1;
        }
#else
        for (task = 0; task < task_count; task++)
            task_local_index[task] = local_task_count++;
#endif
    }
    if ((size_t)MAX((INTEGER)1, local_task_count) > SIZE_MAX/(5*stride)
        || (values = (size_t)MAX((INTEGER)1, local_task_count)*5*stride)
            > SIZE_MAX/sizeof(*task_storage)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree local result size overflow");
        goto cleanup;
    }
    task_storage = calloc(values, sizeof(*task_storage));
    local_tasks = malloc(
        (size_t)MAX((INTEGER)1, local_task_count)*sizeof(*local_tasks));
    task_body_pairs = calloc(
        (size_t)MAX((INTEGER)1, local_task_count),
        sizeof(*task_body_pairs));
    task_cell_pairs = calloc(
        (size_t)MAX((INTEGER)1, local_task_count),
        sizeof(*task_cell_pairs));
    if (getenv("CBALLS_SHEAR_PROFILE") != NULL)
        profiles = calloc((size_t)threads, sizeof(*profiles));
    if (task_storage == NULL || local_tasks == NULL
        || task_body_pairs == NULL || task_cell_pairs == NULL
        || (getenv("CBALLS_SHEAR_PROFILE") != NULL && profiles == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result allocation failed");
        goto cleanup;
    }
    for (INTEGER task = 0; task < task_count; task++)
        if (task_local_index[task] >= 0)
            local_tasks[task_local_index[task]] = task;

#pragma omp parallel
    {
        int thread = 0;
#ifdef OPENMPCODE
        thread = omp_get_thread_num();
#endif
        shear_profile_counters *profile = profiles != NULL
            ? &profiles[thread] : NULL;
        const double started = profile != NULL
            ? shear_profile_wall_time() : 0.0;

#pragma omp for schedule(dynamic,1) nowait
        for (INTEGER local_task = 0;
             local_task < local_task_count; local_task++) {
            const INTEGER task = local_tasks[local_task];
            real *base = task_storage + (size_t)local_task*5*stride;
            shear_sphere_pair_histogram hist = {
                base, base + stride, base + 2*stride, base + 3*stride,
                base + 4*stride, 0, 0, profile
            };
            const INTEGER first_index = task_first[task];
            const INTEGER second_index = task_second[task];
            const INTEGER last_index = task_last[task];

            if (auto_correlation) {
                if (first_index == second_index)
                    shear_sphere_process_auto(
                        &context, &hist, frontier1[first_index]);
                else for (INTEGER other = second_index;
                          other <= last_index; other++)
                        shear_sphere_process_pair(
                            &context, &hist, frontier1[first_index],
                            frontier2[other], TRUE);
            } else {
                for (INTEGER other = second_index;
                     other <= last_index; other++)
                    shear_sphere_process_pair(
                        &context, &hist, frontier1[first_index],
                        frontier2[other], FALSE);
            }
            task_body_pairs[local_task] = hist.body_pairs;
            task_cell_pairs[local_task] = hist.cell_pairs;
            if (profile != NULL) profile->pivots++;
        }
        if (profile != NULL) {
            profile->elapsed_seconds = shear_profile_wall_time() - started;
            profile->walk_seconds = profile->elapsed_seconds;
        }
    }
    shear_print_profiles(profiles, threads);

    for (INTEGER task = 0; task < task_count; task++) {
        const INTEGER local_task = task_local_index[task];
        const real *base;

        if (local_task < 0) continue;
        base = task_storage + (size_t)local_task*5*stride;
        for (int bin = 0; bin < bins; bin++) {
            gd->histShearXiPlusRe[bin] += base[bin];
            gd->histShearXiPlusIm[bin] += base[stride + (size_t)bin];
            gd->histShearXiMinusRe[bin] += base[2*stride + (size_t)bin];
            gd->histShearXiMinusIm[bin] += base[3*stride + (size_t)bin];
            gd->histShearXiWeight[bin] += base[4*stride + (size_t)bin];
        }
        gd->nbbcalc += task_body_pairs[local_task];
        gd->nbccalc += task_cell_pairs[local_task];
    }
    verb_print(cmd->verbose,
               SHEAR_ENGINE_NAME ": spherical dual tree used %"
               INTEGER_FMT " tasks, %" INTEGER_FMT " body pairs and %"
               INTEGER_FMT " accepted node pairs\n",
               task_count, gd->nbbcalc, gd->nbccalc);
    status = SUCCESS;

cleanup:
    free(profiles);
    free(task_cell_pairs);
    free(task_body_pairs);
    free(local_tasks);
    free(task_storage);
    free(rank_work);
    free(task_work);
    free(task_owner);
    free(task_local_index);
    free(task_last);
    free(task_second);
    free(task_first);
    if (!auto_correlation)
        free(frontier2);
    free(frontier1);
    return status;
}

#undef SHEAR_SPHERE_PAIR_MAX_FRONTIER
#undef SHEAR_SPHERE_PAIR_FRONTIER
#undef SHEAR_SPHERE_PAIR_TASK_CHUNK
#undef SHEAR_SPHERE_DUAL_NODE_SPLIT_FACTOR
#endif
