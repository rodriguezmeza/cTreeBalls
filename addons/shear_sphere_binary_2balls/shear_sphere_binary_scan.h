#ifndef CTREEBALLS_SHEAR_SPHERE_BINARY_SCAN_H
#define CTREEBALLS_SHEAR_SPHERE_BINARY_SCAN_H

/* Shared KD/ball-tree scan, included after the spin-2 helper routines. */

#ifndef SHEAR_SPHERE_BINARY_TREE_HEADER
#define SHEAR_SPHERE_BINARY_TREE_HEADER "kdtree_shear_sphere_tree.h"
#endif
#include SHEAR_SPHERE_BINARY_TREE_HEADER

#ifndef SHEAR_SPHERE_BINARY_TREE_FRONTIER
#define SHEAR_SPHERE_BINARY_TREE_FRONTIER kdtree_shear_sphere_frontier
#endif
#ifndef SHEAR_SPHERE_BINARY_TREE_FREE
#define SHEAR_SPHERE_BINARY_TREE_FREE kdtree_shear_sphere_free
#endif

#define KD_SHEAR_PAIR_FRONTIER ((INTEGER)256)
#define KD_SHEAR_PAIR_TASK_CHUNK ((INTEGER)64)
#define KD_SHEAR_SPLIT_FACTOR ((real)0.585)

typedef enum {
    KD_SHEAR_NODE,
    KD_SHEAR_BODY
} kd_shear_ref_kind;

typedef struct {
    kd_shear_ref_kind kind;
    fcfc_balltreeptr tree;
    INTEGER index;
    bodyptr body;
} kd_shear_ref;

typedef struct {
    real *xi_plus_re;
    real *xi_plus_im;
    real *xi_minus_re;
    real *xi_minus_im;
    real *weight;
    INTEGER body_pairs;
    INTEGER cell_pairs;
    shear_profile_counters *profile;
} kd_shear_pair_histogram;

static kd_shear_ref kd_shear_node_ref(fcfc_balltreeptr tree, INTEGER index)
{
    kd_shear_ref ref = {KD_SHEAR_NODE, tree, index, NULL};
    return ref;
}

static kd_shear_ref kd_shear_body_ref(fcfc_balltreeptr tree, bodyptr body)
{
    kd_shear_ref ref = {KD_SHEAR_BODY, tree, -1, body};
    return ref;
}

static fcfc_ballnode *kd_shear_node(kd_shear_ref ref)
{
    return &ref.tree->nodes[ref.index];
}

static const real *kd_shear_position(kd_shear_ref ref)
{
    return ref.kind == KD_SHEAR_NODE
        ? kd_shear_node(ref)->center : Pos(ref.body);
}

static real kd_shear_radius(kd_shear_ref ref)
{
    return ref.kind == KD_SHEAR_NODE
        ? (real)kd_shear_node(ref)->radius : 0.0;
}

static real kd_shear_weight(kd_shear_ref ref)
{
    return ref.kind == KD_SHEAR_NODE
        ? kd_shear_node(ref)->weight : Weight(ref.body);
}

static shear_complex kd_shear_weighted_gamma(kd_shear_ref ref)
{
    if (ref.kind == KD_SHEAR_NODE)
        return shear_make(kd_shear_node(ref)->shear_gamma_re,
                          kd_shear_node(ref)->shear_gamma_im);
    return shear_scale(shear_make(Gamma1(ref.body), Gamma2(ref.body)),
                       Weight(ref.body));
}

static real kd_shear_transport_error(kd_shear_ref ref)
{
    return ref.kind == KD_SHEAR_NODE
        ? kd_shear_node(ref)->shear_transport_error : 0.0;
}

static INTEGER kd_shear_child_count(kd_shear_ref parent)
{
    fcfc_ballnode *node;

    if (parent.kind != KD_SHEAR_NODE) return 0;
    node = kd_shear_node(parent);
    return node->left >= 0 && node->right >= 0
        ? 2 : node->last - node->first + 1;
}

static kd_shear_ref kd_shear_child(kd_shear_ref parent, INTEGER child)
{
    fcfc_ballnode *node = kd_shear_node(parent);

    if (node->left >= 0 && node->right >= 0)
        return kd_shear_node_ref(parent.tree,
                                 child == 0 ? node->left : node->right);
    return kd_shear_body_ref(parent.tree,
                             parent.tree->bptr[node->first + child]);
}

/* Return -2 outside the histogram, -1 when subdivision is required, or a bin. */
static int kd_shear_pair_bin(struct cmdline_data *cmd,
                             struct global_data *gd,
                             kd_shear_pair_histogram *hist,
                             kd_shear_ref first, kd_shear_ref second,
                             bool permit_nodes)
{
    compute_vector first_unit;
    compute_vector second_unit;
    compute_vector delta;
    real distance2;
    real distance;
    real size;
    real lower;
    real upper;
    real tolerance;
    real separation;
    real error;
    int bin;

    if (!shear_load_unit3(kd_shear_position(first), first_unit)
        || !shear_load_unit3(kd_shear_position(second), second_unit))
        return -1;
    DOTPSUBV(distance2, delta, first_unit, second_unit);
    if (!(distance2 > 0.0) || !isfinite(distance2)) return -2;
    distance = rsqrt(distance2);
    size = kd_shear_radius(first) + kd_shear_radius(second);
    if (distance + size <= cmd->rminHist
        || distance - size >= cmd->rangeN)
        return -2;
    if (hist->profile == NULL) {
        bin = shear_radial_bin(cmd, gd, distance);
    } else {
        double radial_started = 0.0;
        bool radial_active = shear_profile_sample_begin(
            &hist->profile->radial, &radial_started);
        bin = shear_radial_bin(cmd, gd, distance);
        shear_profile_sample_end(
            &hist->profile->radial, radial_started, radial_active);
    }
    if (first.kind == KD_SHEAR_BODY && second.kind == KD_SHEAR_BODY)
        return bin < 0 ? -2 : bin;
    if (!permit_nodes || !(cmd->theta > 0.0)
        || cballs_opt_no_two_balls(cmd) || cballs_opt_no_one_ball(cmd)
        || !isfinite(kd_shear_transport_error(first))
        || !isfinite(kd_shear_transport_error(second)))
        return -1;
    if (bin < 0) return -1;
    lower = distance - size;
    upper = distance + size;
    if (scanopt(cmd->options, "dual-node-bin-slop")) {
        real width = cmd->useLogHist
            ? cmd->theta*(cmd->rminHist > 0.0
                ? rlog(10.0)*gd->deltaR
                : rlog(10.0)/(real)cmd->logHistBinsPD)*distance
            : cmd->theta*gd->deltaR;
        if (!(width > 0.0) || size > width) return -1;
    } else if (!shear_interval_within_radial_bin(
                   cmd, gd, lower, upper, bin)) {
        return -1;
    }
    tolerance = MIN(0.5*PI, cmd->theta*PI/9.0);
    if (distance <= size || size/distance > rsin(MAX(0.0, tolerance)))
        return -1;
    separation = 2.0*rasin(MIN(1.0, 0.5*distance));
    error = kd_shear_transport_error(first)
          + kd_shear_transport_error(second)
          + 2.0*(2.0*rasin(MIN(1.0, 0.5*kd_shear_radius(first)))
                 + 2.0*rasin(MIN(1.0, 0.5*kd_shear_radius(second))))
            * separation;
    return isfinite(error) && error <= tolerance ? bin : -1;
}

typedef struct {
    compute_vector unit;
    compute_vector east;
    compute_vector north;
    shear_complex weighted_gamma;
    real weight;
} kd_shear_pair_frame;

static bool kd_shear_pair_phase(const real *pivot_unit,
                                const real *pivot_east,
                                const real *pivot_north,
                                const real *neighbor_unit,
                                shear_complex *phase)
{
    const real projection = shear_dot3(pivot_unit, neighbor_unit);
    const real norm = rsqrt(MAX(0.0, 1.0 - projection*projection));

    if (!(norm > 64.0*DBL_EPSILON) || !isfinite(norm)) return FALSE;
    phase->re = shear_dot3(neighbor_unit, pivot_east)/norm;
    phase->im = shear_dot3(neighbor_unit, pivot_north)/norm;
    return isfinite(phase->re) && isfinite(phase->im);
}

static bool kd_shear_prepare_pair_frame(kd_shear_ref ref,
                                        kd_shear_pair_frame *frame)
{
    if (!shear_load_unit3(kd_shear_position(ref), frame->unit)
        || !shear_spherical_basis(frame->unit, frame->east, frame->north))
        return FALSE;
    frame->weighted_gamma = kd_shear_weighted_gamma(ref);
    frame->weight = kd_shear_weight(ref);
    return TRUE;
}

static void kd_shear_accumulate_oriented_frames(
        kd_shear_pair_histogram *hist, const kd_shear_pair_frame *pivot,
        const kd_shear_pair_frame *neighbor, shear_complex rotation, int bin)
{
    shear_complex phase;
    shear_complex neighbor_gamma;
    shear_complex z2;
    shear_complex z4;
    shear_complex value;

    if (!kd_shear_pair_phase(pivot->unit, pivot->east, pivot->north,
                             neighbor->unit, &phase))
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

static void kd_shear_accumulate_pair(kd_shear_pair_histogram *hist,
                                     kd_shear_ref first,
                                     kd_shear_ref second, int bin,
                                     bool bidirectional)
{
    kd_shear_pair_frame first_frame;
    kd_shear_pair_frame second_frame;
    shear_complex rotation;

    if (!kd_shear_prepare_pair_frame(first, &first_frame)
        || !kd_shear_prepare_pair_frame(second, &second_frame))
        return;
    bool rotation_valid;
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
    kd_shear_accumulate_oriented_frames(
        hist, &first_frame, &second_frame, rotation, bin);
    if (bidirectional)
        kd_shear_accumulate_oriented_frames(
            hist, &second_frame, &first_frame, shear_conj(rotation), bin);
}

static void kd_shear_process_pair(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  kd_shear_pair_histogram *hist,
                                  kd_shear_ref first, kd_shear_ref second,
                                  bool bidirectional)
{
    bool split_first = first.kind == KD_SHEAR_NODE;
    bool split_second = second.kind == KD_SHEAR_NODE;
    INTEGER first_count = 0;
    INTEGER second_count = 0;
    int bin = kd_shear_pair_bin(cmd, gd, hist, first, second, TRUE);

    if (bin == -2) return;
    if (bin >= 0) {
        kd_shear_accumulate_pair(hist, first, second, bin, bidirectional);
        if (first.kind == KD_SHEAR_NODE || second.kind == KD_SHEAR_NODE)
            hist->cell_pairs++;
        else
            hist->body_pairs++;
        return;
    }
    if (!split_first && !split_second) return;
    if (split_first && split_second) {
        const real first_radius = kd_shear_radius(first);
        const real second_radius = kd_shear_radius(second);
        compute_vector delta;
        real distance2;
        real distance;
        real effective_width;

        DOTPSUBV(distance2, delta, kd_shear_position(first),
                 kd_shear_position(second));
        distance = distance2 > 0.0 ? rsqrt(distance2) : 0.0;
        effective_width = cmd->useLogHist
            ? cmd->theta*(cmd->rminHist > 0.0
                ? rlog(10.0)*gd->deltaR
                : rlog(10.0)/(real)cmd->logHistBinsPD)*distance
            : cmd->theta*gd->deltaR;
        if (second_radius > first_radius) {
            split_first = !(second_radius > 2.0*first_radius)
                && first_radius > KD_SHEAR_SPLIT_FACTOR*effective_width;
            split_second = TRUE;
        } else {
            split_first = TRUE;
            split_second = !(first_radius > 2.0*second_radius)
                && second_radius > KD_SHEAR_SPLIT_FACTOR*effective_width;
        }
    }
    if (split_first)
        first_count = kd_shear_child_count(first);
    if (split_second)
        second_count = kd_shear_child_count(second);
    if (split_first && split_second) {
        for (INTEGER i = 0; i < first_count; i++)
            for (INTEGER j = 0; j < second_count; j++)
                kd_shear_process_pair(cmd, gd, hist, kd_shear_child(first, i),
                                      kd_shear_child(second, j), bidirectional);
    } else if (split_first) {
        for (INTEGER i = 0; i < first_count; i++)
            kd_shear_process_pair(cmd, gd, hist, kd_shear_child(first, i),
                                  second, bidirectional);
    } else {
        for (INTEGER j = 0; j < second_count; j++)
            kd_shear_process_pair(cmd, gd, hist, first,
                                  kd_shear_child(second, j), bidirectional);
    }
}

static void kd_shear_process_auto(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  kd_shear_pair_histogram *hist,
                                  kd_shear_ref parent)
{
    const INTEGER count = kd_shear_child_count(parent);

    for (INTEGER first = 0; first < count; first++) {
        kd_shear_ref first_child = kd_shear_child(parent, first);
        if (first_child.kind == KD_SHEAR_NODE)
            kd_shear_process_auto(cmd, gd, hist, first_child);
        for (INTEGER second = first + 1; second < count; second++)
            kd_shear_process_pair(cmd, gd, hist, first_child,
                                  kd_shear_child(parent, second), TRUE);
    }
}

static bool kd_shear_pair_frontier_overlap(
        const struct cmdline_data *cmd, const fcfc_ballnode *first,
        const fcfc_ballnode *second)
{
    compute_vector delta;
    real distance2;
    real distance;
    const real radius = (real)first->radius + (real)second->radius;

    DOTPSUBV(distance2, delta, first->center, second->center);
    if (!(distance2 >= 0.0) || !isfinite(distance2)) return TRUE;
    distance = rsqrt(distance2);
    return distance + radius > cmd->rminHist
        && distance - radius < cmd->rangeN;
}

static long double kd_shear_pair_task_work(
        const struct cmdline_data *cmd, fcfc_balltreeptr first_tree,
        fcfc_balltreeptr second_tree, const INTEGER *frontier1,
        const INTEGER *frontier2, INTEGER first_index,
        INTEGER second_index, INTEGER last_index, bool auto_correlation)
{
    const fcfc_ballnode *first = &first_tree->nodes[frontier1[first_index]];
    const long double first_points =
        (long double)(first->last - first->first + 1);
    long double work = 0.0L;

    if (auto_correlation && first_index == second_index)
        return 0.5L*first_points*MAX(0.0L, first_points - 1.0L);
    for (INTEGER second = second_index; second <= last_index; second++) {
        const fcfc_ballnode *other =
            &second_tree->nodes[frontier2[second]];

        if (kd_shear_pair_frontier_overlap(cmd, first, other))
            work += first_points
                * (long double)(other->last - other->first + 1);
    }
    return work;
}

static int kd_shear_dual_tree_2pcf(
        struct cmdline_data *cmd, struct global_data *gd,
        fcfc_balltreeptr first_tree, fcfc_balltreeptr second_tree,
        bool auto_correlation, int bins)
{
    INTEGER *frontier1 = NULL;
    INTEGER *frontier2 = NULL;
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
    INTEGER *body_pairs = NULL;
    INTEGER *cell_pairs = NULL;
    shear_profile_counters *profiles = NULL;
    const size_t stride = (size_t)bins;
    size_t values;
    int threads = 1;
    int status = FAILURE;

#ifdef OPENMPCODE
    threads = omp_get_max_threads();
#endif
    if (threads < 1) threads = 1;

    if (SHEAR_SPHERE_BINARY_TREE_FRONTIER(
            cmd, first_tree, KD_SHEAR_PAIR_FRONTIER,
            &frontier1, &count1) == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        frontier2 = frontier1;
        count2 = count1;
    } else if (SHEAR_SPHERE_BINARY_TREE_FRONTIER(
                   cmd, second_tree, KD_SHEAR_PAIR_FRONTIER,
                   &frontier2, &count2) == FAILURE) {
        goto cleanup;
    }
    if (count1 <= 0 || count2 <= 0 || stride == 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree task count overflow");
        goto cleanup;
    }
    for (INTEGER first = 0; first < count1; first++) {
        const INTEGER pair_count = auto_correlation
            ? count2 - first - 1 : count2;
        const INTEGER chunks = pair_count > 0
            ? 1 + (pair_count - 1)/KD_SHEAR_PAIR_TASK_CHUNK : 0;
        task_count += chunks + (auto_correlation ? 1 : 0);
    }
    if (task_count <= 0 || stride > SIZE_MAX/5
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
                 SHEAR_ENGINE_NAME ": dual-tree result allocation failed");
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
                    count2 - 1, second + KD_SHEAR_PAIR_TASK_CHUNK - 1);
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
            task_work[task] = kd_shear_pair_task_work(
                cmd, first_tree, second_tree, frontier1, frontier2,
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
        for (task = 0; task < task_count; task++) {
            task_local_index[task] = local_task_count++;
        }
#endif
    }
    if ((size_t)MAX((INTEGER)1, local_task_count)
            > SIZE_MAX/(5*stride)
        || (values = (size_t)MAX((INTEGER)1, local_task_count)*5*stride)
            > SIZE_MAX/sizeof(*task_storage)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result size overflow");
        goto cleanup;
    }
    task_storage = calloc(values, sizeof(*task_storage));
    local_tasks = malloc(
        (size_t)MAX((INTEGER)1, local_task_count)*sizeof(*local_tasks));
    body_pairs = calloc(
        (size_t)MAX((INTEGER)1, local_task_count), sizeof(*body_pairs));
    cell_pairs = calloc(
        (size_t)MAX((INTEGER)1, local_task_count), sizeof(*cell_pairs));
    if (getenv("CBALLS_SHEAR_PROFILE") != NULL)
        profiles = calloc((size_t)threads, sizeof(*profiles));
    if (task_storage == NULL || local_tasks == NULL
        || body_pairs == NULL || cell_pairs == NULL
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
            real *base;

            base = task_storage + (size_t)local_task*5*stride;
            kd_shear_pair_histogram hist = {
                base, base + stride, base + 2*stride, base + 3*stride,
                base + 4*stride, 0, 0, profile
            };
            const INTEGER first_index = task_first[task];
            const INTEGER second_index = task_second[task];
            const INTEGER last_index = task_last[task];
            kd_shear_ref first = kd_shear_node_ref(
                first_tree, frontier1[first_index]);

            if (auto_correlation) {
                if (first_index == second_index)
                    kd_shear_process_auto(cmd, gd, &hist, first);
                else for (INTEGER other = second_index;
                          other <= last_index; other++)
                        kd_shear_process_pair(
                            cmd, gd, &hist, first,
                            kd_shear_node_ref(
                                first_tree, frontier2[other]), TRUE);
            } else {
                for (INTEGER other = second_index;
                     other <= last_index; other++)
                    kd_shear_process_pair(
                        cmd, gd, &hist, first,
                        kd_shear_node_ref(
                            second_tree, frontier2[other]), FALSE);
            }
            body_pairs[local_task] = hist.body_pairs;
            cell_pairs[local_task] = hist.cell_pairs;
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
        gd->nbbcalc += body_pairs[local_task];
        gd->nbccalc += cell_pairs[local_task];
    }
    status = SUCCESS;

cleanup:
    free(profiles);
    free(cell_pairs);
    free(body_pairs);
    free(local_tasks);
    free(task_storage);
    free(rank_work);
    free(task_work);
    free(task_owner);
    free(task_local_index);
    free(task_last);
    free(task_second);
    free(task_first);
    if (!auto_correlation) free(frontier2);
    free(frontier1);
    return status;
}

static void kd_shear_release_trees(fcfc_balltreeptr pivot,
                                   fcfc_balltreeptr first,
                                   fcfc_balltreeptr second)
{
    if (second != NULL && second != first && second != pivot)
        SHEAR_SPHERE_BINARY_TREE_FREE(second);
    if (first != NULL && first != pivot)
        SHEAR_SPHERE_BINARY_TREE_FREE(first);
    SHEAR_SPHERE_BINARY_TREE_FREE(pivot);
}

static bool kd_shear_cell_geometry(shear_pivot_workspace *work,
                                   fcfc_ballnode *node,
                                   real *distance, shear_complex *phase,
                                   int *bin)
{
    compute_vector center_unit;
    compute_vector delta;
    real distance2;
    const real radius = (real)node->radius;
    real lower;
    real upper;
    real angular_tolerance;
    real cell_angle;
    real separation;

    if (!shear_load_unit3(node->center, center_unit)) return FALSE;
    DOTPSUBV(distance2, delta, work->pivot_unit, center_unit);
    if (!(distance2 >= 0.0) || !isfinite(distance2)) return FALSE;
    *distance = rsqrt(distance2);
    *bin = shear_radial_bin_profiled(work, *distance);
    if (*distance > 0.0
        && !shear_spherical_phase(work, node->center, phase))
        return FALSE;
    if (*distance + radius <= work->cmd->rminHist
        || *distance - radius >= work->cmd->rangeN)
        return FALSE;
    if (*bin < 0 || !work->allow_cells || !(work->cmd->theta > 0.0)
        || cballs_opt_no_two_balls(work->cmd)
        || node->weight <= 0.0 || *distance <= radius)
        return FALSE;
    lower = *distance - radius;
    upper = *distance + radius;
    if (!shear_interval_within_radial_bin(
            work->cmd, work->gd, lower, upper, *bin))
        return FALSE;
    angular_tolerance = MIN(0.5*PI,
        work->cmd->theta*PI/(2.0*(real)work->ring_max + 1.0));
    if (radius/(*distance) > rsin(MAX(0.0, angular_tolerance)))
        return FALSE;
    cell_angle = 2.0*rasin(MIN(1.0, 0.5*radius));
    separation = 2.0*rasin(MIN(1.0, 0.5*(*distance)));
    return isfinite(node->shear_transport_error)
        && node->shear_transport_error + 2.0*cell_angle*separation
           <= angular_tolerance;
}

static void kd_shear_accumulate_node_2pcf(shear_pivot_workspace *work,
                                           fcfc_ballnode *node,
                                           int bin, shear_complex phase)
{
    shear_complex rotation;
    shear_complex gamma = shear_make(node->shear_gamma_re,
                                     node->shear_gamma_im);

    if (!shear_transport_rotation_profiled(work, node->center, &rotation))
        return;
    shear_accumulate_2pcf_sample(work, bin, phase,
                                 shear_mul(gamma, rotation), node->weight);
}

static void kd_shear_accumulate_node_3pcf(shear_pivot_workspace *work,
                                           fcfc_ballnode *node,
                                           int bin, shear_complex phase)
{
    shear_complex rotation;
    shear_complex weighted_gamma = shear_make(node->shear_gamma_re,
                                              node->shear_gamma_im);
    shear_complex gamma2 = shear_make(node->shear_gamma2_re,
                                      node->shear_gamma2_im);

    if (!shear_transport_rotation_profiled(work, node->center, &rotation))
        return;
    weighted_gamma = shear_mul(weighted_gamma, rotation);
    gamma2 = shear_mul(gamma2, shear_mul(rotation, rotation));
    shear_accumulate_3pcf_sample(
        work, bin, phase, weighted_gamma, node->weight, gamma2,
        node->shear_gamma_abs2, node->shear_weight2);
}

static void kd_shear_walk(shear_pivot_workspace *work,
                          fcfc_balltreeptr tree, INTEGER node_index,
                          bool run_2pcf, bool run_3pcf)
{
    fcfc_ballnode *node = &tree->nodes[node_index];
    real distance = 0.0;
    shear_complex phase = shear_make(0.0, 0.0);
    int bin = -1;
    bool accept = kd_shear_cell_geometry(work, node, &distance, &phase, &bin);

    if (distance + (real)node->radius <= work->cmd->rminHist
        || distance - (real)node->radius >= work->cmd->rangeN)
        return;
    if (accept) {
        if (run_2pcf) kd_shear_accumulate_node_2pcf(work, node, bin, phase);
        if (run_3pcf) kd_shear_accumulate_node_3pcf(work, node, bin, phase);
        return;
    }
    if (node->left >= 0 && node->right >= 0) {
        kd_shear_walk(work, tree, node->left, run_2pcf, run_3pcf);
        kd_shear_walk(work, tree, node->right, run_2pcf, run_3pcf);
        return;
    }
    for (INTEGER point = node->first; point <= node->last; point++) {
        bodyptr body = tree->bptr[point];
        if (body == work->pivot) continue;
        if (run_2pcf) shear_accumulate_body_2pcf(work, (nodeptr)body);
        if (run_3pcf) shear_accumulate_body_3pcf(work, (nodeptr)body);
    }
}

#ifdef SHEAR_SPHERE_BINARY_FRONTIER_SCHEDULER

/* Split expensive pivot cells first, then retain only conservative neighbor
 * frontier roots.  Exact body pivots and in-root traversal order are kept. */
#define KD_SHEAR_PIVOT_TASK_TARGET ((INTEGER)256)
#define KD_SHEAR_NEIGHBOR_FRONTIER_TARGET ((INTEGER)256)
#define KD_SHEAR_TASK_HISTOGRAM_MEMORY ((size_t)256 << 20)

typedef struct {
    INTEGER pivot_node;
    INTEGER first;
    INTEGER last;
    int owner;
    INTEGER local_index;
    size_t first_offset;
    INTEGER first_count;
    size_t second_offset;
    INTEGER second_count;
    long double estimated_work;
} kd_shear_pivot_task;

typedef struct {
    kd_shear_pivot_task *tasks;
    INTEGER task_count;
    INTEGER *first_nodes;
    INTEGER *second_nodes;
    size_t first_node_count;
    size_t second_node_count;
} kd_shear_frontier_schedule;

static INTEGER kd_shear_node_point_count(const fcfc_ballnode *node)
{
    return node->last - node->first + 1;
}

static bool kd_shear_frontier_nodes_overlap(
        const struct cmdline_data *cmd, const fcfc_ballnode *pivot,
        const fcfc_ballnode *neighbor)
{
    compute_vector delta;
    real distance2;
    real distance;
    const real radius = (real)pivot->radius + (real)neighbor->radius;

    DOTPSUBV(distance2, delta, pivot->center, neighbor->center);
    if (!(distance2 >= 0.0) || !isfinite(distance2)) return TRUE;
    distance = rsqrt(distance2);
    return distance + radius > cmd->rminHist
        && distance - radius < cmd->rangeN;
}

static long double kd_shear_estimate_pivot_task(
        const struct cmdline_data *cmd, fcfc_balltreeptr pivot_tree,
        INTEGER pivot_node, fcfc_balltreeptr first_tree,
        const INTEGER *first_frontier, INTEGER first_count,
        fcfc_balltreeptr second_tree, const INTEGER *second_frontier,
        INTEGER second_count)
{
    const fcfc_ballnode *pivot = &pivot_tree->nodes[pivot_node];
    long double neighbor_points = 0.0L;

    for (INTEGER i = 0; i < first_count; i++) {
        const fcfc_ballnode *neighbor = &first_tree->nodes[first_frontier[i]];

        if (kd_shear_frontier_nodes_overlap(cmd, pivot, neighbor))
            neighbor_points += (long double)kd_shear_node_point_count(neighbor);
    }
    if (second_tree != NULL) {
        for (INTEGER i = 0; i < second_count; i++) {
            const fcfc_ballnode *neighbor =
                &second_tree->nodes[second_frontier[i]];

            if (kd_shear_frontier_nodes_overlap(cmd, pivot, neighbor))
                neighbor_points +=
                    (long double)kd_shear_node_point_count(neighbor);
        }
    }
    return (long double)kd_shear_node_point_count(pivot) * neighbor_points;
}

static kd_shear_pivot_task kd_shear_make_pivot_task(
        const struct cmdline_data *cmd, fcfc_balltreeptr pivot_tree,
        INTEGER pivot_node, fcfc_balltreeptr first_tree,
        const INTEGER *first_frontier, INTEGER first_count,
        fcfc_balltreeptr second_tree, const INTEGER *second_frontier,
        INTEGER second_count)
{
    const fcfc_ballnode *node = &pivot_tree->nodes[pivot_node];
    kd_shear_pivot_task task;

    memset(&task, 0, sizeof(task));
    task.pivot_node = pivot_node;
    task.first = node->first;
    task.last = node->last;
    task.owner = -1;
    task.estimated_work = kd_shear_estimate_pivot_task(
        cmd, pivot_tree, pivot_node, first_tree, first_frontier, first_count,
        second_tree, second_frontier, second_count);
    return task;
}

static void kd_shear_sort_pivot_tasks(kd_shear_pivot_task *tasks,
                                      INTEGER count)
{
    for (INTEGER i = 1; i < count; i++) {
        kd_shear_pivot_task value = tasks[i];
        INTEGER j = i;

        while (j > 0 && tasks[j - 1].first > value.first) {
            tasks[j] = tasks[j - 1];
            j--;
        }
        tasks[j] = value;
    }
}

static int kd_shear_assign_pivot_tasks(
        struct cmdline_data *cmd, kd_shear_pivot_task *tasks,
        INTEGER count, int rank_count)
{
    long double *rank_work;

    rank_count = MAX(1, rank_count);
    if ((size_t)rank_count > SIZE_MAX / sizeof(*rank_work)
        || (rank_work = calloc((size_t)rank_count, sizeof(*rank_work)))
               == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": rank-work allocation failed");
        return FAILURE;
    }
    for (INTEGER task = 0; task < count; task++) tasks[task].owner = -1;
    for (INTEGER assigned = 0; assigned < count; assigned++) {
        INTEGER largest = -1;
        int owner = 0;

        for (INTEGER task = 0; task < count; task++) {
            if (tasks[task].owner < 0
                && (largest < 0
                    || tasks[task].estimated_work
                       > tasks[largest].estimated_work))
                largest = task;
        }
        for (int rank = 1; rank < rank_count; rank++)
            if (rank_work[rank] < rank_work[owner]) owner = rank;
        tasks[largest].owner = owner;
        rank_work[owner] += tasks[largest].estimated_work;
    }
    for (INTEGER task = 0; task < count; task++) {
        INTEGER local_index = 0;

        for (INTEGER previous = 0; previous < task; previous++)
            if (tasks[previous].owner == tasks[task].owner) local_index++;
        tasks[task].local_index = local_index;
    }
    free(rank_work);
    return SUCCESS;
}

static size_t kd_shear_count_task_neighbors(
        const struct cmdline_data *cmd, const fcfc_ballnode *pivot,
        fcfc_balltreeptr neighbor_tree, const INTEGER *frontier,
        INTEGER frontier_count)
{
    size_t count = 0;

    for (INTEGER i = 0; i < frontier_count; i++)
        if (kd_shear_frontier_nodes_overlap(
                cmd, pivot, &neighbor_tree->nodes[frontier[i]]))
            count++;
    return count;
}

static void kd_shear_fill_task_neighbors(
        const struct cmdline_data *cmd, const fcfc_ballnode *pivot,
        fcfc_balltreeptr neighbor_tree, const INTEGER *frontier,
        INTEGER frontier_count, INTEGER *result)
{
    size_t count = 0;

    for (INTEGER i = 0; i < frontier_count; i++)
        if (kd_shear_frontier_nodes_overlap(
                cmd, pivot, &neighbor_tree->nodes[frontier[i]]))
            result[count++] = frontier[i];
}

static void kd_shear_release_frontier_schedule(
        kd_shear_frontier_schedule *schedule)
{
    if (schedule == NULL) return;
    free(schedule->second_nodes);
    free(schedule->first_nodes);
    free(schedule->tasks);
    memset(schedule, 0, sizeof(*schedule));
}

static int kd_shear_build_frontier_schedule(
        struct cmdline_data *cmd, fcfc_balltreeptr pivot_tree,
        fcfc_balltreeptr first_tree, fcfc_balltreeptr second_tree,
        INTEGER target, int rank_count,
        kd_shear_frontier_schedule *schedule)
{
    INTEGER *first_frontier = NULL;
    INTEGER *second_frontier = NULL;
    INTEGER first_frontier_count = 0;
    INTEGER second_frontier_count = 0;
    INTEGER count = 1;
    size_t first_total = 0;
    size_t second_total = 0;
    int status = FAILURE;

    memset(schedule, 0, sizeof(*schedule));
    if (target < 1 || (size_t)target > SIZE_MAX / sizeof(*schedule->tasks)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid pivot-frontier target");
        return FAILURE;
    }
    if (SHEAR_SPHERE_BINARY_TREE_FRONTIER(
            cmd, first_tree, KD_SHEAR_NEIGHBOR_FRONTIER_TARGET,
            &first_frontier, &first_frontier_count) == FAILURE)
        goto cleanup;
    if (second_tree != NULL
        && SHEAR_SPHERE_BINARY_TREE_FRONTIER(
               cmd, second_tree, KD_SHEAR_NEIGHBOR_FRONTIER_TARGET,
               &second_frontier, &second_frontier_count) == FAILURE)
        goto cleanup;

    schedule->tasks = calloc((size_t)target, sizeof(*schedule->tasks));
    if (schedule->tasks == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": pivot-frontier allocation failed");
        goto cleanup;
    }
    schedule->tasks[0] = kd_shear_make_pivot_task(
        cmd, pivot_tree, 0, first_tree, first_frontier,
        first_frontier_count, second_tree, second_frontier,
        second_frontier_count);
    while (count < target) {
        INTEGER best = -1;
        long double largest_work = -1.0L;

        for (INTEGER i = 0; i < count; i++) {
            const fcfc_ballnode *node =
                &pivot_tree->nodes[schedule->tasks[i].pivot_node];

            if (node->left >= 0 && node->right >= 0
                && schedule->tasks[i].estimated_work > largest_work) {
                best = i;
                largest_work = schedule->tasks[i].estimated_work;
            }
        }
        if (best < 0) break;
        {
            const fcfc_ballnode *node =
                &pivot_tree->nodes[schedule->tasks[best].pivot_node];

            schedule->tasks[best] = kd_shear_make_pivot_task(
                cmd, pivot_tree, node->left, first_tree, first_frontier,
                first_frontier_count, second_tree, second_frontier,
                second_frontier_count);
            schedule->tasks[count++] = kd_shear_make_pivot_task(
                cmd, pivot_tree, node->right, first_tree, first_frontier,
                first_frontier_count, second_tree, second_frontier,
                second_frontier_count);
        }
    }
    kd_shear_sort_pivot_tasks(schedule->tasks, count);
    if (kd_shear_assign_pivot_tasks(
            cmd, schedule->tasks, count, rank_count) == FAILURE)
        goto cleanup;

    for (INTEGER task = 0; task < count; task++) {
        const fcfc_ballnode *pivot =
            &pivot_tree->nodes[schedule->tasks[task].pivot_node];
        const size_t first_count = kd_shear_count_task_neighbors(
            cmd, pivot, first_tree, first_frontier, first_frontier_count);
        size_t second_count = 0;

        if (second_tree != NULL)
            second_count = kd_shear_count_task_neighbors(
                cmd, pivot, second_tree, second_frontier,
                second_frontier_count);
        if (first_count > SIZE_MAX - first_total
            || second_count > SIZE_MAX - second_total
            || first_count > (size_t)
#ifdef LONGINT
                LONG_MAX
#else
                INT_MAX
#endif
            || second_count > (size_t)
#ifdef LONGINT
                LONG_MAX
#else
                INT_MAX
#endif
            ) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": frontier node count overflow");
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
         && first_total > SIZE_MAX / sizeof(*schedule->first_nodes))
        || (second_total > 0
            && second_total > SIZE_MAX / sizeof(*schedule->second_nodes))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": frontier allocation size overflow");
        goto cleanup;
    }
    if (first_total > 0)
        schedule->first_nodes = malloc(
            first_total * sizeof(*schedule->first_nodes));
    if (second_total > 0)
        schedule->second_nodes = malloc(
            second_total * sizeof(*schedule->second_nodes));
    if ((first_total > 0 && schedule->first_nodes == NULL)
        || (second_total > 0 && schedule->second_nodes == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": neighbor-frontier allocation failed");
        goto cleanup;
    }
    for (INTEGER task = 0; task < count; task++) {
        const fcfc_ballnode *pivot =
            &pivot_tree->nodes[schedule->tasks[task].pivot_node];

        if (schedule->tasks[task].first_count > 0)
            kd_shear_fill_task_neighbors(
                cmd, pivot, first_tree, first_frontier,
                first_frontier_count,
                schedule->first_nodes + schedule->tasks[task].first_offset);
        if (second_tree != NULL && schedule->tasks[task].second_count > 0)
            kd_shear_fill_task_neighbors(
                cmd, pivot, second_tree, second_frontier,
                second_frontier_count,
                schedule->second_nodes + schedule->tasks[task].second_offset);
    }
    schedule->task_count = count;
    schedule->first_node_count = first_total;
    schedule->second_node_count = second_total;
    status = SUCCESS;

cleanup:
    free(second_frontier);
    free(first_frontier);
    if (status == FAILURE) kd_shear_release_frontier_schedule(schedule);
    return status;
}

static INTEGER kd_shear_pivot_task_target(int threads, int ranks,
                                          size_t accumulator_stride)
{
    size_t workers = (size_t)MAX(1, threads);
    size_t target;
    size_t bytes_per_task;
    size_t memory_target;

    ranks = MAX(1, ranks);
    if (workers > SIZE_MAX / (size_t)ranks)
        workers = SIZE_MAX;
    else
        workers *= (size_t)ranks;
    target = workers > (size_t)KD_SHEAR_PIVOT_TASK_TARGET / 8
        ? (size_t)KD_SHEAR_PIVOT_TASK_TARGET : 8 * workers;
    target = MAX((size_t)64, target);
    target = MIN(target, (size_t)KD_SHEAR_PIVOT_TASK_TARGET);
    if (accumulator_stride == 0
        || accumulator_stride > SIZE_MAX / sizeof(real))
        return 1;
    bytes_per_task = accumulator_stride * sizeof(real);
    memory_target = KD_SHEAR_TASK_HISTOGRAM_MEMORY / bytes_per_task;
    target = MIN(target, MAX((size_t)1, memory_target));
    return (INTEGER)MAX((size_t)1, target);
}

#undef KD_SHEAR_PIVOT_TASK_TARGET
#undef KD_SHEAR_NEIGHBOR_FRONTIER_TARGET
#undef KD_SHEAR_TASK_HISTOGRAM_MEMORY
#endif /* SHEAR_SPHERE_BINARY_FRONTIER_SCHEDULER */

#undef KD_SHEAR_PAIR_FRONTIER
#undef KD_SHEAR_SPLIT_FACTOR
#endif
