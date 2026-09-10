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
                             kd_shear_ref first, kd_shear_ref second,
                             bool permit_nodes)
{
    compute_vector first_unit;
    compute_vector first_east;
    compute_vector first_north;
    compute_vector second_unit;
    compute_vector second_east;
    compute_vector second_north;
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

    if (!shear_unit3(kd_shear_position(first), first_unit)
        || !shear_spherical_basis(first_unit, first_east, first_north)
        || !shear_unit3(kd_shear_position(second), second_unit)
        || !shear_spherical_basis(second_unit, second_east, second_north))
        return -1;
    (void)first_east;
    (void)first_north;
    (void)second_east;
    (void)second_north;
    DOTPSUBV(distance2, delta, first_unit, second_unit);
    if (!(distance2 > 0.0) || !isfinite(distance2)) return -2;
    distance = rsqrt(distance2);
    size = kd_shear_radius(first) + kd_shear_radius(second);
    if (distance + size <= cmd->rminHist
        || distance - size >= cmd->rangeN)
        return -2;
    bin = shear_radial_bin(cmd, gd, distance);
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
    } else if (!(lower > cmd->rminHist && upper < cmd->rangeN)
               || shear_radial_bin(cmd, gd, lower) != bin
               || shear_radial_bin(cmd, gd, upper) != bin) {
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

static bool kd_shear_pair_phase(const real *pivot_unit,
                                const real *pivot_east,
                                const real *pivot_north,
                                const real *neighbor_unit,
                                shear_complex *phase)
{
    compute_vector tangent;
    const real projection = shear_dot3(pivot_unit, neighbor_unit);
    real norm;

    tangent[0] = neighbor_unit[0] - projection*pivot_unit[0];
    tangent[1] = neighbor_unit[1] - projection*pivot_unit[1];
    tangent[2] = neighbor_unit[2] - projection*pivot_unit[2];
    norm = rsqrt(shear_dot3(tangent, tangent));
    if (!(norm > 64.0*DBL_EPSILON) || !isfinite(norm)) return FALSE;
    phase->re = shear_dot3(tangent, pivot_east)/norm;
    phase->im = shear_dot3(tangent, pivot_north)/norm;
    return isfinite(phase->re) && isfinite(phase->im);
}

static void kd_shear_accumulate_oriented(
        kd_shear_pair_histogram *hist, kd_shear_ref pivot,
        kd_shear_ref neighbor, int bin)
{
    compute_vector pivot_unit;
    compute_vector pivot_east;
    compute_vector pivot_north;
    compute_vector neighbor_unit;
    shear_complex phase;
    shear_complex rotation;
    shear_complex pivot_gamma;
    shear_complex neighbor_gamma;
    shear_complex z2;
    shear_complex z4;
    shear_complex value;

    if (!shear_unit3(kd_shear_position(pivot), pivot_unit)
        || !shear_spherical_basis(pivot_unit, pivot_east, pivot_north)
        || !shear_unit3(kd_shear_position(neighbor), neighbor_unit)
        || !kd_shear_pair_phase(pivot_unit, pivot_east, pivot_north,
                                neighbor_unit, &phase)
        || !shear_transport_rotation_to_pivot(
               pivot_unit, pivot_east, pivot_north,
               kd_shear_position(neighbor), &rotation))
        return;
    pivot_gamma = kd_shear_weighted_gamma(pivot);
    neighbor_gamma = shear_mul(kd_shear_weighted_gamma(neighbor), rotation);
    value = shear_mul(pivot_gamma, shear_conj(neighbor_gamma));
    hist->xi_plus_re[bin] += value.re;
    hist->xi_plus_im[bin] += value.im;
    z2 = shear_mul(phase, phase);
    z4 = shear_mul(z2, z2);
    value = shear_mul(shear_mul(pivot_gamma, neighbor_gamma), shear_conj(z4));
    hist->xi_minus_re[bin] += value.re;
    hist->xi_minus_im[bin] += value.im;
    hist->weight[bin] += kd_shear_weight(pivot)*kd_shear_weight(neighbor);
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
    int bin = kd_shear_pair_bin(cmd, gd, first, second, TRUE);

    if (bin == -2) return;
    if (bin >= 0) {
        kd_shear_accumulate_oriented(hist, first, second, bin);
        if (bidirectional)
            kd_shear_accumulate_oriented(hist, second, first, bin);
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

static int kd_shear_dual_tree_2pcf(
        struct cmdline_data *cmd, struct global_data *gd,
        fcfc_balltreeptr first_tree, fcfc_balltreeptr second_tree,
        bool auto_correlation, int bins)
{
    INTEGER *frontier1 = NULL;
    INTEGER *frontier2 = NULL;
    INTEGER count1 = 0;
    INTEGER count2 = 0;
    real *task_storage = NULL;
    INTEGER *body_pairs = NULL;
    INTEGER *cell_pairs = NULL;
    const size_t stride = (size_t)bins;
    size_t values;
    int status = FAILURE;

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
    if (count1 <= 0 || stride == 0 || stride > SIZE_MAX/5
        || (size_t)count1 > SIZE_MAX/(5*stride)
        || (values = (size_t)count1*5*stride)
             > SIZE_MAX/sizeof(*task_storage)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result size overflow");
        goto cleanup;
    }
    task_storage = calloc(values, sizeof(*task_storage));
    body_pairs = calloc((size_t)count1, sizeof(*body_pairs));
    cell_pairs = calloc((size_t)count1, sizeof(*cell_pairs));
    if (task_storage == NULL || body_pairs == NULL || cell_pairs == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result allocation failed");
        goto cleanup;
    }
#pragma omp parallel for schedule(dynamic,1)
    for (INTEGER task = 0; task < count1; task++) {
        real *base = task_storage + (size_t)task*5*stride;
        kd_shear_pair_histogram hist = {
            base, base + stride, base + 2*stride, base + 3*stride,
            base + 4*stride, 0, 0
        };
        kd_shear_ref first = kd_shear_node_ref(first_tree, frontier1[task]);

        if (auto_correlation) {
            kd_shear_process_auto(cmd, gd, &hist, first);
            for (INTEGER other = task + 1; other < count2; other++)
                kd_shear_process_pair(
                    cmd, gd, &hist, first,
                    kd_shear_node_ref(first_tree, frontier2[other]), TRUE);
        } else {
            for (INTEGER other = 0; other < count2; other++)
                kd_shear_process_pair(
                    cmd, gd, &hist, first,
                    kd_shear_node_ref(second_tree, frontier2[other]), FALSE);
        }
        body_pairs[task] = hist.body_pairs;
        cell_pairs[task] = hist.cell_pairs;
    }
    for (INTEGER task = 0; task < count1; task++) {
        const real *base = task_storage + (size_t)task*5*stride;
        for (int bin = 0; bin < bins; bin++) {
            gd->histShearXiPlusRe[bin] += base[bin];
            gd->histShearXiPlusIm[bin] += base[stride + (size_t)bin];
            gd->histShearXiMinusRe[bin] += base[2*stride + (size_t)bin];
            gd->histShearXiMinusIm[bin] += base[3*stride + (size_t)bin];
            gd->histShearXiWeight[bin] += base[4*stride + (size_t)bin];
        }
        gd->nbbcalc += body_pairs[task];
        gd->nbccalc += cell_pairs[task];
    }
    status = SUCCESS;

cleanup:
    free(cell_pairs);
    free(body_pairs);
    free(task_storage);
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

    if (!shear_unit3(node->center, center_unit)) return FALSE;
    DOTPSUBV(distance2, delta, work->pivot_unit, center_unit);
    if (!(distance2 >= 0.0) || !isfinite(distance2)) return FALSE;
    *distance = rsqrt(distance2);
    *bin = shear_radial_bin(work->cmd, work->gd, *distance);
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
    if (!(lower > work->cmd->rminHist && upper < work->cmd->rangeN)
        || shear_radial_bin(work->cmd, work->gd, lower) != *bin
        || shear_radial_bin(work->cmd, work->gd, upper) != *bin)
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

    if (!shear_transport_rotation_to_pivot(
            work->pivot_unit, work->pivot_east, work->pivot_north,
            node->center, &rotation))
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

    if (!shear_transport_rotation_to_pivot(
            work->pivot_unit, work->pivot_east, work->pivot_north,
            node->center, &rotation))
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

#undef KD_SHEAR_PAIR_FRONTIER
#undef KD_SHEAR_SPLIT_FACTOR
#endif
