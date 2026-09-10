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
        nodeptr second, bool permit_cells)
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
    bin = shear_radial_bin(context->cmd, context->gd, distance);
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

static void shear_sphere_accumulate_oriented(
        const shear_sphere_pair_context *context,
        shear_sphere_pair_histogram *hist, nodeptr pivot, nodeptr neighbor,
        int bin)
{
    compute_vector pivot_unit;
    compute_vector pivot_east;
    compute_vector pivot_north;
    compute_vector neighbor_unit;
    compute_vector neighbor_east;
    compute_vector neighbor_north;
    shear_complex phase;
    shear_complex rotation;
    shear_complex pivot_gamma;
    shear_complex neighbor_gamma;
    shear_complex z2;
    shear_complex z4;
    shear_complex value;
    real pivot_radius;
    real neighbor_radius;
    real pivot_angle;
    real neighbor_angle;
    real pair_weight;

    if (bin < 0 || bin >= context->bins
        || !shear_sphere_pair_node_frame(
               context, pivot, pivot_unit, pivot_east, pivot_north,
               &pivot_radius, &pivot_angle)
        || !shear_sphere_pair_node_frame(
               context, neighbor, neighbor_unit, neighbor_east,
               neighbor_north, &neighbor_radius, &neighbor_angle)
        || !shear_sphere_pair_phase(
               pivot_unit, pivot_east, pivot_north, neighbor_unit, &phase)
        || !shear_transport_rotation_to_pivot(
               pivot_unit, pivot_east, pivot_north, neighbor_unit,
               &rotation))
        return;
    (void)neighbor_east;
    (void)neighbor_north;
    (void)pivot_radius;
    (void)neighbor_radius;
    (void)pivot_angle;
    (void)neighbor_angle;
    pivot_gamma = shear_sphere_node_weighted_gamma(pivot);
    neighbor_gamma = shear_mul(
        shear_sphere_node_weighted_gamma(neighbor), rotation);
    pair_weight = shear_sphere_node_weight(pivot)
                * shear_sphere_node_weight(neighbor);
    value = shear_mul(pivot_gamma, shear_conj(neighbor_gamma));
    hist->xi_plus_re[bin] += value.re;
    hist->xi_plus_im[bin] += value.im;
    z2 = shear_mul(phase, phase);
    z4 = shear_mul(z2, z2);
    value = shear_mul(shear_mul(pivot_gamma, neighbor_gamma),
                      shear_conj(z4));
    hist->xi_minus_re[bin] += value.re;
    hist->xi_minus_im[bin] += value.im;
    hist->weight[bin] += pair_weight;
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
    bin = shear_sphere_pair_bin(context, first, second, TRUE);
    if (bin == -2)
        return;
    if (bin >= 0) {
        shear_sphere_accumulate_oriented(context, hist, first, second, bin);
        if (bidirectional)
            shear_sphere_accumulate_oriented(
                context, hist, second, first, bin);
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
    frontier[0] = root;
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
        frontier[selected] = children[0];
        for (int child = 1; child < child_count; child++)
            frontier[count++] = children[child];
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

static int shear_sphere_dual_tree_2pcf(
        struct cmdline_data *cmd, struct global_data *gd,
        nodeptr root1, nodeptr root2, bool auto_correlation, int bins)
{
    shear_sphere_pair_context context;
    nodeptr *frontier1 = NULL;
    nodeptr *frontier2 = NULL;
    INTEGER count1 = 0;
    INTEGER count2 = 0;
    real *task_storage = NULL;
    INTEGER *task_body_pairs = NULL;
    INTEGER *task_cell_pairs = NULL;
    size_t stride;
    size_t values;
    int status = FAILURE;

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
    stride = (size_t)bins;
    if (count1 <= 0 || stride == 0 || stride > SIZE_MAX/5
        || (size_t)count1 > SIZE_MAX/(5*stride)
        || (values = (size_t)count1*(5*stride))
             > SIZE_MAX/sizeof(*task_storage)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result size overflow");
        goto cleanup;
    }
    task_storage = calloc(values, sizeof(*task_storage));
    task_body_pairs = calloc((size_t)count1, sizeof(*task_body_pairs));
    task_cell_pairs = calloc((size_t)count1, sizeof(*task_cell_pairs));
    if (task_storage == NULL || task_body_pairs == NULL
        || task_cell_pairs == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": dual-tree result allocation failed");
        goto cleanup;
    }

#ifdef OPENMPCODE
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (INTEGER task = 0; task < count1; task++) {
        real *base = task_storage + (size_t)task*5*stride;
        shear_sphere_pair_histogram hist;

        hist.xi_plus_re = base;
        hist.xi_plus_im = base + stride;
        hist.xi_minus_re = base + 2*stride;
        hist.xi_minus_im = base + 3*stride;
        hist.weight = base + 4*stride;
        hist.body_pairs = 0;
        hist.cell_pairs = 0;
        if (auto_correlation) {
            shear_sphere_process_auto(&context, &hist, frontier1[task]);
            for (INTEGER other = task + 1; other < count2; other++)
                shear_sphere_process_pair(
                    &context, &hist, frontier1[task], frontier2[other], TRUE);
        } else {
            for (INTEGER other = 0; other < count2; other++)
                shear_sphere_process_pair(
                    &context, &hist, frontier1[task], frontier2[other], FALSE);
        }
        task_body_pairs[task] = hist.body_pairs;
        task_cell_pairs[task] = hist.cell_pairs;
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
        gd->nbbcalc += task_body_pairs[task];
        gd->nbccalc += task_cell_pairs[task];
    }
    verb_print(cmd->verbose,
               SHEAR_ENGINE_NAME ": spherical dual tree used %"
               INTEGER_FMT " tasks, %" INTEGER_FMT " body pairs and %"
               INTEGER_FMT " accepted node pairs\n",
               count1, gd->nbbcalc, gd->nbccalc);
    status = SUCCESS;

cleanup:
    free(task_cell_pairs);
    free(task_body_pairs);
    free(task_storage);
    if (!auto_correlation)
        free(frontier2);
    free(frontier1);
    return status;
}

#undef SHEAR_SPHERE_PAIR_MAX_FRONTIER
#undef SHEAR_SPHERE_PAIR_FRONTIER
#undef SHEAR_SPHERE_DUAL_NODE_SPLIT_FACTOR
#endif
