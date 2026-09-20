#ifndef CTREEBALLS_SHEAR_PIVOT_REUSE_H
#define CTREEBALLS_SHEAR_PIVOT_REUSE_H

/* Partial rings stay in their acceptance frame. Children inherit only sparse
 * unresolved lists; completed branches rotate each ancestor directly once. */
#define SHEAR_REUSE_MEMORY_LIMIT ((size_t)64 << 20)

#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
typedef kd_shear_ref shear_reuse_ref;
static bool shear_reuse_cell(shear_reuse_ref q) { return q.kind == KD_SHEAR_NODE; }
static real shear_reuse_error(shear_reuse_ref q) { return kd_shear_transport_error(q); }
#else
typedef nodeptr shear_reuse_ref;
static bool shear_reuse_cell(shear_reuse_ref q) { return Type(q) == CELL; }
static real shear_reuse_error(shear_reuse_ref q)
{ return Type(q) == CELL ? ShearTransportError(q) : 0.0; }
#endif

typedef struct shear_reuse_level {
    struct shear_reuse_level *parent, *child;
    shear_pivot_workspace work;
    shear_complex *storage;
    real *diagonal;
    unsigned char *active;
    shear_reuse_ref *unresolved[2];
    size_t count[2], capacity[2];
    real radius;
    real pivot_error;
    bool frame_valid, has_rings;
} shear_reuse_level;

typedef struct {
    shear_reuse_level *root;
    real tolerance;
    int nmax;
    size_t ring_count, weight_ring_count, bytes;
    bool failed;
    uint64_t nodes, aggregated, bodies, accepted_cells, accepted_bodies;
    uint64_t inherited, ancestor_merges, pruned, peak_unresolved, max_depth;
} shear_reuse_context;

static int shear_reuse_tolerance(struct cmdline_data *cmd, real *result)
{
    const char *text = getenv("CBALLS_SHEAR_PIVOT_TOL");
    char *end;
    double value = 0.1;

    if (text != NULL) {
        value = strtod(text, &end);
        if (end == text || *end != '\0' || !isfinite(value)
            || value < 0.0 || value > 3.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "CBALLS_SHEAR_PIVOT_TOL must be finite in [0,3] radians");
            return FAILURE;
        }
    }
    *result = (real)value;
    return SUCCESS;
}

static void shear_reuse_free(shear_reuse_context *context)
{
    shear_reuse_level *level = context->root;
    while (level != NULL) {
        shear_reuse_level *next = level->child;
        free(level->unresolved[0]);
        free(level->unresolved[1]);
        free(level->storage);
        free(level->diagonal);
        free(level->active);
        free(level);
        level = next;
    }
    context->root = NULL;
}

static int shear_reuse_allocate_level(shear_reuse_context *context,
        shear_pivot_workspace *base, shear_reuse_level **slot,
        shear_reuse_level *parent)
{
    const size_t legs = base->same_neighbor_catalog ? 1 : 2;
    size_t count;
    shear_reuse_level *level;
    shear_complex *cursor;
    size_t bytes;

    if (*slot != NULL) return SUCCESS;
    if (context->ring_count > SHEAR_REUSE_MEMORY_LIMIT/sizeof(*cursor)
        || context->weight_ring_count > SHEAR_REUSE_MEMORY_LIMIT/sizeof(*cursor)
        || (size_t)base->bins > SHEAR_REUSE_MEMORY_LIMIT/(3*sizeof(*cursor)))
        return FAILURE;
    count = legs*(context->ring_count + context->weight_ring_count)
            + (size_t)3*base->bins;
    if (shear_size_mul(count, sizeof(*cursor), &bytes) == FAILURE
        || bytes > SHEAR_REUSE_MEMORY_LIMIT)
        return FAILURE;
    const size_t total_bytes = bytes + sizeof(*level)
        + (size_t)base->bins*(sizeof(real) + 1);
    if (total_bytes > SHEAR_REUSE_MEMORY_LIMIT
        || context->bytes > SHEAR_REUSE_MEMORY_LIMIT - total_bytes)
        return FAILURE;
    level = calloc(1, sizeof(*level));
    if (level == NULL) return FAILURE;
    *slot = level;
    level->parent = parent;
    level->work = *base;
    level->storage = malloc(bytes);
    level->diagonal = malloc((size_t)base->bins*sizeof(real));
    level->active = malloc((size_t)base->bins);
    if (level->storage == NULL || level->diagonal == NULL
        || level->active == NULL) return FAILURE;
    context->bytes += total_bytes;
    cursor = level->storage;
    level->work.g_ring_first = cursor;
    cursor += context->ring_count;
    level->work.w_ring_first = cursor;
    cursor += context->weight_ring_count;
    if (base->same_neighbor_catalog) {
        level->work.g_ring_second = level->work.g_ring_first;
        level->work.w_ring_second = level->work.w_ring_first;
    } else {
        level->work.g_ring_second = cursor;
        cursor += context->ring_count;
        level->work.w_ring_second = cursor;
        cursor += context->weight_ring_count;
    }
    level->work.diag_g6 = cursor; cursor += base->bins;
    level->work.diag_g2 = cursor; cursor += base->bins;
    level->work.diag_abs2 = cursor;
    level->work.diag_w2 = level->diagonal;
    return SUCCESS;
}

static int shear_reuse_retain(shear_reuse_context *context,
                              shear_reuse_level *level, int leg, shear_reuse_ref q)
{
    if (level->count[leg] == level->capacity[leg]) {
        size_t capacity = MAX((size_t)32, level->capacity[leg]*2);
        size_t extra = (capacity - level->capacity[leg])*sizeof(shear_reuse_ref);
        shear_reuse_ref *nodes;
        if (capacity < level->capacity[leg]
            || extra > SHEAR_REUSE_MEMORY_LIMIT
            || context->bytes > SHEAR_REUSE_MEMORY_LIMIT - extra)
            return FAILURE;
        nodes = realloc(level->unresolved[leg], capacity*sizeof(*nodes));
        if (nodes == NULL) return FAILURE;
        level->unresolved[leg] = nodes;
        level->capacity[leg] = capacity;
        context->bytes += extra;
    }
    level->unresolved[leg][level->count[leg]++] = q;
    return SUCCESS;
}

/* The bearing derivative along a spherical geodesic is bounded by 1/sin(d).
 * Reject both coincident and antipodal caps. The holonomy allowance encloses
 * the two source/destination transport triangles, in spin-2 phase radians.
 * This controls phase perturbations, NOT relative error after cancellations. */
static bool shear_reuse_phase_bound(shear_reuse_context *context,
        shear_reuse_level *level, shear_reuse_ref q, real distance, real radius)
{
    const real a = 2.0*rasin(MIN(1.0, 0.5*level->radius));
    const real b = 2.0*rasin(MIN(1.0, 0.5*radius));
    const real d = 2.0*rasin(MIN(1.0, 0.5*distance));
    const real lower = d - a - b, upper = d + a + b;
    const real qerror = shear_reuse_error(q);
    const real perror = level->pivot_error;
    real sine, triangle, bearing, holonomy;

    if (!(lower > 0.0 && upper < PI)
        || !isfinite(qerror) || !isfinite(perror)
        || perror > context->tolerance/3.0)
        return FALSE;
    sine = MIN(rsin(lower), rsin(upper));
    bearing = (a + b)/sine;
    triangle = rtan(0.5*(a + b))*rtan(0.5*upper);
    if (!(triangle >= 0.0 && triangle < 1.0)) return FALSE;
    holonomy = 8.0*rasin(triangle);
    return level->work.ring_max*bearing + holonomy + qerror
        <= context->tolerance/3.0;
}

static int shear_reuse_neighbor(shear_reuse_context *context,
        shear_reuse_level *level, shear_reuse_ref q, int leg)
{
    shear_pivot_workspace *work = &level->work;
    compute_vector unit, difference;
    real radius = INFINITY, distance2, distance, size;
    int bin = -1;
    bool accept = FALSE;
    bool geometry_valid = FALSE;
    shear_complex phase;
    bool fully_valid = TRUE;
    nodeptr body;

#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
    body = shear_reuse_cell(q) ? NULL : (nodeptr)q.body;
    if (!work->aggregate_pivot && q.kind == KD_SHEAR_BODY
        && q.body == work->pivot) return SUCCESS;
    if (level->frame_valid) {
        geometry_valid = shear_load_unit3(kd_shear_position(q), unit);
        radius = kd_shear_radius(q);
    }
#else
    if (q == NULL || (cballs_opt_read_mask(work->cmd)
                      && Mask(q) == MASK_NODE_MASKED)) return SUCCESS;
    body = shear_reuse_cell(q) ? NULL : q;
    if (!work->aggregate_pivot && q == (nodeptr)work->pivot) return SUCCESS;
    fully_valid = !cballs_opt_read_mask(work->cmd) || Mask(q) == MASK_NODE_VALID;
    if (level->frame_valid) {
        if (Type(q) != CELL) {
            for (int axis = 0; axis < NDIM; axis++) unit[axis] = Pos(q)[axis];
            radius = 0.0;
            geometry_valid = TRUE;
        } else {
            geometry_valid = shear_spherical_node_geometry(
                work->cmd, q, (real)Size(q), unit, &radius);
        }
    }
#endif
    if (geometry_valid) {
        DOTPSUBV(distance2, difference, work->pivot_unit, unit);
        distance = rsqrt(MAX(0.0, distance2));
        size = level->radius + radius;
        if (distance + size <= work->cmd->rminHist
            || distance - size >= work->cmd->rangeN) {
            context->pruned++;
            return SUCCESS;
        }
        if (distance > size && fully_valid) {
            /* Body pairs have no aggregation or radial uncertainty. */
            if (!work->aggregate_pivot && !shear_reuse_cell(q)) {
                bin = shear_radial_bin_profiled(work, distance);
                accept = bin >= 0;
            } else if (size*work->ring_max <= context->tolerance/3.0
                       *distance*rsqrt(MAX(0.0, 1.0 - distance2/4.0))) {
                bin = shear_radial_bin_profiled(work, distance);
                accept = bin >= 0 && shear_interval_within_radial_bin(
                    work->cmd, work->gd, distance - size, distance + size, bin)
                    && shear_reuse_phase_bound(context, level, q, distance, radius);
            }
        }
        if (accept && shear_spherical_phase_unit(work, unit, &phase)) {
            if (shear_reuse_cell(q)) {
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
                kd_shear_accumulate_node_3pcf(work, kd_shear_node(q), bin, phase);
#else
                shear_accumulate_cell_3pcf(work, q, bin, phase);
#endif
                context->accepted_cells++;
            } else {
                shear_complex rotation, gamma;
                const real weight = Weight(body);
                if (!shear_transport_rotation_node_profiled(work, body, &rotation))
                    return SUCCESS;
                gamma = shear_scale(shear_mul(
                    shear_make(Gamma1(body), Gamma2(body)), rotation), weight);
                shear_accumulate_3pcf_sample(work, bin, phase, gamma, weight,
                    shear_mul(gamma, gamma), shear_abs2(gamma), weight*weight);
                context->accepted_bodies++;
            }
            level->active[bin] = 1;
            level->has_rings = TRUE;
            return SUCCESS;
        }
    }
    /* Split the larger side. If pivot extent dominates, postpone this node
     * for the children instead of expanding a useless large neighbor list. */
    if (shear_reuse_cell(q) && level->frame_valid
        && (!work->aggregate_pivot || radius > level->radius)) {
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
        for (INTEGER child = 0; child < kd_shear_child_count(q); child++)
            if (shear_reuse_neighbor(context, level, kd_shear_child(q, child), leg)
                == FAILURE) return FAILURE;
#else
        for (nodeptr child = More(q); child != Next(q); child = Next(child))
            if (shear_reuse_neighbor(context, level, child, leg) == FAILURE)
                return FAILURE;
#endif
        return SUCCESS;
    }
    if (!work->aggregate_pivot) return SUCCESS;
    return shear_reuse_retain(context, level, leg, q);
}

static bool shear_reuse_vector_rotation(const shear_pivot_workspace *target,
        const shear_pivot_workspace *source, shear_complex *rotation)
{
    const real dot = shear_dot3(target->pivot_unit, source->pivot_unit);
    const real denominator = 1.0 + dot;
    const real projection = shear_dot3(source->pivot_east, target->pivot_unit);
    real c, s, norm;
    if (!(denominator > 64.0*DBL_EPSILON)) return FALSE;
    c = shear_dot3(source->pivot_east, target->pivot_east)
        - shear_dot3(target->pivot_east, source->pivot_unit)*projection/denominator;
    s = shear_dot3(source->pivot_east, target->pivot_north)
        - shear_dot3(target->pivot_north, source->pivot_unit)*projection/denominator;
    norm = rsqrt(c*c + s*s);
    if (!(norm > 0.0) || !isfinite(norm)) return FALSE;
    *rotation = shear_make(c/norm, s/norm);
    return TRUE;
}

static int shear_reuse_merge_level(shear_reuse_context *context,
        shear_pivot_workspace *target, const shear_reuse_level *level)
{
    const shear_pivot_workspace *source = &level->work;
    shear_complex u, spin;
    const int max = target->ring_max;
    if (!level->has_rings) return SUCCESS;
    if (target->reuse_anchor == source->reuse_anchor) u = shear_make(1.0, 0.0);
    else if (!shear_reuse_vector_rotation(target, source, &u)) return FAILURE;
    spin = shear_mul(u, u);
    context->inherited++;
    context->ancestor_merges += target->reuse_anchor != source->reuse_anchor;
    for (int bin = 0; bin < target->bins; bin++) {
        shear_complex power = shear_make(1.0, 0.0);
        if (!level->active[bin]) continue;
        for (int order = 0; order <= max; order++) {
            const size_t plus = shear_ring_index(order, bin, max, target->bins);
            const size_t minus = shear_ring_index(-order, bin, max, target->bins);
            const size_t weight = shear_weight_ring_index(order, bin, max);
            const shear_complex positive = shear_mul(spin, power);
            const shear_complex negative = shear_mul(spin, shear_conj(power));
            for (int leg = 0; leg < (target->same_neighbor_catalog ? 1 : 2); leg++) {
                shear_complex *g = leg ? target->g_ring_second : target->g_ring_first;
                shear_complex *w = leg ? target->w_ring_second : target->w_ring_first;
                const shear_complex *sg = leg ? source->g_ring_second : source->g_ring_first;
                const shear_complex *sw = leg ? source->w_ring_second : source->w_ring_first;
                shear_accumulate_complex_product(&g[plus], sg[plus], positive);
                if (order != 0)
                    shear_accumulate_complex_product(&g[minus], sg[minus], negative);
                shear_accumulate_complex_product(&w[weight], sw[weight], power);
            }
            power = shear_mul(power, u);
        }
        shear_accumulate_complex_product(&target->diag_g6[bin],
                                          source->diag_g6[bin], shear_conj(spin));
        shear_accumulate_complex_product(&target->diag_g2[bin],
                                          source->diag_g2[bin], spin);
        shear_accumulate_complex_product(&target->diag_abs2[bin],
                                          source->diag_abs2[bin], shear_conj(spin));
        target->diag_w2[bin] += source->diag_w2[bin];
    }
    return SUCCESS;
}

static int shear_reuse_visit(shear_reuse_context *context,
        shear_pivot_workspace *target, shear_result_accumulator *result,
        shear_reuse_level *level, shear_reuse_ref pivot,
        shear_reuse_ref *first, size_t nfirst,
        shear_reuse_ref *second, size_t nsecond, unsigned depth)
{
    shear_pivot_workspace *work = &level->work;
    INTEGER represented;
    level->pivot_error = shear_reuse_error(pivot);
    work->aggregate_pivot = shear_reuse_cell(pivot);
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
    work->pivot = pivot.body;
    if (work->aggregate_pivot) {
        fcfc_ballnode *node = kd_shear_node(pivot);
        work->reuse_anchor = node;
        work->aggregate_gamma = kd_shear_weighted_gamma(pivot);
        work->aggregate_weight = node->weight;
        level->radius = node->radius;
        level->frame_valid = shear_load_unit3(node->center, work->pivot_unit)
            && shear_spherical_basis(work->pivot_unit,
                                    work->pivot_east, work->pivot_north);
        represented = node->last - node->first + 1;
    } else {
        if (!Update(pivot.body)) return SUCCESS;
        work->reuse_anchor = pivot.body;
        level->radius = 0.0;
        level->frame_valid = shear_prepare_pivot_geometry(work);
        represented = 1;
    }
#else
    if (cballs_opt_read_mask(work->cmd) && Mask(pivot) == MASK_NODE_MASKED)
        return SUCCESS;
    if (Type(pivot) != CELL && !Update(pivot)) return SUCCESS;
    work->pivot = (bodyptr)pivot;
    work->reuse_anchor = pivot;
    represented = work->aggregate_pivot ? Nb(pivot) : 1;
    level->radius = 0.0;
    level->frame_valid = work->aggregate_pivot
        ? shear_spherical_node_frame(work->cmd, pivot, (real)Size(pivot),
            work->pivot_unit, work->pivot_east, work->pivot_north,
            &level->radius, NULL)
        : shear_prepare_pivot_geometry(work);
#endif
    level->has_rings = FALSE;
    level->count[0] = level->count[1] = 0;
    memset(level->active, 0, (size_t)work->bins);
    shear_clear_pivot_workspace(work, context->ring_count,
        context->weight_ring_count, work->bins, TRUE, FALSE,
        !work->same_neighbor_catalog);
    context->nodes++;
    context->max_depth = MAX(context->max_depth, depth);
    for (int leg = 0; leg < (work->same_neighbor_catalog ? 1 : 2); leg++) {
        shear_reuse_ref *list = leg ? second : first;
        const size_t count = leg ? nsecond : nfirst;
        work->g_ring_active = leg ? work->g_ring_second : work->g_ring_first;
        work->w_ring_active = leg ? work->w_ring_second : work->w_ring_first;
        work->collect_first_leg = leg == 0;
        for (size_t index = 0; index < count; index++)
            if (shear_reuse_neighbor(context, level, list[index], leg) == FAILURE)
                return FAILURE;
    }
    context->peak_unresolved = MAX(context->peak_unresolved,
                                   level->count[0] + level->count[1]);
    if (work->aggregate_pivot && (level->count[0] || level->count[1])) {
        if (shear_reuse_allocate_level(context, target,
                                       &level->child, level) == FAILURE)
            return FAILURE;
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
        for (INTEGER index = 0; index < kd_shear_child_count(pivot); index++) {
            shear_reuse_ref child = kd_shear_child(pivot, index);
#else
        for (nodeptr child = More(pivot); child != Next(pivot); child = Next(child)) {
#endif
            if (shear_reuse_visit(context, target, result, level->child, child,
                    level->unresolved[0], level->count[0],
                    level->unresolved[1], level->count[1], depth + 1) == FAILURE)
                return FAILURE;
        }
        return SUCCESS;
    }
    if (!level->frame_valid) return SUCCESS;
    target->pivot = work->pivot;
    target->aggregate_pivot = work->aggregate_pivot;
    target->reuse_anchor = work->reuse_anchor;
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
    if (work->aggregate_pivot) {
        target->aggregate_gamma = work->aggregate_gamma;
        target->aggregate_weight = work->aggregate_weight;
    }
#endif
    memcpy(target->pivot_unit, work->pivot_unit, sizeof(compute_vector));
    memcpy(target->pivot_east, work->pivot_east, sizeof(compute_vector));
    memcpy(target->pivot_north, work->pivot_north, sizeof(compute_vector));
    shear_clear_pivot_workspace(target, context->ring_count,
        context->weight_ring_count, work->bins, TRUE, FALSE,
        !work->same_neighbor_catalog);
    for (shear_reuse_level *ancestor = level; ancestor; ancestor = ancestor->parent)
        if (shear_reuse_merge_level(context, target, ancestor) == FAILURE)
            return FAILURE;
    shear_reduce_pivot_profiled(result, target, context->nmax, FALSE, TRUE);
    context->aggregated += work->aggregate_pivot;
    context->bodies += represented;
    if (target->profile != NULL)
        target->profile->pivots += represented;
    target->aggregate_pivot = FALSE;
    return SUCCESS;
}

static int shear_reuse_task(shear_reuse_context *context,
        shear_pivot_workspace *work, shear_result_accumulator *result,
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
        const kd_shear_pivot_task *task,
        const kd_shear_frontier_schedule *schedule,
        fcfc_balltreeptr pivot, fcfc_balltreeptr first, fcfc_balltreeptr second)
#else
        const shear_native_pivot_task *task,
        const shear_native_frontier_schedule *schedule)
#endif
{
    const double started = work->profile ? shear_profile_wall_time() : 0.0;
    int status;
    if (shear_reuse_allocate_level(context, work, &context->root, NULL) == FAILURE)
        return FAILURE;
#ifdef SHEAR_SPHERE_BINARY_PIVOT_REUSE
    /* The schedule holds node indices; materialize the small root lists once
     * per task. Deeper levels inherit only unresolved typed references. */
    const size_t count = (size_t)task->first_count + (size_t)task->second_count;
    if (count > (SHEAR_REUSE_MEMORY_LIMIT - context->bytes)/sizeof(shear_reuse_ref))
        return FAILURE;
    shear_reuse_ref *roots = count ? malloc(count*sizeof(*roots)) : NULL;
    if (count && roots == NULL) return FAILURE;
    for (INTEGER i = 0; i < task->first_count; i++)
        roots[i] = kd_shear_node_ref(first, schedule->first_nodes[task->first_offset+i]);
    for (INTEGER i = 0; i < task->second_count; i++)
        roots[task->first_count+i] = kd_shear_node_ref(second,
            schedule->second_nodes[task->second_offset+i]);
    context->bytes += count*sizeof(*roots);
    status = shear_reuse_visit(context, work, result, context->root,
        kd_shear_node_ref(pivot, task->pivot_node), roots, task->first_count,
        roots ? roots + task->first_count : NULL, task->second_count, 0);
    context->bytes -= count*sizeof(*roots);
    free(roots);
#else
    status = shear_reuse_visit(context, work, result, context->root,
        task->pivot_node,
        schedule->first_nodes ? schedule->first_nodes + task->first_offset : NULL,
        (size_t)task->first_count,
        schedule->second_nodes ? schedule->second_nodes + task->second_offset : NULL,
        (size_t)task->second_count, 0);
#endif
    if (work->profile)
        work->profile->walk_seconds += shear_profile_wall_time() - started;
    return status;
}

static void shear_reuse_print_profile(const shear_reuse_context *context, int thread)
{
    fprintf(stderr, "SHEAR_REUSE thread=%d nodes=%" PRIu64
            " aggregate_pivots=%" PRIu64 " represented_bodies=%" PRIu64
            " accepted_cells=%" PRIu64 " accepted_bodies=%" PRIu64
            " ring_merges=%" PRIu64 " pruned=%" PRIu64
            " ancestor_merges=%" PRIu64
            " peak_unresolved=%" PRIu64 " max_depth=%" PRIu64
            " scratch_bytes=%zu phase_budget=%.6g\n", thread, context->nodes,
            context->aggregated, context->bodies, context->accepted_cells,
            context->accepted_bodies, context->inherited, context->pruned,
            context->ancestor_merges,
            context->peak_unresolved, context->max_depth,
            context->bytes, context->tolerance);
}

#endif
