/* Shared Fourier moments, self-term removal, transport and pivot traversal.
 * Included under THREEPCFCONVERGENCE with backend specialization in scope. */
typedef struct {
    const unsigned char *selected_pairs;
    INTEGER reuse_pairs, reuse_represented_pairs, reuse_parent_reductions;
    real *field_cos;
    real *field_sin;
    real *distinct_coscos;
    real *distinct_sinsin;
    real *distinct_sincos;
    real *normalization;
    real *normalization_distinct;
    real *neighbor_count; /* Exact 0, 1, or 2-or-more; shared scratch layout. */
    real *window_cos;
    real *window_sin;
    real *window_distinct;
    size_t stride;
    size_t values;
    int orders;
    int window_orders;
    INTEGER body_visits;
    INTEGER accepted_nodes;
    INTEGER pair_tests;
    INTEGER pivot_restarts;
    INTEGER pivot_finishes;
#ifdef DUAL_NODE_PIVOT_PROGRESS
    INTEGER progress_pending;
#endif
    double pivot_transport_seconds;
    double scratch_clear_seconds;
    double multipole_product_seconds;
} dual_node_multipole_scratch;

#ifdef DUAL_NODE_PIVOT_PROGRESS
static inline void dual_node_record_completed_pivots(
        const dual_node_search_context *context,
        dual_node_multipole_scratch *statistics, INTEGER count)
{
    if (context->progress == NULL || !context->progress->enabled) return;
    statistics->progress_pending += count;
    if (statistics->progress_pending >= context->progress->batch) {
        dual_node_publish_pivot_progress(context, statistics->progress_pending);
        statistics->progress_pending = 0;
    }
}
#endif

typedef struct {
    INTEGER *nodes;
    INTEGER count;
    INTEGER capacity;
    bool allocation_failed;
} dual_node_neighbor_frontier;

static bool dual_node_neighbor_frontier_append(
        dual_node_neighbor_frontier *frontier, INTEGER node)
{
    if (frontier->count == frontier->capacity) {
        const INTEGER capacity = frontier->capacity > 0
            ? 2 * frontier->capacity : 64;
        INTEGER *nodes;

        if (capacity < frontier->capacity
            || (size_t)capacity > SIZE_MAX / sizeof(*nodes)) {
            frontier->allocation_failed = TRUE;
            return FALSE;
        }
        nodes = realloc(frontier->nodes,
                        (size_t)capacity * sizeof(*nodes));
        if (nodes == NULL) {
            frontier->allocation_failed = TRUE;
            return FALSE;
        }
        frontier->nodes = nodes;
        frontier->capacity = capacity;
    }
    frontier->nodes[frontier->count++] = node;
    return TRUE;
}

typedef struct {
    const cballs_storage_real *position;
    bodyptr source;
    real radius;
    real field_sum;
    real normalization_sum;
    INTEGER first;
    INTEGER last;
    bool can_split;
#if NDIM == 3
    real position_norm;
    real normal[3];
    real first_axis[3];
    real second_axis[3];
    bool angular_basis_valid;
#endif
} dual_node_multipole_pivot;

static inline void dual_node_multipole_initialize_pivot_geometry(
        dual_node_multipole_pivot *pivot)
{
#if NDIM == 3 && !defined(DUAL_NODE_DISABLE_CACHED_GEOMETRY)
    pivot->position_norm = hypot(
        hypot((real)pivot->position[0], (real)pivot->position[1]),
        (real)pivot->position[2]);
    pivot->angular_basis_valid = cballs_angular_basis(
        pivot->position, pivot->normal,
        pivot->first_axis, pivot->second_axis);
#else
    (void)pivot;
#endif
}

static inline bool dual_node_multipole_angular_phase(
        const dual_node_multipole_pivot *pivot,
        const compute_vector dr, real *cosphi, real *sinphi)
{
#if NDIM == 3 && !defined(DUAL_NODE_DISABLE_CACHED_GEOMETRY)
    real x = 0.0;
    real y = 0.0;
    real norm;

    if (!pivot->angular_basis_valid) return FALSE;
    for (int k = 0; k < 3; k++) {
        x -= dr[k] * pivot->first_axis[k];
        y -= dr[k] * pivot->second_axis[k];
    }
    norm = hypot(x, y);
    if (!(norm > 32.0 * DBL_EPSILON
          * hypot(hypot(dr[0], dr[1]), dr[2])))
        return FALSE;
    *cosphi = x / norm;
    *sinphi = y / norm;
    return isfinite(*cosphi) && isfinite(*sinphi);
#else
    return cballs_angular_phase(pivot->position, dr, cosphi, sinphi);
#endif
}

static inline real dual_node_multipole_angular_extent(
        const dual_node_multipole_pivot *pivot,
        const compute_vector dr, real neighbor_radius)
{
#if NDIM == 3 && !defined(DUAL_NODE_DISABLE_CACHED_GEOMETRY)
    real cross[3];
    real distance;
    real transverse;
    real error;

    if (!pivot->angular_basis_valid
        || !(pivot->position_norm > pivot->radius))
        return 1.0;
    CROSSVP(cross, pivot->normal, dr);
    distance = hypot(hypot(dr[0], dr[1]), dr[2]);
    transverse = hypot(hypot(cross[0], cross[1]), cross[2]);
    error = pivot->radius + neighbor_radius
        + distance * pivot->radius / (pivot->position_norm - pivot->radius);
    return transverse > error ? error / transverse : 1.0;
#else
    return cballs_angular_extent(
        pivot->position, dr, pivot->radius, neighbor_radius);
#endif
}

static void dual_node_initialize_multipole_scratch(
        dual_node_multipole_scratch *, real *, size_t, int, int, size_t);

static inline real dual_node_node_field_distinct_sum(
        const dual_node_search_context *context,
        const fcfc_ballnode *ball_node)
{
    return context->weighted ? ball_node->weighted_kappa_distinct
                             : ball_node->kappa_distinct;
}

static inline real dual_node_node_normalization_distinct_sum(
        const dual_node_search_context *context,
        const fcfc_ballnode *ball_node)
{
    const real count = (real)dual_node_node_count(ball_node);
    return context->weighted ? ball_node->field_weight_distinct : count*(count-1.0);
}

static inline real dual_node_body_normalization(
        const dual_node_search_context *context, bodyptr p)
{
    return context->weighted ? Weight(p) : 1.0;
}

static inline real dual_node_body_weighted_field(
        const dual_node_search_context *context, bodyptr p)
{
    const real field = dual_node_body_field(p);

    return context->weighted ? Weight(p) * field : field;
}

static inline real dual_node_point_normalization(
        const dual_node_search_context *context,
        const fcfc_ballpoint *point)
{
    return context->weighted ? (real)point->weight : 1.0;
}

static inline real dual_node_point_weighted_field(
        const dual_node_search_context *context,
        const fcfc_ballpoint *point)
{
    return context->weighted
        ? (real)point->weighted_kappa : (real)point->kappa;
}

static void dual_node_multipole_clear(dual_node_multipole_scratch *scratch)
{
    memset(scratch->field_cos, 0, scratch->values * sizeof(real));
}

static void dual_node_multipole_clear_profiled(
        const dual_node_search_context *context,
        dual_node_multipole_scratch *scratch)
{
    const double started = context->profile ? dual_node_timer_now() : 0.0;

    dual_node_multipole_clear(scratch);
    if (context->profile)
        scratch->scratch_clear_seconds += dual_node_timer_now() - started;
}

static void dual_node_multipole_add(
        dual_node_multipole_scratch *scratch, int radial_bin,
        real field_sum, real field_distinct_sum,
        real normalization_sum, real normalization_distinct_sum,
        INTEGER neighbors, real cosphi1, real sinphi1)
{
    real cosphi = cosphi1;
    real sinphi = sinphi1;
    size_t index = (size_t)radial_bin;

    /* Arguments and distinct_* slots carry internal ordered distinct products.
     * Singleton bodies supply zero; accepted nodes supply stable tree merges. */
    scratch->normalization_distinct[index] = cballs_distinct_merge(
        scratch->normalization_distinct[index], normalization_distinct_sum,
        scratch->normalization[index], normalization_sum);
    scratch->normalization[index] += normalization_sum;
    scratch->neighbor_count[index] = MIN(2.0,
        scratch->neighbor_count[index] + (real)MIN(neighbors, 2));
    for (int order = 0;
         order < MAX(scratch->orders, scratch->window_orders); order++) {
        const real c = order ? cosphi : 1.0;
        const real sn = order ? sinphi : 0.0;
        if (order < scratch->orders) {
            const real x = cballs_product(field_sum, c);
            const real y = cballs_product(field_sum, sn);
            const real px = scratch->field_cos[index], py = scratch->field_sin[index];
            scratch->distinct_coscos[index] = cballs_distinct_merge(
                scratch->distinct_coscos[index], field_distinct_sum*c*c, px, x);
            scratch->distinct_sinsin[index] = cballs_distinct_merge(
                scratch->distinct_sinsin[index], field_distinct_sum*sn*sn, py, y);
            scratch->distinct_sincos[index] += field_distinct_sum*c*sn
                + cballs_product(px,y) + cballs_product(py,x);
            scratch->field_cos[index] += x;
            scratch->field_sin[index] += y;
        }
        if (order < scratch->window_orders) {
            const real x = normalization_sum*c, y = normalization_sum*sn;
            scratch->window_distinct[index] += normalization_distinct_sum
                + cballs_distinct_complex_cross(scratch->window_cos[index],
                    scratch->window_sin[index], x, y);
            scratch->window_cos[index] += x;
            scratch->window_sin[index] += y;
        }
        if (order) {
            const real next_cos = cosphi*cosphi1-sinphi*sinphi1;
            sinphi = sinphi*cosphi1+cosphi*sinphi1;
            cosphi = next_cos;
        }
        index += scratch->stride;
    }
}

/* Individual neighbors have no internal distinct pairs. Keep this hot path
 * separate from aggregate-node/window accumulation, and expose the disjoint
 * scratch arrays to the compiler. The ordered distinct-product recurrence and
 * checked multiplications are unchanged. */
static inline void dual_node_multipole_add_singleton(
        dual_node_multipole_scratch *scratch, int radial_bin,
        real field, real normalization, real cosphi1, real sinphi1)
{
    if (scratch->window_orders > 0) {
        dual_node_multipole_add(scratch, radial_bin, field, 0.0,
                               normalization, 0.0, 1, cosphi1, sinphi1);
        return;
    }
    const size_t stride = scratch->stride;
    size_t index = (size_t)radial_bin;
    real *restrict field_cos = scratch->field_cos;
    real *restrict field_sin = scratch->field_sin;
    real *restrict distinct_coscos = scratch->distinct_coscos;
    real *restrict distinct_sinsin = scratch->distinct_sinsin;
    real *restrict distinct_sincos = scratch->distinct_sincos;
    real cosphi = cosphi1, sinphi = sinphi1;

    scratch->normalization_distinct[index] = cballs_distinct_merge(
        scratch->normalization_distinct[index], 0.0,
        scratch->normalization[index], normalization);
    scratch->normalization[index] += normalization;
    scratch->neighbor_count[index] = MIN(2.0, scratch->neighbor_count[index]+1.0);
    if (scratch->orders == 0) return;

    /* The monopole phase is exactly (1,0). Its sine components stay zero;
     * the cosine recurrence still detects overflow/underflow of products. */
    distinct_coscos[index] = cballs_distinct_merge(
        distinct_coscos[index], 0.0, field_cos[index], field);
    field_cos[index] += field;
    index += stride;
    for (int order = 1; order < scratch->orders; order++, index += stride) {
        const real px = field_cos[index], py = field_sin[index];
        real x = field*cosphi, y = field*sinphi;
        real xx = px*x, yy = py*y, xy = px*y, yx = py*x;
        /* A nonzero product needs no underflow-to-zero handling. Collect the
         * rare zero cases before checking their operands, retaining the same
         * checked arithmetic (including NaN propagation) on that path. */
        if (x == 0.0 || y == 0.0 || xx == 0.0 || yy == 0.0
            || xy == 0.0 || yx == 0.0) {
            x = cballs_product(field, cosphi);
            y = cballs_product(field, sinphi);
            xx = cballs_product(px,x); yy = cballs_product(py,y);
            xy = cballs_product(px,y); yx = cballs_product(py,x);
        }
        distinct_coscos[index] = distinct_coscos[index]+0.0+2.0*xx;
        distinct_sinsin[index] = distinct_sinsin[index]+0.0+2.0*yy;
        distinct_sincos[index] += xy+yx;
        field_cos[index] += x;
        field_sin[index] += y;
        const real next_cos = cosphi*cosphi1-sinphi*sinphi1;
        sinphi = sinphi*cosphi1+cosphi*sinphi1;
        cosphi = next_cos;
    }
}

static int dual_node_multipole_pair_status_limited(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        const cballs_storage_real *neighbor_position, real neighbor_radius,
        real upper_limit, int *radial_bin,
        real *cosphi, real *sinphi, real *upper_extent)
{
    compute_vector dr;
    real distance;
    const real size = pivot->radius + neighbor_radius;
    const real lower_limit = context->cmd->rminHist;

    distance = dual_node_position_distance(
        context, pivot->position, neighbor_position, dr);
    *upper_extent = MIN(upper_limit, distance + size);
    if (!(distance > 0.0))
        return DUAL_NODE_TRIPLE_SPLIT;
    if (distance + size <= lower_limit
        || distance - size >= upper_limit)
        return DUAL_NODE_TRIPLE_OUTSIDE;
    if (!context->use_two_balls || !(distance > size)
        || size > dual_node_bin_theta_width(context, distance)
        || dual_node_multipole_angular_extent(pivot, dr, neighbor_radius)
            > context->max_angular_ratio)
        return DUAL_NODE_TRIPLE_SPLIT;
    *radial_bin = dual_node_bin_index(context, distance);
    if (*radial_bin < 0 || !(distance < upper_limit))
        return DUAL_NODE_TRIPLE_OUTSIDE;
    if (!dual_node_multipole_angular_phase(
            pivot, dr, cosphi, sinphi))
        return DUAL_NODE_TRIPLE_SPLIT;
    return DUAL_NODE_TRIPLE_ACCEPT;
}

static int dual_node_multipole_pair_status(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        const cballs_storage_real *neighbor_position, real neighbor_radius,
        int *radial_bin, real *cosphi, real *sinphi)
{
    real upper_extent;

    return dual_node_multipole_pair_status_limited(
        context, pivot, neighbor_position, neighbor_radius,
        context->cmd->rangeN, radial_bin, cosphi, sinphi, &upper_extent);
}

static int dual_node_multipole_add_body(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        dual_node_multipole_scratch *scratch,
        const fcfc_ballpoint *neighbor)
{
    compute_vector dr;
    real distance;
    real cosphi;
    real sinphi;
    int radial_bin;

    if (pivot->source != NULL && pivot->source == neighbor->source)
        return SUCCESS;
    if (pivot->radius > 0.0) {
        scratch->pair_tests++;
        const int pair_status = dual_node_multipole_pair_status(
            context, pivot, neighbor->pos, 0.0,
            &radial_bin, &cosphi, &sinphi);
        if (pair_status == DUAL_NODE_TRIPLE_OUTSIDE) return SUCCESS;
        if (pair_status != DUAL_NODE_TRIPLE_ACCEPT) return FAILURE;
    } else {
        distance = dual_node_position_distance(
            context, pivot->position, neighbor->pos, dr);
        if (!(distance > 0.0)) return SUCCESS;
        radial_bin = dual_node_bin_index(context, distance);
        if (radial_bin < 0) return SUCCESS;
        /* Out-of-range neighbors contribute no moment: avoid constructing
         * their angular phase, without changing the bin or phase formulas. */
        if (!dual_node_multipole_angular_phase(
                pivot, dr, &cosphi, &sinphi))
            return SUCCESS;
    }
    {
        const real field = context->weighted
            ? (real)neighbor->weighted_kappa : (real)neighbor->kappa;
        const real normalization = context->weighted
            ? (real)neighbor->weight : 1.0;

        dual_node_multipole_add_singleton(
            scratch, radial_bin, field, normalization, cosphi, sinphi);
    }
    scratch->body_visits++;
    return SUCCESS;
}

static bool dual_node_ranges_overlap(INTEGER first1, INTEGER last1,
                                    INTEGER first2, INTEGER last2)
{
    return first1 <= last2 && first2 <= last1;
}

/*
 * Scan one neighbor tree for a pivot ball.  FAILURE means that the neighbor
 * side can no longer resolve the requested radial/angular accuracy, so the
 * caller must discard this scratch buffer and split the pivot ball.
 */
static int dual_node_multipole_scan_neighbors(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        const fcfc_balltreeptr neighbor_tree, INTEGER neighbor_index,
        bool same_tree, dual_node_multipole_scratch *scratch)
{
    const fcfc_ballnode *neighbor = &neighbor_tree->nodes[neighbor_index];
    int radial_bin;
    real cosphi;
    real sinphi;
    int pair_status;

    if (same_tree && dual_node_ranges_overlap(
            pivot->first, pivot->last, neighbor->first, neighbor->last)) {
        if (neighbor->first >= pivot->first
            && neighbor->last <= pivot->last)
            return SUCCESS;
        if (dual_node_node_is_leaf(neighbor)) {
            for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
                if (i >= pivot->first && i <= pivot->last) continue;
                if (dual_node_multipole_add_body(
                        context, pivot, scratch,
                        &neighbor_tree->packed_points[i]) == FAILURE)
                    return FAILURE;
            }
            return SUCCESS;
        }
        if (dual_node_multipole_scan_neighbors(
                context, pivot, neighbor_tree, neighbor->left,
                same_tree, scratch) == FAILURE)
            return FAILURE;
        return dual_node_multipole_scan_neighbors(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch);
    }

    scratch->pair_tests++;
    pair_status = dual_node_multipole_pair_status(
        context, pivot, neighbor->center, (real)neighbor->radius,
        &radial_bin, &cosphi, &sinphi);
    if (pair_status == DUAL_NODE_TRIPLE_OUTSIDE) return SUCCESS;
    if (pair_status == DUAL_NODE_TRIPLE_ACCEPT) {
        dual_node_multipole_add(
            scratch, radial_bin,
            dual_node_node_field_sum(context, neighbor),
            dual_node_node_field_distinct_sum(context, neighbor),
            dual_node_node_normalization_sum(context, neighbor),
            dual_node_node_normalization_distinct_sum(context, neighbor),
            dual_node_node_count(neighbor), cosphi, sinphi);
        scratch->accepted_nodes++;
        return SUCCESS;
    }

    if (!dual_node_node_is_leaf(neighbor)) {
        if (pivot->can_split
            && pivot->radius > (real)neighbor->radius)
            return FAILURE;
        if (dual_node_multipole_scan_neighbors(
                context, pivot, neighbor_tree, neighbor->left,
                same_tree, scratch) == FAILURE)
            return FAILURE;
        return dual_node_multipole_scan_neighbors(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch);
    }

    for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
        if (same_tree && i >= pivot->first && i <= pivot->last) continue;
        if (dual_node_multipole_add_body(
                context, pivot, scratch,
                &neighbor_tree->packed_points[i]) == FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}

static void dual_node_multipole_finish_pivot_range(
        dual_node_triple_histogram *hist,
        const dual_node_multipole_pivot *pivot,
        const dual_node_multipole_scratch *scratch,
        int first_radial_bin, int radial_bin_limit, int radial_bins)
{
    /* Keep each order's contiguous histogram row hot. The former pair-major
     * loop jumped between all order/component planes for every radial pair.
     * Each SIMD iteration owns one distinct upper/lower symmetric pair. */
    const size_t stride=hist->stride;
    const unsigned char *selected=scratch->selected_pairs;
#define DUAL_NODE_PAIR_PRESENT(n1,n2) \
    (scratch->neighbor_count[n2] > 0.0 \
     && ((n1)!=(n2) || scratch->neighbor_count[n1]>=2.0) \
     && (selected==NULL || selected[(size_t)(n1)*stride+(n2)]))
    for (int n1=first_radial_bin;n1<radial_bin_limit;n1++) {
        if (scratch->neighbor_count[n1]==0.0) continue;
        const size_t row=(size_t)n1*stride;
        for (int n2=n1;n2<=radial_bins;n2++) {
            if (!DUAL_NODE_PAIR_PRESENT(n1,n2)) continue;
            real denominator=n1==n2 ? scratch->normalization_distinct[n1]
                : cballs_product(scratch->normalization[n1],scratch->normalization[n2]);
            hist->normalization[row+n2]+=cballs_product(pivot->normalization_sum,denominator);
            if (n1!=n2) hist->normalization[(size_t)n2*stride+n1]+=cballs_product(pivot->normalization_sum,denominator);
        }
        for (int order=0;order<scratch->window_orders;order++) {
            const size_t offset=(size_t)order*scratch->stride;
            const real x=scratch->window_cos[offset+n1],y=scratch->window_sin[offset+n1];
            real *restrict re=hist->window_re+(size_t)order*hist->plane;
            real *restrict im=hist->window_im+(size_t)order*hist->plane;
#pragma omp simd
            for (int n2=n1;n2<=radial_bins;n2++) {
                if (!DUAL_NODE_PAIR_PRESENT(n1,n2)) continue;
                const real a=n1==n2 ? scratch->window_distinct[offset+n1]
                    : cballs_product(x,scratch->window_cos[offset+n2])+cballs_product(y,scratch->window_sin[offset+n2]);
                /* A complex diagonal is real by definition, including under FMA. */
                const real b=n1==n2 ? 0.0
                    : cballs_product(y,scratch->window_cos[offset+n2])-cballs_product(x,scratch->window_sin[offset+n2]);
                re[row+n2]+=cballs_product(pivot->normalization_sum,a);
                im[row+n2]+=cballs_product(pivot->normalization_sum,b);
                if (n1!=n2) {
                    re[(size_t)n2*stride+n1]+=cballs_product(pivot->normalization_sum,a);
                    im[(size_t)n2*stride+n1]-=cballs_product(pivot->normalization_sum,b);
                }
            }
        }
        real *restrict monopole=dual_node_zeta_component(hist,DUAL_NODE_ZETA_COS,0);
#pragma omp simd
        for (int n2=n1;n2<=radial_bins;n2++) {
            if (!DUAL_NODE_PAIR_PRESENT(n1,n2)) continue;
            real value=n1==n2 ? scratch->distinct_coscos[n1]
                : cballs_product(scratch->field_cos[n1],scratch->field_cos[n2]);
            monopole[row+n2]+=cballs_product(pivot->field_sum,value);
            if (n1!=n2) monopole[(size_t)n2*stride+n1]+=cballs_product(pivot->field_sum,value);
        }
        for (int order=1;order<scratch->orders;order++) {
            const size_t offset=(size_t)order*scratch->stride;
            const real x=scratch->field_cos[offset+n1],y=scratch->field_sin[offset+n1];
            real *restrict cc=dual_node_zeta_component(hist,DUAL_NODE_ZETA_COS,order);
            real *restrict ss=dual_node_zeta_component(hist,DUAL_NODE_ZETA_SIN,order);
            real *restrict sc=dual_node_zeta_component(hist,DUAL_NODE_ZETA_SINCOS,order);
            real *restrict cs=dual_node_zeta_component(hist,DUAL_NODE_ZETA_COSSIN,order);
#pragma omp simd
            for (int n2=n1;n2<=radial_bins;n2++) {
                if (!DUAL_NODE_PAIR_PRESENT(n1,n2)) continue;
                real a=n1==n2 ? 0.0 : cballs_product(x,scratch->field_cos[offset+n2]);
                real b=n1==n2 ? 0.0 : cballs_product(y,scratch->field_sin[offset+n2]);
                real c=n1==n2 ? 0.0 : cballs_product(y,scratch->field_cos[offset+n2]);
                real d=n1==n2 ? 0.0 : cballs_product(x,scratch->field_sin[offset+n2]);
                if (n1==n2) {
                    a=scratch->distinct_coscos[offset+n1];b=scratch->distinct_sinsin[offset+n1];
                    c=scratch->distinct_sincos[offset+n1];d=scratch->distinct_sincos[offset+n1];
                }
                cc[row+n2]+=cballs_product(pivot->field_sum,a);ss[row+n2]+=cballs_product(pivot->field_sum,b);
                sc[row+n2]+=cballs_product(pivot->field_sum,c);cs[row+n2]+=cballs_product(pivot->field_sum,d);
                if (n1!=n2) {
                    const size_t column=(size_t)n2*stride+n1;
                    cc[column]+=cballs_product(pivot->field_sum,a);ss[column]+=cballs_product(pivot->field_sum,b);
                    sc[column]+=cballs_product(pivot->field_sum,d);cs[column]+=cballs_product(pivot->field_sum,c);
                }
            }
        }
    }
#undef DUAL_NODE_PAIR_PRESENT
}

static void dual_node_multipole_finish_pivot(
        const dual_node_search_context *context,
        dual_node_triple_histogram *hist,
        const dual_node_multipole_pivot *pivot,
        dual_node_multipole_scratch *scratch, int radial_bins)
{
    const double started = context->profile ? dual_node_timer_now() : 0.0;

    dual_node_multipole_finish_pivot_range(
        hist, pivot, scratch, 1, radial_bins + 1, radial_bins);
    if (context->profile)
        scratch->multipole_product_seconds += dual_node_timer_now() - started;
    scratch->pivot_finishes++;
#ifdef DUAL_NODE_PIVOT_PROGRESS
    dual_node_record_completed_pivots(context, scratch, 1);
#endif
}

static void dual_node_multipole_finish_pivot_range_profiled(
        const dual_node_search_context *context,
        dual_node_triple_histogram *hist,
        const dual_node_multipole_pivot *pivot,
        dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics,
        int first_radial_bin, int radial_bin_limit, int radial_bins)
{
    const double started = context->profile ? dual_node_timer_now() : 0.0;

    dual_node_multipole_finish_pivot_range(
        hist, pivot, scratch, first_radial_bin,
        radial_bin_limit, radial_bins);
    if (context->profile)
        statistics->multipole_product_seconds += dual_node_timer_now() - started;
}

static void dual_node_multipole_process_body(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        dual_node_multipole_scratch *scratch,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballpoint *pivot_point =
        &pivot_tree->packed_points[pivot_index];
    dual_node_multipole_pivot pivot;

    pivot.position = pivot_point->pos;
    pivot.source = pivot_point->source;
    pivot.radius = 0.0;
    pivot.field_sum = dual_node_point_weighted_field(context, pivot_point);
    pivot.normalization_sum =
        dual_node_point_normalization(context, pivot_point);
    pivot.first = pivot_index;
    pivot.last = pivot_index;
    pivot.can_split = FALSE;
    dual_node_multipole_initialize_pivot_geometry(&pivot);
    dual_node_multipole_clear_profiled(context, scratch);
    if (dual_node_multipole_scan_neighbors(
            context, &pivot, neighbor_tree, 0,
            same_tree, scratch) == SUCCESS)
        dual_node_multipole_finish_pivot(
            context, hist, &pivot, scratch, context->cmd->sizeHistN);
}

#ifdef DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE
static void dual_node_multipole_process_body_pivots(
        const dual_node_search_context *, const fcfc_balltreeptr, INTEGER,
        const fcfc_balltreeptr, bool, dual_node_multipole_scratch *,
        dual_node_triple_histogram *);

#ifdef DUAL_NODE_PERSISTENT_NEIGHBOR_FRONTIER
#ifndef DUAL_NODE_ADAPTIVE_FRONTIER_BASE
#define DUAL_NODE_ADAPTIVE_FRONTIER_BASE ((INTEGER)32)
#endif
#ifndef DUAL_NODE_ADAPTIVE_FRONTIER_MAX
#define DUAL_NODE_ADAPTIVE_FRONTIER_MAX ((INTEGER)256)
#endif
static INTEGER dual_node_adaptive_frontier_limit(
        const dual_node_search_context *context,
        const fcfc_ballnode *pivot_node,
        const fcfc_balltreeptr neighbor_tree)
{
    const real pivot_count = (real)dual_node_node_count(pivot_node);
    const real radial_fraction = MIN(
        1.0, 0.25 * rsqr(context->cmd->rangeN));
    const real expected_neighbors = MAX(
        1.0, radial_fraction * (real)neighbor_tree->npoint);
    const real pivot_scale = MAX(1.0, rlog2(pivot_count + 1.0));
    const INTEGER limit = DUAL_NODE_ADAPTIVE_FRONTIER_BASE
        + (INTEGER)(4.0 * rsqrt(expected_neighbors)
                    + 2.0 * pivot_scale);

    return MIN(DUAL_NODE_ADAPTIVE_FRONTIER_MAX,
               MAX(DUAL_NODE_ADAPTIVE_FRONTIER_BASE, limit));
}

static bool dual_node_neighbor_frontier_append_bounded(
        dual_node_neighbor_frontier *frontier, INTEGER node, INTEGER limit)
{
    return frontier->count < limit
        && dual_node_neighbor_frontier_append(frontier, node);
}

static bool dual_node_neighbor_frontier_refine_node(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        const fcfc_balltreeptr neighbor_tree, INTEGER neighbor_index,
        bool same_tree, dual_node_neighbor_frontier *frontier, INTEGER limit)
{
    const fcfc_ballnode *neighbor = &neighbor_tree->nodes[neighbor_index];
    const bool overlap = same_tree && dual_node_ranges_overlap(
        pivot->first, pivot->last, neighbor->first, neighbor->last);
    compute_vector dr;
    const real distance = dual_node_position_distance(
        context, pivot->position, neighbor->center, dr);
    const real size = pivot->radius + (real)neighbor->radius;
    const bool outside = !overlap
        && (distance + size <= context->cmd->rminHist
            || distance - size >= context->cmd->rangeN);
    const real split_radius = MAX(
        pivot->radius, 0.125 * context->cmd->rangeN);

    if (outside) return TRUE;
    if (dual_node_node_is_leaf(neighbor)
        || (!overlap && (real)neighbor->radius <= split_radius))
        return dual_node_neighbor_frontier_append_bounded(
            frontier, neighbor_index, limit);
    /* Preserve a coarse candidate when the bounded sparse list cannot hold
     * both refinements.  Descendant pivots may split it after pruning more
     * of the surrounding volume. */
    if (frontier->count + 1 >= limit)
        return dual_node_neighbor_frontier_append_bounded(
            frontier, neighbor_index, limit);
    if (!dual_node_neighbor_frontier_refine_node(
            context, pivot, neighbor_tree, neighbor->left,
            same_tree, frontier, limit - 1))
        return FALSE;
    return dual_node_neighbor_frontier_refine_node(
        context, pivot, neighbor_tree, neighbor->right,
        same_tree, frontier, limit);
}

static bool dual_node_neighbor_frontier_refine(
        const dual_node_search_context *context,
        const fcfc_ballnode *pivot_node,
        const fcfc_balltreeptr neighbor_tree,
        const INTEGER *parent_nodes, INTEGER parent_count,
        bool same_tree, dual_node_neighbor_frontier *frontier, INTEGER limit)
{
    dual_node_multipole_pivot pivot;

    pivot.position = pivot_node->center;
    pivot.source = NULL;
    pivot.radius = (real)pivot_node->radius;
    pivot.field_sum = 0.0;
    pivot.normalization_sum = 0.0;
    pivot.first = pivot_node->first;
    pivot.last = pivot_node->last;
    pivot.can_split = !dual_node_node_is_leaf(pivot_node);
    dual_node_multipole_initialize_pivot_geometry(&pivot);
    frontier->count = 0;
    frontier->allocation_failed = FALSE;
    for (INTEGER i = 0; i < parent_count; i++) {
        /* Reserve one slot for each untouched disjoint parent subtree. */
        const INTEGER local_limit =
            limit - (parent_count - i - 1);

        if (!dual_node_neighbor_frontier_refine_node(
                context, &pivot, neighbor_tree, parent_nodes[i],
                same_tree, frontier, local_limit))
            return FALSE;
    }
    return TRUE;
}

static void dual_node_multipole_process_body_frontier(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        const dual_node_neighbor_frontier *frontier,
        dual_node_multipole_scratch *scratch,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballpoint *pivot_point =
        &pivot_tree->packed_points[pivot_index];
    dual_node_multipole_pivot pivot;

    pivot.position = pivot_point->pos;
    pivot.source = pivot_point->source;
    pivot.radius = 0.0;
    pivot.field_sum = dual_node_point_weighted_field(context, pivot_point);
    pivot.normalization_sum =
        dual_node_point_normalization(context, pivot_point);
    pivot.first = pivot_index;
    pivot.last = pivot_index;
    pivot.can_split = FALSE;
    dual_node_multipole_initialize_pivot_geometry(&pivot);
    dual_node_multipole_clear_profiled(context, scratch);
    for (INTEGER i = 0; i < frontier->count; i++) {
        if (dual_node_multipole_scan_neighbors(
                context, &pivot, neighbor_tree, frontier->nodes[i],
                same_tree, scratch) == FAILURE)
            return;
    }
    dual_node_multipole_finish_pivot(
        context, hist, &pivot, scratch, context->cmd->sizeHistN);
}

static void dual_node_multipole_process_body_pivots_frontier(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        const INTEGER *parent_nodes, INTEGER parent_count,
        dual_node_neighbor_frontier *levels, int depth,
        dual_node_multipole_scratch *scratch,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node = &pivot_tree->nodes[pivot_node_index];
    dual_node_neighbor_frontier *frontier = &levels[depth];
    const INTEGER frontier_limit = MAX(
        parent_count,
        dual_node_adaptive_frontier_limit(
            context, pivot_node, neighbor_tree));

    if (!dual_node_neighbor_frontier_refine(
            context, pivot_node, neighbor_tree,
            parent_nodes, parent_count, same_tree, frontier, frontier_limit)) {
        scratch->pivot_restarts += dual_node_node_count(pivot_node);
        dual_node_multipole_process_body_pivots(
            context, pivot_tree, pivot_node_index, neighbor_tree,
            same_tree, scratch, hist);
        return;
    }
    if (!dual_node_node_is_leaf(pivot_node)) {
        dual_node_multipole_process_body_pivots_frontier(
            context, pivot_tree, pivot_node->left, neighbor_tree,
            same_tree, frontier->nodes, frontier->count,
            levels, depth + 1, scratch, hist);
        dual_node_multipole_process_body_pivots_frontier(
            context, pivot_tree, pivot_node->right, neighbor_tree,
            same_tree, frontier->nodes, frontier->count,
            levels, depth + 1, scratch, hist);
        return;
    }
    for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
        dual_node_multipole_process_body_frontier(
            context, pivot_tree, i, neighbor_tree,
            same_tree, frontier, scratch, hist);
}
#endif

static void dual_node_multipole_process_body_pivots(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        dual_node_multipole_scratch *scratch,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node =
        &pivot_tree->nodes[pivot_node_index];

    if (!dual_node_node_is_leaf(pivot_node)) {
        dual_node_multipole_process_body_pivots(
            context, pivot_tree, pivot_node->left,
            neighbor_tree, same_tree, scratch, hist);
        dual_node_multipole_process_body_pivots(
            context, pivot_tree, pivot_node->right,
            neighbor_tree, same_tree, scratch, hist);
        return;
    }
    for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
        dual_node_multipole_process_body(
            context, pivot_tree, i, neighbor_tree,
            same_tree, scratch, hist);
}
#endif

static void dual_node_multipole_process_pivots(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        dual_node_multipole_scratch *scratch,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node =
        &pivot_tree->nodes[pivot_node_index];

    if (!context->use_two_balls
        || (same_tree
            && 2.0 * (real)pivot_node->radius > context->cmd->rminHist)) {
        if (!dual_node_node_is_leaf(pivot_node)) {
            dual_node_multipole_process_pivots(
                context, pivot_tree, pivot_node->left,
                neighbor_tree, same_tree, scratch, hist);
            dual_node_multipole_process_pivots(
                context, pivot_tree, pivot_node->right,
                neighbor_tree, same_tree, scratch, hist);
        } else {
            for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
                dual_node_multipole_process_body(
                    context, pivot_tree, i, neighbor_tree,
                    same_tree, scratch, hist);
        }
        return;
    }

    {
        dual_node_multipole_pivot pivot;

        pivot.position = pivot_node->center;
        pivot.source = NULL;
        pivot.radius = (real)pivot_node->radius;
        pivot.field_sum = dual_node_node_field_sum(context, pivot_node);
        pivot.normalization_sum =
            dual_node_node_normalization_sum(context, pivot_node);
        pivot.first = pivot_node->first;
        pivot.last = pivot_node->last;
        pivot.can_split = !dual_node_node_is_leaf(pivot_node);
        dual_node_multipole_initialize_pivot_geometry(&pivot);
        dual_node_multipole_clear_profiled(context, scratch);
        if (dual_node_multipole_scan_neighbors(
                context, &pivot, neighbor_tree, 0,
                same_tree, scratch) == SUCCESS) {
            dual_node_multipole_finish_pivot(
                context, hist, &pivot, scratch, context->cmd->sizeHistN);
            return;
        }
        scratch->pivot_restarts++;
    }

    if (!dual_node_node_is_leaf(pivot_node)) {
        dual_node_multipole_process_pivots(
            context, pivot_tree, pivot_node->left,
            neighbor_tree, same_tree, scratch, hist);
        dual_node_multipole_process_pivots(
            context, pivot_tree, pivot_node->right,
            neighbor_tree, same_tree, scratch, hist);
    } else {
        for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++)
            dual_node_multipole_process_body(
                context, pivot_tree, i, neighbor_tree,
                same_tree, scratch, hist);
    }
}

static real dual_node_radial_upper_edge(
        const dual_node_search_context *context, int radial_bin)
{
    const struct cmdline_data *cmd = context->cmd;

    if (radial_bin >= cmd->sizeHistN) return cmd->rangeN;
    if (radial_bin <= 0) return cmd->rminHist;
    if (!cmd->useLogHist)
        return cmd->rminHist + (real)radial_bin * context->gd->deltaR;
    if (cmd->rminHist > 0.0)
        return cmd->rminHist
             * rpow(10.0, (real)radial_bin * context->gd->deltaR);
    return rpow(10.0,
        ((real)radial_bin - (real)cmd->sizeHistN) / cmd->logHistBinsPD
        + rlog10(cmd->rangeN));
}

static int dual_node_unresolved_radial_bin(
        const dual_node_search_context *context,
        real upper_extent, int radial_bin_limit)
{
    int radial_bin;

    if (!(upper_extent > context->cmd->rminHist)) return 0;
    if (upper_extent >= context->cmd->rangeN)
        return MIN(context->cmd->sizeHistN, radial_bin_limit - 1);
    radial_bin = dual_node_bin_index(context, upper_extent);
    if (radial_bin < 1) radial_bin = radial_bin_limit - 1;
    return MIN(radial_bin, radial_bin_limit - 1);
}

static void dual_node_multipole_mark_unresolved(
        const dual_node_search_context *context,
        real upper_extent, real *max_unresolved_extent)
{
    if (upper_extent > context->cmd->rminHist)
        *max_unresolved_extent = MAX(*max_unresolved_extent, upper_extent);
}

static void dual_node_multipole_clear_radial_through(
        const dual_node_search_context *context,
        dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics, int radial_bin)
{
    const double started = context->profile ? dual_node_timer_now() : 0.0;
    const size_t count = (size_t)radial_bin + 1;
    real *order_arrays[5] = {
        scratch->field_cos, scratch->field_sin,
        scratch->distinct_coscos, scratch->distinct_sinsin,
        scratch->distinct_sincos
    };

    for (size_t array = 0; array < 5; array++)
        for (int order = 0; order < scratch->orders; order++)
            memset(order_arrays[array] + (size_t)order * scratch->stride,
                   0, count * sizeof(real));
    memset(scratch->normalization, 0, count * sizeof(real));
    memset(scratch->normalization_distinct, 0, count * sizeof(real));
    memset(scratch->neighbor_count, 0, count * sizeof(real));
    for (int order = 0; order < scratch->window_orders; order++) {
        memset(scratch->window_cos + (size_t)order * scratch->stride,
               0, count * sizeof(real));
        memset(scratch->window_sin + (size_t)order * scratch->stride,
               0, count * sizeof(real));
        memset(scratch->window_distinct + (size_t)order * scratch->stride,
               0, count * sizeof(real));
    }
    if (context->profile)
        statistics->scratch_clear_seconds += dual_node_timer_now() - started;
}

static void dual_node_multipole_copy_values(
        dual_node_multipole_scratch *destination,
        const dual_node_multipole_scratch *source)
{
    memcpy(destination->field_cos, source->field_cos,
           source->values * sizeof(real));
}

static void dual_node_multipole_scan_body_partial(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics,
        const fcfc_ballpoint *neighbor, real active_upper,
        real *max_unresolved_extent)
{
    compute_vector dr;
    real distance;
    real cosphi;
    real sinphi;
    real upper_extent;
    int radial_bin;

    if (pivot->source != NULL && pivot->source == neighbor->source) return;
    if (pivot->radius > 0.0) {
        const int pair_status = dual_node_multipole_pair_status_limited(
            context, pivot, neighbor->pos, 0.0, active_upper,
            &radial_bin, &cosphi, &sinphi, &upper_extent);

        statistics->pair_tests++;
        if (pair_status == DUAL_NODE_TRIPLE_OUTSIDE) return;
        if (pair_status != DUAL_NODE_TRIPLE_ACCEPT) {
            dual_node_multipole_mark_unresolved(
                context, upper_extent, max_unresolved_extent);
            return;
        }
    } else {
        distance = dual_node_position_distance(
            context, pivot->position, neighbor->pos, dr);
        if (!(distance > 0.0)
            || !dual_node_multipole_angular_phase(
                pivot, dr, &cosphi, &sinphi)
            || !(distance < active_upper))
            return;
        radial_bin = dual_node_bin_index(context, distance);
        if (radial_bin < 0) return;
    }
    {
        const real field = dual_node_point_weighted_field(context, neighbor);
        const real normalization =
            dual_node_point_normalization(context, neighbor);

        dual_node_multipole_add(
            scratch, radial_bin, field, 0.0,
            normalization, 0.0,
            1, cosphi, sinphi);
    }
    statistics->body_visits++;
}

static void dual_node_multipole_scan_neighbors_partial(
        const dual_node_search_context *context,
        const dual_node_multipole_pivot *pivot,
        const fcfc_balltreeptr neighbor_tree, INTEGER neighbor_index,
        bool same_tree, dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics,
        real active_upper, real *max_unresolved_extent,
        dual_node_neighbor_frontier *unresolved)
{
    const fcfc_ballnode *neighbor = &neighbor_tree->nodes[neighbor_index];
    int radial_bin;
    real cosphi;
    real sinphi;
    real upper_extent;
    int pair_status;

    if (same_tree && dual_node_ranges_overlap(
            pivot->first, pivot->last, neighbor->first, neighbor->last)) {
        if (neighbor->first >= pivot->first
            && neighbor->last <= pivot->last) {
            dual_node_multipole_mark_unresolved(
                context, MIN(active_upper, 2.0 * pivot->radius),
                max_unresolved_extent);
            if (unresolved != NULL)
                dual_node_neighbor_frontier_append(
                    unresolved, neighbor_index);
            return;
        }
        if (dual_node_node_is_leaf(neighbor)) {
            for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
                if (i >= pivot->first && i <= pivot->last) continue;
                dual_node_multipole_scan_body_partial(
                    context, pivot, scratch, statistics,
                    &neighbor_tree->packed_points[i], active_upper,
                    max_unresolved_extent);
            }
            return;
        }
        dual_node_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->left,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent, unresolved);
        if (unresolved != NULL && unresolved->allocation_failed) return;
        dual_node_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent, unresolved);
        return;
    }

    statistics->pair_tests++;
    pair_status = dual_node_multipole_pair_status_limited(
        context, pivot, neighbor->center, (real)neighbor->radius,
        active_upper, &radial_bin, &cosphi, &sinphi, &upper_extent);
    if (pair_status == DUAL_NODE_TRIPLE_OUTSIDE) return;
    if (pair_status == DUAL_NODE_TRIPLE_ACCEPT) {
        dual_node_multipole_add(
            scratch, radial_bin,
            dual_node_node_field_sum(context, neighbor),
            dual_node_node_field_distinct_sum(context, neighbor),
            dual_node_node_normalization_sum(context, neighbor),
            dual_node_node_normalization_distinct_sum(context, neighbor),
            dual_node_node_count(neighbor), cosphi, sinphi);
        statistics->accepted_nodes++;
        /* Acceptance is provisional until every neighbor has established the
         * completed radial range. If a different node leaves this bin
         * unresolved, the child must re-evaluate this node in its own basis. */
        if (unresolved != NULL)
            dual_node_neighbor_frontier_append(unresolved, neighbor_index);
        return;
    }

    if (!pivot->can_split && pivot->radius > 0.0 && unresolved != NULL) {
        dual_node_multipole_mark_unresolved(
            context, upper_extent, max_unresolved_extent);
        dual_node_neighbor_frontier_append(unresolved, neighbor_index);
        return;
    }

    if (!dual_node_node_is_leaf(neighbor)
        && (!pivot->can_split
            || (real)neighbor->radius > pivot->radius)) {
        dual_node_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->left,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent, unresolved);
        if (unresolved != NULL && unresolved->allocation_failed) return;
        dual_node_multipole_scan_neighbors_partial(
            context, pivot, neighbor_tree, neighbor->right,
            same_tree, scratch, statistics,
            active_upper, max_unresolved_extent, unresolved);
        return;
    }
    if (pivot->can_split) {
        dual_node_multipole_mark_unresolved(
            context, upper_extent, max_unresolved_extent);
        if (unresolved != NULL)
            dual_node_neighbor_frontier_append(unresolved, neighbor_index);
        return;
    }

    for (INTEGER i = neighbor->first; i <= neighbor->last; i++) {
        if (same_tree && i >= pivot->first && i <= pivot->last) continue;
        dual_node_multipole_scan_body_partial(
            context, pivot, scratch, statistics,
            &neighbor_tree->packed_points[i], active_upper,
            max_unresolved_extent);
    }
}

/* Completed bins are expressed in their parent's tangent basis. Transport
 * them before mixing them with moments evaluated at a descendant pivot. */
static void dual_node_transport_moments(
        const dual_node_search_context *context,
        dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics,
        const cballs_storage_real *from, const cballs_storage_real *to)
{
    const double started = context->profile ? dual_node_timer_now() : 0.0;
#if NDIM == 3
    real n0[3], a0[3], b0[3], n1[3], a1[3], b1[3], transported[3];
    if (!cballs_angular_basis(from, n0, a0, b0)
        || !cballs_angular_basis(to, n1, a1, b1)) goto transport_done;
    real dot, along, c, s;
    DOTVP(dot, n0, n1);
    if (!(1.0 + dot > 32.0*DBL_EPSILON)) goto transport_done;
    DOTVP(along, a0, n1);
    for (int k = 0; k < 3; k++)
        transported[k] = a0[k] - along*(n0[k]+n1[k])/(1.0+dot);
    DOTVP(c, transported, a1);
    DOTVP(s, transported, b1);
    const real norm = hypot(c, s);
    if (!(norm > 0.0)) goto transport_done;
    c /= norm; s /= norm;
    real cm = 1.0, sm = 0.0;
    for (int m = 0; m < MAX(scratch->orders, scratch->window_orders); m++) {
        for (size_t b = 1; b < scratch->stride; b++) {
            const size_t k = (size_t)m*scratch->stride+b;
            if (m < scratch->orders) {
                const real x = scratch->field_cos[k], y = scratch->field_sin[k];
                const real cc = scratch->distinct_coscos[k], ss = scratch->distinct_sinsin[k];
                const real cs = scratch->distinct_sincos[k];
                scratch->field_cos[k] = cm*x-sm*y;
                scratch->field_sin[k] = sm*x+cm*y;
                scratch->distinct_coscos[k] = cm*cm*cc+sm*sm*ss-2*cm*sm*cs;
                scratch->distinct_sinsin[k] = sm*sm*cc+cm*cm*ss+2*cm*sm*cs;
                scratch->distinct_sincos[k] = cm*sm*(cc-ss)+(cm*cm-sm*sm)*cs;
            }
            if (m < scratch->window_orders) {
                const real x = scratch->window_cos[k], y = scratch->window_sin[k];
                scratch->window_cos[k] = cm*x-sm*y;
                scratch->window_sin[k] = sm*x+cm*y;
            }
        }
        const real next = cm*c-sm*s;
        sm = sm*c+cm*s; cm = next;
    }
#else
    (void)scratch; (void)from; (void)to;
#endif
transport_done:
    if (context->profile)
        statistics->pivot_transport_seconds += dual_node_timer_now() - started;
}

static int dual_node_balltree_depth(
        const fcfc_balltreeptr tree, INTEGER node_index)
{
    const fcfc_ballnode *ball_node = &tree->nodes[node_index];

    if (dual_node_node_is_leaf(ball_node)) return 1;
    return 1 + MAX(dual_node_balltree_depth(tree, ball_node->left),
                   dual_node_balltree_depth(tree, ball_node->right));
}

/*
 * Complete radial bins from large to small.  Before splitting a pivot, retain
 * its completed high-radius moments and clear only the unresolved low bins;
 * each child then inherits the completed work, as in dual-node's LogMultipole
 * traversal.
 */
static int dual_node_multipole_process_pivots_partial(
        const dual_node_search_context *context,
        const fcfc_balltreeptr pivot_tree, INTEGER pivot_node_index,
        const fcfc_balltreeptr neighbor_tree, bool same_tree,
        dual_node_multipole_scratch *scratch,
        dual_node_multipole_scratch *statistics,
        real *scratch_levels, size_t scratch_values_per_level,
        int depth, const INTEGER *parent_nodes, INTEGER parent_count,
        dual_node_neighbor_frontier *frontier_levels,
        real active_upper, int radial_bin_limit,
        dual_node_triple_histogram *hist)
{
    const fcfc_ballnode *pivot_node = &pivot_tree->nodes[pivot_node_index];
    dual_node_neighbor_frontier *frontier = &frontier_levels[depth];
    dual_node_multipole_pivot pivot;
    real max_unresolved_extent = 0.0;
    int unresolved_bin;

    pivot.position = pivot_node->center;
    pivot.source = NULL;
    pivot.radius = (real)pivot_node->radius;
    pivot.field_sum = dual_node_node_field_sum(context, pivot_node);
    pivot.normalization_sum =
        dual_node_node_normalization_sum(context, pivot_node);
    pivot.first = pivot_node->first;
    pivot.last = pivot_node->last;
    pivot.can_split = !dual_node_node_is_leaf(pivot_node);
    dual_node_multipole_initialize_pivot_geometry(&pivot);

    frontier->count = 0;
    frontier->allocation_failed = FALSE;
    for (INTEGER i = 0; i < parent_count; i++) {
        dual_node_multipole_scan_neighbors_partial(
            context, &pivot, neighbor_tree, parent_nodes[i], same_tree,
            scratch, statistics, active_upper, &max_unresolved_extent,
            frontier);
        if (frontier->allocation_failed) return FAILURE;
    }
    unresolved_bin = dual_node_unresolved_radial_bin(
        context, max_unresolved_extent, radial_bin_limit);

    if (unresolved_bin + 1 < radial_bin_limit) {
        dual_node_multipole_finish_pivot_range_profiled(
            context, hist, &pivot, scratch, statistics, unresolved_bin + 1,
            radial_bin_limit, context->cmd->sizeHistN);
        statistics->pivot_finishes++;
    }
    if (unresolved_bin == 0) {
#ifdef DUAL_NODE_PIVOT_PROGRESS
        /* Outer rings may finish at several ancestors. Count the represented
         * active pivots only once, when no radial work remains in this node. */
        dual_node_record_completed_pivots(
            context, statistics, dual_node_node_count(pivot_node));
#endif
        return SUCCESS;
    }

    statistics->pivot_restarts++;
    dual_node_multipole_clear_radial_through(
        context, scratch, statistics, unresolved_bin);
    active_upper = dual_node_radial_upper_edge(context, unresolved_bin);
    radial_bin_limit = unresolved_bin + 1;

    if (!dual_node_node_is_leaf(pivot_node)) {
        const INTEGER children[2] = {pivot_node->left, pivot_node->right};

        for (int child_index = 0; child_index < 2; child_index++) {
            dual_node_multipole_scratch child_scratch;
            real *child_base = scratch_levels
                + (size_t)(depth + 1) * scratch_values_per_level;

            dual_node_initialize_multipole_scratch(
                &child_scratch, child_base, scratch->stride,
                scratch->orders, scratch->window_orders, scratch_values_per_level);
            dual_node_multipole_copy_values(&child_scratch, scratch);
            dual_node_transport_moments(
                context, &child_scratch, statistics, pivot.position,
                pivot_tree->nodes[children[child_index]].center);
            if (dual_node_multipole_process_pivots_partial(
                context, pivot_tree, children[child_index],
                neighbor_tree, same_tree, &child_scratch, statistics,
                scratch_levels, scratch_values_per_level,
                depth + 1, frontier->nodes, frontier->count,
                frontier_levels, active_upper, radial_bin_limit, hist)
                    == FAILURE)
                return FAILURE;
        }
    } else {
        for (INTEGER i = pivot_node->first; i <= pivot_node->last; i++) {
            dual_node_multipole_scratch body_scratch;
            dual_node_multipole_pivot body_pivot;
            const fcfc_ballpoint *pivot_point =
                &pivot_tree->packed_points[i];
            real body_unresolved_extent = 0.0;
            real *body_base = scratch_levels
                + (size_t)(depth + 1) * scratch_values_per_level;

            dual_node_initialize_multipole_scratch(
                &body_scratch, body_base, scratch->stride,
                scratch->orders, scratch->window_orders, scratch_values_per_level);
            dual_node_multipole_copy_values(&body_scratch, scratch);
            dual_node_transport_moments(
                context, &body_scratch, statistics,
                pivot.position, pivot_point->pos);
            body_pivot.position = pivot_point->pos;
            body_pivot.source = pivot_point->source;
            body_pivot.radius = 0.0;
            body_pivot.field_sum =
                dual_node_point_weighted_field(context, pivot_point);
            body_pivot.normalization_sum =
                dual_node_point_normalization(context, pivot_point);
            body_pivot.first = i;
            body_pivot.last = i;
            body_pivot.can_split = FALSE;
            dual_node_multipole_initialize_pivot_geometry(&body_pivot);
            for (INTEGER candidate = 0;
                 candidate < frontier->count; candidate++)
                dual_node_multipole_scan_neighbors_partial(
                    context, &body_pivot, neighbor_tree,
                    frontier->nodes[candidate], same_tree,
                    &body_scratch, statistics, active_upper,
                    &body_unresolved_extent, NULL);
            dual_node_multipole_finish_pivot_range_profiled(
                context, hist, &body_pivot, &body_scratch, statistics,
                1, radial_bin_limit, context->cmd->sizeHistN);
            statistics->pivot_finishes++;
#ifdef DUAL_NODE_PIVOT_PROGRESS
            dual_node_record_completed_pivots(context, statistics, 1);
#endif
        }
    }
    return SUCCESS;
}

static int dual_node_allocate_multipole_scratch(
        struct cmdline_data *cmd, INTEGER task_count,
        size_t stride, int orders, int levels, real **storage,
        size_t *values_per_level, size_t *values_per_task)
{
    const size_t order_arrays = 5;
    size_t order_values;
    size_t values;
    size_t tasks;
    const int window_orders = dual_node_window_orders(cmd);

    *storage = NULL;
    *values_per_level = 0;
    *values_per_task = 0;
    if (task_count <= 0 || orders <= 0 || levels <= 0 || window_orders < 0
        || (size_t)orders > SIZE_MAX / order_arrays
        || (size_t)orders * order_arrays > SIZE_MAX / stride)
        goto invalid_size;
    order_values = order_arrays * (size_t)orders * stride;
    if (stride > (SIZE_MAX - order_values) / 3) goto invalid_size;
    *values_per_level = order_values + 3 * stride;
    if ((size_t)window_orders > (SIZE_MAX - *values_per_level) / 3 / stride)
        goto invalid_size;
    *values_per_level += 3 * (size_t)window_orders * stride;
    if ((size_t)levels > SIZE_MAX / *values_per_level)
        goto invalid_size;
    *values_per_task = (size_t)levels * *values_per_level;
    tasks = (size_t)task_count;
    if (tasks > SIZE_MAX / *values_per_task
        || (values = tasks * *values_per_task) > SIZE_MAX / sizeof(**storage))
        goto invalid_size;
    *storage = calloc(values, sizeof(**storage));
    if (*storage == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: multipole scratch allocation failed", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    return SUCCESS;

invalid_size:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: multipole scratch size overflow", DUAL_NODE_METHOD_NAME);
    return FAILURE;
}

static void dual_node_initialize_multipole_scratch(
        dual_node_multipole_scratch *scratch, real *base,
        size_t stride, int orders, int window_orders, size_t values)
{
    const size_t order_values = (size_t)orders * stride;

    scratch->selected_pairs = NULL;
    scratch->reuse_pairs = scratch->reuse_represented_pairs = 0;
    scratch->reuse_parent_reductions = 0;
    scratch->field_cos = base;
    scratch->field_sin = scratch->field_cos + order_values;
    scratch->distinct_coscos = scratch->field_sin + order_values;
    scratch->distinct_sinsin = scratch->distinct_coscos + order_values;
    scratch->distinct_sincos = scratch->distinct_sinsin + order_values;
    scratch->normalization = scratch->distinct_sincos + order_values;
    scratch->normalization_distinct = scratch->normalization + stride;
    scratch->neighbor_count = scratch->normalization_distinct + stride;
    scratch->window_cos = window_orders ? scratch->neighbor_count + stride : NULL;
    scratch->window_sin = window_orders
        ? scratch->window_cos + (size_t)window_orders * stride : NULL;
    scratch->window_distinct = window_orders
        ? scratch->window_sin + (size_t)window_orders * stride : NULL;
    scratch->stride = stride;
    scratch->values = values;
    scratch->orders = orders;
    scratch->window_orders = window_orders;
    scratch->body_visits = 0;
    scratch->accepted_nodes = 0;
    scratch->pair_tests = 0;
    scratch->pivot_restarts = 0;
    scratch->pivot_finishes = 0;
#ifdef DUAL_NODE_PIVOT_PROGRESS
    scratch->progress_pending = 0;
#endif
    scratch->pivot_transport_seconds = 0.0;
    scratch->scratch_clear_seconds = 0.0;
    scratch->multipole_product_seconds = 0.0;
}
