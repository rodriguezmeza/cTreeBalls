/* Opt-in scalar hierarchy. Accepted neighbor moments are frozen in their
 * original coordinate phase: acceptance bounds that phase against EVERY
 * descendant pivot's own native basis. No successive transport errors accrue.
 * A radial pair is published once, when both bins have no unresolved source.
 * Included after dual_node_multipole.h. */
static int dual_node_reuse_control(struct cmdline_data *cmd, const char *name,
                                   real fallback, real maximum, real *value)
{
    const char *text = getenv(name);
    char *end = NULL;
    double parsed;
    *value = fallback;
    if (text == NULL) return SUCCESS;
    parsed = strtod(text, &end);
    if (end == text || *end != '\0' || !isfinite(parsed)
        || parsed < 0.0 || parsed > maximum) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s must be finite in [0,%g]", name, (double)maximum);
        return FAILURE;
    }
    *value = (real)parsed;
    return SUCCESS;
}

static int dual_node_reuse_prepare(dual_node_search_context *context,
                                   struct cmdline_data *cmd,
                                   struct global_data *gd, bool run_3pcf)
{
    const bool requested = run_3pcf && scanopt(cmd->options, "scalar-pivot-reuse");
    int status = SUCCESS;
    real phase = 0.1, bin = 0.0;
    char controls[160];
    if (requested) {
        status = dual_node_reuse_control(cmd, "CBALLS_SCALAR_PIVOT_TOL",
                                         0.1, 3.0, &phase);
        if (status == SUCCESS)
            status = dual_node_reuse_control(cmd, "CBALLS_SCALAR_BIN_THETA",
                                             0.0, 1.0, &bin);
    }
    if (dual_node_distributed_consensus(cmd, status,
            "scalar pivot-reuse controls") == FAILURE) return FAILURE;
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    snprintf(controls, sizeof(controls), "%d:%a:%a", requested,
             (double)phase, (double)bin);
    if (cballs_mpi_agree_text(cmd, controls,
            "scalar pivot-reuse controls") == FAILURE) return FAILURE;
#else
    (void)controls;
#endif
    context->reuse_enabled = requested && phase > 0.0 && cmd->theta > 0.0
        && !cmd->usePeriodic && !cballs_opt_smooth_pivot(cmd)
        && !cballs_opt_no_one_ball(cmd) && !cballs_opt_no_two_balls(cmd);
    context->reuse_bin_theta = bin;
    context->reuse_max_relative_width = expm1(log(10.0)*gd->deltaR);
    const int order = MAX(cmd->mChebyshev, dual_node_window_orders(cmd) - 1);
    /* Each leg may change phase by phase/(2*order); the corresponding pair
     * product's highest harmonic then changes phase by at most phase. */
    context->reuse_phase_ratio = sin(MIN(0.5*PI, phase/(2.0*MAX(order, 1))));
    gd->scalarReuseEnabled = context->reuse_enabled;
    gd->scalarReusePhaseBudget = requested ? phase : 0.0;
    gd->scalarReuseBinTheta = requested ? bin : 0.0;
    return SUCCESS;
}

typedef struct {
    dual_node_multipole_scratch moments;
    dual_node_neighbor_frontier pending;
    unsigned char *done, *selected, *unresolved;
} dual_node_reuse_level;

typedef struct {
    const dual_node_search_context *context;
    fcfc_balltreeptr pivots, neighbors;
    bool same_tree;
    dual_node_reuse_level *levels;
    real *edges;
    size_t stride, plane;
    dual_node_multipole_scratch *statistics;
    dual_node_triple_histogram *hist;
} dual_node_reuse_work;

/* Bounds the CHANGE in the native angular coordinates, not merely geodesic
 * transport. For u = reference - dot(reference,n)*n, |du| <= 2|dn|,
 * |d(normalize(u))| <= 2|du|/|u|, and |d(a x n)| <= |da|+|dn|.
 * Reject caps crossing angular_basis's least-aligned Cartesian-axis switch.
 * The normalized-vector inequality gives |dn| <= 2*radius/|position|. */
static real dual_node_reuse_basis_error(const dual_node_multipole_pivot *pivot)
{
#if NDIM == 3
    if (!pivot->angular_basis_valid || !(pivot->position_norm > pivot->radius))
        return HUGE_VAL;
    const real dn = 2.0*pivot->radius/pivot->position_norm;
    int axis=0;
    for (int k=1;k<3;k++)
        if (fabs(pivot->normal[k]) < fabs(pivot->normal[axis])) axis=k;
    if (pivot->radius > 0.0)
        for (int k=0;k<3;k++)
            if (k!=axis && fabs(pivot->normal[axis])+dn >= fabs(pivot->normal[k])-dn)
                return HUGE_VAL;
    const real dot = pivot->normal[axis];
    const real h = sqrt(MAX(0.0, 1.0-dot*dot));
    if (!(h > 0.0)) return HUGE_VAL;
    const real da = MIN(2.0, 4.0*dn/h);
    const real db = MIN(2.0, da+dn);
    return hypot(da, db);
#else
    (void)pivot;
    return 0.0;
#endif
}

static void dual_node_reuse_mark(dual_node_reuse_work *work,
                                 dual_node_reuse_level *level,
                                 real lower, real upper)
{
    /* Include both sides of numerical bin boundaries. Empty lower extension
     * of the legacy rmin=0 log grid is harmlessly treated as unresolved. */
    size_t lo=1,hi=work->stride;
    while(lo<hi) {
        const size_t mid=lo+(hi-lo)/2;
        if (work->edges[mid]<lower) lo=mid+1; else hi=mid;
    }
    const size_t first=lo;
    lo=1;hi=work->stride;
    while(lo<hi) {
        const size_t mid=lo+(hi-lo)/2;
        if (work->edges[mid-1]<=upper) lo=mid+1; else hi=mid;
    }
    if (first<lo) memset(level->unresolved+first,1,lo-first);
}

/* A negative reference encodes one packed body; nonnegative is a tree node. */
static int dual_node_reuse_scan(dual_node_reuse_work *work,
        dual_node_reuse_level *level, const dual_node_multipole_pivot *pivot,
        real basis_error, INTEGER ref)
{
    const dual_node_search_context *context = work->context;
    if (ref>=0 && work->neighbors->nodes[ref].first==work->neighbors->nodes[ref].last)
        ref=-work->neighbors->nodes[ref].first-1;
    const bool body = ref < 0;
    const INTEGER point_index = body ? -ref-1 : 0;
    const fcfc_ballnode *node = body ? NULL : &work->neighbors->nodes[ref];
    const fcfc_ballpoint *point = body ? &work->neighbors->packed_points[point_index] : NULL;
    const cballs_storage_real *position = body ? point->pos : node->center;
    const real radius = body ? 0.0 : (real)node->radius;
    const real size = pivot->radius + radius;
    compute_vector dr;
    const real distance = dual_node_position_distance(context, pivot->position, position, dr);
    const real pad = 64.0*DBL_EPSILON*(distance+size+context->cmd->rangeN);
    const real lower = MAX(0.0, distance-size-pad), upper = distance+size+pad;
    const bool overlap = work->same_tree && (body
        ? point_index >= pivot->first && point_index <= pivot->last
        : dual_node_ranges_overlap(pivot->first, pivot->last, node->first, node->last));
    const bool single_pivot = pivot->source != NULL || pivot->first == pivot->last;
    int bin;
    real c, s;
    work->statistics->pair_tests++;
    if (upper <= context->cmd->rminHist || lower >= context->cmd->rangeN) return SUCCESS;
    if (body && single_pivot && (overlap || (pivot->source != NULL && pivot->source == point->source)))
        return SUCCESS;
    /* No source can be accepted across a pivot's native-basis switch.
     * Split that pivot first. Refining its neighbors here would permanently
     * replace reusable cells by enormous body lists for all descendants. */
    if (!single_pivot && !isfinite(basis_error)) {
        dual_node_reuse_mark(work,level,lower,upper);
        return dual_node_neighbor_frontier_append(&level->pending,ref) ? SUCCESS : FAILURE;
    }
    /* A source wholly inside the pivot cap remains one pending tree node.
     * Expanding it into bodies here would destroy hierarchical frontier reuse. */
    if (overlap && !single_pivot && !body
        && node->first>=pivot->first && node->last<=pivot->last) {
        dual_node_reuse_mark(work,level,lower,upper);
        return dual_node_neighbor_frontier_append(&level->pending,ref) ? SUCCESS : FAILURE;
    }
    if (body && pivot->radius == 0.0 && single_pivot) {
        /* Exact terminal classification, including same-position/axial skips. */
        if (!(distance > 0.0) || (bin = dual_node_bin_index(context, distance)) < 1
            || !dual_node_multipole_angular_phase(pivot, dr, &c, &s)) return SUCCESS;
        const real field = dual_node_point_weighted_field(context, point);
        const real norm = dual_node_point_normalization(context, point);
        dual_node_multipole_add(&level->moments, bin, field, field*field,
                                norm, norm*norm, 1, c, s);
        work->statistics->body_visits++;
        return SUCCESS;
    }
    bool accept = !overlap && isfinite(basis_error) && distance>size
        && lower > context->cmd->rminHist && upper < context->cmd->rangeN;
    /* A cheap upper bound on any containing bin's width avoids logarithms
     * for large rejected cells. rmin=0 retains its legacy extended first bin. */
    if (accept && (!context->cmd->useLogHist || context->cmd->rminHist>0.0)) {
        const real width_bound=context->cmd->useLogHist
            ? distance*context->reuse_max_relative_width : context->gd->deltaR;
        if (size+pad > MAX(0.5,context->reuse_bin_theta)*width_bound) accept=FALSE;
    }
    bin=-1;
    if (accept) {
        bin=dual_node_bin_index(context,distance);
        accept=bin>0;
    }
    if (accept) {
        const real width = work->edges[bin]-work->edges[bin-1];
        const bool slack=context->reuse_bin_theta>0.0
            && size+pad<=context->reuse_bin_theta*width;
        if (slack && (!context->cmd->useLogHist || context->cmd->rminHist>0.0)) {
            /* Strict outer bounds already establish valid bins. */
        } else {
            const int low_bin=dual_node_bin_index(context,lower);
            const int high_bin=dual_node_bin_index(context,upper);
            accept=(low_bin==bin && high_bin==bin)
                || (slack && low_bin>0 && high_bin>0);
        }
    }
    if (accept) {
        real transverse;
#if NDIM == 3
        real x = 0.0, y = 0.0;
        for (int k=0;k<3;k++) {
            x -= dr[k]*pivot->first_axis[k]; y -= dr[k]*pivot->second_axis[k];
        }
        transverse = hypot(x,y);
#else
        transverse = distance;
#endif
        /* Projection in the descendant's orthonormal frame is contractive.
         * Perturb its separation by <=size and its basis by <=basis_error. */
        const real error = size + distance*basis_error + pad;
        accept = transverse > error && error <= context->reuse_phase_ratio*transverse
            && dual_node_multipole_angular_phase(pivot, dr, &c, &s);
    }
    if (accept) {
        const real field = body ? dual_node_point_weighted_field(context, point)
                                : dual_node_node_field_sum(context, node);
        const real norm = body ? dual_node_point_normalization(context, point)
                               : dual_node_node_normalization_sum(context, node);
        dual_node_multipole_add(&level->moments, bin, field,
            body ? field*field : dual_node_node_field_sq_sum(context,node), norm,
            body ? norm*norm : dual_node_node_normalization_sq_sum(context,node),
            body ? 1 : dual_node_node_count(node), c,s);
        work->statistics->accepted_nodes++;
        return SUCCESS;
    }
    /* Refine the larger neighbor, or resolve all neighbors at an exact pivot.
     * Self overlap always descends: an aggregate never contains its own leg. */
    if (!body && (single_pivot || overlap || radius > pivot->radius)) {
        if (!dual_node_node_is_leaf(node)) {
            if (dual_node_reuse_scan(work,level,pivot,basis_error,node->left) == FAILURE)
                return FAILURE;
            return dual_node_reuse_scan(work,level,pivot,basis_error,node->right);
        }
        for (INTEGER i=node->first;i<=node->last;i++)
            if (dual_node_reuse_scan(work,level,pivot,basis_error,-i-1) == FAILURE)
                return FAILURE;
        return SUCCESS;
    }
    dual_node_reuse_mark(work,level,lower,upper);
    return dual_node_neighbor_frontier_append(&level->pending,ref) ? SUCCESS : FAILURE;
}

static int dual_node_reuse_visit(dual_node_reuse_work *work, INTEGER ref,
                                 int depth, const INTEGER *pending, INTEGER count)
{
    dual_node_reuse_level *level = &work->levels[depth];
    const bool body = ref < 0;
    const INTEGER index = body ? -ref-1 : 0;
    const fcfc_ballnode *node = body ? NULL : &work->pivots->nodes[ref];
    const fcfc_ballpoint *point = body ? &work->pivots->packed_points[index] : NULL;
    dual_node_multipole_pivot pivot;
    const INTEGER represented = body ? 1 : dual_node_node_count(node);
    /* Canonicalize singleton nodes, including cache entries with source=NULL. */
    if (!body && node->first == node->last)
        return dual_node_reuse_visit(work,-node->first-1,depth,pending,count);
    pivot.position = body ? point->pos : node->center;
    pivot.source = body ? point->source : NULL;
    pivot.radius = body ? 0.0 : (real)node->radius;
    pivot.field_sum = body ? dual_node_point_weighted_field(work->context,point)
                          : dual_node_node_field_sum(work->context,node);
    pivot.normalization_sum = body ? dual_node_point_normalization(work->context,point)
                                  : dual_node_node_normalization_sum(work->context,node);
    pivot.first = body ? index : node->first; pivot.last = body ? index : node->last;
    pivot.can_split = !body;
    dual_node_multipole_initialize_pivot_geometry(&pivot);
#if NDIM == 3 && defined(DUAL_NODE_DISABLE_CACHED_GEOMETRY)
    pivot.position_norm=hypot(hypot((real)pivot.position[0],(real)pivot.position[1]),(real)pivot.position[2]);
    pivot.angular_basis_valid=cballs_angular_basis(pivot.position,pivot.normal,pivot.first_axis,pivot.second_axis);
#endif
    if (depth > 0) {
        dual_node_reuse_level *parent=&work->levels[depth-1];
        dual_node_multipole_copy_values(&level->moments,&parent->moments);
        memcpy(level->done,parent->done,work->plane);
    }
    memset(level->unresolved,0,work->stride);
    memset(level->selected,0,work->plane);
    level->pending.count=0; level->pending.allocation_failed=FALSE;
    const real basis_error = dual_node_reuse_basis_error(&pivot);
    for (INTEGER i=0;i<count;i++)
        if (dual_node_reuse_scan(work,level,&pivot,basis_error,pending[i]) == FAILURE)
            return FAILURE;
    INTEGER combinations=0;
    bool complete=TRUE;
    for (size_t i=1;i<work->stride;i++) for (size_t j=i;j<work->stride;j++) {
        const size_t pair=i*work->stride+j;
        if (level->done[pair]) continue;
        if (level->unresolved[i] || level->unresolved[j]) { complete=FALSE; continue; }
        level->done[pair]=level->selected[pair]=1;
        if (level->moments.neighbor_count[i] > 0 && level->moments.neighbor_count[j] > 0
            && (i!=j || level->moments.neighbor_count[i]>=2)) combinations++;
    }
    if (combinations) {
        level->moments.selected_pairs=level->selected;
        dual_node_multipole_finish_pivot_range_profiled(work->context,work->hist,
            &pivot,&level->moments,work->statistics,1,(int)work->stride,(int)work->stride-1);
        work->statistics->reuse_pairs+=combinations;
        work->statistics->reuse_represented_pairs+=combinations*represented;
        work->statistics->reuse_parent_reductions+=represented>1;
        work->statistics->pivot_finishes++;
    }
    if (complete) {
#ifdef DUAL_NODE_PIVOT_PROGRESS
        dual_node_record_completed_pivots(work->context,work->statistics,represented);
#endif
        return SUCCESS;
    }
    /* Body scans always resolve every candidate. */
    if (body) return FAILURE;
    work->statistics->pivot_restarts++;
    if (!dual_node_node_is_leaf(node)) {
        if (dual_node_reuse_visit(work,node->left,depth+1,level->pending.nodes,level->pending.count) == FAILURE)
            return FAILURE;
        return dual_node_reuse_visit(work,node->right,depth+1,level->pending.nodes,level->pending.count);
    }
    for (INTEGER i=node->first;i<=node->last;i++)
        if (dual_node_reuse_visit(work,-i-1,depth+1,level->pending.nodes,level->pending.count) == FAILURE)
            return FAILURE;
    return SUCCESS;
}

static int dual_node_reuse_task(const dual_node_search_context *context,
        fcfc_balltreeptr pivots, INTEGER root, fcfc_balltreeptr neighbors,
        bool same_tree, dual_node_multipole_scratch *statistics,
        dual_node_triple_histogram *hist)
{
    const size_t stride=(size_t)context->cmd->sizeHistN+1;
    const int levels=dual_node_balltree_depth(pivots,root)+1;
    real *storage=NULL;
    unsigned char *masks=NULL;
    size_t per_level=0, total=0;
    int status=FAILURE;
    dual_node_reuse_work work={context,pivots,neighbors,same_tree,NULL,NULL,
                              stride,0,statistics,hist};
    if (stride > SIZE_MAX/stride) return FAILURE;
    work.plane=stride*stride;
    if (work.plane > (SIZE_MAX-stride)/2
        || (size_t)levels > SIZE_MAX/(2*work.plane+stride)) return FAILURE;
    /* Bound live DFS storage per worker; never allocate depth scratch for all
     * frontier tasks. Allocation failure is reported collectively after join. */
    const size_t limit=(size_t)64*1024*1024;
    if (statistics->values > limit/sizeof(real)/(size_t)levels
        || (2*work.plane+stride) > limit/(size_t)levels
        || (size_t)levels*(statistics->values*sizeof(real)+2*work.plane+stride) > limit)
        return FAILURE;
    struct cmdline_data allocation_cmd=*context->cmd;
    if (dual_node_allocate_multipole_scratch(&allocation_cmd,1,stride,
            context->cmd->mChebyshev+1,levels,&storage,&per_level,&total) == FAILURE)
        return FAILURE;
    work.levels=calloc((size_t)levels,sizeof(*work.levels));
    work.edges=malloc(stride*sizeof(*work.edges));
    masks=calloc((size_t)levels,2*work.plane+stride);
    if (!work.levels || !work.edges || !masks) goto done;
    for (size_t b=0;b<stride;b++) work.edges[b]=dual_node_radial_upper_edge(context,(int)b);
    for (int d=0;d<levels;d++) {
        dual_node_reuse_level *level=&work.levels[d];
        dual_node_initialize_multipole_scratch(&level->moments,storage+(size_t)d*per_level,
            stride,context->cmd->mChebyshev+1,dual_node_window_orders(context->cmd),per_level);
        level->done=masks+(size_t)d*(2*work.plane+stride);
        level->selected=level->done+work.plane;
        level->unresolved=level->selected+work.plane;
    }
    const INTEGER neighbor_root=0;
    status=dual_node_reuse_visit(&work,root,0,&neighbor_root,1);
done:
    if (work.levels) for (int d=0;d<levels;d++) free(work.levels[d].pending.nodes);
    free(masks); free(work.edges); free(work.levels); free(storage);
    return status;
}
