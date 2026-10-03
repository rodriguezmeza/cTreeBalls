/* Persistent 3PCF traversal; uses lya_forest_cells.h. */
/* status: -1 outside the physical domain, 0 unresolved, 1 usable leg bins. */
typedef struct {
    size_t node;
    REAL dr[3], radius, lo[3], hi[3];
    int radial, polar, status, approximate;
    REAL rmin, rmax, tmin, tmax;
    size_t valid;
} lya_cell_leg;

/* Bounded direct-mapped cache: collisions only cause recomputation. Keys are
 * immutable node IDs within this search; no cached pointer survives a search. */
typedef struct {
    size_t pivot_plus_one, neighbor;
    lya_cell_leg geometry;
} lya_cell_cache_entry;
#define LYA_CELL_PAIR_CACHE_SIZE 8192

/* Leaf tiles need one cache line including their tag, rather than a full
 * interval record. Reconstruct interval fields only for hierarchy traversal. */
typedef struct {
    REAL dr[3], radius, theta;
    int radial, polar, status;
} lya_cell_pixel_geometry;

typedef struct {
    lya_cell_cache_entry *pair_cache;
    lya_cell_leg *legs, *point_cache;
    lya_cell_pixel_geometry *pixel_cache;
    size_t *point_tags, *pixel_tags;
    size_t task_begin;
    size_t cached_pivot;
    size_t count, capacity, radial_bytes;
    REAL radial_scale, angular_scale, polar_scale;
    unsigned long long pivot_aggregates;
    unsigned long long geometry_evaluations, cache_hits, pruned_nodes, leaf_evaluations, pair_cache_hits;
} lya_cell_workspace;

/* Intersect box bounds with a forest-aligned capsule viewed from the pivot
 * bounding sphere. The capsule is [first,last] + ball(tube); the pivot sphere
 * adds its radius. This remains an enclosure for non-collinear forests. */
static void lya_cell_capsule_bounds(const lya_cells *tree,const lya_cell *p,
                                     const lya_cell *q,REAL *rlo,REAL *rhi,
                                     REAL *cone_cos,REAL *cone_sin)
{
    *cone_cos=-1;*cone_sin=0;
    if(!isfinite(q->tube)) return;
    const REAL *a=Pos(tree->points[q->begin]),*b=Pos(tree->points[q->end-1]);
    REAL scale=0,d0[3],d1[3],axis[3],a2=0,b2=0,c2=0;
    for(int k=0;k<3;k++) {
        scale=MAX(scale,MAX(MAX(fabs(p->lo[k]),fabs(p->hi[k])),MAX(fabs(a[k]),fabs(b[k]))));
        d0[k]=a[k]-p->center[k];d1[k]=b[k]-p->center[k];
        axis[k]=q->center[k]-p->center[k];
        a2+=d0[k]*d0[k];b2+=d1[k]*d1[k];c2+=axis[k]*axis[k];
    }
    if(!(scale<sqrt(DBL_MAX)/32)) return;
    REAL pad=4096*DBL_EPSILON*scale+16*sqrt(DBL_MIN);
    REAL distance=lya_cell_segment_distance(p->center,a,b);
    REAL error=p->radius+q->tube+pad;
    REAL low=MAX(0.,distance-error),high=MAX(sqrt(a2),sqrt(b2))+error;
    if(!isfinite(high)) return;
    *rlo=MAX(*rlo,low);*rhi=MIN(*rhi,high);
    /* A positive endpoint cosine places the whole segment in the same
     * hemisphere. A ball error expands its angular cone by asin(error/d). */
    if(!(distance>error) || !(a2>DBL_MIN) || !(b2>DBL_MIN) || !(c2>DBL_MIN)) return;
    REAL c0=0,c1=0;
    for(int k=0;k<3;k++) {c0+=axis[k]*d0[k];c1+=axis[k]*d1[k];}
    REAL cosine=MIN(c0/(sqrt(c2)*sqrt(a2)),c1/(sqrt(c2)*sqrt(b2)))-4096*DBL_EPSILON;
    if(!(cosine>0)) return;
    cosine=MIN(1.,cosine);
    REAL sine=sqrt(MAX(0.,1-cosine*cosine)),extra=MIN(1.,error/distance);
    REAL ec=sqrt(MAX(0.,1-extra*extra));
    REAL cc=cosine*ec-sine*extra-4096*DBL_EPSILON;
    if(!(cc>0)) return;
    *cone_cos=cc;*cone_sin=MIN(1.,sqrt(MAX(0.,1-cc*cc))+4096*DBL_EPSILON);
}

static lya_cell_leg lya_cell_geometry(const lya_cells *tree,struct cmdline_data *cmd,
                                       size_t ip,size_t iq)
{
    const lya_cell *p=&tree->nodes[ip],*q=&tree->nodes[iq];
    lya_cell_leg g={0};g.node=iq;
    REAL d2=0;
    for(int k=0;k<3;k++) {g.dr[k]=q->center[k]-p->center[k];d2+=g.dr[k]*g.dr[k];}
    g.radius=sqrt(d2);
    if (p->left==SIZE_MAX && q->left==SIZE_MAX) {
        if (!(g.radius>0) || !(g.radius<cmd->lya3RMax)) {g.status=-1;return g;}
        REAL dot;DOTVP(dot,g.dr,p->sight);
        g.radial=lya_bin_positive(g.radius,cmd->lya3RMax,cmd->lya3RBins);
        g.tmin=g.tmax=acos(lya_clamp(dot/g.radius,-1.,1.));
        g.polar=lya_bin_theta(g.tmin,cmd->lya3ThetaBins);
        for(int k=0;k<3;k++) g.lo[k]=g.hi[k]=g.dr[k]/g.radius;
        g.rmin=g.rmax=g.radius;g.valid=1;g.status=1;return g;
    }
#ifdef __FAST_MATH__
    return g;
#endif
    if (sizeof(REAL)!=sizeof(double)) return g;
    REAL dl[3],dh[3],rlo2=0,rhi2=0;
    for(int k=0;k<3;k++) {
        /* Subtraction is rounded in the reference too. The absolute endpoint
         * padding accounts for cancellation and later norm evaluation. */
        REAL scale=MAX(MAX(fabs(p->lo[k]),fabs(p->hi[k])),MAX(fabs(q->lo[k]),fabs(q->hi[k])));
        REAL pad=512*DBL_EPSILON*scale+DBL_MIN;
        dl[k]=q->lo[k]-p->hi[k]-pad;dh[k]=q->hi[k]-p->lo[k]+pad;
        REAL near=dl[k]>0?dl[k]:(dh[k]<0?-dh[k]:0);
        REAL far=MAX(fabs(dl[k]),fabs(dh[k]));rlo2+=near*near;rhi2+=far*far;
    }
    REAL rlo=sqrt(rlo2)*(1-512*DBL_EPSILON),rhi=sqrt(rhi2)*(1+512*DBL_EPSILON);
    if (!isfinite(rlo) || !isfinite(rhi)) return g;
    REAL cone_cos,cone_sin;
    lya_cell_capsule_bounds(tree,p,q,&rlo,&rhi,&cone_cos,&cone_sin);
    if(!(rlo<=rhi)) return g;
    if (rlo>=cmd->lya3RMax) {g.status=-1;return g;}
    if (rlo<16*sqrt(DBL_MIN) || rhi>sqrt(DBL_MAX)/16 || !(rhi<cmd->lya3RMax)) return g;
    g.radial=lya_cell_bin(rlo,rhi,g.radius,cmd->lya3RMax,cmd->lya3RBins,cmd->lya3RadialSlop,0,&g.approximate);
    if(g.radial<0) return g;
    for(int k=0;k<3;k++) {
        REAL a=dl[k]/rlo,b=dl[k]/rhi,c=dh[k]/rlo,d=dh[k]/rhi;
        g.lo[k]=MAX(-1.,MIN(MIN(a,b),MIN(c,d))-512*DBL_EPSILON);
        g.hi[k]=MIN(1.,MAX(MAX(a,b),MAX(c,d))+512*DBL_EPSILON);
    }
    if(cone_cos>0 && g.radius>0) for(int k=0;k<3;k++) {
        REAL axis=lya_clamp(g.dr[k]/g.radius,-1.,1.);
        REAL transverse=sqrt(MAX(0.,1-axis*axis));
        REAL low=-axis>=cone_cos?-1.:axis*cone_cos-transverse*cone_sin;
        REAL high=axis>=cone_cos?1.:axis*cone_cos+transverse*cone_sin;
        g.lo[k]=MAX(g.lo[k],low-4096*DBL_EPSILON);
        g.hi[k]=MIN(g.hi[k],high+4096*DBL_EPSILON);
    }
    REAL lo,hi,dot;lya_cell_dot_bounds(g.lo,g.hi,p->loslo,p->loshi,&lo,&hi);
    DOTVP(dot,g.dr,p->sight);
    g.polar=lya_cell_bin(acos(hi),acos(lo),acos(lya_clamp(dot/g.radius,-1.,1.)),LYA_PI,cmd->lya3ThetaBins,cmd->lya3PolarSlop,1,&g.approximate);
    if(g.polar<0) return g;
    g.status=1;return g;
}

/* Reject with a padded lower bound before visiting any descendant. Equality
 * is outside the strict r < Rmax domain. Unsupported arithmetic profiles fall
 * through to the reference leaf geometry, just as aggregation does. */
static int lya_cell_outside(const lya_cells *tree,const lya_cell *p,const lya_cell *q,REAL maximum)
{
#ifdef __FAST_MATH__
    return 0;
#else
    if(sizeof(REAL)!=sizeof(double)) return 0;
    REAL sum=0;
    for(int k=0;k<3;k++) {
        REAL scale=MAX(MAX(fabs(p->lo[k]),fabs(p->hi[k])),MAX(fabs(q->lo[k]),fabs(q->hi[k])));
        REAL pad=512*DBL_EPSILON*scale+DBL_MIN;
        REAL near=MAX(0.,MAX(q->lo[k]-p->hi[k],p->lo[k]-q->hi[k])-pad);
        sum+=near*near;
    }
    REAL low=sqrt(sum)*(1-512*DBL_EPSILON);
    if(isfinite(low) && low>=maximum) return 1;
    /* Capsule distance often rejects oblique forests whose boxes overlap. No
     * angular work is necessary for this coarse reachability test. */
    if(isfinite(q->tube)) {
        const REAL *a=Pos(tree->points[q->begin]),*b=Pos(tree->points[q->end-1]);
        REAL scale=0;
        for(int k=0;k<3;k++) scale=MAX(scale,MAX(MAX(fabs(p->lo[k]),fabs(p->hi[k])),MAX(fabs(a[k]),fabs(b[k]))));
        if(scale<sqrt(DBL_MAX)/32) {
            REAL pad=4096*DBL_EPSILON*scale+16*sqrt(DBL_MIN);
            REAL d=lya_cell_segment_distance(p->center,a,b)-p->radius-q->tube-pad;
            if(isfinite(d) && d>=maximum) return 1;
        }
    }
    return 0;
#endif
}

/* Direct-indexed pixel cache for the active pivot task. Unlike the bounded
 * node hash, these slots cannot evict one another while its children split.
 * Singleton tasks reuse the original dense node cache without extra storage. */
static const lya_cell_leg *lya_cell_leaf_cached(const lya_cells *tree,struct cmdline_data *cmd,
                                                 size_t ip,size_t iq,lya_cell_workspace *ws)
{
    lya_cell_leg *g=&ws->point_cache[iq];size_t *tag=&ws->point_tags[iq];
    if(*tag==ip+1) {ws->cache_hits++;return g;}
    ws->geometry_evaluations++;ws->leaf_evaluations++;
    *g=lya_cell_geometry(tree,cmd,ip,iq);
    *tag=ip+1;return g;
}

static const lya_cell_pixel_geometry *lya_cell_pixel_cached(const lya_cells *tree,
        struct cmdline_data *cmd,size_t ip,size_t iq,lya_cell_workspace *ws)
{
    size_t slot=(tree->nodes[ip].begin-ws->task_begin)*tree->count+tree->nodes[iq].begin;
    lya_cell_pixel_geometry *g=&ws->pixel_cache[slot];
    if(ws->pixel_tags[slot]==ip+1) {ws->cache_hits++;return g;}
    ws->geometry_evaluations++;ws->leaf_evaluations++;
    lya_cell_leg leg=lya_cell_geometry(tree,cmd,ip,iq);
    for(int k=0;k<3;k++) g->dr[k]=leg.dr[k];
    g->radius=leg.radius;g->theta=leg.tmin;
    g->radial=leg.radial;g->polar=leg.polar;g->status=leg.status;
    ws->pixel_tags[slot]=ip+1;return g;
}

static lya_cell_leg lya_cell_get_geometry(const lya_cells *,struct cmdline_data *,
                                           size_t,size_t,lya_cell_workspace *);

/* Point-pivot cache: measure reachable leaf geometry once, then reuse bottom-up
 * extrema on the persistent hierarchy. This avoids O(triples) acos/norm work
 * on unresolved nodes and bounds directions more tightly than box division. */
static lya_cell_leg lya_cell_cache_build(const lya_cells *tree,struct cmdline_data *cmd,
                                         size_t ip,size_t iq,lya_cell_workspace *ws)
{
    const lya_cell *q=&tree->nodes[iq];
    lya_cell_leg g={0};g.node=iq;
    if(q->left!=SIZE_MAX && lya_cell_outside(tree,&tree->nodes[ip],q,cmd->lya3RMax)) {
        g.status=-1;ws->pruned_nodes++;
    } else if(q->left==SIZE_MAX) {
        ws->leaf_evaluations++;
        g=lya_cell_geometry(tree,cmd,ip,iq);
    } else {
        lya_cell_leg a=lya_cell_get_geometry(tree,cmd,ip,q->left,ws);
        lya_cell_leg b=lya_cell_get_geometry(tree,cmd,ip,q->right,ws);
        g.valid=a.valid+b.valid;
        if(!g.valid) g.status=-1;
        else if(g.valid==q->end-q->begin) {
            g.rmin=MIN(a.rmin,b.rmin);g.rmax=MAX(a.rmax,b.rmax);
            g.tmin=MIN(a.tmin,b.tmin);g.tmax=MAX(a.tmax,b.tmax);
            REAL d2=0,dot;
            for(int k=0;k<3;k++) {
                g.dr[k]=q->center[k]-tree->nodes[ip].center[k];d2+=g.dr[k]*g.dr[k];
                g.lo[k]=MIN(a.lo[k],b.lo[k]);g.hi[k]=MAX(a.hi[k],b.hi[k]);
            }
            g.radius=sqrt(d2);DOTVP(dot,g.dr,tree->nodes[ip].sight);
            g.radial=lya_cell_bin(g.rmin,g.rmax,g.radius,cmd->lya3RMax,cmd->lya3RBins,cmd->lya3RadialSlop,0,&g.approximate);
            g.polar=lya_cell_bin(g.tmin,g.tmax,acos(lya_clamp(dot/g.radius,-1.,1.)),LYA_PI,cmd->lya3ThetaBins,cmd->lya3PolarSlop,1,&g.approximate);
            if(g.radial>=0 && g.polar>=0 && g.rmin>=16*sqrt(DBL_MIN) && g.rmax<=sqrt(DBL_MAX)/16) g.status=1;
        }
    }
    return g;
}

static lya_cell_leg lya_cell_get_geometry(const lya_cells *tree,struct cmdline_data *cmd,
                                           size_t ip,size_t iq,lya_cell_workspace *ws)
{
    if(tree->nodes[ip].left==SIZE_MAX && tree->nodes[iq].left==SIZE_MAX) {
        if(ws->cached_pivot==ip) return *lya_cell_leaf_cached(tree,cmd,ip,iq,ws);
        const lya_cell_pixel_geometry *pixel=lya_cell_pixel_cached(tree,cmd,ip,iq,ws);
        lya_cell_leg g={0};g.node=iq;g.status=pixel->status;
        if(g.status==1) {
            g.valid=1;g.radius=g.rmin=g.rmax=pixel->radius;
            g.radial=pixel->radial;g.polar=pixel->polar;g.tmin=g.tmax=pixel->theta;
            for(int k=0;k<3;k++) {g.dr[k]=pixel->dr[k];g.lo[k]=g.hi[k]=g.dr[k]/g.radius;}
        }
        return g;
    }
    const int dense=ws->cached_pivot==ip;
    lya_cell_cache_entry *entry=NULL;
    if(dense) {
        if(ws->point_tags[iq]==ip+1) {ws->cache_hits++;return ws->point_cache[iq];}
    } else {
        /* Unsigned wrap is intentional, independent of pointer layout. */
        uint64_t key=(uint64_t)ip*UINT64_C(0x9e3779b97f4a7c15)+(uint64_t)iq;
        key^=key>>30;key*=UINT64_C(0xbf58476d1ce4e5b9);key^=key>>27;
        entry=&ws->pair_cache[key&(LYA_CELL_PAIR_CACHE_SIZE-1)];
        if(entry->pivot_plus_one==ip+1 && entry->neighbor==iq) {
            ws->cache_hits++;ws->pair_cache_hits++;return entry->geometry;
        }
    }
    ws->geometry_evaluations++;
    lya_cell_leg g=tree->nodes[ip].left==SIZE_MAX
        ? lya_cell_cache_build(tree,cmd,ip,iq,ws) : lya_cell_geometry(tree,cmd,ip,iq);
    if(dense) {ws->point_cache[iq]=g;ws->point_tags[iq]=ip+1;}
    else {entry->geometry=g;entry->pivot_plus_one=ip+1;entry->neighbor=iq;}
    return g;
}

static int lya_cell_append(lya_cell_workspace *ws,lya_cell_leg leg,
                            const lya_cells *tree,struct cmdline_data *cmd,ErrorMsg err)
{
    if(ws->count==ws->capacity) {
        size_t cap=ws->capacity?ws->capacity*2:128,bytes,total;
        if(cap<ws->capacity || !cballs_size_mul(cap,sizeof(leg),&bytes)
            || !cballs_size_mul(bytes,2,&total)
            || !cballs_size_add(total,ws->radial_bytes,&total)
            || !cballs_size_mul(total,(size_t)MAX(1,cmd->numthreads),&total)
            || !cballs_size_add(total,tree->plan,&total)) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha node frontier size overflow");return FAILURE;
        }
        if(cballs_memory_preflight(total,"Ly-alpha histograms, forest nodes and frontiers",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
        lya_cell_leg *next=realloc(ws->legs,bytes);
        if(!next) {snprintf(err,_ERRORMSGSIZE_,"Ly-alpha node frontier allocation failed");return FAILURE;}
        ws->legs=next;ws->capacity=cap;
    }
    ws->legs[ws->count++]=leg;return SUCCESS;
}

static int lya_cell_frontier(const lya_cells *tree,struct cmdline_data *cmd,
                              size_t ip,size_t iq,lya_cell_workspace *ws,ErrorMsg err)
{
    const lya_cell *p=&tree->nodes[ip],*q=&tree->nodes[iq];
    if(p->forest==q->forest) return SUCCESS;
    lya_cell_leg g=lya_cell_get_geometry(tree,cmd,ip,iq,ws);
    if(g.status<0) return SUCCESS;
    if(g.status==0 && q->left!=SIZE_MAX && (p->left==SIZE_MAX || q->radius>=p->radius)) {
        if(lya_cell_frontier(tree,cmd,ip,q->left,ws,err)==FAILURE) return FAILURE;
        return lya_cell_frontier(tree,cmd,ip,q->right,ws,err);
    }
    return lya_cell_append(ws,g,tree,cmd,err);
}

static uintmax_t lya_cell_count_max(void)
{return (((uintmax_t)1 << (sizeof(INTEGER)*CHAR_BIT-1))-1);}

static int lya_cell_count(lya_worker_hist *w,size_t np,size_t nq,size_t nr,ErrorMsg err)
{
    size_t count;
    if(!cballs_size_mul(np,nq,&count) || !cballs_size_mul(count,nr,&count)
        || (uintmax_t)count>(lya_cell_count_max()-(uintmax_t)w->ordered_triplet_count)/2) {
        snprintf(err,_ERRORMSGSIZE_,"Ly-alpha ordered triplet count overflow");return FAILURE;
    }
    w->ordered_triplet_count+=(INTEGER)(2*count);return SUCCESS;
}

/* Small unresolved products are cheaper as cached direct tiles than as
 * repeated node-bound calculations. Every coordinate uses reference pixels. */
static lya_neighbor lya_cell_cached_pixel(const lya_cells *tree,struct cmdline_data *cmd,
                                           lya_cell_workspace *ws,size_t ip,size_t index)
{
    lya_neighbor a={0};int radial,polar;
    if(ws->cached_pivot==ip) {
        /* A singleton frontier refines every mixed leg until all descendants
         * are valid. Its accepted nodes (and any subsequent children) already
         * have their leaf geometry cached for this pivot. Avoid billions of
         * redundant tag checks in the direct fallback. Multi-pivot tasks below
         * still fill their dedicated slots on demand. */
        const lya_cell_leg *g=&ws->point_cache[tree->leaves[index]];
        ws->cache_hits++;
        if(g->status!=1) return a;
        for(int k=0;k<3;k++) a.displacement[k]=g->dr[k];
        a.radius=g->radius;radial=g->radial;polar=g->polar;
    } else {
        const lya_cell_pixel_geometry *g=lya_cell_pixel_cached(tree,cmd,ip,tree->leaves[index],ws);
        if(g->status!=1) return a;
        for(int k=0;k<3;k++) a.displacement[k]=g->dr[k];
        a.radius=g->radius;radial=g->radial;polar=g->polar;
    }
    bodyptr q=tree->points[index];a.weight=Weight(q);a.weighted_delta=Weight(q)*Kappa(q);
    a.first_index=(((size_t)radial*cmd->lya3RBins*cmd->lya3ThetaBins+polar)*cmd->lya3ThetaBins)*cmd->lya3MuBins;
    a.second_index=((size_t)radial*cmd->lya3ThetaBins*cmd->lya3ThetaBins+polar)*cmd->lya3MuBins;
    return a;
}

static int lya_cell_direct_tile(const lya_cells *tree,struct cmdline_data *cmd,
                                 const lya_cell *p,const lya_cell *q,const lya_cell *r,
                                 lya_worker_hist *w,lya_cell_workspace *ws,ErrorMsg err)
{
    lya_neighbor aa[8],bb[8];
    for(size_t i=p->begin;i<p->end;i++) {
        bodyptr pivot=tree->points[i];
        size_t na=0,nb=0;
        for(size_t j=q->begin;j<q->end;j++) {
            lya_neighbor a=lya_cell_cached_pixel(tree,cmd,ws,tree->leaves[i],j);if(a.radius>0) aa[na++]=a;
        }
        for(size_t j=r->begin;j<r->end;j++) {
            lya_neighbor b=lya_cell_cached_pixel(tree,cmd,ws,tree->leaves[i],j);if(b.radius>0) bb[nb++]=b;
        }
        if(lya_cell_count(w,1,na,nb,err)==FAILURE) return FAILURE;
        w->direct_pairs+=na*nb;
        for(size_t j=0;j<na;j++) {
            const lya_neighbor *a=&aa[j];REAL pw=Weight(pivot)*a->weight,pn=(Weight(pivot)*Kappa(pivot))*a->weighted_delta;
            REAL cosine[8];
#pragma omp simd
            for(size_t k=0;k<nb;k++) {REAL dot;DOTVP(dot,a->displacement,bb[k].displacement);cosine[k]=dot/(a->radius*bb[k].radius);}
            for(size_t k=0;k<nb;k++) {
                int m=lya_bin_mu(cosine[k],cmd->lya3MuBins);
                lya_deposit3(w,a->first_index+bb[k].second_index+m,bb[k].first_index+a->second_index+m,pn*bb[k].weighted_delta,pw*bb[k].weight);
            }
        }
    }
    return SUCCESS;
}

/* Estimate bin uncertainty contributed by each splittable cell. The pivot
 * affects both legs and its LOS; resolved leg bins need only mu refinement.
 * This is a traversal heuristic, never an aggregation acceptance criterion. */
static int lya_cell_split_choice(const lya_cell *p,const lya_cell *q,const lya_cell *r,
                                  lya_cell_leg a,lya_cell_leg b,const lya_cell_workspace *ws)
{
    REAL ia=a.radius>DBL_MIN?1/a.radius:1/DBL_MIN;
    REAL ib=b.radius>DBL_MIN?1/b.radius:1/DBL_MIN;
    REAL wa=ws->angular_scale*ia+(a.status==0?ws->radial_scale+ws->polar_scale*ia:0);
    REAL wb=ws->angular_scale*ib+(b.status==0?ws->radial_scale+ws->polar_scale*ib:0);
    REAL sp=p->left==SIZE_MAX?-1:p->radius*(wa+wb)
        +(a.status==0 || b.status==0?p->los_spread*ws->polar_scale:0);
    REAL sq=q->left==SIZE_MAX?-1:q->radius*wa;
    REAL sr=r->left==SIZE_MAX?-1:r->radius*wb;
    /* Always choose an actual branch, even when extreme geometry overflows
     * the heuristic score. Acceptance still uses the guarded bounds above. */
    int split=-1;REAL best=-1;
    if(p->left!=SIZE_MAX) {split=0;best=isfinite(sp)?sp:DBL_MAX;}
    if(q->left!=SIZE_MAX && (split<0 || sq>best || !isfinite(sq))) {
        split=1;best=isfinite(sq)?sq:DBL_MAX;
    }
    if(r->left!=SIZE_MAX && (split<0 || sr>best || !isfinite(sr))) split=2;
    return split;
}

static int lya_cell_triple(const lya_cells *tree,struct cmdline_data *cmd,size_t ip,
                            lya_cell_leg a,lya_cell_leg b,lya_worker_hist *w,
                            lya_cell_workspace *ws,ErrorMsg err)
{
    if(a.status<0 || b.status<0) return SUCCESS;
    const lya_cell *p=&tree->nodes[ip],*q=&tree->nodes[a.node],*r=&tree->nodes[b.node];
    const int leaf=p->left==SIZE_MAX && q->left==SIZE_MAX && r->left==SIZE_MAX;
    int bm=-1,approximate=a.approximate||b.approximate;
    if(a.status==1 && b.status==1) {
        REAL dot;DOTVP(dot,a.dr,b.dr);
        REAL mu=dot/(a.radius*b.radius);
        if(leaf) bm=lya_bin_mu(mu,cmd->lya3MuBins);
        else if(sizeof(REAL)==sizeof(double)) {
#ifndef __FAST_MATH__
            REAL lo,hi;lya_cell_dot_bounds(a.lo,a.hi,b.lo,b.hi,&lo,&hi);
            bm=lya_cell_bin(lo,hi,mu,2.,cmd->lya3MuBins,cmd->lya3MuSlop,2,&approximate);
#endif
        }
    }
    if(bm>=0) {
        REAL d,n;
        long double ld=p->weight*q->weight*r->weight,ln=p->field*q->field*r->field;
        if(leaf) {
            /* Deliberately match the independent direct kernel's operation order. */
            bodyptr pp=tree->points[p->begin],qq=tree->points[q->begin],rr=tree->points[r->begin];
            d=(Weight(pp)*Weight(qq))*Weight(rr);
            n=((Weight(pp)*Kappa(pp))*(Weight(qq)*Kappa(qq)))*(Weight(rr)*Kappa(rr));
        } else {d=(REAL)ld;n=(REAL)ln;}
        if(leaf || (p->safe && q->safe && r->safe && isfinite(d) && isfinite(n)
                     && (d>=DBL_MIN || ld==0))) {
            size_t np=p->end-p->begin,nq=q->end-q->begin,nr=r->end-r->begin;
            if(lya_cell_count(w,np,nq,nr,err)==FAILURE) return FAILURE;
            size_t forward=lya_index3(a.radial,b.radial,a.polar,b.polar,bm,cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3MuBins);
            size_t reverse=lya_index3(b.radial,a.radial,b.polar,a.polar,bm,cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3MuBins);
            lya_deposit3(w,forward,reverse,n,d);
            unsigned long long count=(unsigned long long)np*nq*nr;
            if(leaf) w->direct_pairs++; else {w->aggregated_pairs+=count;w->segment_accepts++;}
            if(approximate) w->approximate_pairs+=count;
            if(np>1) ws->pivot_aggregates++;
            return SUCCESS;
        }
    }
    if(p->end-p->begin<=8 && q->end-q->begin<=8 && r->end-r->begin<=8)
        return lya_cell_direct_tile(tree,cmd,p,q,r,w,ws,err);
    /* Refine only this unresolved Cartesian product; already accepted triples
     * are not revisited when the pivot splits. Children partition each node. */
    int split=lya_cell_split_choice(p,q,r,a,b,ws);
    if(split==0 && p->left!=SIZE_MAX) {
        size_t children[2]={p->left,p->right};
        for(int k=0;k<2;k++) {
            lya_cell_leg aa=lya_cell_get_geometry(tree,cmd,children[k],a.node,ws),bb=lya_cell_get_geometry(tree,cmd,children[k],b.node,ws);
            if(lya_cell_triple(tree,cmd,children[k],aa,bb,w,ws,err)==FAILURE) return FAILURE;
        }
    } else if(split==1 && q->left!=SIZE_MAX) {
        lya_cell_leg aa=lya_cell_get_geometry(tree,cmd,ip,q->left,ws);
        if(lya_cell_triple(tree,cmd,ip,aa,b,w,ws,err)==FAILURE) return FAILURE;
        aa=lya_cell_get_geometry(tree,cmd,ip,q->right,ws);
        if(lya_cell_triple(tree,cmd,ip,aa,b,w,ws,err)==FAILURE) return FAILURE;
    } else if(r->left!=SIZE_MAX) {
        lya_cell_leg bb=lya_cell_get_geometry(tree,cmd,ip,r->left,ws);
        if(lya_cell_triple(tree,cmd,ip,a,bb,w,ws,err)==FAILURE) return FAILURE;
        bb=lya_cell_get_geometry(tree,cmd,ip,r->right,ws);
        if(lya_cell_triple(tree,cmd,ip,a,bb,w,ws,err)==FAILURE) return FAILURE;
    }
    return SUCCESS;
}

#include "lya_radial_moments.h"

static int lya_cells_run(struct cmdline_data *cmd,struct global_data *gd,
                          bodyptr base,INTEGER count,INTEGER first,INTEGER last,
                          nodeptr root,int compute2,size_t bins2,size_t bins3,size_t plan,
                          REAL *num2,REAL *den2,REAL *num3,REAL *den3,
                          INTEGER *visits,INTEGER *pairs,INTEGER *triplets,
                          unsigned long long *aggregated,unsigned long long *direct,
                          unsigned long long *accepts,unsigned long long *approximate)
{
    lya_cells tree;
    double build_start=CPUTIME;
    if(lya_cells_build(&tree,cmd,base,count,first,last,plan,cmd->error_message)==FAILURE) return FAILURE;
    double forest_build_cpu=CPUTIME-build_start;
    if(compute2 && cmd->lya2Kernel==1) {
        verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
            "Ly-alpha 2PCF forest build: nodes=%zu build_CPU=%g shared_with_3pcf=1\n",
            tree.nodes_count,forest_build_cpu);
        if(lya_pairs_run(&tree,cmd,gd,base,first,last,bins2,tree.plan,num2,den2,pairs,visits)==FAILURE) {
            lya_cells_free(&tree);return FAILURE;
        }
        compute2=0; /* Pair products are complete; preserve the shared hierarchy. */

    }
    /* Sparse pivot groups cannot amortize persistent per-pivot geometry and
     * mixed-forest frontier setup. Preserve the exact segment path in that
     * regime. This fixed geometric decision is independent of threads/ranks;
     * explicitly requested geometry slop retains the cell traversal. */
    if(cmd->lya3Kernel==5 && cmd->lya3MuSlop==0 && cmd->lya3RadialSlop==0
        && cmd->lya3PolarSlop==0 && tree.tasks_count>tree.count/4) {
        lya_cells_free(&tree);
        gd->lyaHierarchyPixelFallback=TRUE;
        verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
            "Ly-alpha hierarchy: sparse pivot groups; exact pixel-segment fallback\n");
        return 2; /* private dispatch result, distinct from SUCCESS/FAILURE */
    }
    build_start=CPUTIME; /* Cache setup excludes the already measured pair pass. */
    size_t cache_bytes,cache_total,pixel_slots,pixel_bytes;
    if(!cballs_size_mul(tree.max_pivots>1?tree.max_pivots:0,tree.count,&pixel_slots)
        || !cballs_size_mul(pixel_slots,sizeof(lya_cell_pixel_geometry)+sizeof(size_t),&pixel_bytes)
        || !cballs_size_mul(tree.nodes_count,sizeof(lya_cell_leg)+sizeof(size_t),&cache_bytes)
        || !cballs_size_add(cache_bytes,pixel_bytes,&cache_bytes)
        || !cballs_size_add(cache_bytes,LYA_CELL_PAIR_CACHE_SIZE*sizeof(lya_cell_cache_entry),&cache_bytes)
        || !cballs_size_mul(cache_bytes,(size_t)MAX(1,cmd->numthreads),&cache_total)
        || !cballs_size_add(tree.plan,cache_total,&tree.plan)) {
        lya_cells_free(&tree);snprintf(cmd->error_message,_ERRORMSGSIZE_,"Ly-alpha geometry cache dimensions overflow");return FAILURE;
    }
    if(cballs_memory_preflight(tree.plan,"Ly-alpha persistent nodes, histograms and geometry caches",cmd->error_message,_ERRORMSGSIZE_)==FAILURE) {lya_cells_free(&tree);return FAILURE;}
    double build_cpu=forest_build_cpu+CPUTIME-build_start;
    int failed=0;ErrorMsg error="";
    size_t block_size=cmd->lya3PivotBlock?(size_t)cmd->lya3PivotBlock:8;
    size_t blocks=tree.tasks_count/block_size+(tree.tasks_count%block_size!=0);
    const size_t first_block=lya_parallel_first(cmd),block_stride=lya_parallel_stride(cmd);
    unsigned long long pivot_aggregates=0;
    unsigned long long geometry_evaluations=0,cache_hits=0,pruned_nodes=0,leaf_evaluations=0,pair_cache_hits=0;
#pragma omp parallel
    {
        lya_worker_hist w;
        lya_cell_workspace ws={0};
        lya_radial_workspace radial={0};
        ErrorMsg local_error="";
        int ready=lya_worker_init(cmd,&w,bins2,bins3,compute2,1,tree.plan,local_error)==SUCCESS;
        if(ready && (cballs_calloc_checked((void**)&ws.point_cache,tree.nodes_count,sizeof(lya_cell_leg),
            "Ly-alpha point geometry cache",local_error,_ERRORMSGSIZE_)==FAILURE
            || cballs_calloc_checked((void**)&ws.point_tags,tree.nodes_count,sizeof(size_t),
            "Ly-alpha geometry cache tags",local_error,_ERRORMSGSIZE_)==FAILURE
            || cballs_calloc_checked((void**)&ws.pair_cache,LYA_CELL_PAIR_CACHE_SIZE,sizeof(lya_cell_cache_entry),
            "Ly-alpha pivot subdivision cache",local_error,_ERRORMSGSIZE_)==FAILURE
            || (pixel_slots && (cballs_calloc_checked((void**)&ws.pixel_cache,pixel_slots,sizeof(lya_cell_pixel_geometry),
            "Ly-alpha active pivot pixel cache",local_error,_ERRORMSGSIZE_)==FAILURE
            || cballs_calloc_checked((void**)&ws.pixel_tags,pixel_slots,sizeof(size_t),
            "Ly-alpha active pivot pixel tags",local_error,_ERRORMSGSIZE_)==FAILURE)))) {lya_worker_free(&w);ready=0;}
        ws.cached_pivot=SIZE_MAX;
        ws.radial_scale=cmd->lya3RBins/cmd->lya3RMax;
        ws.polar_scale=cmd->lya3ThetaBins/LYA_PI;
        ws.angular_scale=.5*cmd->lya3MuBins;
        int worker_failed=!ready;
#pragma omp for schedule(static,1) ordered
        for(size_t block=first_block;block<blocks;block+=block_stride) {
            w.accepted_visits=w.pair_count=w.ordered_triplet_count=0;
            size_t begin=block*block_size,end=MIN(tree.tasks_count,begin+block_size);
            for(size_t task=begin;!worker_failed && task<end;task++) {
                size_t ip=tree.tasks[task];const lya_cell *p=&tree.nodes[ip];
                /* Combined estimators retain the independent exact 2PCF path. */
                if(compute2) for(size_t j=p->begin;j<p->end;j++)
                    if(lya_walktree(cmd,tree.points[j],root,hypot(cmd->lya2RpMax,cmd->lya2RtMax),1,0,&w,local_error)==FAILURE) {worker_failed=1;break;}
                ws.count=0;ws.cached_pivot=p->left==SIZE_MAX?ip:SIZE_MAX;ws.task_begin=p->begin;
                for(size_t j=0;!worker_failed && j<tree.roots_count;j++)
                    if(lya_cell_frontier(&tree,cmd,ip,tree.roots[j],&ws,local_error)==FAILURE) worker_failed=1;
                if(cmd->lya3Kernel==5) {
                    if(!worker_failed && lya_radial_run(&tree,cmd,ip,&w,&ws,&radial,local_error)==FAILURE) worker_failed=1;
                } else for(size_t j=0;!worker_failed && j<ws.count;j++) for(size_t k=j+1;!worker_failed && k<ws.count;k++) {
                    if(tree.nodes[ws.legs[j].node].forest==tree.nodes[ws.legs[k].node].forest) continue;
                    if(lya_cell_triple(&tree,cmd,ip,ws.legs[j],ws.legs[k],&w,&ws,local_error)==FAILURE) worker_failed=1;
                }
            }
#pragma omp ordered
            {
                if(!worker_failed && !failed) {
                    if((uintmax_t)w.ordered_triplet_count>lya_cell_count_max()-(uintmax_t)*triplets) {
                        worker_failed=1;snprintf(local_error,sizeof(local_error),"Ly-alpha global triplet count overflow");
                    } else {
                        lya_commit_worker(&w,num2,den2,num3,den3,compute2,1);
                        *visits+=w.accepted_visits;*pairs+=w.pair_count;*triplets+=w.ordered_triplet_count;
                    }
                }
                if(worker_failed && !failed) {failed=1;snprintf(error,sizeof(error),"%s",local_error);}
            }
        }
#pragma omp critical(lya_cells_finish)
        {
            if(worker_failed && !failed) {failed=1;snprintf(error,sizeof(error),"%s",local_error);}
            *aggregated+=w.aggregated_pairs;*direct+=w.direct_pairs;*accepts+=w.segment_accepts;*approximate+=w.approximate_pairs;
            gd->lyaHierarchyCounts[0]+=radial.built;gd->lyaHierarchyCounts[1]+=radial.visited;
            gd->lyaHierarchyCounts[2]+=radial.accepted;gd->lyaHierarchyCounts[3]+=radial.represented;
            pivot_aggregates+=ws.pivot_aggregates;
            geometry_evaluations+=ws.geometry_evaluations;cache_hits+=ws.cache_hits;pruned_nodes+=ws.pruned_nodes;
            leaf_evaluations+=ws.leaf_evaluations;pair_cache_hits+=ws.pair_cache_hits;
        }
        free(radial.nodes);free(ws.legs);free(ws.point_cache);free(ws.point_tags);free(ws.pair_cache);free(ws.pixel_cache);free(ws.pixel_tags);if(ready) lya_worker_free(&w);
    }
    verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
        "Persistent forest cells: forests=%zu nodes=%zu pivot_tasks=%zu pivot_aggregates=%llu build_CPU=%g task_block=%zu\n",
        tree.roots_count,tree.nodes_count,tree.tasks_count,pivot_aggregates,build_cpu,block_size);
    verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
        "Forest geometry: evaluations=%llu cache_hits=%llu pruned_nodes=%llu leaf_evaluations=%llu pair_cache_hits=%llu\n",
        geometry_evaluations,cache_hits,pruned_nodes,leaf_evaluations,pair_cache_hits);
    lya_cells_free(&tree);
    if(failed) {snprintf(cmd->error_message,_ERRORMSGSIZE_,"%s",error);return FAILURE;}
    return SUCCESS;
}
