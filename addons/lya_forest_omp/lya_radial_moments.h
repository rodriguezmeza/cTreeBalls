/* Exact hard-bin moment hierarchy. A node collects accepted neighbor cells
 * sharing one (radial, polar) bin. Forest-ID enclosures certify disjointness;
 * no subtraction of a large same-forest auto product is needed. Angular
 * enclosures cover every source pixel and every pivot in the current task.
 * Included after lya_cell_triple. Finite Legendre truncation is NOT used. */
typedef struct {
    size_t begin,end,left,right,count;
    INTEGER forest_lo,forest_hi;
    uint64_t forests[4];
    REAL lo[3],hi[3];
    long double weight,field;
    int radial,polar,valid,safe,approximate;
} lya_radial_moment;

typedef struct {size_t leg,bin,node;REAL azimuth;} lya_radial_entry;

typedef struct {
    lya_radial_moment *nodes;
    lya_radial_entry *entries;
    size_t count,capacity;
    uint64_t built,visited,accepted,represented;
} lya_radial_workspace;

/* Angular sorting improves the direction enclosure. A hashed forest mask is
 * used only as a conservative disjointness certificate: collisions force
 * refinement and can never admit a same-forest pair. ID ranges are an exact
 * additional certificate for larger disjoint groups. */
static int lya_radial_leg_compare(const void *aa,const void *bb)
{
    const lya_radial_entry *a=aa,*b=bb;
    if(a->bin!=b->bin) return a->bin<b->bin?-1:1;
    if(a->azimuth!=b->azimuth) return a->azimuth<b->azimuth?-1:1;
    return a->node<b->node?-1:(a->node>b->node);
}

static size_t lya_radial_build(const lya_cells *tree,const lya_cell_workspace *ws,
                                lya_radial_workspace *rw,size_t begin,size_t end)
{
    size_t index=rw->count++;
    lya_radial_moment *s=&rw->nodes[index];
    s->begin=begin;s->end=end;s->left=s->right=SIZE_MAX;
    if(end-begin==1) {
        const lya_cell_leg *g=&ws->legs[rw->entries[begin].leg];const lya_cell *q=&tree->nodes[g->node];
        s->count=q->end-q->begin;s->weight=q->weight;s->field=q->field;
        s->forest_lo=s->forest_hi=q->forest;s->safe=q->safe;
        uint64_t key=(uint64_t)q->forest*UINT64_C(0x9e3779b97f4a7c15);key^=key>>33;
        unsigned bit=(unsigned)(key&255);
        memset(s->forests,0,sizeof(s->forests));s->forests[bit/64]=UINT64_C(1)<<(bit%64);
        s->valid=g->status==1;s->radial=g->radial;s->polar=g->polar;s->approximate=g->approximate;
        for(int k=0;k<3;k++) {s->lo[k]=g->lo[k];s->hi[k]=g->hi[k];}
    } else {
        size_t mid=begin+(end-begin)/2;
        s->left=lya_radial_build(tree,ws,rw,begin,mid);
        s->right=lya_radial_build(tree,ws,rw,mid,end);
        const lya_radial_moment *a=&rw->nodes[s->left],*b=&rw->nodes[s->right];
        s->count=a->count+b->count;s->weight=a->weight+b->weight;s->field=a->field+b->field;
        s->forest_lo=MIN(a->forest_lo,b->forest_lo);s->forest_hi=MAX(a->forest_hi,b->forest_hi);
        s->valid=a->valid&&b->valid&&a->radial==b->radial&&a->polar==b->polar;
        s->radial=a->radial;s->polar=a->polar;s->safe=a->safe&&b->safe;
        s->approximate=a->approximate||b->approximate;
        for(int k=0;k<4;k++) s->forests[k]=a->forests[k]|b->forests[k];
        for(int k=0;k<3;k++) {s->lo[k]=MIN(a->lo[k],b->lo[k]);s->hi[k]=MAX(a->hi[k],b->hi[k]);}
    }
    return index;
}

static int lya_radial_product(const lya_cells *tree,struct cmdline_data *cmd,
                               size_t ip,size_t ia,size_t ib,lya_worker_hist *w,
                               lya_cell_workspace *ws,lya_radial_workspace *rw,ErrorMsg err)
{
    const lya_radial_moment *a=&rw->nodes[ia],*b=&rw->nodes[ib];
    rw->visited++;
    if(a->forest_lo==a->forest_hi && b->forest_lo==a->forest_lo && b->forest_hi==a->forest_lo) return SUCCESS;
    if(ia==ib) {
        if(a->left==SIZE_MAX) return SUCCESS;
        if(lya_radial_product(tree,cmd,ip,a->left,a->left,w,ws,rw,err)==FAILURE
           || lya_radial_product(tree,cmd,ip,a->left,a->right,w,ws,rw,err)==FAILURE) return FAILURE;
        return lya_radial_product(tree,cmd,ip,a->right,a->right,w,ws,rw,err);
    }
    const lya_cell *p=&tree->nodes[ip];
    int bounded=sizeof(REAL)==sizeof(double);
#ifdef __FAST_MATH__
    bounded=0;
#endif
    uint64_t common=0;for(int k=0;k<4;k++) common|=a->forests[k]&b->forests[k];
    if(bounded && a->valid && b->valid && a->safe && b->safe && p->safe
        && (common==0 || a->forest_hi<b->forest_lo || b->forest_hi<a->forest_lo)) {
        REAL lo,hi;lya_cell_dot_bounds(a->lo,a->hi,b->lo,b->hi,&lo,&hi);
        int approximate=a->approximate||b->approximate;
        int mu=isfinite(lo)&&isfinite(hi)
            ? lya_cell_bin(lo,hi,.5*lo+.5*hi,2.,cmd->lya3MuBins,cmd->lya3MuSlop,2,&approximate) : -1;
        long double ld=p->weight*a->weight*b->weight,ln=p->field*a->field*b->field;
        if(mu>=0 && isfinite((REAL)ld) && isfinite((REAL)ln) && ((REAL)ld>=DBL_MIN || ld==0)) {
            size_t np=p->end-p->begin;
            if(lya_cell_count(w,np,a->count,b->count,err)==FAILURE) return FAILURE;
            size_t ab=lya_index3(a->radial,b->radial,a->polar,b->polar,mu,cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3MuBins);
            size_t ba=lya_index3(b->radial,a->radial,b->polar,a->polar,mu,cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3MuBins);
            lya_deposit3(w,ab,ba,(REAL)ln,(REAL)ld);
            uint64_t represented=(uint64_t)np*a->count*b->count;
            w->aggregated_pairs+=represented;w->segment_accepts++;
            if(approximate) w->approximate_pairs+=represented;
            if(np>1) ws->pivot_aggregates++;
            rw->accepted++;rw->represented+=represented;
            return SUCCESS;
        }
    }
    if(a->left==SIZE_MAX && b->left==SIZE_MAX)
        return lya_cell_triple(tree,cmd,ip,ws->legs[rw->entries[a->begin].leg],ws->legs[rw->entries[b->begin].leg],w,ws,err);
    /* Split an unresolved radial combination, retaining every earlier deposit.
     * Prefer mixed-bin nodes, then the larger frontier, independently of w. */
    int split_a=a->left!=SIZE_MAX && (b->left==SIZE_MAX
        || (!a->valid && b->valid) || (a->valid==b->valid && a->end-a->begin>=b->end-b->begin));
    if(split_a) {
        if(lya_radial_product(tree,cmd,ip,a->left,ib,w,ws,rw,err)==FAILURE) return FAILURE;
        return lya_radial_product(tree,cmd,ip,a->right,ib,w,ws,rw,err);
    }
    if(lya_radial_product(tree,cmd,ip,ia,b->left,w,ws,rw,err)==FAILURE) return FAILURE;
    return lya_radial_product(tree,cmd,ip,ia,b->right,w,ws,rw,err);
}

static int lya_radial_run(const lya_cells *tree,struct cmdline_data *cmd,size_t ip,
                           lya_worker_hist *w,lya_cell_workspace *ws,
                           lya_radial_workspace *rw,ErrorMsg err)
{
    if(ws->count<2) return SUCCESS;
    if(ws->capacity>rw->capacity) {
        size_t bytes,per_worker,total;
        /* The shared plan reserves the fixed geometry caches; frontier memory
         * and possible realloc overlap are additional, for every worker. */
        if(!cballs_size_mul(ws->capacity,2*sizeof(lya_radial_moment)+sizeof(lya_radial_entry),&bytes)
            || !cballs_size_mul(ws->capacity,2*sizeof(lya_cell_leg)+4*sizeof(lya_radial_moment)+2*sizeof(lya_radial_entry),&per_worker)
            || !cballs_size_mul(per_worker,(size_t)MAX(1,cmd->numthreads),&total)
            || !cballs_size_add(total,tree->plan,&total)) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha radial moment scratch overflow");return FAILURE;
        }
        if(cballs_memory_preflight(total,"Ly-alpha radial moment hierarchy",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
        lya_radial_moment *next=realloc(rw->nodes,bytes);
        if(!next) {snprintf(err,_ERRORMSGSIZE_,"Ly-alpha radial moment allocation failed");return FAILURE;}
        rw->nodes=next;rw->capacity=ws->capacity;ws->radial_bytes=bytes;
        rw->entries=(lya_radial_entry*)(next+2*ws->capacity);
    }
    for(size_t j=0;j<ws->count;j++) {
        const lya_cell_leg *g=&ws->legs[j];
        rw->entries[j]=(lya_radial_entry){j,g->status==1?(size_t)g->radial*cmd->lya3ThetaBins+g->polar:SIZE_MAX,
            g->node,isfinite(g->dr[0])&&isfinite(g->dr[1])?atan2(g->dr[1],g->dr[0]):0};
    }
    qsort(rw->entries,ws->count,sizeof(*rw->entries),lya_radial_leg_compare);
    rw->count=0;lya_radial_build(tree,ws,rw,0,ws->count);rw->built+=rw->count;
    return lya_radial_product(tree,cmd,ip,0,0,w,ws,rw,err);
}
