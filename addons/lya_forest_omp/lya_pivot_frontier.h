/* Forest-preserving pivot smoothing and an immutable octree scan frontier.
 * Smoothing is an explicit geometry approximation; radius zero never groups
 * pixels. No catalog body, field, LOS, Update flag or neighbor is modified.
 * Included privately by search_lya_forest_omp.c. */
typedef struct {
    bodyptr pivot;
    bodyptr *members;
    INTEGER idlo, idhi;
    REAL weight, field, radius;
    INTEGER count;
} lya_pivot_group;

typedef struct {
    lya_pivot_group *groups;
    bodyptr *members;
    size_t *order, *offsets;
    size_t count, ordered, blocks, active, plan;
    REAL max_radius;
} lya_pivot_frontier;

static void lya_pivot_frontier_free(lya_pivot_frontier *f)
{
    free(f->groups); free(f->members); free(f->order); free(f->offsets);
    memset(f,0,sizeof(*f));
}

static int lya_pivot_compare(const void *aa,const void *bb)
{
    bodyptr a=*(bodyptr const *)aa,b=*(bodyptr const *)bb;
    if(LyaForestId(a)!=LyaForestId(b)) return LyaForestId(a)<LyaForestId(b)?-1:1;
    if(LyaDistance(a)!=LyaDistance(b)) return LyaDistance(a)<LyaDistance(b)?-1:1;
    return Id(a)<Id(b)?-1:(Id(a)>Id(b));
}

static REAL lya_pivot_distance(bodyptr a,bodyptr b)
{
    return hypot(hypot((REAL)Pos(a)[0]-Pos(b)[0],
                       (REAL)Pos(a)[1]-Pos(b)[1]),
                       (REAL)Pos(a)[2]-Pos(b)[2]);
}

static void lya_pivot_frontier_flush(lya_pivot_frontier *f)
{
    if(f->ordered!=f->offsets[f->blocks])
        f->offsets[++f->blocks]=f->ordered;
}

static void lya_pivot_frontier_append(lya_pivot_frontier *f,size_t group,size_t block)
{
    f->order[f->ordered++]=group;
    if(f->ordered-f->offsets[f->blocks]>=block) lya_pivot_frontier_flush(f);
}

static void lya_pivot_frontier_walk(lya_pivot_frontier *f,nodeptr q,
                                    bodyptr base,const size_t *mapping,
                                    unsigned depth,unsigned level,size_t block)
{
    if(Type(q)==CELL) {
        for(nodeptr child=More(q);child!=Next(q);child=Next(child))
            lya_pivot_frontier_walk(f,child,base,mapping,depth+1,level,block);
    } else {
        size_t group=mapping[(bodyptr)q-base];
        if(group!=SIZE_MAX) lya_pivot_frontier_append(f,group,block);
    }
    /* Bound task size within a spatial cell. Shallow leaves form complete
     * cells too. The partition is independent of the OpenMP team size. */
    if(depth==level || (Type(q)!=CELL && depth<level)) lya_pivot_frontier_flush(f);
}

static int lya_pivot_frontier_build(lya_pivot_frontier *f,
                                    struct cmdline_data *cmd,bodyptr base,
                                    INTEGER count,INTEGER first,INTEGER last,
                                    nodeptr root,size_t block,int compute3,
                                    size_t plan,ErrorMsg err)
{
    bodyptr *points=NULL;
    size_t *mapping=NULL,n,bytes,slots,all;
    memset(f,0,sizeof(*f));
    if(count<0 || (uintmax_t)count>SIZE_MAX || first<0 || last<first || last>count)
        goto overflow;
    n=(size_t)count;
    /* Exact spatial tasks need only index arrays, never smoothing records or
     * a sortable member buffer. Keep pair-only setup small and linear. */
    const size_t point_bytes=2*sizeof(size_t)+(cmd->lyaPivotRadius>0
                              ?sizeof(lya_pivot_group)+sizeof(bodyptr):0);
    if(!cballs_size_add(n,1,&slots)
        || !cballs_size_mul(n,point_bytes,&bytes)
        || !cballs_size_mul(slots,sizeof(size_t),&all)
        || !cballs_size_add(bytes,all,&bytes)
        || !cballs_size_add(plan,bytes,&f->plan)) goto overflow;
    if(cballs_memory_preflight(f->plan,"Ly-alpha pivot frontier, sorting and histograms",err,_ERRORMSGSIZE_)==FAILURE)
        return FAILURE;
#define LYA_PIVOT_ALLOC(p,len) \
    if(cballs_calloc_checked((void**)&(p),(len),sizeof(*(p)),"Ly-alpha pivot frontier",err,_ERRORMSGSIZE_)==FAILURE) goto fail
    LYA_PIVOT_ALLOC(f->order,n);
    LYA_PIVOT_ALLOC(f->offsets,slots);
    if(cmd->lyaPivotRadius>0) {
        LYA_PIVOT_ALLOC(f->groups,n);
        LYA_PIVOT_ALLOC(points,n);
    }
    if(cmd->lyaScanLevel) {
        LYA_PIVOT_ALLOC(mapping,n);
        for(size_t j=0;j<n;j++) mapping[j]=SIZE_MAX;
    }
#undef LYA_PIVOT_ALLOC
    for(INTEGER j=first;j<last;j++)
        if(Update(base+j)!=FALSE && Mask(base+j)==MASK_NODE_VALID) {
            if(points) points[f->active]=base+j;
            else if(mapping) mapping[j]=(size_t)j;
            else lya_pivot_frontier_append(f,(size_t)j,block);
            f->active++;
        }

    if(cmd->lyaPivotRadius>0) {
        /* Multiplicity-weighted counts remain signed INTEGER values. Check a
         * conservative full-catalog bound before the accelerated arithmetic. */
        const uintmax_t maximum=((uintmax_t)1<<(sizeof(INTEGER)*CHAR_BIT-1))-1;
        size_t bound;
        if(!cballs_size_mul(f->active,n,&bound)
            || (compute3 && !cballs_size_mul(bound,n,&bound))
            || (uintmax_t)bound>maximum) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha smoothed pivot count bound exceeds INTEGER; reduce the catalog or disable lyaPivotRadius");
            goto fail;
        }
        qsort(points,f->active,sizeof(*points),lya_pivot_compare);
    }
    for(size_t begin=0;points && begin<f->active;) {
        bodyptr p=points[begin];size_t end=begin+1;
        long double weight=Weight(p),field=(REAL)(Weight(p)*Kappa(p));
        long double absolute_field=fabsl(field);
        REAL radius=0;
        if(cmd->lyaPivotRadius>0 && isfinite((REAL)field)) {
            while(end<f->active && end-begin<(size_t)cmd->lyaPivotMax) {
                bodyptr q=points[end];
                if(LyaForestId(q)!=LyaForestId(p)) break;
                REAL distance=lya_pivot_distance(p,q), qfield=Weight(q)*Kappa(q);
                if(!isfinite(distance) || distance>cmd->lyaPivotRadius || !isfinite(qfield)) break;
                long double next_weight=weight+Weight(q),next_field=field+qfield;
                long double next_absolute=absolute_field+fabsl(qfield);
                if(!isfinite((REAL)next_weight) || !isfinite((REAL)next_field)
                    || !isfinite((REAL)next_absolute)) break;
                weight=next_weight;field=next_field;absolute_field=next_absolute;
                radius=MAX(radius,distance);end++;
            }
        }
        size_t id=f->count++;
        INTEGER idlo=Id(p),idhi=Id(p);
        for(size_t j=begin+1;j<end;j++) {idlo=MIN(idlo,Id(points[j]));idhi=MAX(idhi,Id(points[j]));}
        f->groups[id]=(lya_pivot_group){p,points+begin,idlo,idhi,(REAL)weight,(REAL)field,radius,(INTEGER)(end-begin)};
        f->max_radius=MAX(f->max_radius,radius);
        if(mapping) mapping[p-base]=id;
        begin=end;
    }
    if(!points) f->count=f->active;
    if(mapping) lya_pivot_frontier_walk(f,root,base,mapping,0,(unsigned)cmd->lyaScanLevel,block);
    else if(points) for(size_t j=0;j<f->count;j++) lya_pivot_frontier_append(f,j,block);
    lya_pivot_frontier_flush(f);
    if(f->ordered!=f->count) {
        snprintf(err,_ERRORMSGSIZE_,"Ly-alpha pivot frontier did not cover every representative");goto fail;
    }
    f->members=points;free(mapping);return SUCCESS;
overflow:
    snprintf(err,_ERRORMSGSIZE_,"Ly-alpha pivot frontier dimensions overflow");
fail:
    free(points);free(mapping);lya_pivot_frontier_free(f);return FAILURE;
}
