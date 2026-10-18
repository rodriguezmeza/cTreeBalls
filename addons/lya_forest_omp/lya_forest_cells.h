/* Private persistent forest hierarchy shared by 2PCF and 3PCF.
 * Cells contain one forest only; geometric bounds use all actual pixels.
 * Included once by search_lya_forest_omp.c, after estimator helpers. */
typedef struct {
    size_t begin, end, left, right, pivots;
    INTEGER forest, idlo, idhi;
    REAL chilo, chihi;
    REAL lo[3], hi[3], center[3], loslo[3], loshi[3], sight[3], radius, tube, los_spread;
    long double weight, field;
    int safe;
} lya_cell;

typedef struct {
    bodyptr *points;
    lya_cell *nodes;
    size_t *roots, *tasks, *leaves;
    size_t count, nodes_count, roots_count, tasks_count, plan, max_pivots;
    double sort_wall,nodes_wall;
    size_t sort_scratch_bytes;
    int sort_threads,build_threads;
} lya_cells;

static int lya_build_team_size(void)
{
#ifdef OPENMPCODE
    return omp_get_num_threads();
#else
    return 1;
#endif
}

static int lya_cell_point_compare(const void *aa, const void *bb)
{
    bodyptr a=*(bodyptr const *)aa, b=*(bodyptr const *)bb;
    if (LyaForestId(a)!=LyaForestId(b)) return LyaForestId(a)<LyaForestId(b)?-1:1;
    if (LyaDistance(a)!=LyaDistance(b)) return LyaDistance(a)<LyaDistance(b)?-1:1;
    return Id(a)<Id(b)?-1:(Id(a)>Id(b));
}

/* Distance to a finite segment, used only on the guarded finite scale below.
 * Capsule radii and distance bounds receive absolute-coordinate padding. */
static REAL lya_cell_segment_distance(const REAL *x,const REAL *a,const REAL *b)
{
    REAL v[3],d[3],vv=0,dv=0;
    for(int k=0;k<3;k++) {v[k]=b[k]-a[k];d[k]=x[k]-a[k];vv+=v[k]*v[k];dv+=d[k]*v[k];}
    REAL t=vv>DBL_MIN?lya_clamp(dv/vv,0.,1.):0.,r2=0;
    for(int k=0;k<3;k++) {REAL r=d[k]-t*v[k];r2+=r*r;}
    return sqrt(r2);
}

/* A binary subtree with n pixels occupies exactly 2*n-1 preorder slots.
 * Reserving those slots before workers run removes the shared node counter. */
static size_t lya_cell_build(lya_cells *tree,size_t begin,size_t end,size_t index,
                              bodyptr base,INTEGER first,INTEGER last)
{
    lya_cell *s=&tree->nodes[index];
    s->begin=begin;s->end=end;s->left=s->right=SIZE_MAX;
    s->forest=LyaForestId(tree->points[begin]);s->safe=1;
    if (end-begin==1) {
        bodyptr p=tree->points[begin];tree->leaves[begin]=index;
        REAL v=Weight(p),q=v*Kappa(p);
        s->weight=v;s->field=q;
        s->chilo=s->chihi=LyaDistance(p);s->idlo=s->idhi=Id(p);
        s->pivots=(p-base>=first && p-base<last);
        s->safe=(v==0 || (v>=1e-50 && v<=1e50))
             && (q==0 || (fabs(q)>=1e-50 && fabs(q)<=1e50));
        for(int k=0;k<3;k++) {
            s->lo[k]=s->hi[k]=s->center[k]=Pos(p)[k];
            s->loslo[k]=s->loshi[k]=s->sight[k]=LyaLOS(p)[k];
        }
    } else {
        size_t mid=begin+(end-begin)/2;
        s->left=lya_cell_build(tree,begin,mid,index+1,base,first,last);
        s->right=lya_cell_build(tree,mid,end,index+2*(mid-begin),base,first,last);
        const lya_cell *a=&tree->nodes[s->left],*b=&tree->nodes[s->right];
        s->weight=a->weight+b->weight;s->field=a->field+b->field;
        s->safe=a->safe&&b->safe;s->pivots=a->pivots+b->pivots;
        s->chilo=a->chilo;s->chihi=b->chihi;
        s->idlo=MIN(a->idlo,b->idlo);s->idhi=MAX(a->idhi,b->idhi);
        REAL r2=0,l2=0;
        for(int k=0;k<3;k++) {
            s->lo[k]=MIN(a->lo[k],b->lo[k]);s->hi[k]=MAX(a->hi[k],b->hi[k]);
            s->center[k]=.5*s->lo[k]+.5*s->hi[k];
            REAL d=MAX(fabs(s->center[k]-s->lo[k]),fabs(s->hi[k]-s->center[k]));r2+=d*d;
            s->loslo[k]=MIN(a->loslo[k],b->loslo[k]);s->loshi[k]=MAX(a->loshi[k],b->loshi[k]);
            s->sight[k]=.5*s->loslo[k]+.5*s->loshi[k];l2+=s->sight[k]*s->sight[k];
        }
        s->radius=sqrt(r2);
        REAL scale=0,los2=0;
        for(int k=0;k<3;k++) {
            scale=MAX(scale,MAX(fabs(s->lo[k]),fabs(s->hi[k])));
            REAL d=s->loshi[k]-s->loslo[k];los2+=d*d;
        }
        s->los_spread=sqrt(los2);
        s->tube=INFINITY;
#ifndef __FAST_MATH__
        if(sizeof(REAL)==sizeof(double) && scale<sqrt(DBL_MAX)/32) {
            /* The segment need not follow an ideal observer ray. Every actual
             * point is enclosed, so bent/reversed/degenerate forests are safe. */
            REAL tube=0;
            for(size_t j=begin;j<end;j++)
                tube=MAX(tube,lya_cell_segment_distance(Pos(tree->points[j]),
                    Pos(tree->points[begin]),Pos(tree->points[end-1])));
            s->tube=tube+4096*DBL_EPSILON*scale+16*sqrt(DBL_MIN);
        }
#endif
        for(int k=0;k<3;k++) s->sight[k]=l2>0?s->sight[k]/sqrt(l2):0;
    }
    return index;
}

static void lya_cell_tasks(lya_cells *tree,size_t node,size_t cap,REAL max_radius)
{
    lya_cell *s=&tree->nodes[node];
    if (!s->pivots) return;
    if (s->pivots==s->end-s->begin && s->pivots<=cap && s->radius<=max_radius) tree->tasks[tree->tasks_count++]=node;
    else {lya_cell_tasks(tree,s->left,cap,max_radius);lya_cell_tasks(tree,s->right,cap,max_radius);}
}

static void lya_cells_free(lya_cells *t)
{free(t->points);free(t->nodes);free(t->roots);free(t->tasks);free(t->leaves);memset(t,0,sizeof(*t));}

/* Co-rank an output boundary in a stable merge. This lets workers merge
 * independent output slices even at the final, whole-catalog merge level. */
static size_t lya_cells_merge_rank(bodyptr const *a,size_t na,bodyptr const *b,size_t nb,size_t k)
{
    size_t lo=k>nb?k-nb:0,hi=MIN(k,na);
    while(lo<hi) {
        size_t i=lo+(hi-lo)/2,j=k-i;
        if(j>0 && i<na && lya_cell_point_compare(b+j-1,a+i)>=0) lo=i+1;
        else hi=i;
    }
    return lo;
}

/* Sort disjoint fixed chunks and merge disjoint runs. The comparator includes
 * the immutable original ID, so serial and parallel layouts are identical.
 * The pointer scratch is temporary and checked against the peak build plan. */
static int lya_cells_sort(lya_cells *t,int threads,ErrorMsg err)
{
    const size_t grain=32768,count=t->count;
    double start=lya_progress_wall_time();
    t->sort_threads=1;
    if(threads<=1 || count<2*grain) {
        qsort(t->points,count,sizeof(bodyptr),lya_cell_point_compare);
    } else {
        bodyptr *scratch=NULL,*src=t->points,*dst;
        size_t bytes,peak;
        if(!cballs_size_mul(count,sizeof(*scratch),&bytes) || !cballs_size_add(t->plan,bytes,&peak)) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha parallel sort dimensions overflow");return FAILURE;
        }
        if(cballs_memory_preflight(peak,"Ly-alpha parallel sort peak workspace",err,_ERRORMSGSIZE_)==FAILURE
            || cballs_calloc_checked((void**)&scratch,count,sizeof(*scratch),
                "Ly-alpha parallel sort pointers",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
        t->sort_scratch_bytes=bytes;dst=scratch;
#pragma omp parallel num_threads(threads) shared(src,dst)
        {
#pragma omp single
            t->sort_threads=lya_build_team_size();
#pragma omp for schedule(static)
            for(size_t chunk=0;chunk<(count+grain-1)/grain;chunk++) {
                size_t begin=chunk*grain;
                qsort(t->points+begin,MIN(grain,count-begin),sizeof(bodyptr),lya_cell_point_compare);
            }
            for(size_t width=grain;width<count;width=width>count/2?count:2*width) {
                /* Pointer allocation bounds make 2*width safe. */
                size_t stride=2*width;
#pragma omp for schedule(static)
                for(size_t chunk=0;chunk<(count+grain-1)/grain;chunk++) {
                    size_t out=chunk*grain,begin=out/stride*stride;
                    size_t mid=MIN(begin+width,count),end=MIN(begin+stride,count);
                    size_t offset=out-begin;
                    size_t left=lya_cells_merge_rank(src+begin,mid-begin,src+mid,end-mid,offset);
                    size_t i=begin+left,j=mid+offset-left,stop=MIN(out+grain,count);
                    for(size_t k=out;k<stop;k++) {
                        if(i<mid && (j==end || lya_cell_point_compare(src+i,src+j)<=0)) dst[k]=src[i++];
                        else dst[k]=src[j++];
                    }
                }
#pragma omp single
                {bodyptr *swap=src;src=dst;dst=swap;}
            }
            if(src!=t->points) {
#pragma omp for schedule(static)
                for(size_t i=0;i<count;i++) t->points[i]=src[i];
            }
        }
        free(scratch);
    }
    t->sort_wall=lya_progress_wall_time()-start;
    return SUCCESS;
}

static int lya_cells_build(lya_cells *t,struct cmdline_data *cmd,bodyptr base,
                            INTEGER count,INTEGER first,INTEGER last,size_t plan,ErrorMsg err)
{
    size_t bytes,nodes,capacity;
    memset(t,0,sizeof(*t));
    if (count<0 || (uintmax_t)count>SIZE_MAX) goto overflow;
    capacity=(size_t)count;
    if (!cballs_size_mul(capacity,2,&nodes)
        || !cballs_size_mul(capacity,2*sizeof(lya_cell)+sizeof(bodyptr)+3*sizeof(size_t),&bytes)
        || !cballs_size_add(plan,bytes,&t->plan)) goto overflow;
    if (cballs_memory_preflight(t->plan,"Ly-alpha persistent forest nodes and histograms",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
    if (cballs_calloc_checked((void**)&t->points,capacity,sizeof(bodyptr),"Ly-alpha forest points",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&t->nodes,nodes,sizeof(lya_cell),"Ly-alpha forest nodes",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&t->roots,capacity,sizeof(size_t),"Ly-alpha forest roots",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&t->tasks,capacity,sizeof(size_t),"Ly-alpha pivot cells",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&t->leaves,capacity,sizeof(size_t),"Ly-alpha leaf indices",err,_ERRORMSGSIZE_)==FAILURE) {
        lya_cells_free(t);return FAILURE;
    }
    for(INTEGER i=0;i<count;i++) if(Update(base+i)!=FALSE && Mask(base+i)==MASK_NODE_VALID) t->points[t->count++]=base+i;
    int threads=lya_forest_is_rprt_method(cmd->searchMethod) && !scanopt(cmd->options,"lya-serial-build")
        ?MAX(1,cmd->numthreads):1;
    if(lya_cells_sort(t,threads,err)==FAILURE) {lya_cells_free(t);return FAILURE;}
    double nodes_start=lya_progress_wall_time();
    /* Discover forests and reserve disjoint preorder ranges serially. The
     * expensive bounds, capsules and moments are then independent by forest. */
    for(size_t begin=0;begin<t->count;) {
        size_t end=begin+1;
        while(end<t->count && LyaForestId(t->points[end])==LyaForestId(t->points[begin])) end++;
        size_t root=t->nodes_count;t->nodes_count+=2*(end-begin)-1;
        t->roots[t->roots_count++]=root;
        t->nodes[root].begin=begin;t->nodes[root].end=end;
        begin=end;
    }
    t->build_threads=1;
#pragma omp parallel num_threads(threads) if(threads>1 && t->count>=32768 && t->roots_count>1)
    {
#pragma omp single
        t->build_threads=lya_build_team_size();
#pragma omp for schedule(dynamic,1)
        for(size_t row=0;row<t->roots_count;row++) {
            size_t root=t->roots[row],begin=t->nodes[root].begin,end=t->nodes[root].end;
            lya_cell_build(t,begin,end,root,base,first,last);
        }
    }
    /* This 2PCF-only method visits forest roots directly. Its unused 3PCF
     * pivot frontier need not be traversed and materialized after the build. */
    if(!lya_forest_is_rprt_method(cmd->searchMethod))
    for(size_t row=0;row<t->roots_count;row++)
        lya_cell_tasks(t,t->roots[row],cmd->lya3Kernel>=4?(size_t)cmd->lya3PivotCellMax:1,
            cmd->lya3RBins>0?.125*(cmd->lya3RMax/cmd->lya3RBins):0);
    t->nodes_wall=lya_progress_wall_time()-nodes_start;
    for(size_t i=0;i<t->tasks_count;i++) {
        const lya_cell *p=&t->nodes[t->tasks[i]];
        t->max_pivots=MAX(t->max_pivots,p->end-p->begin);
    }
    return SUCCESS;
overflow:
    snprintf(err,_ERRORMSGSIZE_,"Ly-alpha persistent forest dimensions overflow");return FAILURE;
}

/* Enclosing interval for a dot product. Input directions include their own
 * normalization error; pad again for the three multiplications and additions. */
static void lya_cell_dot_bounds(const REAL *al,const REAL *ah,const REAL *bl,
                                 const REAL *bh,REAL *lo,REAL *hi)
{
    *lo=*hi=0;
    for(int k=0;k<3;k++) {
        REAL a=al[k]*bl[k],b=al[k]*bh[k],c=ah[k]*bl[k],d=ah[k]*bh[k];
        *lo+=MIN(MIN(a,b),MIN(c,d));*hi+=MAX(MAX(a,b),MAX(c,d));
    }
    *lo=MAX(-1.,*lo-512*DBL_EPSILON);*hi=MIN(1.,*hi+512*DBL_EPSILON);
}

static int lya_cell_bin(REAL lo,REAL hi,REAL center,REAL maximum,int bins,
                         REAL slop,int kind,int *approximate)
{
    int low=kind==0?lya_bin_positive(lo,maximum,bins):kind==1?lya_bin_theta(lo,bins):lya_bin_mu(lo,bins);
    int high=kind==0?lya_bin_positive(hi,maximum,bins):kind==1?lya_bin_theta(hi,bins):lya_bin_mu(hi,bins);
    if (low>=0 && low==high) return low;
    if (slop<=0 || !isfinite(center)) return -1;
    int bin=kind==0?lya_bin_positive(center,maximum,bins):kind==1?lya_bin_theta(center,bins):lya_bin_mu(center,bins);
    if (bin<0) return -1;
    REAL width=maximum/bins,lower=(kind==2?-1.:0.)+bin*width;
    if (lo<lower-slop*width || hi>lower+(1+slop)*width) return -1;
    *approximate=1;return bin;
}

