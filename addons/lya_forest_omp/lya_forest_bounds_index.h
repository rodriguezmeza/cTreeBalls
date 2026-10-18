#ifndef LYA_FOREST_BOUNDS_INDEX_H
#define LYA_FOREST_BOUNDS_INDEX_H

/* Broad phase for lya2Kernel=1. Median partitions use forest-box centres,
 * but rejection ALWAYS uses unions of complete forest pixel bounds. The
 * existing radial trees and exact pair ownership/binning remain unchanged.
 * No globals, catalog mutation, geometry approximation or pair cache. */
#define LYA_FOREST_KD_LEAF 8
#define LYA_FOREST_KD_MIN 32

typedef struct {
    size_t begin, end, left, right, last_row;
    lya_cell bounds;
} lya_forest_kd_node;

typedef struct {
    size_t *rows;
    lya_forest_kd_node *nodes;
    size_t count, capacity, words, summaries, plan;
    int enabled;
} lya_forest_kd;

typedef struct { uint64_t *bits, *summary; } lya_forest_kd_bits;

static void lya_forest_kd_free(lya_forest_kd *t)
{ free(t->rows); free(t->nodes); memset(t,0,sizeof(*t)); }

static int lya_forest_kd_before(const lya_cells *t,size_t a,size_t b,int axis)
{
    REAL x=t->nodes[t->roots[a]].center[axis],y=t->nodes[t->roots[b]].center[axis];
    return x<y || (x==y && a<b);
}

/* Deterministic median-of-three QuickSelect; ties use immutable forest rows.
 * Each iteration removes the pivot, including for coincident box centres. */
static void lya_forest_kd_select(lya_forest_kd *kd,const lya_cells *t,
                                  size_t begin,size_t end,size_t target,int axis)
{
    while(end-begin>1) {
        size_t a=kd->rows[begin],b=kd->rows[begin+(end-begin)/2],c=kd->rows[end-1];
        if(lya_forest_kd_before(t,b,a,axis)) {size_t x=a;a=b;b=x;}
        if(lya_forest_kd_before(t,c,b,axis)) {size_t x=b;b=c;c=x;}
        if(lya_forest_kd_before(t,b,a,axis)) b=a;
        size_t at=begin;
        while(kd->rows[at]!=b) at++;
        kd->rows[at]=kd->rows[end-1];kd->rows[end-1]=b;
        size_t store=begin;
        for(size_t j=begin;j+1<end;j++) if(lya_forest_kd_before(t,kd->rows[j],b,axis)) {
            size_t x=kd->rows[store];kd->rows[store++]=kd->rows[j];kd->rows[j]=x;
        }
        kd->rows[end-1]=kd->rows[store];kd->rows[store]=b;
        if(store==target) return;
        if(target<store) end=store; else begin=store+1;
    }
}

static size_t lya_forest_kd_build_node(lya_forest_kd *kd,const lya_cells *t,
                                        size_t begin,size_t end)
{
    size_t id=kd->count++;
    lya_forest_kd_node *n=kd->nodes+id;
    n->begin=begin;n->end=end;n->left=n->right=SIZE_MAX;n->last_row=0;
    n->bounds.chilo=INFINITY;n->bounds.chihi=-INFINITY;
    REAL lo[3]={INFINITY,INFINITY,INFINITY},hi[3]={-INFINITY,-INFINITY,-INFINITY};
    for(int k=0;k<3;k++) {
        n->bounds.lo[k]=n->bounds.loslo[k]=INFINITY;
        n->bounds.hi[k]=n->bounds.loshi[k]=-INFINITY;
    }
    for(size_t j=begin;j<end;j++) {
        size_t row=kd->rows[j];const lya_cell *f=t->nodes+t->roots[row];
        n->last_row=MAX(n->last_row,row);
        n->bounds.chilo=MIN(n->bounds.chilo,f->chilo);
        n->bounds.chihi=MAX(n->bounds.chihi,f->chihi);
        for(int k=0;k<3;k++) {
            n->bounds.lo[k]=MIN(n->bounds.lo[k],f->lo[k]);
            n->bounds.hi[k]=MAX(n->bounds.hi[k],f->hi[k]);
            n->bounds.loslo[k]=MIN(n->bounds.loslo[k],f->loslo[k]);
            n->bounds.loshi[k]=MAX(n->bounds.loshi[k],f->loshi[k]);
            lo[k]=MIN(lo[k],f->center[k]);hi[k]=MAX(hi[k],f->center[k]);
        }
    }
    if(end-begin>LYA_FOREST_KD_LEAF) {
        int axis=0;
        for(int k=1;k<3;k++) if(.5*hi[k]-.5*lo[k]>.5*hi[axis]-.5*lo[axis]) axis=k;
        size_t mid=begin+(end-begin)/2;
        lya_forest_kd_select(kd,t,begin,end,mid,axis);
        n->left=lya_forest_kd_build_node(kd,t,begin,mid);
        n->right=lya_forest_kd_build_node(kd,t,mid,end);
    }
    return id;
}

static int lya_forest_kd_build(lya_forest_kd *kd,const lya_cells *t,
                                struct cmdline_data *cmd,size_t plan,ErrorMsg error)
{
    size_t leaves=1,target,bytes,extra,scratch;
    memset(kd,0,sizeof(*kd));kd->plan=plan;
    /* Match the existing exact interval-arithmetic coverage. Other profiles
     * retain the original forest-row traversal, without speculative pruning. */
    if(t->roots_count<LYA_FOREST_KD_MIN || sizeof(REAL)!=sizeof(double)
        || sizeof(real)!=sizeof(double) || sizeof(cballs_storage_real)!=sizeof(double)) return SUCCESS;
#ifdef __FAST_MATH__
    return SUCCESS;
#endif
    target=t->roots_count/LYA_FOREST_KD_LEAF+(t->roots_count%LYA_FOREST_KD_LEAF!=0);
    while(leaves<target) {if(leaves>SIZE_MAX/2) goto overflow;leaves*=2;}
    if(!cballs_size_add(leaves,leaves-1,&kd->capacity)) goto overflow;
    kd->words=t->roots_count/64+(t->roots_count%64!=0);
    kd->summaries=kd->words/64+(kd->words%64!=0);
    if(!cballs_size_mul(kd->capacity,sizeof(*kd->nodes),&bytes)
        || !cballs_size_mul(t->roots_count,sizeof(*kd->rows),&extra)
        || !cballs_size_add(bytes,extra,&bytes)
        || !cballs_size_add(kd->words,kd->summaries,&scratch)
        || !cballs_size_mul(scratch,sizeof(uint64_t),&scratch)
        || !cballs_size_mul(scratch,(size_t)MAX(1,cmd->numthreads),&scratch)
        || !cballs_size_add(bytes,scratch,&bytes)
        || !cballs_size_add(plan,bytes,&kd->plan)) goto overflow;
    if(cballs_memory_preflight(kd->plan,"Ly-alpha forest-bounds index and pair workspace",error,_ERRORMSGSIZE_)==FAILURE)
        return FAILURE;
    if(cballs_calloc_checked((void**)&kd->nodes,kd->capacity,sizeof(*kd->nodes),
            "Ly-alpha forest-bounds nodes",error,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&kd->rows,t->roots_count,sizeof(*kd->rows),
            "Ly-alpha forest-bounds rows",error,_ERRORMSGSIZE_)==FAILURE) {
        lya_forest_kd_free(kd);return FAILURE;
    }
    for(size_t i=0;i<t->roots_count;i++) kd->rows[i]=i;
    lya_forest_kd_build_node(kd,t,0,t->roots_count);
    kd->enabled=1;
    return SUCCESS;
overflow:
    snprintf(error,_ERRORMSGSIZE_,"Ly-alpha forest-bounds index dimensions overflow");
    return FAILURE;
}

static int lya_forest_kd_bits_init(lya_forest_kd_bits *w,const lya_forest_kd *kd,ErrorMsg error)
{
    memset(w,0,sizeof(*w));
    if(!kd->enabled) return SUCCESS;
    if(cballs_calloc_checked((void**)&w->bits,kd->words,sizeof(*w->bits),
            "Ly-alpha forest candidate bits",error,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&w->summary,kd->summaries,sizeof(*w->summary),
            "Ly-alpha forest candidate summary",error,_ERRORMSGSIZE_)==FAILURE) {
        free(w->bits);free(w->summary);memset(w,0,sizeof(*w));return FAILURE;
    }
    return SUCCESS;
}

static unsigned lya_forest_kd_first_bit(uint64_t bits)
{
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_ctzll((unsigned long long)bits);
#else
    unsigned n=0;while(!(bits&1)) {bits>>=1;n++;}return n;
#endif
}

/* Conservative rejection only. Bounds have the same absolute-coordinate and
 * observer-angle padding as the existing pair traversal. Whole long forests,
 * bent rays and pixels far from their forest centroid are all enclosed. */
static int lya_forest_kd_outside(const lya_cell *p,const lya_cell *q,
                                  struct cmdline_data *cmd,REAL cutoff)
{
    REAL upper;
    if(lya_pair_box_outside(p,q,cutoff,&upper)) return 1;
    lya_pair_angle a=lya_pair_angles(p,q);
    REAL gap=MAX(0.,MAX(p->chilo-q->chihi,q->chilo-p->chihi));
    REAL sumlo=p->chilo+q->chilo,sumhi=p->chihi+q->chihi;
    REAL pad=4096*DBL_EPSILON*sumhi+16*sqrt(DBL_MIN);
    REAL rp=MAX(0.,gap*a.cmin-pad),rt=MAX(0.,sumlo*a.smin-pad);
    return isfinite(rp) && isfinite(rt) && isfinite(pad)
        && (rp>=cmd->lya2RpMax || rt>=cmd->lya2RtMax);
}

static void lya_forest_kd_query(const lya_forest_kd *kd,const lya_cells *t,
                                 const lya_cell *p,size_t row,size_t id,
                                 struct cmdline_data *cmd,REAL cutoff,
                                 lya_forest_kd_bits *w,lya_pair_work *work)
{
    const lya_forest_kd_node *n=kd->nodes+id;
    if(n->last_row<=row) return;
    work->forest_nodes++;
    if(lya_forest_kd_outside(p,&n->bounds,cmd,cutoff)) {work->forest_pruned++;return;}
    if(n->left!=SIZE_MAX) {
        lya_forest_kd_query(kd,t,p,row,n->left,cmd,cutoff,w,work);
        lya_forest_kd_query(kd,t,p,row,n->right,cmd,cutoff,w,work);
    } else for(size_t j=n->begin;j<n->end;j++) {
        size_t candidate=kd->rows[j];
        if(candidate<=row) continue;
        size_t word=candidate/64;
        w->bits[word]|=UINT64_C(1)<<(candidate%64);
        w->summary[word/64]|=UINT64_C(1)<<(word%64);
        work->forest_candidates++;
    }
}
#endif
