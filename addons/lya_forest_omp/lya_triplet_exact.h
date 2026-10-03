/* Included by search_lya_forest_omp.c after the independent direct kernel.
 * Zero slop requires every opening cosine to have the same bin. Positive
 * mu slop explicitly permits bounded bin leakage. Sorting changes summation order. */
#include "lya_neighbor_sort.h"
static int lya_neighbor_compare(const void *aa, const void *bb)
{
    const lya_neighbor *a=aa, *b=bb;
    if (a->forest_id != b->forest_id) return a->forest_id < b->forest_id ? -1 : 1;
    if (a->leg_bin != b->leg_bin) return a->leg_bin < b->leg_bin ? -1 : 1;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    return a->ordinal < b->ordinal ? -1 : (a->ordinal > b->ordinal);
}

static void lya_deposit3(lya_worker_hist *w, size_t a, size_t b, REAL n, REAL d)
{
    if (d <= 0) return;
    if (w->den3[a] == 0) w->touched3[w->touched3_count++]=a;
    w->num3[a]+=n; w->den3[a]+=d;
    if (w->den3[b] == 0) w->touched3[w->touched3_count++]=b;
    w->num3[b]+=n; w->den3[b]+=d;
}

static size_t lya_segment_build(lya_worker_hist *w, size_t begin, size_t end)
{
    size_t index=w->segment_count++;
    lya_segment *s=&w->segments[index];
    s->begin=begin; s->end=end; s->left=s->right=SIZE_MAX;
    if (end-begin > 8) {
        size_t middle=begin+(end-begin)/2;
        s->left=lya_segment_build(w,begin,middle);
        s->right=lya_segment_build(w,middle,end);
        const lya_segment *a=&w->segments[s->left],*b=&w->segments[s->right];
        s->weight=a->weight+b->weight;s->weighted_delta=a->weighted_delta+b->weighted_delta;
        for(int k=0;k<3;k++) {s->lo[k]=MIN(a->lo[k],b->lo[k]);s->hi[k]=MAX(a->hi[k],b->hi[k]);}
    } else {
    s->weight=s->weighted_delta=0;
    for (int k=0;k<3;k++) {s->lo[k]=1; s->hi[k]=-1;}
    for (size_t j=begin;j<end;j++) {
        const lya_neighbor *q=&w->neighbors[j];
        s->weight+=q->weight; s->weighted_delta+=q->weighted_delta;
        for (int k=0;k<3;k++) {
            REAL u=q->displacement[k]/q->radius;
            /* Disable interval acceptance for unsafe exponent ranges. The
             * reference dot/(r1*r2) can underflow/overflow there. */
            if (!isfinite(u) || q->radius < 16*sqrt(DBL_MIN)
                || q->radius > sqrt(DBL_MAX)/16) {
                s->lo[k]=-INFINITY; s->hi[k]=INFINITY;
            } else {s->lo[k]=MIN(s->lo[k],u); s->hi[k]=MAX(s->hi[k],u);}
        }
    }
    }
    return index;
}

static int lya_segment_mu(const lya_worker_hist *w, const lya_segment *a,
                           const lya_segment *b, int bins, REAL slop, int *approximate)
{
#ifdef __FAST_MATH__
    return -1;
#endif
    if (sizeof(REAL)!=sizeof(double)) return -1;
    REAL lo=0,hi=0;
    for (int k=0;k<3;k++) {
        REAL x=a->lo[k]*b->lo[k], y=a->lo[k]*b->hi[k];
        REAL z=a->hi[k]*b->lo[k], t=a->hi[k]*b->hi[k];
        lo+=MIN(MIN(x,y),MIN(z,t)); hi+=MAX(MAX(x,y),MAX(z,t));
    }
    if (!isfinite(lo) || !isfinite(hi)) return -1;
    /* Unit components are bounded by 1+roundoff. 512 epsilon exceeds the
     * absolute error of normalization, six endpoint products, three sums,
     * and the reference dot/(radius product), including fused operations.
     * Edge-straddling intervals always descend; do not round toward centers. */
    int low=lya_bin_mu(lo-512*DBL_EPSILON,bins);
    int high=lya_bin_mu(hi+512*DBL_EPSILON,bins);
    *approximate=0;
    if (low==high) return low;
    if (slop <= 0) return -1;
    /* Geometric midpoint of segment endpoints; never use signed field weights
     * to locate a representative. The enclosing interval remains authoritative. */
    REAL ca[3],cb[3],aa=0,bb=0,dot=0;
    for (int k=0;k<3;k++) {
        ca[k]=.5*w->neighbors[a->begin].displacement[k]+.5*w->neighbors[a->end-1].displacement[k];
        cb[k]=.5*w->neighbors[b->begin].displacement[k]+.5*w->neighbors[b->end-1].displacement[k];
        aa+=ca[k]*ca[k];bb+=cb[k]*cb[k];dot+=ca[k]*cb[k];
    }
    if (!(aa>0) || !(bb>0) || !isfinite(aa) || !isfinite(bb)) return -1;
    REAL mu=dot/(sqrt(aa)*sqrt(bb));
    if (!isfinite(mu)) return -1;
    int bin=lya_bin_mu(mu,bins);
    REAL width=2.0/bins, lower=-1+bin*width, upper=lower+width;
    lo=MAX(-1.0,lo-512*DBL_EPSILON);hi=MIN(1.0,hi+512*DBL_EPSILON);
    if (lo < lower-slop*width || hi > upper+slop*width) return -1;
    *approximate=1;
    return bin;
}

static void lya_direct_ranges(struct cmdline_data *cmd, bodyptr p,
                              lya_worker_hist *w, size_t ab, size_t ae,
                              size_t bb, size_t be)
{
    for (size_t i=ab;i<ae;i++) {
        const lya_neighbor *q=&w->neighbors[i];
        REAL pw=lya_pivot_weight(w,p)*q->weight, pn=lya_pivot_field(w,p)*q->weighted_delta;
        /* Separate vectorizable geometry from dependent histogram scatters.
         * Preserve the reference multiply/add/divide order at mu boundaries. */
        for (size_t base=bb;base<be;base+=32) {
            size_t n=MIN((size_t)32,be-base);
            REAL cosine[32];
#pragma omp simd
            for (size_t k=0;k<n;k++) {
                const lya_neighbor *r=&w->neighbors[base+k];
                REAL dot; DOTVP(dot,q->displacement,r->displacement);
                cosine[k]=dot/(q->radius*r->radius);
            }
            for (size_t k=0;k<n;k++) {
                const lya_neighbor *r=&w->neighbors[base+k];
                if (q->forest_id==r->forest_id) continue;
                int bm=lya_bin_mu(cosine[k],cmd->lya3MuBins);
                w->ordered_triplet_count+=2; w->direct_pairs++;
                lya_deposit3(w,q->first_index+r->second_index+bm,
                             r->first_index+q->second_index+bm,
                             pn*r->weighted_delta,pw*r->weight);
            }
        }
    }
}

/* Both segments have one radial/polar bin and distinct forests. Batch the
 * remaining exact mu decisions into a small local histogram, then publish
 * each occupied mu bin once. Large mu grids use the generic tiled fallback. */
static void lya_direct_segment_ranges(struct cmdline_data *cmd, bodyptr p,
                                      lya_worker_hist *w, const lya_segment *a,
                                      const lya_segment *b)
{
    int bins=cmd->lya3MuBins;
    if (bins>64 || (a->end-a->begin)*(b->end-b->begin)<8) {lya_direct_ranges(cmd,p,w,a->begin,a->end,b->begin,b->end);return;}
    REAL numerator[64],denominator[64];
    memset(numerator,0,(size_t)bins*sizeof(REAL));
    memset(denominator,0,(size_t)bins*sizeof(REAL));
    const lya_neighbor *qa=&w->neighbors[a->begin], *rb=&w->neighbors[b->begin];
    for (size_t i=a->begin;i<a->end;i++) {
        const lya_neighbor *q=&w->neighbors[i];
        REAL pw=lya_pivot_weight(w,p)*q->weight,pn=lya_pivot_field(w,p)*q->weighted_delta;
        REAL cosine[8];
        size_t count=b->end-b->begin; /* binary segment leaves have <=8 pixels */
#pragma omp simd
        for (size_t k=0;k<count;k++) {
            const lya_neighbor *r=&w->neighbors[b->begin+k];
            REAL dot;DOTVP(dot,q->displacement,r->displacement);
            cosine[k]=dot/(q->radius*r->radius);
        }
        for (size_t k=0;k<count;k++) {
            const lya_neighbor *r=&w->neighbors[b->begin+k];
            REAL d=pw*r->weight;
            if (d<=0) continue;
            int bin=lya_bin_mu(cosine[k],bins);
            numerator[bin]+=pn*r->weighted_delta;denominator[bin]+=d;
        }
    }
    size_t count=(a->end-a->begin)*(b->end-b->begin);
    w->ordered_triplet_count+=(INTEGER)(2*count);w->direct_pairs+=count;
    for (int m=0;m<bins;m++) if (denominator[m]>0)
        lya_deposit3(w,qa->first_index+rb->second_index+m,
                     rb->first_index+qa->second_index+m,numerator[m],denominator[m]);
}

static void lya_segment_pair(struct cmdline_data *cmd, bodyptr p,
                             lya_worker_hist *w, size_t ia, size_t ib)
{
    const lya_segment *a=&w->segments[ia], *b=&w->segments[ib];
    size_t na=a->end-a->begin, nb=b->end-b->begin;
    if ((na > 1 || nb > 1) && w->aggregation_safe) {
        int approximate=0;
        int bm=lya_segment_mu(w,a,b,cmd->lya3MuBins,cmd->lya3MuSlop,&approximate);
        if (bm >= 0) {
            long double d=(long double)lya_pivot_weight(w,p)*a->weight*b->weight;
            long double n=(w->pivot_multiplicity?(long double)w->pivot_field:(long double)Weight(p)*Kappa(p))*a->weighted_delta*b->weighted_delta;
            /* Preserve leaf behavior for underflow or overflowing products. */
            if (isfinite((double)d) && isfinite((double)n) && d >= DBL_MIN) {
                const lya_neighbor *q=&w->neighbors[a->begin], *r=&w->neighbors[b->begin];
                lya_deposit3(w,q->first_index+r->second_index+bm,
                             r->first_index+q->second_index+bm,(REAL)n,(REAL)d);
                w->ordered_triplet_count+=(INTEGER)(2*na*nb);
                w->aggregated_pairs+=(unsigned long long)na*nb;
                if (approximate) w->approximate_pairs+=(unsigned long long)na*nb;
                w->segment_accepts++;
                return;
            }
        }
    }
    if (a->left==SIZE_MAX && b->left==SIZE_MAX) {
        lya_direct_segment_ranges(cmd,p,w,a,b);
    } else if (a->left!=SIZE_MAX && (b->left==SIZE_MAX || na>=nb)) {
        lya_segment_pair(cmd,p,w,a->left,ib); lya_segment_pair(cmd,p,w,a->right,ib);
    } else {
        lya_segment_pair(cmd,p,w,ia,b->left); lya_segment_pair(cmd,p,w,ia,b->right);
    }
}

#include "lya_los_moments.h"

static int lya_accumulate_segments(struct cmdline_data *cmd, bodyptr p,
                                   lya_worker_hist *w, ErrorMsg err)
{
    size_t count=w->neighbor_count, bytes, scratch, total, roots=0, forests=0;
    if (count<2) return SUCCESS;
    if (cmd->lya3Kernel==2) {
        for (size_t i=0;i<count;i++) lya_direct_ranges(cmd,p,w,i,i+1,i+1,count);
        return SUCCESS;
    }
    if (count > w->segment_capacity) {
        /* Include neighbor capacity and all workers in the scratch preflight.
         * realloc can transiently hold both old and new storage: reserve 2x. */
        if (!cballs_size_mul(w->neighbor_capacity,2*sizeof(lya_segment)+sizeof(size_t),&bytes)
            || !cballs_size_mul(w->neighbor_capacity,sizeof(lya_neighbor),&scratch)
            || !cballs_size_mul(bytes,2,&bytes)
            || !cballs_size_add(bytes,scratch,&total)
            || !cballs_size_add(total,w->los_moment_bytes,&total)
            || !cballs_size_mul(total,w->scratch_workers,&total)
            || !cballs_size_add(total,w->histogram_plan_bytes,&total)) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha segment scratch dimensions overflow"); return FAILURE;
        }
        if (cballs_memory_preflight(total,"Ly-alpha segment scratch",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
        lya_segment *nodes=realloc(w->segments,2*w->neighbor_capacity*sizeof(*nodes));
        if (!nodes) {snprintf(err,_ERRORMSGSIZE_,"Ly-alpha segment allocation failed");return FAILURE;}
        w->segments=nodes;
        size_t *r=realloc(w->segment_roots,w->neighbor_capacity*sizeof(*r));
        if (!r) {snprintf(err,_ERRORMSGSIZE_,"Ly-alpha segment root allocation failed");return FAILURE;}
        w->segment_roots=r; w->segment_capacity=w->neighbor_capacity;
    }
    REAL pivot_weight=lya_pivot_weight(w,p),pivot_field=fabs(lya_pivot_field(w,p));
    w->aggregation_safe=(pivot_weight==0 || (pivot_weight>=1e-50 && pivot_weight<=1e50))
        && (pivot_field==0 || (pivot_field>=1e-50 && pivot_field<=1e50));
    for (size_t i=0;i<count;i++) {
        REAL d=fabs(w->neighbors[i].weighted_delta), v=w->neighbors[i].weight;
        if ((v!=0 && (v<1e-50 || v>1e50)) || (d!=0 && (d<1e-50 || d>1e50))) w->aggregation_safe=0;
    }
    if(!w->neighbors_sorted) lya_sort_neighbors(w->neighbors,count);
    w->segment_count=0;
    for (size_t begin=0;begin<count;) {
        if(begin==0 || w->neighbors[begin].forest_id!=w->neighbors[begin-1].forest_id) forests++;
        size_t end=begin+1;
        while (end<count && w->neighbors[end].forest_id==w->neighbors[begin].forest_id
               && w->neighbors[end].first_index==w->neighbors[begin].first_index) end++;
        w->segment_roots[roots++]=lya_segment_build(w,begin,end); begin=end;
    }
    /* Repeated forest membership across many leg bins forces the mixed
     * hierarchy to refine almost every product. Keep the compact per-forest
     * loop when the frontier averages more than three roots per forest.
     * This work estimate depends only on geometry, never thread/rank count. */
    if(cmd->lya3Kernel==5 && lya_forest_is_los_tree_method(cmd->searchMethod)
        && roots>=32 && (roots-1)/forests<3 && w->aggregation_safe && cmd->lya3MuSlop==0)
        return lya_los_moments_run(cmd,p,w,roots,err);
    for (size_t i=0;i<roots;i++) for (size_t j=i+1;j<roots;j++) {
        size_t a=w->segment_roots[i], b=w->segment_roots[j];
        if (w->neighbors[w->segments[a].begin].forest_id
            ==w->neighbors[w->segments[b].begin].forest_id) continue;
        lya_segment_pair(cmd,p,w,a,b);
    }
    return SUCCESS;
}
