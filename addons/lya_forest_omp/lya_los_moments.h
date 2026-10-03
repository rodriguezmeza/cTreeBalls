/* A second level above the per-forest segments. It reuses radial/polar
 * products across distinct forests without subtracting same-forest totals.
 * Used only by the zero-slop LOS kernel-5 fallback; leaves retain the exact
 * segment/pixel implementation. Bounds are independent of signed delta. */
static int lya_los_moment_compare(const void *aa,const void *bb)
{
    const lya_los_moment_entry *a=aa,*b=bb;
    if(a->bin!=b->bin) return a->bin<b->bin?-1:1;
    if(a->angle!=b->angle) return a->angle<b->angle?-1:1;
    return a->source<b->source?-1:a->source>b->source;
}

static size_t lya_los_moment_build(lya_worker_hist *w,size_t begin,size_t end)
{
    size_t id=w->los_moment_count++;
    lya_los_moment *s=&w->los_moments[id];
    s->left=s->right=SIZE_MAX;
    if(end-begin==1) {
        s->source=w->los_entries[begin].source;
        s->bounds=w->segments[s->source];
        const lya_neighbor *q=&w->neighbors[s->bounds.begin];
        s->count=s->bounds.end-s->bounds.begin;
        s->bin=q->leg_bin;s->first=q->first_index;s->second=q->second_index;
        s->forest_lo=s->forest_hi=q->forest_id;
        uint64_t key=(uint64_t)q->forest_id*UINT64_C(0x9e3779b97f4a7c15);key^=key>>33;
        unsigned bit=(unsigned)(key&255);
        memset(s->forests,0,sizeof(s->forests));s->forests[bit/64]=UINT64_C(1)<<(bit%64);
    } else {
        size_t mid=begin+(end-begin)/2;
        s->left=lya_los_moment_build(w,begin,mid);s->right=lya_los_moment_build(w,mid,end);
        const lya_los_moment *a=&w->los_moments[s->left],*b=&w->los_moments[s->right];
        s->source=SIZE_MAX;s->count=a->count+b->count;
        s->bin=a->bin==b->bin?a->bin:SIZE_MAX;s->first=a->first;s->second=a->second;
        s->forest_lo=MIN(a->forest_lo,b->forest_lo);s->forest_hi=MAX(a->forest_hi,b->forest_hi);
        s->bounds.weight=a->bounds.weight+b->bounds.weight;
        s->bounds.weighted_delta=a->bounds.weighted_delta+b->bounds.weighted_delta;
        for(int k=0;k<3;k++) {s->bounds.lo[k]=MIN(a->bounds.lo[k],b->bounds.lo[k]);s->bounds.hi[k]=MAX(a->bounds.hi[k],b->bounds.hi[k]);}
        for(int k=0;k<4;k++) s->forests[k]=a->forests[k]|b->forests[k];
    }
    return id;
}

static void lya_los_moment_product(struct cmdline_data *cmd,bodyptr p,lya_worker_hist *w,size_t ia,size_t ib)
{
    const lya_los_moment *a=&w->los_moments[ia],*b=&w->los_moments[ib];
    w->los_moment_counts[1]++;
    if(a->forest_lo==a->forest_hi && b->forest_lo==a->forest_lo && b->forest_hi==a->forest_lo) return;
    if(ia==ib) {
        if(a->left==SIZE_MAX) return;
        lya_los_moment_product(cmd,p,w,a->left,a->left);
        lya_los_moment_product(cmd,p,w,a->left,a->right);
        lya_los_moment_product(cmd,p,w,a->right,a->right);return;
    }
    uint64_t common=0;for(int k=0;k<4;k++) common|=a->forests[k]&b->forests[k];
    if(w->aggregation_safe && a->bin!=SIZE_MAX && b->bin!=SIZE_MAX
        && (!common || a->forest_hi<b->forest_lo || b->forest_hi<a->forest_lo)) {
        int unused=0,mu=lya_segment_mu(w,&a->bounds,&b->bounds,cmd->lya3MuBins,0,&unused);
        if(mu>=0) {
            long double d=(long double)lya_pivot_weight(w,p)*a->bounds.weight*b->bounds.weight;
            long double n=(long double)lya_pivot_field(w,p)*a->bounds.weighted_delta*b->bounds.weighted_delta;
            if(isfinite((REAL)d) && isfinite((REAL)n) && d>=DBL_MIN) {
                lya_deposit3(w,a->first+b->second+mu,b->first+a->second+mu,(REAL)n,(REAL)d);
                size_t count=a->count*b->count;
                w->ordered_triplet_count+=(INTEGER)(2*count);
                w->aggregated_pairs+=count;w->segment_accepts++;
                w->los_moment_counts[2]++;w->los_moment_counts[3]+=count;return;
            }
        }
    }
    if(a->left==SIZE_MAX && b->left==SIZE_MAX) {
        lya_segment_pair(cmd,p,w,a->source,b->source);return;
    }
    int split_a=a->left!=SIZE_MAX && (b->left==SIZE_MAX
        || (a->bin==SIZE_MAX && b->bin!=SIZE_MAX)
        || ((a->bin==SIZE_MAX)==(b->bin==SIZE_MAX) && a->count>=b->count));
    if(split_a) {
        lya_los_moment_product(cmd,p,w,a->left,ib);lya_los_moment_product(cmd,p,w,a->right,ib);
    } else {
        lya_los_moment_product(cmd,p,w,ia,b->left);lya_los_moment_product(cmd,p,w,ia,b->right);
    }
}

static int lya_los_moments_run(struct cmdline_data *cmd,bodyptr p,lya_worker_hist *w,size_t roots,ErrorMsg err)
{
    if(roots>w->los_moment_capacity) {
        size_t bytes,segments,total,per_worker;
        if(!cballs_size_mul(roots,2*sizeof(lya_los_moment)+sizeof(lya_los_moment_entry),&bytes)
            || !cballs_size_mul(w->segment_capacity,2*sizeof(lya_segment)+sizeof(size_t),&segments)
            || !cballs_size_mul(w->neighbor_capacity,sizeof(lya_neighbor),&per_worker)
            || !cballs_size_add(per_worker,segments,&per_worker)
            || !cballs_size_add(per_worker,bytes,&per_worker)
            || !cballs_size_add(per_worker,w->los_moment_bytes,&per_worker)
            || !cballs_size_mul(per_worker,w->scratch_workers,&total)
            || !cballs_size_add(total,w->histogram_plan_bytes,&total)) {
            snprintf(err,_ERRORMSGSIZE_,"Ly-alpha LOS radial moments scratch overflow");return FAILURE;
        }
        if(cballs_memory_preflight(total,"Ly-alpha LOS radial moments",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
        lya_los_moment *next=realloc(w->los_moments,bytes);
        if(!next) {snprintf(err,_ERRORMSGSIZE_,"Ly-alpha LOS radial moments allocation failed");return FAILURE;}
        w->los_moments=next;w->los_moment_bytes=bytes;w->los_moment_capacity=roots;
        w->los_entries=(lya_los_moment_entry*)(next+2*roots);
    }
    /* Polar bins are about the pivot LOS, so sort azimuth in that same
     * frame. World-axis azimuth can interleave opposite sides of a polar
     * ring and destroy the useful direction bounds. This basis controls
     * ordering only; all acceptance bounds still use stored displacements. */
    REAL u[3]={0,0,0},v[3],norm=0;
    int axis=0;for(int k=1;k<3;k++) if(fabs(LyaLOS(p)[k])<fabs(LyaLOS(p)[axis])) axis=k;
    int a=(axis+1)%3,b=(axis+2)%3;
    u[a]=LyaLOS(p)[b];u[b]=-LyaLOS(p)[a];
    for(int k=0;k<3;k++) norm+=u[k]*u[k];
    norm=sqrt(norm);for(int k=0;k<3;k++) u[k]/=norm;
    for(int k=0;k<3;k++) v[k]=LyaLOS(p)[(k+1)%3]*u[(k+2)%3]-LyaLOS(p)[(k+2)%3]*u[(k+1)%3];
    for(size_t j=0;j<roots;j++) {
        size_t id=w->segment_roots[j];const lya_segment *s=&w->segments[id];
        const lya_neighbor *q=&w->neighbors[s->begin];
        REAL x=0,y=0;
        for(int k=0;k<3;k++) {x+=q->displacement[k]*u[k];y+=q->displacement[k]*v[k];}
        w->los_entries[j]=(lya_los_moment_entry){id,q->leg_bin,atan2(y,x)};
    }
    qsort(w->los_entries,roots,sizeof(*w->los_entries),lya_los_moment_compare);
    w->los_moment_count=0;lya_los_moment_build(w,0,roots);w->los_moment_counts[0]+=w->los_moment_count;
    lya_los_moment_product(cmd,p,w,0,0);return SUCCESS;
}
