/* A genuine 2D quadtree partition of the requested (rp,rt) BIN PLANE.
 * Pixel coordinates in this plane depend on BOTH endpoints. We therefore
 * project immutable forest-segment products into conservative rectangles,
 * without ever materializing an N-pixel squared pair catalog. Unresolved
 * products refine their radial hierarchies and reuse exact geometric certificates.
 * Default binning is exact. Optional per-axis slop permits a whole product
 * only inside a certified expanded bin; forest exclusions and outer cuts
 * remain exact. Signed-field sums never define geometric representatives.
 */
#ifndef LYA_RPRT_QUADTREE_H
#define LYA_RPRT_QUADTREE_H

typedef struct {
    int p0,p1,t0,t1; /* half-open bin INDEX ranges, so arbitrary grid sizes work */
    REAL pcut,tcut;
    size_t children[4];
} lya_rprt_quad_node;

typedef struct { REAL chi,weight,field; } lya_rprt_pixel;

typedef struct {
    lya_rprt_quad_node *nodes;
    lya_rprt_pixel *pixels;
    size_t *bin_leaves; /* optional O(1) leaf directory for expanded-bin queries */
    size_t count,capacity,plan;
    REAL pscale,tscale,pwidth,twidth,popen,topen;
    double pixel_wall,bin_wall;
    int pixel_threads,bin_threads;
} lya_rprt_quad;

/* Exact compact subtree size, in logarithmic time. At depth k there are
 * min(2^k,np)*min(2^k,nt) retained bin rectangles. Subtract the singleton
 * rectangles that already terminated at earlier levels. This reserves the
 * SAME compact preorder slots as the serial builder, without shared writes. */
static size_t lya_rprt_quad_node_count(size_t np,size_t nt)
{
    size_t count=1,scale=1,old_leaves=(np==1 && nt==1);
    while(scale<np || scale<nt) {
        scale*=2;
        count+=MIN(scale,np)*MIN(scale,nt)-old_leaves;
        size_t sp=scale>=np?np:(2*scale>np?2*scale-np:0);
        size_t st=scale>=nt?nt:(2*scale>nt?2*scale-nt:0);
        old_leaves=sp*st;
    }
    return count;
}

static size_t lya_rprt_quad_node_build(lya_rprt_quad *q,struct cmdline_data *cmd,size_t id,
                                         int p0,int p1,int t0,int t1,int tasks)
{
    lya_rprt_quad_node *n=q->nodes+id;
    n->p0=p0;n->p1=p1;n->t0=t0;n->t1=t1;
    for(int j=0;j<4;j++) n->children[j]=SIZE_MAX;
    if(p1-p0==1 && t1-t0==1) {
        if(q->bin_leaves) q->bin_leaves[(size_t)p0*cmd->lya2RtBins+t0]=id;
        return 1;
    }
    int pm=p0+(p1-p0)/2,tm=t0+(t1-t0)/2;
    n->pcut=cmd->lya2RpMax*((REAL)pm/cmd->lya2RpBins);
    n->tcut=cmd->lya2RtMax*((REAL)tm/cmd->lya2RtBins);
    int pe[3]={p0,pm,p1},te[3]={t0,tm,t1};
    /* A singleton dimension has one child interval, not two empty copies. */
    int np=p1-p0>1?2:1,nt=t1-t0>1?2:1;
    if(np==1) pe[1]=p1;
    if(nt==1) te[1]=t1;
    size_t used[4]={0},next=id+1;
    for(int i=0;i<np;i++) for(int j=0;j<nt;j++) {
            int slot=2*i+j,pa=pe[i],pb=pe[i+1],ta=te[j],tb=te[j+1];
            size_t child=next,leaves=(size_t)(pb-pa)*(tb-ta);
            n->children[slot]=child;next+=lya_rprt_quad_node_count(pb-pa,tb-ta);
            if(tasks && leaves>=1024) {
#pragma omp task shared(used) firstprivate(q,cmd,child,pa,pb,ta,tb,slot,tasks)
                used[slot]=lya_rprt_quad_node_build(q,cmd,child,pa,pb,ta,tb,tasks);
            } else used[slot]=lya_rprt_quad_node_build(q,cmd,child,pa,pb,ta,tb,0);
    }
    if(tasks) {
#pragma omp taskwait
    }
    return 1+used[0]+used[1]+used[2]+used[3];
}

static int lya_rprt_quad_build(lya_rprt_quad *q,const lya_cells *t,struct cmdline_data *cmd,
                                 size_t bins,size_t plan,ErrorMsg error)
{
    size_t twice,bytes,pixel_bytes,leaf_bytes;
    size_t leaf_count=(cmd->lya2RpSlop>0 || cmd->lya2RtSlop>0)?bins:0;
    memset(q,0,sizeof(*q));q->plan=plan;
    q->pscale=cmd->lya2RpBins/cmd->lya2RpMax;q->tscale=cmd->lya2RtBins/cmd->lya2RtMax;
    q->pwidth=cmd->lya2RpMax/cmd->lya2RpBins;q->twidth=cmd->lya2RtMax/cmd->lya2RtBins;
    q->popen=q->pscale/(1+2*cmd->lya2RpSlop);q->topen=q->tscale/(1+2*cmd->lya2RtSlop);
    /* Each internal node has at least two children: at most 2*leaves-1.
     * This includes non-square/non-power-of-two and singleton grids. */
    if(!bins || !cballs_size_mul(bins,2,&twice)
        || !cballs_size_mul(twice-1,sizeof(*q->nodes),&bytes)
        || !cballs_size_mul(t->count,sizeof(*q->pixels),&pixel_bytes)
        || !cballs_size_mul(leaf_count,sizeof(*q->bin_leaves),&leaf_bytes)
        || !cballs_size_add(plan,bytes,&q->plan)
        || !cballs_size_add(q->plan,pixel_bytes,&q->plan)
        || !cballs_size_add(q->plan,leaf_bytes,&q->plan)) {
        snprintf(error,_ERRORMSGSIZE_,"Ly-alpha rp-rt quadtree dimensions overflow");return FAILURE;
    }
    q->capacity=twice-1;
    if(cballs_memory_preflight(q->plan,"Ly-alpha rp-rt bin quadtree and pair workspace",error,_ERRORMSGSIZE_)==FAILURE)
        return FAILURE;
    if(cballs_calloc_checked((void**)&q->nodes,q->capacity,sizeof(*q->nodes),
            "Ly-alpha rp-rt bin quadtree",error,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
    if(cballs_calloc_checked((void**)&q->pixels,t->count,sizeof(*q->pixels),
            "Ly-alpha rp-rt compact pixels",error,_ERRORMSGSIZE_)==FAILURE) {
        free(q->nodes);q->nodes=NULL;return FAILURE;
    }
    if(leaf_count && cballs_calloc_checked((void**)&q->bin_leaves,leaf_count,sizeof(*q->bin_leaves),
            "Ly-alpha rp-rt slop leaf directory",error,_ERRORMSGSIZE_)==FAILURE) {
        free(q->pixels);free(q->nodes);q->pixels=NULL;q->nodes=NULL;return FAILURE;
    }
    int threads=scanopt(cmd->options,"lya-serial-build")?1:MAX(1,cmd->numthreads);
    q->pixel_threads=q->bin_threads=1;
    double start=lya_progress_wall_time();
#pragma omp parallel num_threads(threads) if(threads>1 && t->count>=32768)
    {
#pragma omp single
        q->pixel_threads=lya_build_team_size();
#pragma omp for schedule(static)
        for(size_t i=0;i<t->count;i++) {
            bodyptr p=t->points[i];
            q->pixels[i]=(lya_rprt_pixel){LyaDistance(p),Weight(p),Weight(p)*Kappa(p)};
        }
    }
    q->pixel_wall=lya_progress_wall_time()-start;start=lya_progress_wall_time();
    /* A normal 50x50 bin plane is too small to amortize task/team overhead.
     * Large grids initialize their physical Rp/Rt cuts concurrently. */
#pragma omp parallel num_threads(threads) if(threads>1 && bins>=16384)
    {
#pragma omp single
        {
            q->bin_threads=lya_build_team_size();
            q->count=lya_rprt_quad_node_build(q,cmd,0,0,cmd->lya2RpBins,0,cmd->lya2RtBins,q->bin_threads>1);
        }
    }
    q->bin_wall=lya_progress_wall_time()-start;
    return SUCCESS;
}

/* Find the smallest quadtree rectangle containing the ENTIRE projected bin
 * range. A straddling rectangle remains at the parent; it is never assigned
 * to its centre's quadrant. The returned node is immutable and shared. */
static size_t lya_rprt_quad_locate(const lya_rprt_quad *q,size_t id,REAL plo,REAL phi,REAL tlo,REAL thi,
                                     lya_pair_work *work)
{
    for(;;) {
        work->quad_visits++;
        const lya_rprt_quad_node *n=q->nodes+id;
        int p=0,t=0;
        if(n->p1-n->p0>1) {
            if(phi<n->pcut) p=0;
            else if(plo>=n->pcut) p=1;
            else return id; /* the projected cell straddles a quadrant edge */
        }
        if(n->t1-n->t0>1) {
            if(thi<n->tcut) t=0;
            else if(tlo>=n->tcut) t=1;
            else return id;
        }
        size_t next=n->children[2*p+t];
        if(next==SIZE_MAX) return id;
        id=next;work->quad_descents++;
    }
}

/* Positive-slop queries know the proposed bin from the enclosure midpoint.
 * Compare physical interval ends directly against that expanded leaf. Avoid
 * extra endpoint bin divisions on this explicitly approximate path; an axis
 * with zero slop still uses reference endpoint bin agreement. */
static int lya_rprt_slop_bin(REAL lo,REAL hi,REAL center,REAL maximum,int bins,REAL slop)
{
    if(slop==0) {
        int unused=0;return lya_cell_bin(lo,hi,center,maximum,bins,0.,0,&unused);
    }
    int bin=lya_bin_positive(center,maximum,bins);
    if(bin<0) return -1;
    REAL width=maximum/bins,lower=bin*width;
    return lo>=lower-slop*width && hi<=lower+(1+slop)*width?bin:-1;
}

/* Stream the hot fields in radial order. The common case has certified
 * ownership and discovery-sphere containment. An interval certificate still
 * decides each rp/rt bin; edge ambiguity uses the original pixel arithmetic.
 * No LOS direction, coordinate or signed-field centroid is approximated. */
static int lya_rprt_quad_tile(const lya_cells *t,const lya_rprt_quad *quad,
                                struct cmdline_data *cmd,const lya_cell *p,const lya_cell *q,
                                const lya_pair_angle *angle,bodyptr base,INTEGER first,
                                INTEGER last,REAL cutoff,lya_worker_hist *w,
                                lya_pair_work *work,ErrorMsg error)
{
    work->quad_tiles++;
    const REAL ps=quad->pscale,ts=quad->tscale;
    const size_t nq=q->end-q->begin;
    for(size_t i=p->begin;i<p->end;i++) {
        const lya_rprt_pixel a=quad->pixels[i];
        int bp[8],bt[8];
#pragma omp simd
        for(size_t j=0;j<nq;j++) {
            const lya_rprt_pixel b=quad->pixels[q->begin+j];
            REAL diff=fabs(a.chi-b.chi),sum=a.chi+b.chi;
            REAL pad=4096*DBL_EPSILON*sum+16*sqrt(DBL_MIN);
            REAL ph_raw=diff*angle->cmax+pad,th_raw=sum*angle->smax+pad;
            REAL pl=MAX(0.,diff*angle->cmin-pad)*ps,ph=ph_raw*ps;
            REAL tl=MAX(0.,sum*angle->smin-pad)*ts,th=th_raw*ts;
            bp[j]=bt[j]=-1;
            if(ph_raw<cmd->lya2RpMax && th_raw<cmd->lya2RtMax && pl>=0 && ph<cmd->lya2RpBins && pl<=ph
                && tl>=0 && th<cmd->lya2RtBins && tl<=th) {
                int ip=(int)pl,it=(int)tl;
                if(ph<ip+1. && th<it+1.) {bp[j]=ip;bt[j]=it;}
            }
        }
        for(size_t j=0;j<nq;j++) {
            if(bp[j]>=0) {
                const lya_rprt_pixel b=quad->pixels[q->begin+j];
                if(lya_pair_add_count(w,1,1,error)==FAILURE) return FAILURE;
                lya_pair_deposit(w,(size_t)bp[j]*cmd->lya2RtBins+bt[j],a.field*b.field,a.weight*b.weight);
                work->leaf_pairs++;work->certified_leaf_pairs++;work->quad_cached_pairs++;
            } else {
                lya_cell one={0},two={0};one.begin=i;one.end=i+1;
                two.begin=q->begin+j;two.end=two.begin+1;
                if(lya_pair_tile(t,cmd,&one,&two,base,first,last,cutoff,NULL,1,w,work,error)==FAILURE) return FAILURE;
            }
        }
    }
    return SUCCESS;
}

/* Scan products of two compact radial balls together. With radial centres
 * m1,m2 and half-widths h1,h2, the distance error is at most h1+h2 (Jarvis
 * 2004, eqs. 62-63). Project both supports into rp/rt using an enclosing LOS
 * angle, and require the complete rectangle to fit one requested bin.
 * In contrast to a recursive singleton walk, these tiny balls and moments
 * stream from the packed cache. Unresolved blocks retain pixel arithmetic. */
static int lya_rprt_quad_balls(const lya_cells *t,const lya_rprt_quad *quad,
                                struct cmdline_data *cmd,const lya_cell *p,const lya_cell *q,
                                const lya_pair_angle *angle,bodyptr base,INTEGER first,
                                INTEGER last,REAL cutoff,lya_worker_hist *w,
                                lya_pair_work *work,ErrorMsg error)
{
    size_t np=p->end-p->begin,nq=q->end-q->begin;
    int bp[16],bt[16];
    struct {REAL lo,hi,field,weight;} pb[4],qb[4];
    int nbp=(int)((np+1)/2),nbq=(int)((nq+1)/2);
    REAL ps=quad->pscale,ts=quad->tscale;
    work->quad_ball_tests+=(unsigned long long)nbp*nbq;
    /* Build each tiny ball and its moments once, then scan both arrays.
     * Bounds and moments share the same disjoint radial pixel partition. */
    for(int side=0;side<2;side++) {
        size_t begin=side?q->begin:p->begin,end=side?q->end:p->end;
        int count=side?nbq:nbp;
        for(int i=0;i<count;i++) {
            const lya_rprt_pixel a=quad->pixels[begin+2*i];
            const lya_rprt_pixel b=quad->pixels[MIN(begin+2*i+1,end-1)];
            REAL field=a.field,weight=a.weight;
            if(begin+2*i+1<end) {field+=b.field;weight+=b.weight;}
            if(side) {qb[i].lo=a.chi;qb[i].hi=b.chi;qb[i].field=field;qb[i].weight=weight;}
            else {pb[i].lo=a.chi;pb[i].hi=b.chi;pb[i].field=field;pb[i].weight=weight;}
        }
    }
    for(int i=0;i<nbp;i++) {
#pragma omp simd
        for(int j=0;j<nbq;j++) {
            REAL gap=MAX(0.,MAX(pb[i].lo-qb[j].hi,qb[j].lo-pb[i].hi));
            REAL span=MAX(fabs(pb[i].lo-qb[j].hi),fabs(pb[i].hi-qb[j].lo));
            REAL sumlo=pb[i].lo+qb[j].lo,sumhi=pb[i].hi+qb[j].hi;
            REAL pad=4096*DBL_EPSILON*sumhi+16*sqrt(DBL_MIN);
            REAL pl=MAX(0.,gap*angle->cmin-pad)*ps,ph=(span*angle->cmax+pad)*ps;
            REAL tl=MAX(0.,sumlo*angle->smin-pad)*ts,th=(sumhi*angle->smax+pad)*ts;
            int k=i*nbq+j;bp[k]=bt[k]=-1;
            if(pl>=0 && ph<cmd->lya2RpBins && pl<=ph && tl>=0 && th<cmd->lya2RtBins && tl<=th) {
                int ip=(int)pl,it=(int)tl;
                if(ph<ip+1. && th<it+1.) {bp[k]=ip;bt[k]=it;}
            }
        }
    }
    for(int i=0;i<nbp;i++) for(int j=0;j<nbq;j++) {
        size_t pbegin=p->begin+2*i,qbegin=q->begin+2*j;
        size_t pn=MIN((size_t)2,p->end-pbegin),qn=MIN((size_t)2,q->end-qbegin);
        int k=i*nbq+j;
        if(bp[k]>=0) {
            REAL numerator=pb[i].field*qb[j].field,denominator=pb[i].weight*qb[j].weight;
            if(isfinite(numerator) && isfinite(denominator) && (denominator>=DBL_MIN || denominator==0)) {
                if(lya_pair_add_count(w,pn,qn,error)==FAILURE) return FAILURE;
                lya_pair_deposit(w,(size_t)bp[k]*cmd->lya2RtBins+bt[k],numerator,denominator);
                work->aggregates++;work->aggregated_pairs+=(unsigned long long)pn*qn;
                work->quad_single_bins++;
                work->quad_ball_products++;work->quad_ball_pairs+=(unsigned long long)pn*qn;
                continue;
            }
        }
        lya_cell one={0},two={0};one.begin=pbegin;one.end=pbegin+pn;two.begin=qbegin;two.end=qbegin+qn;
        if(lya_rprt_quad_tile(t,quad,cmd,&one,&two,angle,base,first,last,cutoff,w,work,error)==FAILURE) return FAILURE;
    }
    return SUCCESS;
}

static int lya_rprt_quad_visit(const lya_cells *t,const lya_rprt_quad *quad,
                                struct cmdline_data *cmd,size_t ip,size_t iq,
                                const lya_pair_angle *inherited,size_t region,int discovery_inside,
                                bodyptr base,INTEGER first,INTEGER last,REAL cutoff,
                                lya_worker_hist *w,lya_pair_work *work,ErrorMsg error)
{
    const lya_cell *p=t->nodes+ip,*q=t->nodes+iq;
    size_t np=p->end-p->begin,nq=q->end-q->begin;
    work->nodes++;
    if(p->forest==q->forest || (!p->pivots && !q->pivots)) return SUCCESS;
    int bounded=sizeof(REAL)==sizeof(double) && sizeof(real)==sizeof(double)
        && sizeof(cballs_storage_real)==sizeof(double);
#ifdef __FAST_MATH__
    bounded=0;
#endif
    REAL upper=0;
    if(!discovery_inside) {
        if(lya_pair_box_outside(p,q,cutoff,&upper) && bounded) {work->pruned++;return SUCCESS;}
        discovery_inside=bounded && isfinite(upper) && upper<cutoff;
    }
    lya_pair_angle a;
    if(inherited) {a=*inherited;work->angle_reuses++;}
    else {a=lya_pair_angles(p,q);work->angle_evaluations++;}
    REAL gap=MAX(0.,MAX(p->chilo-q->chihi,q->chilo-p->chihi));
    REAL span=MAX(fabs(p->chilo-q->chihi),fabs(p->chihi-q->chilo));
    REAL sumlo=p->chilo+q->chilo,sumhi=p->chihi+q->chihi;
    REAL pad=4096*DBL_EPSILON*sumhi+16*sqrt(DBL_MIN);
    REAL plo=MAX(0.,gap*a.cmin-pad),phi=span*a.cmax+pad;
    REAL tlo=MAX(0.,sumlo*a.smin-pad),thi=sumhi*a.smax+pad;
    bounded=bounded && isfinite(phi) && isfinite(thi) && isfinite(upper);
    if(bounded && (plo>=cmd->lya2RpMax || tlo>=cmd->lya2RtMax)) {work->pruned++;return SUCCESS;}
    int inside=bounded && phi<cmd->lya2RpMax && thi<cmd->lya2RtMax;
    /* A rectangle at least one bin wide cannot be certified by any leaf.
     * Skip bin-tree work until a whole-product certificate is possible. The
     * small outward margin only admits extra work near rounded bin edges. */
    if(inside && discovery_inside && p->safe && q->safe
        && phi-plo<=quad->pwidth*(1+64*DBL_EPSILON)
        && thi-tlo<=quad->twidth*(1+64*DBL_EPSILON)) {
        /* Physical-coordinate quadrant comparisons classify rectangles without
         * four divisions/floors at every product. Only a leaf certificate
         * checks the reference arithmetic, guarding rounded grid edges. */
        {
            region=lya_rprt_quad_locate(quad,region,plo,phi,tlo,thi,work);
            const lya_rprt_quad_node *n=quad->nodes+region;
            int owned=(p->pivots==np && q->pivots==nq)
                || (p->pivots==np && p->idhi<q->idlo)
                || (q->pivots==nq && q->idhi<p->idlo);
            if(inside && n->p1-n->p0==1 && n->t1-n->t0==1
                && owned && p->safe && q->safe && discovery_inside) {
                if(lya_bin_positive(plo,cmd->lya2RpMax,cmd->lya2RpBins)!=n->p0
                    || lya_bin_positive(phi,cmd->lya2RpMax,cmd->lya2RpBins)!=n->p0
                    || lya_bin_positive(tlo,cmd->lya2RtMax,cmd->lya2RtBins)!=n->t0
                    || lya_bin_positive(thi,cmd->lya2RtMax,cmd->lya2RtBins)!=n->t0)
                    goto subdivide;
                long double numerator=p->field*q->field,denominator=p->weight*q->weight;
                if(isfinite((REAL)numerator) && isfinite((REAL)denominator)
                    && ((REAL)denominator>=DBL_MIN || denominator==0)) {
                    if(lya_pair_add_count(w,np,nq,error)==FAILURE) return FAILURE;
                    lya_pair_deposit(w,(size_t)n->p0*cmd->lya2RtBins+n->t0,
                        (REAL)numerator,(REAL)denominator);
                    work->aggregates++;work->aggregated_pairs+=(unsigned long long)np*nq;
                    work->quad_single_bins++;return SUCCESS;
                }
            }
        }
    }
subdivide:
    /* Slop changes only bin placement. A product must still lie completely
     * inside the original histogram rectangle and discovery sphere. Each
     * target bin is selected from the midpoint of its geometric enclosure,
     * never from signed-field/weight centroids. Zero slop in either axis uses
     * the reference bin arithmetic for both ends of that entire interval.
     * Expanded bins can cross an inherited exact region. The leaf directory
     * accesses the target by reference bin indices, independent of that region
     * or rounded physical quadrant cuts. Individual pixels stay exact. */
    if((cmd->lya2RpSlop>0 || cmd->lya2RtSlop>0) && (np>1 || nq>1)
        && inside && discovery_inside && p->safe && q->safe
        && phi-plo<=(1+2*cmd->lya2RpSlop)*quad->pwidth
        && thi-tlo<=(1+2*cmd->lya2RtSlop)*quad->twidth
        && ((p->pivots==np && q->pivots==nq)
            || (p->pivots==np && p->idhi<q->idlo)
            || (q->pivots==nq && q->idhi<p->idlo))) {
        REAL pc=.5*plo+.5*phi,tc=.5*tlo+.5*thi;
        int bp=lya_rprt_slop_bin(plo,phi,pc,cmd->lya2RpMax,cmd->lya2RpBins,cmd->lya2RpSlop);
        int bt=lya_rprt_slop_bin(tlo,thi,tc,cmd->lya2RtMax,cmd->lya2RtBins,cmd->lya2RtSlop);
        if(bp>=0 && bt>=0) {
            size_t target=quad->bin_leaves[(size_t)bp*cmd->lya2RtBins+bt];
            const lya_rprt_quad_node *n=quad->nodes+target;
            work->quad_visits++;
            if(n->p1-n->p0==1 && n->t1-n->t0==1 && n->p0==bp && n->t0==bt) {
                long double numerator=p->field*q->field,denominator=p->weight*q->weight;
                if(isfinite((REAL)numerator) && isfinite((REAL)denominator)
                    && ((REAL)denominator>=DBL_MIN || denominator==0)) {
                    if(lya_pair_add_count(w,np,nq,error)==FAILURE) return FAILURE;
                    lya_pair_deposit(w,(size_t)bp*cmd->lya2RtBins+bt,
                        (REAL)numerator,(REAL)denominator);
                    work->aggregates++;work->aggregated_pairs+=(unsigned long long)np*nq;
                    /* Counts products admitted through the slop certificate,
                     * not the number of individual pairs that changed bins. */
                    work->approximate_pairs+=(unsigned long long)np*nq;
                    work->quad_slop_products++;
                    return SUCCESS;
                }
            }
        }
    }
    if(np<=8 && nq<=8 && bounded && discovery_inside
        && ((p->pivots==np && q->pivots==nq)
            || (p->pivots==np && p->idhi<q->idlo)
            || (q->pivots==nq && q->idhi<p->idlo)))
        return inside && p->safe && q->safe && np>1 && nq>1
            /* Sampling gate only schedules work; fallback tiles are exact.
             * Sparse pixels rarely amortize another whole-product test. */
            && .5*((p->chihi-p->chilo)/(np-1)+(q->chihi-q->chilo)/(nq-1))
                *MAX(a.cmax*quad->popen,a.smax*quad->topen)<.9
            ?lya_rprt_quad_balls(t,quad,cmd,p,q,&a,base,first,last,cutoff,w,work,error)
            :lya_rprt_quad_tile(t,quad,cmd,p,q,&a,base,first,last,cutoff,w,work,error);
    if(np<=8 && nq<=8)
        return lya_pair_tile(t,cmd,p,q,base,first,last,cutoff,bounded?&a:NULL,
                             discovery_inside,w,work,error);
    /* Refine the implicit pair cloud in both source coordinates. Each child
     * is projected afresh into the same 2D quadtree; no pixel pair is stored,
     * revisited, dropped or counted in two bin-plane regions. */
    size_t pp[2]={ip,SIZE_MAX},qq[2]={iq,SIZE_MAX};
    int pn=1,qn=1;
    int splitp=np>8 && p->left!=SIZE_MAX,splitq=nq>8 && q->left!=SIZE_MAX;
    if(splitp && splitq && bounded) {
        /* Two-ball sizes in bin units. For fixed half-angle, a radial ball
         * of radius h displaces rp by at most h*c and rt by at most h*s.
         * Angular enclosures and physical cuts are still checked above;
         * this choice changes traversal only, never geometric acceptance.
         * As in dual node, keep a small partner intact, or split both when
         * both radii are large compared with the permitted bin footprint. */
        REAL scale=MAX(a.cmax*quad->popen,a.smax*quad->topen);
        REAL sp=.5*(p->chihi-p->chilo)*scale,sq=.5*(q->chihi-q->chilo)*scale;
        if(sp>2*sq) splitq=0;
        else if(sq>2*sp) splitp=0;
        else if(sp>=sq) splitq=sq>.2925;
        else splitp=sp>.2925;
    }
    if(splitp) {pp[0]=p->left;pp[1]=p->right;pn=2;}
    if(splitq) {qq[0]=q->left;qq[1]=q->right;qn=2;}
    work->quad_splits++;
    if(pn==2 && qn==2) work->quad_split_both++;else work->quad_split_one++;
    /* An ancestor LOS enclosure covers every descendant, even for bent
     * forests. Tighten it only when uncertainty is appreciable in bin units.
     * Reuse is an exact bound, never a common-ray assumption. */
    const lya_pair_angle *reuse=bounded
        && span*(a.cmax-a.cmin)*quad->pscale
           +sumhi*(a.smax-a.smin)*quad->tscale<.125?&a:NULL;
    for(int i=0;i<pn;i++) for(int j=0;j<qn;j++) {
        if(lya_rprt_quad_visit(t,quad,cmd,pp[i],qq[j],reuse,region,discovery_inside,base,first,last,cutoff,w,work,error)==FAILURE)
            return FAILURE;
    }
    return SUCCESS;
}
#endif
