/* Unordered cross-forest cell pairs, independently controlled in (rp,rt).
 * N = (sum w*delta)_p (sum w*delta)_q; D = (sum w)_p (sum w)_q.
 * No signed-field centroid, ideal-ray assumption, or third-side cut.
 * Each root pair and descendant Cartesian product is visited once. Geometry
 * reuse therefore passes immutable angular bounds down the recursion, rather
 * than allocating an N^2 pixel cache with no repeated-pair reuse opportunity. */
typedef struct { REAL cmin,cmax,smin,smax; } lya_pair_angle;
typedef struct {
    unsigned long long nodes, pruned, aggregates, aggregated_pairs, approximate_pairs;
    unsigned long long leaf_pairs, angle_evaluations, angle_reuses, certified_leaf_pairs, window_skips, scan_products;
    unsigned long long range_aggregates,range_pairs;
} lya_pair_work;

static uintmax_t lya_pair_count_max(void)
{ return ((uintmax_t)1 << (sizeof(INTEGER)*CHAR_BIT-1))-1; }

static int lya_pair_add_count(lya_worker_hist *w,size_t np,size_t nq,ErrorMsg err)
{
    size_t count;
    if (!cballs_size_mul(np,nq,&count)
        || (uintmax_t)count>lya_pair_count_max()-(uintmax_t)w->pair_count) {
        snprintf(err,_ERRORMSGSIZE_,"Ly-alpha pair count overflow");return FAILURE;
    }
    w->pair_count+=(INTEGER)count;return SUCCESS;
}

static lya_pair_angle lya_pair_angles(const lya_cell *p,const lya_cell *q)
{
    REAL lo,hi;
    lya_cell_dot_bounds(p->loslo,p->loshi,q->loslo,q->loshi,&lo,&hi);
    lya_pair_angle a={sqrt(MAX(0.,.5*(1+lo))),sqrt(MAX(0.,.5*(1+hi))),
                      sqrt(MAX(0.,.5*(1-hi))),sqrt(MAX(0.,.5*(1-lo)))};
    return a;
}

static int lya_pair_same_los(const lya_cell *p,const lya_cell *child)
{
    for(int k=0;k<3;k++)
        if(fabs(p->loslo[k]-child->loslo[k])>32*DBL_EPSILON
            || fabs(p->loshi[k]-child->loshi[k])>32*DBL_EPSILON) return 0;
    return 1;
}

/* Cheap box rejection before angular geometry or any descendant pixel work.
 * Padded Cartesian bounds also protect the discovery-sphere convention of
 * the reference walker at floating-point corners of the rp/rt rectangle. */
static int lya_pair_box_outside(const lya_cell *p,const lya_cell *q,REAL cutoff,
                                  REAL *distance_upper)
{
    REAL d2=0,u2=0,scale=0;
    for(int k=0;k<3;k++) {
        REAL lo=q->lo[k]-p->hi[k],hi=q->hi[k]-p->lo[k];
        REAL d=MAX(0.,MAX(lo,-hi)),u=MAX(fabs(lo),fabs(hi));
        d2+=d*d;u2+=u*u;
        scale=MAX(scale,MAX(MAX(fabs(p->lo[k]),fabs(p->hi[k])),MAX(fabs(q->lo[k]),fabs(q->hi[k]))));
    }
    REAL pad=4096*DBL_EPSILON*scale+16*sqrt(DBL_MIN);
    *distance_upper=sqrt(u2)+pad;
    return sqrt(d2)-pad>=cutoff;
}

static void lya_pair_deposit(lya_worker_hist *w,size_t bin,REAL n,REAL d)
{
    if(d<=0) return;
    if(w->den2[bin]==0) w->touched2[w->touched2_count++]=bin;
    w->num2[bin]+=n;w->den2[bin]+=d;
}

/* Small tiles keep the reference arithmetic/order of each pair, with SIMD
 * geometry separated from histogram writes (which can share a bin). */
static int lya_pair_tile(const lya_cells *t,struct cmdline_data *cmd,
                           const lya_cell *p,const lya_cell *q,bodyptr base,
                           INTEGER first,INTEGER last,REAL cutoff,const lya_pair_angle *angle,int discovery_inside,
                           lya_worker_hist *w,lya_pair_work *work,ErrorMsg err)
{
    for(size_t i=p->begin;i<p->end;i++) {
        bodyptr a=t->points[i];
        int parallel[8],transverse[8],certified[8];
        const size_t n=q->end-q->begin;
        const REAL pscale=cmd->lya2RpBins/cmd->lya2RpMax;
        const REAL tscale=cmd->lya2RtBins/cmd->lya2RtMax;
#pragma omp simd
        for(size_t j=0;j<n;j++) {
            bodyptr b=t->points[q->begin+j];
            bodyptr pivot=Id(a)<Id(b)?a:b;
            parallel[j]=transverse[j]=-1;certified[j]=0;
            if(Id(a)==Id(b) || pivot-base<first || pivot-base>=last) continue;
            REAL diff=fabs(LyaDistance(a)-LyaDistance(b)),sum=LyaDistance(a)+LyaDistance(b);
            /* Reuse the inherited observer-angle enclosure for each pixel in
             * the tile. Guard both edges with padding; ambiguous pixels retain
             * the original dot/sqrt/bin arithmetic below. No slop is used. */
            if(angle) {
                REAL pad=4096*DBL_EPSILON*sum+16*sqrt(DBL_MIN);
                REAL pl=MAX(0.,diff*angle->cmin-pad),ph=diff*angle->cmax+pad;
                REAL tl=MAX(0.,sum*angle->smin-pad),th=sum*angle->smax+pad;
                if(ph<cmd->lya2RpMax && th<cmd->lya2RtMax) {
                    REAL plb=pl*pscale,phb=ph*pscale,tlb=tl*tscale,thb=th*tscale;
                    /* Bound floating-to-int conversions even for extreme
                     * accepted parameter values (tiny maxima / huge scales). */
                    if(phb>=0 && phb<cmd->lya2RpBins && thb>=0 && thb<cmd->lya2RtBins
                        && plb>=0 && plb<=phb && tlb>=0 && tlb<=thb
                        && (int)plb==(int)phb && (int)tlb==(int)thb) {
                        int bp=(int)plb,bt=(int)tlb;
                        if(!discovery_inside) {
                            REAL d2;compute_vector dr;DOTPSUBV(d2,dr,Pos(a),Pos(b));
                            if(rsqrt(d2)>=cutoff) continue;
                        }
                        parallel[j]=bp;transverse[j]=bt;certified[j]=1;continue;
                    }
                }
            }
            REAL d2,cosine;compute_vector dr;
            DOTPSUBV(d2,dr,Pos(a),Pos(b));
            if(rsqrt(d2)>=cutoff) continue;
            DOTVP(cosine,LyaLOS(a),LyaLOS(b));
            cosine=lya_clamp(cosine,-1.,1.);
            REAL rp=diff*rsqrt(MAX(0.,.5*(1+cosine)));
            REAL rt=sum*rsqrt(MAX(0.,.5*(1-cosine)));
            parallel[j]=lya_bin_positive(rp,cmd->lya2RpMax,cmd->lya2RpBins);
            transverse[j]=lya_bin_positive(rt,cmd->lya2RtMax,cmd->lya2RtBins);
        }
        work->leaf_pairs+=n;
        for(size_t j=0;j<n;j++) {
            int bp=parallel[j],bt=transverse[j];
            work->certified_leaf_pairs+=certified[j];
            if(bp<0 || bt<0) continue;
            if(lya_pair_add_count(w,1,1,err)==FAILURE) return FAILURE;
            bodyptr b=t->points[q->begin+j];
            lya_pair_deposit(w,(size_t)bp*cmd->lya2RtBins+bt,
                (Weight(a)*Kappa(a))*(Weight(b)*Kappa(b)),Weight(a)*Weight(b));
        }
    }
    return SUCCESS;
}

/* Binary searches in the persistent distance ordering. Inclusive outward
 * windows only remove pairs proved outside rp/rt before any pixel geometry. */
static size_t lya_pair_lower(const lya_cells *t,size_t lo,size_t hi,REAL chi,int strict)
{
    while(lo<hi) {
        size_t mid=lo+(hi-lo)/2;REAL value=LyaDistance(t->points[mid]);
        if(value<chi || (strict && value==chi)) lo=mid+1;else hi=mid;
    }
    return lo;
}

/* Query immutable dyadic moments without subtracting large prefix sums. The
 * bounding box covers every actual pixel, including bent sightlines. */
static void lya_pair_range_moments(const lya_cells *t,size_t node,size_t begin,
                                    size_t end,lya_cell *sum)
{
    const lya_cell *q=&t->nodes[node];
    if(q->end<=begin || q->begin>=end) return;
    if(q->begin>=begin && q->end<=end) {
        sum->weight+=q->weight;sum->field+=q->field;
        for(int k=0;k<3;k++) {sum->lo[k]=MIN(sum->lo[k],q->lo[k]);sum->hi[k]=MAX(sum->hi[k],q->hi[k]);}
        return;
    }
    lya_pair_range_moments(t,q->left,begin,end,sum);
    lya_pair_range_moments(t,q->right,begin,end,sum);
}

static int lya_pair_range_bin(const lya_cells *t,const struct cmdline_data *cmd,
                               REAL chi,size_t begin,size_t end,const lya_pair_angle *angle)
{
    REAL low=LyaDistance(t->points[begin]),high=LyaDistance(t->points[end-1]);
    REAL dlo=MAX(0.,MAX(low-chi,chi-high)),dhi=MAX(fabs(low-chi),fabs(high-chi));
    REAL pad=8192*DBL_EPSILON*(chi+high)+32*sqrt(DBL_MIN);
    REAL pl=MAX(0.,dlo*angle->cmin-pad),ph=dhi*angle->cmax+pad;
    REAL tl=MAX(0.,(chi+low)*angle->smin-pad),th=(chi+high)*angle->smax+pad;
    if(!isfinite(ph)||!isfinite(th)||ph>=cmd->lya2RpMax||th>=cmd->lya2RtMax) return -1;
    int bp=lya_bin_positive(pl,cmd->lya2RpMax,cmd->lya2RpBins);
    int bt=lya_bin_positive(tl,cmd->lya2RtMax,cmd->lya2RtBins);
    if(bp<0 || bt<0 || bp!=lya_bin_positive(ph,cmd->lya2RpMax,cmd->lya2RpBins)
        || bt!=lya_bin_positive(th,cmd->lya2RtMax,cmd->lya2RtBins)) return -1;
    /* Histogram dimensions were preflighted; return the axes via a size_t at
     * the caller instead of multiplying int bin counts here. */
    return bp;
}

static int lya_pair_ranges(const lya_cells *t,struct cmdline_data *cmd,
                            const lya_cell *p,const lya_cell *q,size_t qi,
                            size_t pi,size_t begin,size_t end,const lya_pair_angle *angle,
                            bodyptr base,INTEGER first,INTEGER last,REAL cutoff,lya_worker_hist *w,lya_pair_work *work,ErrorMsg err)
{
    bodyptr a=t->points[pi];REAL chi=LyaDistance(a);
    lya_cell pixel={0};pixel.begin=pi;pixel.end=pi+1;
    for(int k=0;k<3;k++) pixel.lo[k]=pixel.hi[k]=Pos(a)[k];
    for(size_t j=begin;j<end;) {
        int bp=-1;size_t stop=j;
        if(end-j>=8) bp=lya_pair_range_bin(t,cmd,chi,j,j+8,angle);
        if(bp>=0) {
            /* Exponential search limits overhead for short bin runs. The
             * enclosure only widens when the range grows, so certification
             * cannot resume after its first failure. */
            size_t good=j+8,bad=end+1,step=8;
            while(good<end) {
                size_t trial=good+MIN(step,end-good);
                if(lya_pair_range_bin(t,cmd,chi,j,trial,angle)<0) {bad=trial;break;}
                good=trial;step=MIN(step,SIZE_MAX/2)*2;
            }
            if(bad<=end) while(bad-good>1) {
                size_t mid=good+(bad-good)/2;
                if(lya_pair_range_bin(t,cmd,chi,j,mid,angle)>=0) good=mid;else bad=mid;
            }
            stop=good;
            lya_cell sum={0};for(int k=0;k<3;k++) {sum.lo[k]=INFINITY;sum.hi[k]=-INFINITY;}
            lya_pair_range_moments(t,qi,j,stop,&sum);
            REAL upper;lya_pair_box_outside(&pixel,&sum,cutoff,&upper);
            long double ld=(long double)Weight(a)*sum.weight;
            long double ln=(long double)(Weight(a)*Kappa(a))*sum.field;
            if(upper<cutoff && isfinite((REAL)ld) && isfinite((REAL)ln) && ((REAL)ld>=DBL_MIN||ld==0)) {
                REAL pad=8192*DBL_EPSILON*(chi+LyaDistance(t->points[stop-1]))+32*sqrt(DBL_MIN);
                REAL tl=MAX(0.,(chi+LyaDistance(t->points[j]))*angle->smin-pad);
                int bt=lya_bin_positive(tl,cmd->lya2RtMax,cmd->lya2RtBins);
                if(lya_pair_add_count(w,1,stop-j,err)==FAILURE) return FAILURE;
                lya_pair_deposit(w,(size_t)bp*cmd->lya2RtBins+bt,(REAL)ln,(REAL)ld);
                work->range_aggregates++;work->range_pairs+=stop-j;
                j=stop;continue;
            }
        }
        /* Ambiguous bins, geometry or exponent ranges retain exact pixels. */
        lya_cell tile={0};tile.begin=j;tile.end=MIN(j+8,end);
        if(lya_pair_tile(t,cmd,&pixel,&tile,base,first,last,cutoff,angle,0,w,work,err)==FAILURE) return FAILURE;
        j=tile.end;
    }
    (void)p;(void)q;
    return SUCCESS;
}

static int lya_pair_scan(const lya_cells *t,struct cmdline_data *cmd,
                           const lya_cell *p,const lya_cell *q,size_t qi,lya_pair_angle *angle,
                           bodyptr base,INTEGER first,INTEGER last,REAL cutoff,
                           lya_worker_hist *w,lya_pair_work *work,ErrorMsg err)
{
    work->scan_products++;
    /* Sparse regular sampling cannot amortize a moment query. Check once per
     * forest product for any eight-pixel run that could fit the bins. The
     * parallel absolute value can fold at the pivot, hence the factor 1/2.
     * This is a scheduling heuristic only; skipping it retains exact tiles. */
    int ranges=0;
    if(p->safe && q->safe && p->pivots==p->end-p->begin && q->pivots==q->end-q->begin)
        for(size_t j=q->begin;j+7<q->end;j++) {
            REAL span=LyaDistance(t->points[j+7])-LyaDistance(t->points[j]);
            if(.5*span*angle->cmin<cmd->lya2RpMax/cmd->lya2RpBins
                && span*angle->smin<cmd->lya2RtMax/cmd->lya2RtBins) {ranges=1;break;}
        }
    for(size_t i=p->begin;i<p->end;i++) {
        REAL chi=LyaDistance(t->points[i]);
        REAL margin=8192*DBL_EPSILON*(chi+q->chihi)+32*sqrt(DBL_MIN);
        REAL span=angle->cmin>0?cmd->lya2RpMax/angle->cmin:INFINITY;
        REAL low=chi-span-margin,high=chi+span+margin;
        if(angle->smin>0) high=MIN(high,cmd->lya2RtMax/angle->smin-chi+margin);
        size_t begin=lya_pair_lower(t,q->begin,q->end,low,0);
        size_t end=lya_pair_lower(t,begin,q->end,high,1);
        work->window_skips+=(q->end-q->begin)-(end-begin);
        if(ranges) {
            if(lya_pair_ranges(t,cmd,p,q,qi,i,begin,end,angle,base,first,last,cutoff,w,work,err)==FAILURE) return FAILURE;
            continue;
        }
        lya_cell pixel={0};pixel.begin=i;pixel.end=i+1;
        for(size_t j=begin;j<end;j+=8) {
            lya_cell tile={0};tile.begin=j;tile.end=MIN(j+8,end);
            if(lya_pair_tile(t,cmd,&pixel,&tile,base,first,last,cutoff,angle,0,w,work,err)==FAILURE) return FAILURE;
        }
    }
    return SUCCESS;
}

static int lya_pair_visit(const lya_cells *t,struct cmdline_data *cmd,size_t ip,size_t iq,
                            const lya_pair_angle *inherited,bodyptr base,INTEGER first,
                            INTEGER last,REAL cutoff,lya_worker_hist *w,
                            lya_pair_work *work,ErrorMsg err)
{
    const lya_cell *p=&t->nodes[ip],*q=&t->nodes[iq];
    const size_t np=p->end-p->begin,nq=q->end-q->begin;
    work->nodes++;
    if(!p->pivots && !q->pivots) return SUCCESS;
    /* For profiles not covered by double-precision interval padding, descend
     * directly to the retained arithmetic. No approximate acceptance there. */
    int bounded=sizeof(REAL)==sizeof(double);
#ifdef __FAST_MATH__
    bounded=0;
#endif
    REAL upper;
    if(lya_pair_box_outside(p,q,cutoff,&upper) && bounded) {work->pruned++;return SUCCESS;}
    lya_pair_angle a;
    if(inherited) {a=*inherited;work->angle_reuses++;}
    else {a=lya_pair_angles(p,q);work->angle_evaluations++;}
    REAL difflo=MAX(0.,MAX(p->chilo-q->chihi,q->chilo-p->chihi));
    REAL diffhi=MAX(fabs(p->chilo-q->chihi),fabs(p->chihi-q->chilo));
    REAL sumlo=p->chilo+q->chilo,sumhi=p->chihi+q->chihi;
    REAL pad=4096*DBL_EPSILON*sumhi+16*sqrt(DBL_MIN);
    REAL plo=MAX(0.,difflo*a.cmin-pad),phi=diffhi*a.cmax+pad;
    REAL tlo=MAX(0.,sumlo*a.smin-pad),thi=sumhi*a.smax+pad;
    bounded=bounded && isfinite(phi) && isfinite(thi) && isfinite(upper);
    if(bounded && (plo>=cmd->lya2RpMax || tlo>=cmd->lya2RtMax)) {work->pruned++;return SUCCESS;}
    /* Whole-cell ownership is valid if all potential lower-ID pivots are in
     * the requested range. Mixed ownership always reaches exact leaf checks. */
    int owned=(p->pivots==np && q->pivots==nq)
        || (p->pivots==np && p->idhi<q->idlo)
        || (q->pivots==nq && q->idhi<p->idlo);
    if(bounded && owned && p->safe && q->safe
        && phi<cmd->lya2RpMax && thi<cmd->lya2RtMax && upper<cutoff) {
        REAL cp=.5*p->chilo+.5*p->chihi,cq=.5*q->chilo+.5*q->chihi,cosine=0;
        for(int k=0;k<3;k++) cosine+=p->sight[k]*q->sight[k];
        cosine=lya_clamp(cosine,-1.,1.);
        REAL rp=fabs(cp-cq)*sqrt(MAX(0.,.5*(1+cosine)));
        REAL rt=(cp+cq)*sqrt(MAX(0.,.5*(1-cosine)));
        int approx=0;
        int bp=lya_cell_bin(plo,phi,rp,cmd->lya2RpMax,cmd->lya2RpBins,cmd->lya2RpSlop,0,&approx);
        int bt=lya_cell_bin(tlo,thi,rt,cmd->lya2RtMax,cmd->lya2RtBins,cmd->lya2RtSlop,0,&approx);
        if(bp>=0 && bt>=0) {
            if(lya_pair_add_count(w,np,nq,err)==FAILURE) return FAILURE;
            lya_pair_deposit(w,(size_t)bp*cmd->lya2RtBins+bt,
                (REAL)(p->field*q->field),(REAL)(p->weight*q->weight));
            work->aggregates++;work->aggregated_pairs+=(unsigned long long)np*nq;
            if(approx) work->approximate_pairs+=(unsigned long long)np*nq;
            return SUCCESS;
        }
    }
    /* For sparsely sampled, narrow forests and fine bins, subdivision cannot
     * often sum more than one pixel. Scan bounded radial windows instead of
     * spending most time constructing tiny unresolved cell products. This
     * is an exact scheduling choice, independent of the slop contract. */
    REAL bin_scale=a.cmax*cmd->lya2RpBins/cmd->lya2RpMax
                  +a.smax*cmd->lya2RtBins/cmd->lya2RtMax;
    if(bounded && np>=16 && nq>=16 && a.cmax-a.cmin+a.smax-a.smin<1e-5
        && MIN((p->chihi-p->chilo)/np,(q->chihi-q->chilo)/nq)*bin_scale>.25)
        return lya_pair_scan(t,cmd,p,q,iq,&a,base,first,last,cutoff,w,work,err);
    if(np<=8 && nq<=8) return lya_pair_tile(t,cmd,p,q,base,first,last,cutoff,
        bounded?&a:NULL,upper<cutoff,w,work,err);
    /* Split by uncertainty in the two estimator coordinates, including the
     * observer-LOS extent. Both nodes can split; no fixed pivot bottleneck. */
    REAL radial_scale=a.cmax*cmd->lya2RpBins/cmd->lya2RpMax
                     +a.smax*cmd->lya2RtBins/cmd->lya2RtMax;
    REAL angular_scale=diffhi*cmd->lya2RpBins/cmd->lya2RpMax
                      +sumhi*cmd->lya2RtBins/cmd->lya2RtMax;
    REAL ps=(p->chihi-p->chilo)*radial_scale+p->los_spread*angular_scale;
    REAL qs=(q->chihi-q->chilo)*radial_scale+q->los_spread*angular_scale;
    int split_p=p->left!=SIZE_MAX && (q->left==SIZE_MAX || ps>qs || (ps==qs && np>=nq));
    const lya_cell *parent=split_p?p:q;
    size_t children[2]={parent->left,parent->right};
    for(int k=0;k<2;k++) {
        const lya_cell *child=&t->nodes[children[k]];
        /* An ancestor enclosure remains valid after a split. Reuse it when
         * LOS extrema differ only at roundoff scale; recompute for real tightening. */
        const lya_pair_angle *reuse=lya_pair_same_los(parent,child)?&a:NULL;
        if(lya_pair_visit(t,cmd,split_p?children[k]:ip,split_p?iq:children[k],reuse,
                base,first,last,cutoff,w,work,err)==FAILURE) return FAILURE;
    }
    return SUCCESS;
}

static int lya_pairs_run(const lya_cells *t,struct cmdline_data *cmd,struct global_data *gd,
                           bodyptr base,INTEGER first,INTEGER last,size_t bins,size_t plan,
                           REAL *num,REAL *den,INTEGER *pairs,INTEGER *visits)
{
    int failed=0;ErrorMsg error="";lya_pair_work total={0};
    double start=CPUTIME;
    REAL cutoff=hypot(cmd->lya2RpMax,cmd->lya2RtMax);
    const size_t first_row=lya_parallel_first(cmd),row_stride=lya_parallel_stride(cmd);
#pragma omp parallel
    {
        lya_worker_hist w;lya_pair_work work={0};ErrorMsg local_error="";
        int ready=lya_worker_init(cmd,&w,bins,0,1,0,plan,local_error)==SUCCESS;
        int worker_failed=!ready;
        /* Rows have uneven numbers of forest partners. Dynamic claims balance
         * them; ordered row commits make histogram summation reproducible. */
#pragma omp for schedule(dynamic,1) ordered
        for(size_t row=first_row;row<t->roots_count;row+=row_stride) {
            w.pair_count=0;
            for(size_t j=row+1;!worker_failed && j<t->roots_count;j++)
                if(lya_pair_visit(t,cmd,t->roots[row],t->roots[j],NULL,base,first,last,
                                  cutoff,&w,&work,local_error)==FAILURE) worker_failed=1;
#pragma omp ordered
            {
                if(!failed && !worker_failed) {
                    if((uintmax_t)w.pair_count>lya_pair_count_max()-(uintmax_t)*pairs) {
                        worker_failed=1;snprintf(local_error,sizeof(local_error),"Ly-alpha global pair count overflow");
                    } else {lya_commit_worker(&w,num,den,NULL,NULL,1,0);*pairs+=w.pair_count;}
                }
                if(worker_failed && !failed) {failed=1;snprintf(error,sizeof(error),"%s",local_error);}
            }
        }
#pragma omp critical(lya_pair_work_counts)
        {
            if(worker_failed && !failed) {failed=1;snprintf(error,sizeof(error),"%s",local_error);}
            total.nodes+=work.nodes;total.pruned+=work.pruned;total.aggregates+=work.aggregates;
            total.aggregated_pairs+=work.aggregated_pairs;total.approximate_pairs+=work.approximate_pairs;
            total.leaf_pairs+=work.leaf_pairs;total.angle_evaluations+=work.angle_evaluations;total.angle_reuses+=work.angle_reuses;total.certified_leaf_pairs+=work.certified_leaf_pairs;
            total.range_aggregates+=work.range_aggregates;total.range_pairs+=work.range_pairs;
            total.window_skips+=work.window_skips;total.scan_products+=work.scan_products;
        }
        if(ready) lya_worker_free(&w);
    }
    gd->lyaHierarchyCounts[4]+=total.range_aggregates;gd->lyaHierarchyCounts[5]+=total.range_pairs;
    /* Work count: direct candidate pairs, excluding aggregated products. */
    if(total.leaf_pairs<=lya_pair_count_max()-(uintmax_t)*visits) *visits+=(INTEGER)total.leaf_pairs;
    else {failed=1;snprintf(error,sizeof(error),"Ly-alpha pair work count overflow");}
    verb_print_normal_info(cmd->verbose,cmd->verbose_log,gd->outlog,
        "Ly-alpha 2PCF cells: nodes=%llu pruned=%llu aggregates=%llu aggregated_pairs=%llu approximate_pairs=%llu leaf_pairs=%llu angle_evaluations=%llu angle_reuses=%llu certified_leaf_pairs=%llu window_skips=%llu scan_products=%llu search_CPU=%g\n",
        total.nodes,total.pruned,total.aggregates,total.aggregated_pairs,total.approximate_pairs,
        total.leaf_pairs,total.angle_evaluations,total.angle_reuses,total.certified_leaf_pairs,total.window_skips,total.scan_products,CPUTIME-start);
    if(failed) snprintf(cmd->error_message,_ERRORMSGSIZE_,"%s",error);
    return failed?FAILURE:SUCCESS;
}
