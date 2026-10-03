/* Explicitly anisotropic Legendre moments, with exact radial/polar bins.
 * Real basis H_lm=sqrt(4*pi/(2*l+1))*Y_lm, hence dot(H_l(u),H_l(v))=P_l(u.v).
 * Accumulate each forest against the prefix of OTHER forests, publishing both
 * leg orders. This implements same-forest subtraction without catastrophic
 * cancellation of a dominant forest's auto product. */
static void lya_real_harmonics(const lya_neighbor *q, int limit, REAL *h)
{
    REAL z=lya_clamp(q->displacement[2]/q->radius,-1,1);
    REAL transverse=hypot(q->displacement[0],q->displacement[1]);
    REAL sin_theta=transverse/q->radius;
    REAL cp=transverse>0 ? q->displacement[0]/transverse : 1;
    REAL sp=transverse>0 ? q->displacement[1]/transverse : 0;
    REAL cm=1,sm=0,diagonal=1;
    for (int m=0;m<=limit;m++) {
        if (m) {
            REAL next=cm*cp-sm*sp; sm=sm*cp+cm*sp; cm=next;
            diagonal*=-sqrt((2.*m-1)/(2.*m))*sin_theta;
        }
        REAL prev=0,current=diagonal;
        for (int l=m;l<=limit;l++) {
            if (l>m) {
                REAL next=((2.*l-1)*z*current-sqrt((double)((l-1)*(l-1)-m*m))*prev)
                          /sqrt((double)(l*l-m*m));
                prev=current; current=next;
            }
            if (!m) h[l*l]=current;
            else { h[l*l+2*m-1]=sqrt(2.)*current*cm; h[l*l+2*m]=sqrt(2.)*current*sm; }
        }
    }
}

static int lya_multipole_scratch(struct cmdline_data *cmd, lya_worker_hist *w, ErrorMsg err)
{
    if (w->moments) return SUCCESS;
    size_t b,h,count,bytes,total,lists,pairs;
    h=(size_t)(cmd->lya3LMax+1)*(cmd->lya3LMax+1);
    if (!cballs_size_mul((size_t)cmd->lya3RBins,cmd->lya3ThetaBins,&b)
        || !cballs_size_mul(b,6*h,&count) || !cballs_size_add(count,3*h,&count)
        || !cballs_size_mul(count,sizeof(REAL),&bytes)
        || !cballs_size_mul(b,5*sizeof(size_t)+2*sizeof(unsigned char),&lists)
        || !cballs_size_mul(b,b,&pairs) || !cballs_size_add(lists,pairs,&lists)
        || !cballs_size_add(bytes,lists,&total)
        || !cballs_size_mul(total,w->scratch_workers,&total)
        || !cballs_size_add(total,w->histogram_plan_bytes,&total)) {
        snprintf(err,_ERRORMSGSIZE_,"Ly-alpha multipole scratch dimensions overflow");return FAILURE;
    }
    size_t neighbor_bytes,peak;
    if(!cballs_size_mul(w->neighbor_capacity,sizeof(lya_neighbor),&neighbor_bytes)
        || !cballs_size_mul(neighbor_bytes,w->scratch_workers,&neighbor_bytes)
        || !cballs_size_add(total,neighbor_bytes,&peak)) {
        snprintf(err,_ERRORMSGSIZE_,"Ly-alpha multipole scratch dimensions overflow");return FAILURE;
    }
    if (cballs_memory_preflight(peak,"Ly-alpha multipole scratch",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
    if (cballs_calloc_checked((void**)&w->moments,count,sizeof(REAL),"Ly-alpha moments",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&w->moment_bins,5*b,sizeof(size_t),"Ly-alpha moment bins",err,_ERRORMSGSIZE_)==FAILURE
        || cballs_calloc_checked((void**)&w->moment_seen,2*b+pairs,1,"Ly-alpha moment occupancy",err,_ERRORMSGSIZE_)==FAILURE) return FAILURE;
    w->moment_leg_bins=b;
    /* Cache invariant recurrence square roots and radial output strides once
     * per worker. Preserve the reference recurrence's arithmetic ordering. */
    REAL *den=w->moments+6*b*h+h,*back=den+h;
    for(int m=0;m<=cmd->lya3LMax;m++) {
        if(m) back[m*m+m]=-sqrt((2.*m-1)/(2.*m));
        for(int l=m+1;l<=cmd->lya3LMax;l++) {
            den[l*l+m]=sqrt((double)(l*l-m*m));
            back[l*l+m]=sqrt((double)((l-1)*(l-1)-m*m));
        }
    }
    size_t *first=w->moment_bins+3*b,*second=first+b;
    for(size_t a=0;a<b;a++) {
        first[a]=lya_index3(a/cmd->lya3ThetaBins,0,a%cmd->lya3ThetaBins,0,0,
            cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3LMax+1);
        second[a]=lya_index3(0,a/cmd->lya3ThetaBins,0,a%cmd->lya3ThetaBins,0,
            cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3LMax+1);
    }
    /* Neighbor-growth preflight will now include this persistent scratch. */
    w->histogram_plan_bytes=total;
    return SUCCESS;
}

static int lya_accumulate_multipoles_prefix(struct cmdline_data *cmd, bodyptr p,
                                     lya_worker_hist *w, ErrorMsg err)
{
    if (w->neighbor_count<2) return SUCCESS;
    if (lya_multipole_scratch(cmd,w,err)==FAILURE) return FAILURE;
    const int limit=cmd->lya3LMax, orders=limit+1;
    const size_t b=(size_t)cmd->lya3RBins*cmd->lya3ThetaBins, h=(size_t)orders*orders;
    REAL *tn=w->moments, *td=tn+b*h, *fn=td+b*h, *fd=fn+b*h, *basis=fd+b*h;
    size_t *active=w->moment_bins, *current=active+b, active_count=0, previous_count=0;
    memset(tn,0,2*b*h*sizeof(REAL)); memset(w->moment_seen,0,b);
    /* Raw Legendre denominators can be signed; touched state cannot be
     * inferred from denominator==0 as it is for positive hard-bin sums. */
    if (w->touched3_count==0) {
        size_t count=b*b*orders;
        for (size_t k=0;k<count;k++) w->touched3[w->touched3_count++]=k;
    }
    qsort(w->neighbors,w->neighbor_count,sizeof(*w->neighbors),lya_neighbor_compare);
    w->neighbors_sorted=TRUE;
    for (size_t begin=0;begin<w->neighbor_count;) {
        size_t end=begin, current_count=0;
        memset(fn,0,2*b*h*sizeof(REAL));
        while (end<w->neighbor_count && w->neighbors[end].forest_id==w->neighbors[begin].forest_id) {
            const lya_neighbor *q=&w->neighbors[end];
            size_t bin=q->leg_bin;
            if (end==begin || bin!=w->neighbors[end-1].leg_bin) current[current_count++]=bin;
            lya_real_harmonics(q,limit,basis);
            for (size_t k=0;k<h;k++) {fn[bin*h+k]+=q->weighted_delta*basis[k]; fd[bin*h+k]+=q->weight*basis[k];}
            end++;
        }
        for (size_t ia=0;ia<current_count;ia++) for (size_t ib=0;ib<active_count;ib++) {
            size_t a=current[ia], c=active[ib];
            size_t ab=lya_index3(a/cmd->lya3ThetaBins,c/cmd->lya3ThetaBins,
                a%cmd->lya3ThetaBins,c%cmd->lya3ThetaBins,0,cmd->lya3RBins,cmd->lya3ThetaBins,orders);
            size_t ba=lya_index3(c/cmd->lya3ThetaBins,a/cmd->lya3ThetaBins,
                c%cmd->lya3ThetaBins,a%cmd->lya3ThetaBins,0,cmd->lya3RBins,cmd->lya3ThetaBins,orders);
            for (int l=0;l<=limit;l++) {
                REAL n=0,d=0;
                for (int k=l*l;k<(l+1)*(l+1);k++) {n+=fn[a*h+k]*tn[c*h+k];d+=fd[a*h+k]*td[c*h+k];}
                n*=(Weight(p)*Kappa(p)); d*=Weight(p);
                w->num3[ab+l]+=n;w->num3[ba+l]+=n;
                w->den3[ab+l]+=d;w->den3[ba+l]+=d;
            }
        }
        for (size_t j=0;j<current_count;j++) {
            size_t a=current[j];
            if (!w->moment_seen[a]) {active[active_count++]=a; w->moment_seen[a]=1;}
            for (size_t k=0;k<h;k++) {tn[a*h+k]+=fn[a*h+k]; td[a*h+k]+=fd[a*h+k];}
        }
        w->ordered_triplet_count+=(INTEGER)(2*previous_count*(end-begin));
        previous_count+=end-begin; begin=end;
    }
    return SUCCESS;
}

/* Two-level, forest-disjoint moment hierarchy. The cost pass chooses a forest
 * block size from the ACTUAL radial/polar occupancy. Each leaf is one complete
 * forest; a block is contracted internally, then its moments are reused against
 * completed blocks. No total-minus-self subtraction or geometry approximation
 * is involved. Kernel 1 retains the independent legacy prefix implementation. */
static void lya_multipole_harmonics(const lya_neighbor *q,int limit,REAL *h,
                                    const REAL *den,const REAL *back)
{
    REAL z=lya_clamp(q->displacement[2]/q->radius,-1,1);
    REAL transverse=hypot(q->displacement[0],q->displacement[1]);
    REAL st=transverse/q->radius;
    REAL cp=transverse>0?q->displacement[0]/transverse:1;
    REAL sp=transverse>0?q->displacement[1]/transverse:0;
    REAL cm=1,sm=0,diagonal=1;
    for(int m=0;m<=limit;m++) {
        if(m) {REAL next=cm*cp-sm*sp;sm=sm*cp+cm*sp;cm=next;diagonal*=back[m*m+m]*st;}
        REAL prev=0,current=diagonal;
        for(int l=m;l<=limit;l++) {
            if(l>m) {size_t k=l*l+m;REAL next=((2.*l-1)*z*current-back[k]*prev)/den[k];prev=current;current=next;}
            if(!m) h[l*l]=current;
            else {h[l*l+2*m-1]=sqrt(2.)*current*cm;h[l*l+2*m]=sqrt(2.)*current*sm;}
        }
    }
}

/* Sparse clearing and precomputed output offsets keep the innermost loop a
 * pair of contiguous harmonic dot products. Signed moments need an explicit
 * touched flag: a zero denominator does not imply an untouched output bin. */
static void lya_multipole_cross(struct cmdline_data *cmd,bodyptr p,lya_worker_hist *w,
                                REAL *left,size_t *lb,size_t nl,
                                REAL *right,size_t *rb,size_t nr)
{
    size_t h=(size_t)(cmd->lya3LMax+1)*(cmd->lya3LMax+1);
    size_t b=(size_t)cmd->lya3RBins*cmd->lya3ThetaBins;
    REAL pw=Weight(p),pn=pw*Kappa(p);
    size_t *first=w->moment_bins+3*b,*second=first+b;
    for(size_t i=0;i<nl;i++) for(size_t j=0;j<nr;j++) {
        size_t a=lb[i],c=rb[j],ab=first[a]+second[c],ba=first[c]+second[a];
        REAL *an=left+a*h,*ad=an+b*h,*bn=right+c*h,*bd=bn+b*h;
        size_t orders=(size_t)cmd->lya3LMax+1;
        unsigned char *output_seen=w->moment_seen+2*b;
        if(!output_seen[ab/orders]) {
            output_seen[ab/orders]=1;
            for(size_t l=0;l<orders;l++) w->touched3[w->touched3_count++]=ab+l;
        }
        if(!output_seen[ba/orders]) {
            output_seen[ba/orders]=1;
            for(size_t l=0;l<orders;l++) w->touched3[w->touched3_count++]=ba+l;
        }
        for(int l=0;l<=cmd->lya3LMax;l++) {
            REAL n=0,d=0;
            for(int k=l*l;k<(l+1)*(l+1);k++) {n+=an[k]*bn[k];d+=ad[k]*bd[k];}
            n*=pn;d*=pw;
            w->num3[ab+l]+=n;w->num3[ba+l]+=n;
            w->den3[ab+l]+=d;w->den3[ba+l]+=d;
        }
    }
    w->multipole_counts[3]+=(uint64_t)nl*nr;
}

static void lya_multipole_merge(REAL *to,REAL *from,size_t *dest,size_t *nd,
                                 size_t *src,size_t ns,unsigned char *seen,
                                 size_t b,size_t h)
{
    for(size_t j=0;j<ns;j++) {
        size_t a=src[j],offset=a*h;
        if(!seen[a]) {
            seen[a]=1;dest[(*nd)++]=a;
            memcpy(to+offset,from+offset,h*sizeof(REAL));
            memcpy(to+b*h+offset,from+b*h+offset,h*sizeof(REAL));
        } else for(size_t k=0;k<h;k++) {
            to[offset+k]+=from[offset+k];to[b*h+offset+k]+=from[b*h+offset+k];
        }
    }
}

static size_t lya_multipole_block(struct cmdline_data *cmd,lya_worker_hist *w,size_t b)
{
    /* Estimating every candidate requires only integer occupancy updates,
     * independent of Lmax. A single prefix pass remains best for many sparse
     * inputs, and avoids introducing a hierarchy that cannot save products. */
    const size_t candidates[]={1,4,16,64};
    size_t best=1;long double best_cost=0,prefix_cost=0;
    for(size_t candidate=0;candidate<4;candidate++) {
        size_t block=candidates[candidate],total=0,group=0,forests=0;
        long double cost=0;
        unsigned char *seen=w->moment_seen,*group_seen=seen+b;
        size_t *bins=w->moment_bins;
        memset(seen,0,2*b);
        for(size_t begin=0;begin<w->neighbor_count;) {
            size_t end=begin,current=0;
            while(end<w->neighbor_count && w->neighbors[end].forest_id==w->neighbors[begin].forest_id) {
                size_t bin=w->neighbors[end].leg_bin;
                if(end==begin || bin!=w->neighbors[end-1].leg_bin) current++;
                end++;
            }
            cost+=(long double)current*group;
            for(size_t j=begin;j<end;j++) {
                size_t a=w->neighbors[j].leg_bin;
                if(!group_seen[a]) {group_seen[a]=1;bins[group++]=a;}
            }
            if(++forests%block==0 || end==w->neighbor_count) {
                cost+=(long double)group*total;
                for(size_t j=0;j<group;j++) {size_t a=bins[j];if(!seen[a]) {seen[a]=1;total++;}group_seen[a]=0;}
                group=0;
            }
            begin=end;
        }
        if(!candidate) {prefix_cost=best_cost=cost;}
        else if(cost<best_cost) {best_cost=cost;best=block;}
    }
    w->multipole_counts[4]+=(uint64_t)prefix_cost;
    /* Require a useful reduction to pay for the extra block merge. Kernel 2
     * forces four-forest groups for regression/qualification, even if slower. */
    if(cmd->lya3Kernel==2) return 4;
    return best_cost<.85L*prefix_cost?best:1;
}

static int lya_accumulate_multipoles(struct cmdline_data *cmd,bodyptr p,
                                     lya_worker_hist *w,ErrorMsg err)
{
    if(cmd->lya3Kernel==1) return lya_accumulate_multipoles_prefix(cmd,p,w,err);
    if(w->neighbor_count<2) return SUCCESS;
    if(lya_multipole_scratch(cmd,w,err)==FAILURE) return FAILURE;
    size_t b=(size_t)cmd->lya3RBins*cmd->lya3ThetaBins;
    size_t orders=cmd->lya3LMax+1,h=orders*orders;
    REAL *total=w->moments,*forest=total+2*b*h,*group=forest+2*b*h;
    REAL *basis=group+2*b*h,*den=basis+h,*back=den+h;
    size_t *active=w->moment_bins,*current=active+b,*group_bins=current+b;
    unsigned char *seen=w->moment_seen,*group_seen=seen+b;
    lya_sort_neighbors(w->neighbors,w->neighbor_count);
    w->neighbors_sorted=TRUE;
    size_t block=lya_multipole_block(cmd,w,b),active_count=0,group_count=0;
    size_t forests=0,previous_count=0;
    w->multipole_counts[block==1?1:0]++;
    memset(seen,0,2*b);
    for(size_t begin=0;begin<w->neighbor_count;) {
        size_t end=begin,current_count=0;
        while(end<w->neighbor_count && w->neighbors[end].forest_id==w->neighbors[begin].forest_id) {
            const lya_neighbor *q=&w->neighbors[end];size_t bin=q->leg_bin;
            if(end==begin || bin!=w->neighbors[end-1].leg_bin) {
                current[current_count++]=bin;
                memset(forest+bin*h,0,h*sizeof(REAL));memset(forest+b*h+bin*h,0,h*sizeof(REAL));
            }
            lya_multipole_harmonics(q,cmd->lya3LMax,basis,den,back);
            for(size_t k=0;k<h;k++) {forest[bin*h+k]+=q->weighted_delta*basis[k];forest[b*h+bin*h+k]+=q->weight*basis[k];}
            end++;
        }
        if(block==1) {
            lya_multipole_cross(cmd,p,w,forest,current,current_count,total,active,active_count);
            lya_multipole_merge(total,forest,active,&active_count,current,current_count,seen,b,h);
        } else {
            lya_multipole_cross(cmd,p,w,forest,current,current_count,group,group_bins,group_count);
            lya_multipole_merge(group,forest,group_bins,&group_count,current,current_count,group_seen,b,h);
            if((forests+1)%block==0 || end==w->neighbor_count) {
                lya_multipole_cross(cmd,p,w,group,group_bins,group_count,total,active,active_count);
                lya_multipole_merge(total,group,active,&active_count,group_bins,group_count,seen,b,h);
                for(size_t j=0;j<group_count;j++) group_seen[group_bins[j]]=0;
                group_count=0;w->multipole_counts[2]++;
            }
        }
        w->ordered_triplet_count+=(INTEGER)(2*previous_count*(end-begin));
        previous_count+=end-begin;begin=end;forests++;
    }
    return SUCCESS;
}

static REAL lya_window_coefficient(int l, REAL lower, REAL upper)
{
    if (!l) return (upper-lower)/2;
    REAL ends[2]={lower,upper}, primitive[2];
    for (int k=0;k<2;k++) {
        REAL prev=1,current=ends[k], minus=l==1 ? 1 : 0;
        for (int n=1;n<=l;n++) {
            REAL next=((2*n+1)*ends[k]*current-n*prev)/(n+1);
            if (n==l-1) minus=current;
            prev=current;current=next;
        }
        primitive[k]=(current-minus)/2;
    }
    return primitive[1]-primitive[0];
}

static int lya_write_multipoles(struct cmdline_data *cmd, struct global_data *gd,
                                 const REAL *num,const REAL *den,INTEGER count)
{
    char rawpath[MAXLENGTHOFFILES],path[MAXLENGTHOFFILES];
    if (format_checked(rawpath,sizeof(rawpath),"Ly-alpha moments","%s_lya_multipoles%s",gd->fpfnamehistZetaMFileName,EXTFILES)
        || format_checked(path,sizeof(path),"Ly-alpha reconstruction","%s_lya5d_multipole%s",gd->fpfnamehistZetaMFileName,EXTFILES)) return FAILURE;
    size_t cells=(size_t)cmd->lya3RBins*cmd->lya3RBins*cmd->lya3ThetaBins*cmd->lya3ThetaBins*(cmd->lya3LMax+1);
    /* dimensions were checked before allocating these arrays */
    for (size_t k=0;k<cells;k++) if (!isfinite(num[k]) || !isfinite(den[k])) {
        snprintf(cmd->error_message,_ERRORMSGSIZE_,"Ly-alpha multipole raw sums overflowed; rescale catalog weights/field");
        return FAILURE;
    }
    FILE *raw=fopen(rawpath,"w"),*out=fopen(path,"w");
    if (!raw || !out) { if(raw)fclose(raw); if(out)fclose(out);snprintf(cmd->error_message,_ERRORMSGSIZE_,"cannot open Ly-alpha multipole outputs: %s",strerror(errno));return FAILURE;}
    fprintf(raw,"# Exact anisotropic Legendre raw sums; lmax=%d; distinct-forest ordered triplets: %" INTEGER_FMT "\n",cmd->lya3LMax,count);
    fprintf(raw,"# columns: b1 b2 t1 t2 ell numerator denominator (signed moments; no ratio)\n");
    fprintf(out,"# APPROXIMATE top-hat reconstruction from lmax=%d; radial and polar bins exact\n",cmd->lya3LMax);
    if(cmd->lya3MuMode) fprintf(out,"# Exact companion: _lya5d; this file remains the finite-L diagnostic\n");
    fprintf(out,"# denominator<=0 gives NaN zeta; raw signed reconstruction is retained, never clipped\n");
    fprintf(out,"# columns: b1 b2 t1 t2 bmu r1 r2 theta1 theta2 mu zeta numerator denominator\n");
    for (int b1=0;b1<cmd->lya3RBins;b1++) for(int b2=0;b2<cmd->lya3RBins;b2++)
      for(int t1=0;t1<cmd->lya3ThetaBins;t1++) for(int t2=0;t2<cmd->lya3ThetaBins;t2++) {
        size_t index=lya_index3(b1,b2,t1,t2,0,cmd->lya3RBins,cmd->lya3ThetaBins,cmd->lya3LMax+1);
        for(int l=0;l<=cmd->lya3LMax;l++) fprintf(raw,"%d %d %d %d %d %.17g %.17g\n",b1,b2,t1,t2,l,num[index+l],den[index+l]);
        for(int m=0;m<cmd->lya3MuBins;m++) {
            REAL n=0,d=0,lo=-1+2.*m/cmd->lya3MuBins,hi=-1+2.*(m+1)/cmd->lya3MuBins;
            for(int l=0;l<=cmd->lya3LMax;l++) {REAL c=lya_window_coefficient(l,lo,hi);n+=c*num[index+l];d+=c*den[index+l];}
            fprintf(out,"%d %d %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                b1,b2,t1,t2,m,(b1+.5)*cmd->lya3RMax/cmd->lya3RBins,(b2+.5)*cmd->lya3RMax/cmd->lya3RBins,
                (t1+.5)*LYA_PI/cmd->lya3ThetaBins,(t2+.5)*LYA_PI/cmd->lya3ThetaBins,(lo+hi)/2,d>0?n/d:NAN,n,d);
        }
      }
    int failed=ferror(raw)||ferror(out);
    if (fclose(raw)) failed=1;
    if (fclose(out)) failed=1;
    if (failed) snprintf(cmd->error_message,_ERRORMSGSIZE_,"failed writing Ly-alpha multipole outputs");
    return failed?FAILURE:SUCCESS;
}
