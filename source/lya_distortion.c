/* Continuum distortion, Lya2pcf (2507.00129), equations 2.6--2.7.
 * P = I - sum_a u_a (p*u_a)^T, with weighted constant/slope basis u.
 * Expanding P_f M_B P_g^T gives an identity term, two single-projection
 * terms and a rank <=4 product. This is algebraically the four-index sum;
 * inner/model pixels are never restricted by the outer/data pair window.
 */
#include "globaldefs.h"
#include "lya_distortion.h"
#include <float.h>
#include <limits.h>

typedef struct {
    double chi, los[3], w, p, lambda, u[2];
    int64_t forest;
    size_t ordinal;
} dm_pixel;
typedef struct {
    size_t first, end;
    int rank, active;
    double lo[3], hi[3], chi_lo, chi_hi;
} dm_forest;
typedef struct {
    double *matrix, *den, *a, *b, *ma, *mb;
    uint32_t *touched, *ai, *bi, *mai, *mbi;
    unsigned char *seen, *as, *bs, *mas, *mbs;
    size_t touched_count, ac, bc, mac, mbc;
    uint64_t counts[4];
    int failed;
} dm_worker;

static int dm_error(char *error, size_t size, const char *message)
{ if(error && size) snprintf(error,size,"Ly-alpha distortion: %s",message); return 1; }

int cballs_lya_distortion_available(void)
{
#if (defined(LYAFORESTOMP) || defined(LYAFORESTMPI)) && NDIM == 3 && !defined(__FAST_MATH__)
    return sizeof(REAL)==sizeof(double);
#else
    return 0;
#endif
}

static int dm_shape(const cballs_lya_distortion_config *c,
                    size_t *a, size_t *b, size_t *cells)
{
    return c && c->rp_bins>0 && c->rt_bins>0 && c->model_rp_bins>0 && c->model_rt_bins>0
        && c->projection_order>=-1 && c->projection_order<=1 && c->threads>=1 && c->threads<=1024
        && isfinite(c->rp_max) && c->rp_max>0 && isfinite(c->rt_max) && c->rt_max>0
        && isfinite(c->model_rp_max) && c->model_rp_max>0
        && isfinite(c->model_rt_max) && c->model_rt_max>0
        && cballs_size_mul(c->rp_bins,c->rt_bins,a) && *a<=INT_MAX
        && cballs_size_mul(c->model_rp_bins,c->model_rt_bins,b) && *b<=INT_MAX
        && cballs_size_mul(*a,*b,cells) && *cells<=UINT32_MAX;
}

int cballs_lya_distortion_plan(size_t n,const cballs_lya_distortion_config *c,
                               size_t *bytes,char *error,size_t size)
{
    size_t a,b,cells,work,total,extra,axes;
    if(!bytes || !n || !dm_shape(c,&a,&b,&cells))
        return dm_error(error,size,"invalid catalog, grids, projection order, threads or dimensions");
    /* Private matrix, marks, uint32 touched indices; scratch and row weights.
     * Reserve both raw/normalised caller matrices and all input channels even
     * when aliases/absent channels reduce actual storage. Forest count <= n. */
    if(!cballs_size_mul(cells,13,&work) || !cballs_size_add(a,b,&axes)
        || !cballs_size_mul(axes,58,&extra) || !cballs_size_add(work,extra,&work)
        || !cballs_size_mul(a,8,&extra) || !cballs_size_add(work,extra,&work)
        || !cballs_size_add(work,sizeof(dm_worker),&work)
        || !cballs_size_mul(work,c->threads,&total)
        || !cballs_size_mul(cells,16,&extra) || !cballs_size_add(total,extra,&total)
        || !cballs_size_mul(n,sizeof(dm_pixel)+sizeof(dm_forest)+56,&extra)
        || !cballs_size_add(total,extra,&total)
        || !cballs_size_mul(a,32,&extra) || !cballs_size_add(total,extra,&total)
        || !cballs_size_mul(axes,16,&extra) || !cballs_size_add(total,extra,&total))
        return dm_error(error,size,"workspace dimensions overflow");
    *bytes=total;
    if(!c->max_bytes || total>c->max_bytes)
        return dm_error(error,size,"workspace exceeds max_memory_mb; reduce bins/threads or increase the budget");
    if(cballs_memory_preflight(total,"Ly-alpha distortion workspace",error,size)==FAILURE) return 1;
    return 0;
}

static int dm_compare(const void *aa,const void *bb)
{
    const dm_pixel *a=aa,*b=bb;
    if(a->forest!=b->forest) return a->forest<b->forest?-1:1;
    return a->ordinal<b->ordinal?-1:a->ordinal>b->ordinal;
}

static int dm_bin(double value,double maximum,int bins)
{
    if(!(value>=0) || value>=maximum) return -1;
    int bin=(int)(value/maximum*(double)bins);
    return bin>=0 && bin<bins?bin:-1;
}

static void dm_bins(const dm_pixel *p,const dm_pixel *q,
                     const cballs_lya_distortion_config *c,int *a,int *b)
{
    double mu=p->los[0]*q->los[0]+p->los[1]*q->los[1]+p->los[2]*q->los[2];
    mu=fmax(-1.,fmin(1.,mu));
    double rp=fabs(p->chi-q->chi)*sqrt(fmax(0.,.5*(1.+mu)));
    double rt=(p->chi+q->chi)*sqrt(fmax(0.,.5*(1.-mu)));
    int ap=dm_bin(rp,c->rp_max,c->rp_bins),at=dm_bin(rt,c->rt_max,c->rt_bins);
    int bp=dm_bin(rp,c->model_rp_max,c->model_rp_bins),bt=dm_bin(rt,c->model_rt_max,c->model_rt_bins);
    *a=ap<0 || at<0?-1:ap*c->rt_bins+at;
    *b=bp<0 || bt<0?-1:bp*c->model_rt_bins+bt;
}

static int dm_prepare(dm_pixel *p,dm_forest *f,int order,uint64_t *degenerate)
{
    double maximum=0;
    f->active=0; f->rank=0;
    f->chi_lo=INFINITY;f->chi_hi=0;
    for(int k=0;k<3;k++) {f->lo[k]=INFINITY;f->hi[k]=-INFINITY;}
    for(size_t i=f->first;i<f->end;i++) {
        maximum=fmax(maximum,p[i].p);f->active|=p[i].w>0;
        f->chi_lo=fmin(f->chi_lo,p[i].chi);f->chi_hi=fmax(f->chi_hi,p[i].chi);
        for(int k=0;k<3;k++) {f->lo[k]=fmin(f->lo[k],p[i].los[k]);f->hi[k]=fmax(f->hi[k],p[i].los[k]);}
    }
    if(order<0 || !f->active) return 0;
    if(!(maximum>0)) return 1;
    size_t reference=f->first;
    while(p[reference].p==0) reference++;
    double origin=p[reference].lambda;
    long double sum=0,mean=0,variance=0;
    for(size_t i=f->first;i<f->end;i++) {
        if(p[i].p>0 && p[i].p/maximum==0) return 1;
        p[i].p/=maximum;
        sum+=p[i].p;
        mean+=(long double)p[i].p*((long double)p[i].lambda-origin);
    }
    mean/=sum;
    if(!isfinite(sum) || !isfinite(mean)) return 1;
    double scale=0;
    for(size_t i=f->first;i<f->end;i++) {
        p[i].u[1]=(double)(((long double)p[i].lambda-origin)-mean);
        if(p[i].p>0) scale=fmax(scale,fabs(p[i].u[1]));
    }
    if(!isfinite(scale)) return 1;
    if(scale>0 && isfinite(scale)) for(size_t i=f->first;i<f->end;i++) {
        p[i].u[1]/=scale;
        variance+=(long double)p[i].p*p[i].u[1]*p[i].u[1];
    }
    f->rank=1;
    if(order==1 && variance>0 && isfinite(variance)) f->rank=2;
    else if(order==1) (*degenerate)++;
    double constant=1/sqrt((double)sum),slope=f->rank==2?1/sqrt((double)variance):0;
    for(size_t i=f->first;i<f->end;i++) {
        p[i].u[0]=constant;
        p[i].u[1]=f->rank==2?p[i].u[1]*slope:0;
        if(!isfinite(p[i].u[1]) || !isfinite(p[i].u[0])) return 1;
    }
    return 0;
}

static int dm_prune(const dm_forest *a,const dm_forest *b,
                     const cballs_lya_distortion_config *c)
{
    if(!a->active || !b->active) return 1;
    double lo=0,hi=0;
    for(int k=0;k<3;k++) {
        double x=a->lo[k]*b->lo[k],y=a->lo[k]*b->hi[k];
        double z=a->hi[k]*b->lo[k],t=a->hi[k]*b->hi[k];
        lo+=fmin(fmin(x,y),fmin(z,t));hi+=fmax(fmax(x,y),fmax(z,t));
    }
    lo=fmax(-1.,lo-1024*DBL_EPSILON);hi=fmin(1.,hi+1024*DBL_EPSILON);
    double gap=fmax(0.,fmax(a->chi_lo-b->chi_hi,b->chi_lo-a->chi_hi));
    double rp=gap*sqrt(fmax(0.,.5*(1.+lo)));
    double rt=(a->chi_lo+b->chi_lo)*sqrt(fmax(0.,.5*(1.-hi)));
    double guard=4096*DBL_EPSILON*(a->chi_hi+b->chi_hi+c->rp_max+c->rt_max+1);
    return isfinite(guard) && (rp>c->rp_max+guard || rt>c->rt_max+guard);
}

static void dm_add(dm_worker *w,size_t index,double value)
{
    if(value==0) return;
    if(!w->seen[index]) {w->seen[index]=1;w->touched[w->touched_count++]=(uint32_t)index;}
    w->matrix[index]+=value;
}

static void dm_mark(uint32_t *indices,unsigned char *seen,size_t *count,int bin)
{ if(!seen[bin]) {seen[bin]=1;indices[(*count)++]=(uint32_t)bin;} }

static void dm_clear_legs(dm_worker *w)
{
    for(size_t j=0;j<w->ac;j++) {size_t a=w->ai[j];w->as[a]=0;w->a[2*a]=w->a[2*a+1]=0;}
    for(size_t j=0;j<w->bc;j++) {size_t b=w->bi[j];w->bs[b]=0;w->b[2*b]=w->b[2*b+1]=0;}
    w->ac=w->bc=0;
}

static void dm_pair(dm_worker *w,const dm_pixel *p,const dm_forest *f,
                     const dm_forest *g,const cballs_lya_distortion_config *c,size_t nb)
{
    int has_data=0;
    for(size_t j=g->first;j<g->end;j++) {
        for(size_t i=f->first;i<f->end;i++) {
            int a,b;dm_bins(p+i,p+j,c,&a,&b);
            if(a>=0 && p[i].w>0 && p[j].w>0) {
                double weight=p[i].w*p[j].w;
                if(!isfinite(weight)) {w->failed=1;return;}
                has_data=1;w->counts[3]++;w->den[a]+=weight;
                if(b>=0) dm_add(w,(size_t)a*nb+b,weight);
                dm_mark(w->mai,w->mas,&w->mac,a);
                for(int x=0;x<f->rank;x++) for(int y=0;y<g->rank;y++)
                    w->ma[4*(size_t)a+2*x+y]+=weight*p[i].u[x]*p[j].u[y];
            }
            if(b>=0 && p[i].p>0 && p[j].p>0) {
                dm_mark(w->mbi,w->mbs,&w->mbc,b);
                for(int x=0;x<f->rank;x++) for(int y=0;y<g->rank;y++)
                    w->mb[4*(size_t)b+2*x+y]+=p[i].p*p[j].p*p[i].u[x]*p[j].u[y];
            }
            if(a>=0 && p[i].w>0) {
                dm_mark(w->ai,w->as,&w->ac,a);
                for(int x=0;x<f->rank;x++) w->a[2*(size_t)a+x]+=p[i].w*p[i].u[x];
            }
            if(b>=0 && p[i].p>0) {
                dm_mark(w->bi,w->bs,&w->bc,b);
                for(int x=0;x<f->rank;x++) w->b[2*(size_t)b+x]+=p[i].p*p[i].u[x];
            }
        }
        if(p[j].w>0) for(size_t ia=0;ia<w->ac;ia++) for(size_t ib=0;ib<w->bc;ib++) {
            size_t a=w->ai[ia],b=w->bi[ib];double dot=0;
            for(int x=0;x<f->rank;x++) dot+=w->a[2*a+x]*w->b[2*b+x];
            dm_add(w,a*nb+b,-p[j].w*dot);
        }
        dm_clear_legs(w);
    }
    if(has_data) {
        w->counts[2]++;
        for(size_t i=f->first;i<f->end;i++) {
            for(size_t j=g->first;j<g->end;j++) {
                int a,b;dm_bins(p+i,p+j,c,&a,&b);
                if(a>=0 && p[j].w>0) {
                    dm_mark(w->ai,w->as,&w->ac,a);
                    for(int y=0;y<g->rank;y++) w->a[2*(size_t)a+y]+=p[j].w*p[j].u[y];
                }
                if(b>=0 && p[j].p>0) {
                    dm_mark(w->bi,w->bs,&w->bc,b);
                    for(int y=0;y<g->rank;y++) w->b[2*(size_t)b+y]+=p[j].p*p[j].u[y];
                }
            }
            if(p[i].w>0) for(size_t ia=0;ia<w->ac;ia++) for(size_t ib=0;ib<w->bc;ib++) {
                size_t a=w->ai[ia],b=w->bi[ib];double dot=0;
                for(int y=0;y<g->rank;y++) dot+=w->a[2*a+y]*w->b[2*b+y];
                dm_add(w,a*nb+b,-p[i].w*dot);
            }
            dm_clear_legs(w);
        }
        for(size_t ia=0;ia<w->mac;ia++) for(size_t ib=0;ib<w->mbc;ib++) {
            size_t a=w->mai[ia],b=w->mbi[ib];double dot=0;
            for(int x=0;x<f->rank;x++) for(int y=0;y<g->rank;y++)
                dot+=w->ma[4*a+2*x+y]*w->mb[4*b+2*x+y];
            dm_add(w,a*nb+b,dot);
        }
    }
    for(size_t j=0;j<w->mac;j++) {size_t a=w->mai[j];w->mas[a]=0;memset(w->ma+4*a,0,4*sizeof(double));}
    for(size_t j=0;j<w->mbc;j++) {size_t b=w->mbi[j];w->mbs[b]=0;memset(w->mb+4*b,0,4*sizeof(double));}
    w->mac=w->mbc=0;
}

static void dm_free(dm_worker *w)
{
    free(w->matrix);free(w->den);free(w->a);free(w->b);free(w->ma);free(w->mb);
    free(w->touched);free(w->ai);free(w->bi);free(w->mai);free(w->mbi);
    free(w->seen);free(w->as);free(w->bs);free(w->mas);free(w->mbs);
}

int cballs_lya_distortion_compute(const double *pos,const double *weights,
    const int64_t *ids,const double *lambda,const double *pw,size_t n,
    const cballs_lya_distortion_config *c,double *raw,size_t raw_length,
    double *den,size_t den_length,uint64_t stats[8],char *error,size_t size)
{
    size_t na,nb,cells,plan,nf=0,pairs;
    dm_pixel *p=NULL;dm_forest *forests=NULL;dm_worker *workers=NULL;
    int status=1;
    uint64_t degenerate=0;
    if(!cballs_lya_distortion_available()) return dm_error(error,size,"requires a strict double 3D Ly-alpha build");
    if(cballs_lya_distortion_plan(n,c,&plan,error,size)) return 1;
    if(!dm_shape(c,&na,&nb,&cells) || !pos || !weights || !ids || !raw || !den || !stats
        || raw_length!=cells || den_length!=na || (c->projection_order==1 && !lambda))
        return dm_error(error,size,"missing arrays, log_lambda, or incompatible output lengths");
#define DM_ALLOC(pointer,count) do { \
    if(cballs_calloc_checked((void**)&(pointer),(count),sizeof(*(pointer)), \
        "Ly-alpha distortion",error,size)==FAILURE) goto cleanup; \
} while(0)
    DM_ALLOC(p,n);DM_ALLOC(forests,n);DM_ALLOC(workers,(size_t)c->threads);
    for(size_t i=0;i<n;i++) {
        double r2=0;
        for(int k=0;k<3;k++) {if(!isfinite(pos[3*i+k])) {dm_error(error,size,"non-finite position");goto cleanup;}r2+=pos[3*i+k]*pos[3*i+k];}
        p[i].chi=sqrt(r2);p[i].w=weights[i];p[i].p=pw?pw[i]:weights[i];
        p[i].lambda=lambda?lambda[i]:0;p[i].forest=ids[i];p[i].ordinal=i;
        if(!isfinite(p[i].chi) || p[i].chi<=0 || !isfinite(p[i].w) || p[i].w<0
            || !isfinite(p[i].p) || p[i].p<0 || !isfinite(p[i].lambda)) {
            dm_error(error,size,"invalid observer distance, weight or projection coordinate");goto cleanup;
        }
        for(int k=0;k<3;k++) p[i].los[k]=pos[3*i+k]/p[i].chi;
    }
    qsort(p,n,sizeof(*p),dm_compare);
    for(size_t first=0;first<n;) {
        size_t end=first+1;while(end<n && p[end].forest==p[first].forest) end++;
        dm_forest *f=forests+nf++;f->first=first;f->end=end;
        if(dm_prepare(p,f,c->projection_order,&degenerate)) {
            dm_error(error,size,"zero projection-weight sum, unrepresentable weight range or non-finite projection basis");goto cleanup;
        }
        first=end;
    }
    if(!cballs_size_mul(nf,nf-1,&pairs)) {dm_error(error,size,"forest-pair count overflow");goto cleanup;}
    memset(raw,0,cells*sizeof(double));memset(den,0,na*sizeof(double));memset(stats,0,8*sizeof(uint64_t));
    stats[0]=n;stats[1]=nf;stats[2]=pairs/2;stats[6]=degenerate;stats[7]=plan;
    for(int t=0;t<c->threads;t++) {
        dm_worker *w=workers+t;
        DM_ALLOC(w->matrix,cells);DM_ALLOC(w->seen,cells);DM_ALLOC(w->touched,cells);DM_ALLOC(w->den,na);
        DM_ALLOC(w->a,2*na);DM_ALLOC(w->b,2*nb);DM_ALLOC(w->ma,4*na);DM_ALLOC(w->mb,4*nb);
        DM_ALLOC(w->ai,na);DM_ALLOC(w->bi,nb);DM_ALLOC(w->mai,na);DM_ALLOC(w->mbi,nb);
        DM_ALLOC(w->as,na);DM_ALLOC(w->bs,nb);DM_ALLOC(w->mas,na);DM_ALLOC(w->mbs,nb);
    }
    status=0;
    /* One fixed primary-forest task, committed in forest order. Scheduling and
     * worker assignment cannot change floating-point addition order. */
#ifdef OPENMPCODE
#pragma omp parallel num_threads(c->threads)
#endif
    {
        int thread=0;
#ifdef OPENMPCODE
        thread=omp_get_thread_num();
#pragma omp for ordered schedule(dynamic,1)
#endif
        for(size_t f=0;f<nf;f++) {
            dm_worker *w=workers+thread;
            memset(w->counts,0,sizeof(w->counts));
            for(size_t g=f+1;g<nf && !w->failed;g++) {
                if(dm_prune(forests+f,forests+g,c)) {w->counts[1]++;continue;}
                dm_pair(w,p,forests+f,forests+g,c,nb);
            }
#ifdef OPENMPCODE
#pragma omp ordered
#endif
            {
                if(w->failed) {status=1;dm_error(error,size,"non-finite pair product");}
                for(size_t j=0;j<w->touched_count;j++) {
                    size_t k=w->touched[j];raw[k]+=w->matrix[k];
                    w->matrix[k]=0;w->seen[k]=0;
                }
                w->touched_count=0;
                for(size_t a=0;a<na;a++) {den[a]+=w->den[a];w->den[a]=0;}
                stats[3]+=w->counts[1];stats[4]+=w->counts[2];stats[5]+=w->counts[3];
            }
        }
    }
    if(!status) {
        for(size_t i=0;i<cells;i++) if(!isfinite(raw[i])) {status=dm_error(error,size,"matrix raw sums overflowed");break;}
        for(size_t i=0;i<na;i++) if(!isfinite(den[i])) {status=dm_error(error,size,"row weights overflowed");break;}
    }
cleanup:
    if(workers) for(int t=0;t<c->threads;t++) dm_free(workers+t);
    free(workers);free(forests);free(p);
    return status;
#undef DM_ALLOC
}
