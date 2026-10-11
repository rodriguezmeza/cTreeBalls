#ifndef CBALLS_LEGENDRE_WINDOW_SOLVER_H
#define CBALLS_LEGENDRE_WINDOW_SOLVER_H
/* Pure, bounded Legendre product solve. Coefficients multiply P_l(mu).
 * M_ij=(2i+1) sum_k (R_k/R_0) (i j k;0 0 0)^2.
 * Window modes through 2*solve_lmax are mandatory; higher modes cannot couple
 * within this finite signal space. Padding tests signal truncation separately.
 * Condition diagnostics use the standard Legendre basis, infinity norm.
 */
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#define CB_SURVEY_MAX_SOLVE 32
#define CB_SURVEY_MAX_WINDOW (2*CB_SURVEY_MAX_SOLVE)
enum { CB_SURVEY_EMPTY=0, CB_SURVEY_VALID=1, CB_SURVEY_SINGULAR=2,
       CB_SURVEY_ILL_CONDITIONED=3, CB_SURVEY_RESIDUAL=4,
       CB_SURVEY_NONFINITE=5 };

static inline double cb_survey_gaunt(int a,int b,int c)
{
    int sum=a+b+c, g=sum/2;
    if ((sum&1) || a+b<c || a+c<b || b+c<a) return 0.;
    return exp(2.*lgamma(g+1.)-2.*lgamma(g-a+1.)-2.*lgamma(g-b+1.)
        -2.*lgamma(g-c+1.)+lgamma(sum-2*a+1.)+lgamma(sum-2*b+1.)
        +lgamma(sum-2*c+1.)-lgamma(sum+2.));
}

/* x must hold solve_lmax+1 entries. Invalid outputs are NaN, never a valid zero.
 * rcond=0 for empty/singular support; residual=NaN when no solve was possible.
 * Return codes are per-bin scientific validity, not fatal application errors.
 */
static inline int cb_survey_solve(const double *num,const double *random,
    int solve_lmax,int window_lmax,double rcond_min,double *x,
    double *rcond,double *residual)
{
    const int n=solve_lmax+1;
    double matrix[(CB_SURVEY_MAX_SOLVE+1)*(CB_SURVEY_MAX_SOLVE+1)];
    double work[(CB_SURVEY_MAX_SOLVE+1)*(CB_SURVEY_MAX_SOLVE+1)];
    double inverse[(CB_SURVEY_MAX_SOLVE+1)*(CB_SURVEY_MAX_SOLVE+1)];
    double rhs[CB_SURVEY_MAX_SOLVE+1];
    double anorm=0., invnorm=0.;
    *rcond=0.; *residual=NAN;
    if (solve_lmax<0 || solve_lmax>CB_SURVEY_MAX_SOLVE || window_lmax<2*solve_lmax
        || window_lmax>CB_SURVEY_MAX_WINDOW || !isfinite(rcond_min)
        || rcond_min<=0. || rcond_min>1.) return CB_SURVEY_NONFINITE;
    for(int i=0;i<n;i++) x[i]=NAN;
    for(int k=0;k<=window_lmax;k++) if(!isfinite(random[k])) return CB_SURVEY_NONFINITE;
    for(int i=0;i<n;i++) if(!isfinite(num[i])) return CB_SURVEY_NONFINITE;
    if(random[0]==0.) return CB_SURVEY_EMPTY;
    if(random[0]<0.) return CB_SURVEY_NONFINITE;
    for(int i=0;i<n;i++) {
        rhs[i]=num[i]/random[0];
        if(!isfinite(rhs[i]) || (num[i]!=0. && rhs[i]==0.)) return CB_SURVEY_NONFINITE;
        double rowsum=0.;
        for(int j=0;j<n;j++) {
            long double value=0.;
            for(int k=abs(i-j);k<=i+j;k+=2)
                value+=(long double)random[k]/random[0]*cb_survey_gaunt(i,j,k);
            matrix[i*n+j]=work[i*n+j]=(double)((2*i+1)*value);
            if(!isfinite(matrix[i*n+j])) return CB_SURVEY_NONFINITE;
            rowsum+=fabs(matrix[i*n+j]);
            inverse[i*n+j]=i==j ? 1. : 0.;
        }
        if(rowsum>anorm) anorm=rowsum;
    }
    if(!isfinite(anorm) || anorm==0.) return CB_SURVEY_NONFINITE;
    for(int k=0;k<n;k++) {
        int pivot=k;
        for(int i=k+1;i<n;i++) if(fabs(work[i*n+k])>fabs(work[pivot*n+k])) pivot=i;
        if(fabs(work[pivot*n+k])<=DBL_EPSILON*n*anorm) return CB_SURVEY_SINGULAR;
        for(int j=0;j<n;j++) {
            double v=work[k*n+j];work[k*n+j]=work[pivot*n+j];work[pivot*n+j]=v;
            v=inverse[k*n+j];inverse[k*n+j]=inverse[pivot*n+j];inverse[pivot*n+j]=v;
        }
        double d=work[k*n+k];
        for(int j=0;j<n;j++) {work[k*n+j]/=d;inverse[k*n+j]/=d;}
        for(int i=0;i<n;i++) if(i!=k) {
            double f=work[i*n+k];
            for(int j=0;j<n;j++) {work[i*n+j]-=f*work[k*n+j];inverse[i*n+j]-=f*inverse[k*n+j];}
        }
    }
    for(int i=0;i<n;i++) {
        double rowsum=0.;
        for(int j=0;j<n;j++) rowsum+=fabs(inverse[i*n+j]);
        if(rowsum>invnorm) invnorm=rowsum;
    }
    *rcond=(double)(1.L/((long double)anorm*invnorm));
    if(!isfinite(*rcond) || *rcond<=0.) {*rcond=0.;return CB_SURVEY_SINGULAR;}
    if(*rcond<rcond_min) return CB_SURVEY_ILL_CONDITIONED;
    long double xnorm=0., bnorm=0., err=0.;
    for(int i=0;i<n;i++) {
        long double v=0.;
        for(int j=0;j<n;j++) v+=(long double)inverse[i*n+j]*rhs[j];
        x[i]=(double)v;
        if(!isfinite(x[i]) || (v!=0. && x[i]==0.)) goto nonfinite;
        if(fabsl(v)>xnorm) xnorm=fabsl(v);
        if(fabs(rhs[i])>bnorm) bnorm=fabs(rhs[i]);
    }
    for(int i=0;i<n;i++) {
        long double v=-rhs[i];
        for(int j=0;j<n;j++) v+=(long double)matrix[i*n+j]*x[j];
        if(fabsl(v)>err) err=fabsl(v);
    }
    long double scale=(long double)anorm*xnorm+bnorm;
    *residual=scale>0. ? (double)(err/scale) : 0.;
    if(!isfinite(*residual) || *residual>128.*n*DBL_EPSILON) {
        for(int i=0;i<n;i++) x[i]=NAN;
        return CB_SURVEY_RESIDUAL;
    }
    return CB_SURVEY_VALID;
nonfinite:
    for(int i=0;i<n;i++) x[i]=NAN;
    return CB_SURVEY_NONFINITE;
}
#endif
