/* Ordered distinct products. A merge combines disjoint sets A and B:
 * D(A+B)=D(A)+D(B)+S(A)conj(S(B))+S(B)conj(S(A)).
 * The diagonal is real; exchanging different radial bins conjugates the
 * complex product. No dominant auto product is formed and subtracted.
 * Inputs and results use the production double-precision accumulator type.
 */
#ifndef CBALLS_DISTINCT_PRODUCTS_H
#define CBALLS_DISTINCT_PRODUCTS_H
#include <math.h>
#include <float.h>
static inline double cballs_product(double a, double b)
{
    double p = a*b;
    /* Preserve representable subnormals, but never silently lose a nonzero
     * intermediate to zero. NaN is an internal error marker, checked before
     * publication; it is never a replacement scientific result. */
    return (a != 0.0 && b != 0.0 && p == 0.0) ? NAN : p;
}
/* Zero denominators denote an empty/all-zero-weight bin only if the
 * numerator is also zero. A nonzero quotient lost to zero is a range error. */
static inline double cballs_checked_ratio(double numerator, double denominator)
{
    if (denominator == 0.0) return numerator == 0.0 ? 0.0 : NAN;
    double value = numerator/denominator;
    return value == 0.0 && numerator != 0.0 ? NAN : value;
}
static inline double cballs_distinct_merge(double da, double db,
                                           double sa, double sb)
{
    return da + db + 2.0*cballs_product(sa, sb);
}
static inline double cballs_distinct_complex_cross(double ar, double ai,
                                                    double br, double bi)
{
    return 2.0*(cballs_product(ar, br)+cballs_product(ai, bi));
}
/* Compensated sum with a separate binary exponent: finite input factors
 * need not have representable intermediate products. */
typedef struct { double sum, correction; int exponent; } cballs_scaled_sum;
static inline void cballs_scaled_add(cballs_scaled_sum *s, const double *v, int n)
{
    double x=1.0; int exponent=0;
    for (int k=0;k<n;k++) {
        if (!isfinite(v[k])) { s->sum=NAN; return; }
        if (v[k]==0.0) return;
        int e; x*=frexp(v[k],&e); exponent+=e;
    }
    int e; x=frexp(x,&e); exponent+=e;
    if (s->sum==0.0 && s->correction==0.0) s->exponent=exponent;
    if (exponent>s->exponent) {
        s->sum=scalbn(s->sum,s->exponent-exponent);
        s->correction=scalbn(s->correction,s->exponent-exponent);
        s->exponent=exponent;
    }
    x=scalbn(x,exponent-s->exponent);
    const double total=s->sum+x;
    s->correction += fabs(s->sum)>=fabs(x) ? (s->sum-total)+x : (x-total)+s->sum;
    s->sum=total;
}
static inline double cballs_scaled_ratio(cballs_scaled_sum n, cballs_scaled_sum d)
{
    const double numerator=n.sum+n.correction, denominator=d.sum+d.correction;
    if (denominator==0.0) return numerator==0.0 ? 0.0 : NAN;
    double value=scalbn(numerator/denominator,n.exponent-d.exponent);
    return value==0.0 && numerator!=0.0 ? NAN : value;
}
#endif
