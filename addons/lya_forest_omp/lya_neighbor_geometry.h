/* Compare with cosine edges instead of evaluating acos for each neighbor.
 * Only the double, strict-math profile enables this path, with <=64 polar
 * bins. The outward boundary guard covers cos/edge construction and the
 * reference acos/divide/multiply rounding; ambiguous values use the retained
 * reference arithmetic. This is a bin lookup, not an angular approximation. */
static int lya_polar_bin(REAL cosine,int bins,const REAL *edges)
{
    int lower=0,upper=bins;
    if(!isfinite(cosine)) return lya_bin_theta(racos(cosine),bins);
    while(upper-lower>1) {
        int middle=lower+(upper-lower)/2;
        REAL delta=cosine-edges[middle];
        if(fabs(delta)<=256*DBL_EPSILON)
            return lya_bin_theta(racos(cosine),bins);
        if(delta>0) upper=middle;else lower=middle;
    }
    return lower;
}
