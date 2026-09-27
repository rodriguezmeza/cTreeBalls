#include <math.h>
#include <float.h>
#include <assert.h>
typedef double REAL;
#define racos acos
#define LYA_PI 3.141592653589793238462643383279502884
static int lya_bin_theta(REAL theta,int bins)
{
    if(theta<0)theta=0;if(theta>LYA_PI)theta=LYA_PI;
    int bin=(int)(theta/LYA_PI*bins);return bin==bins?bins-1:bin;
}
#include "../addons/lya_forest_omp/lya_neighbor_geometry.h"
static void check(double value,int bins,const double *edges)
{
    if(value>=-1 && value<=1)
        assert(lya_polar_bin(value,bins,edges)==lya_bin_theta(acos(value),bins));
}
int main(void)
{
    for(int bins=1;bins<=64;bins++) {
        double edges[65];
        for(int k=1;k<bins;k++)edges[k]=cos(LYA_PI*k/bins);
        for(int k=0;k<=32768;k++)check(-1.+2.*k/32768,bins,edges);
        for(int k=1;k<bins;k++) {
            check(edges[k],bins,edges);
            double left=edges[k],right=edges[k];
            for(int j=0;j<16;j++) {
                left=nextafter(left,-INFINITY);right=nextafter(right,INFINITY);
                check(left,bins,edges);check(right,bins,edges);
            }
            for(int j=-1024;j<=1024;j++)check(edges[k]+j*DBL_EPSILON,bins,edges);
        }
    }
    return 0;
}
