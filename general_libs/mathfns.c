//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include <math.h>
#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"

#include <stdlib.h>


real rsqr(real x)
{
    return (x * x);
}

real rqbe(real x)
{
    return (x * x * x);
}

real rdex(real x)
{
    return (rexp(M_LN10 * x));
}

#if defined(SINGLEPREC)


float fcbrt(float x)
{
    return ((float) cbrt((double) x));
}

#endif

#ifndef USEGSL
void xsrandom(long idum)
{
#ifdef RANDOMZENOTYPE
	srandom(idum);					                    // Init randoms
#endif
}
#endif

#ifdef USEGSL
double xrandom(double xl, double xh)
{
    double ran2 = gsl_rng_uniform (r_gsl);
    return ( xl + (xh - xl) * ( ran2 ) );

}
#else
double xrandom(double xl, double xh)
{
#ifdef RANDOMZENOTYPE
    // use srandom(seed) to init
    return (xl + (xh - xl) * ((double) random()) / 2147483647.0);
#else
// idum is a long int, must not be changed between calls
    return ( xl + (xh - xl) * ((double) ran2(&idum)) ); // idum < 0
                                                        // period 2x10^18
#endif
//    return ( xl + (xh - xl) * ((double) rand()) / ((double)(RAND_MAX+1)) );
//    return ( xl + (xh - xl) * ((double) ran0(&idum)) ); // idum any int
                                                        // period 2.1x10^9
//    return ( xl + (xh - xl) * ((double) ran1(&idum)) ); // idum < 0
// Random generator by D. Knuth:
//    return ( xl + (xh - xl) * ((double) ran3(&idum)) ); // idum < 0
// ran4 not used, in case it is needed look for better option
//    return ( xl + (xh - xl) * ((double) ran4(&idum)) ); // idum >< 0
}
#endif

double grandom(double mean, double sdev)
{
    double v1, v2, s;

    do {
        v1 = xrandom(-1.0, 1.0);
        v2 = xrandom(-1.0, 1.0);
        s = v1*v1 + v2*v2;
    } while (s >= 1.0);
    return (mean + sdev * v1 * sqrt(-2.0 * log(s) / s));
}

void pickshell(real vec[], int ndim, real rad)
{
    real rsq, rscale;
    int i;

    do {
        rsq = 0.0;
        for (i = 0; i < ndim; i++) {
            vec[i] = xrandom(-1.0, 1.0);
            rsq = rsq + vec[i] * vec[i];
        }
    } while (rsq > 1.0);
    rscale = rad / rsqrt(rsq);
    for (i = 0; i < ndim; i++)
        vec[i] = vec[i] * rscale;
}

void pickball(real vec[], int ndim, real rad) 
{
    real rsq;
    int i;

    do {
        rsq = 0.0;
        for (i = 0; i < ndim; i++) {
            vec[i] = xrandom(-1.0, 1.0);
            rsq = rsq + vec[i] * vec[i];
        }
    } while (rsq > 1.0);
    for (i = 0; i < ndim; i++)
        vec[i] = vec[i] * rad;
}

void pickbox(real vec[], int ndim, real size) 
{
    int i;

    for (i = 0; i < ndim; i++)
        vec[i] = xrandom(- size, size);
}

#if NDIM == 2

void VRand(vector p)
{
  real s;

  s = 2. * PI * xrandom(0.,1.);
  p[0] = rcos (s);
  p[1] = rsin (s);
}

#elif NDIM == 3

void VRand(vector p)
{
  real s, x, y;

  s = 2.;
  while (s > 1.) {
    x = xrandom(-1.,1.);
    y = xrandom(-1.,1.);
    s = rsqr(x) + rsqr(y);
  }
  p[2] = 1. - 2. * s;
  s = 2. * rsqrt(1. - s);
  p[0] = s * x;
  p[1] = s * y;
}

#endif

// Entero mas cercano (Aumentar el rango de enteros mas que 15 digitos)
int nint(real x) 
{
if (x<0.) return (-(int) (0.5 -x));
else return((int) (0.5+x));
}

// Version previa 'X' estaba en minuscula ...
void bessel(real X,realptr BI0,realptr BI1,realptr BK0,realptr BK1)
{
//   SUBROUTINE TO COMPUTE THE VALUE OF THE MODIFIED BESSEL FUNCTIONS 
//   I0,I1,K0 AND K1 AT X.

	real T,SQX,CT1,CT2,T2,HX,HX2,A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,
			B6,B7,B8,B9,C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,D5,D6,D7,D8,D9,
			E1,E2,E3,E4,E5,E6,E7,F1,F2,F3,F4,F5,F6,F7,G1,G2,G3,G4,G5,
			G6,H1,H2,H3,H4,H5,H6,H7,BI0tmp,BI1tmp,BK0tmp,BK1tmp,TI,HXI;
//	real X;				// Version previa 'X' estaba activada...

	T=X/3.75;
	SQX=rsqrt(X);
	CT1=SQX*rexp(-X);
	CT2=SQX*rexp(X);
	T2=T*T;
	HX=X/2.;
	HX2=HX*HX;

	A1=3.5156229;
	A2=3.0899424;
	A3=1.2067492;
	A4=0.2659732;
	A5=0.0360768;
	A6=0.0045813;
 
	B1=0.39894228;
	B2=0.01328592;
	B3=0.00225319;
	B4=-0.00157565;
	B5=0.00916281;
	B6=-0.02057706;
	B7=0.02635537;
	B8=-0.01647633;
	B9=0.00392377;
 
	C1=0.87890594;
	C2=0.51498869;
	C3=0.15084934;
	C4=0.02658733;
	C5=0.00301532;
	C6=0.00032411;
 
	D1=0.39894228;
	D2=-0.03988024;
	D3=-0.00362018;
	D4=0.00163801;
	D5=-0.01031555;
	D6=0.02282967;
	D7=-0.02895312;
	D8=0.01787654;
	D9=-0.00420059;
 
	E1=-0.57721566;
	E2=0.42278420;
	E3=0.23069756;
	E4=0.03488590;
	E5=0.00262698;
	E6=0.00010750;
	E7=0.00000740;
 
	F1=1.25331414;
	F2=-0.07832358;
	F3=0.02189568;
	F4=-0.01062446;
	F5=0.00587872;
	F6=-0.00251540;
	F7=0.00053208;
 
	G1=0.15443144;
	G2=-0.67278579;
	G3=-0.18156897;
	G4=-0.01919402;
	G5=-0.00110404;
	G6=-0.00004686;
 
	H1=1.25331414;
	H2=0.23498619;
	H3=-0.03655620;
	H4=0.01504268;
	H5=-0.00780353;
	H6=0.00325614;
	H7=-0.00068245;
 
	if (X <= 3.75) {
		BI0tmp=1.+T2*(A1+T2*(A2+T2*(A3+T2*(A4+T2*(A5+T2*A6)))));
		*BI0=BI0tmp;
		BI1tmp=0.5+T2*(C1+T2*(C2+T2*(C3+T2*(C4+T2*(C5+T2*C6)))));
		BI1tmp=BI1tmp*X;
		*BI1=BI1tmp;
	} else {
		TI=1./T;
		BI0tmp=B1+TI*(B2+TI*(B3+TI*(B4+TI*(B5+TI*(B6+TI*
				(B7+TI*(B8+TI*B9)))))));
		BI0tmp=BI0tmp/CT1;
		*BI0=BI0tmp;
		BI1tmp=D1+TI*(D2+TI*(D3+TI*(D4+TI*(D5+TI*(D6+TI*
			(D7+TI*(D8+TI*D9)))))));
		BI1tmp=BI1tmp/CT1;
		*BI1=BI1tmp;
	}
 
	if (X <= 2.) {						// Version previa 'X' estaba en minuscula ...
		*BK0=-rlog(HX)*BI0tmp+E1+HX2*
			(E2+HX2*(E3+HX2*(E4+HX2*(E5+HX2*(E6+HX2*E7)))));
		BK1tmp=X*rlog(HX)*BI1tmp+1.+HX2*(G1+HX2*(G2+HX2*(G3+HX2*
			(G4+HX2*(G5+HX2*G6)))));
		*BK1=BK1tmp/X;
	} else {
		HXI=2./X;
		BK0tmp=F1+HXI*(F2+HXI*(F3+HXI*(F4+HXI*(F5+HXI*(F6+HXI*F7)))));
		*BK0=BK0tmp/CT2;
		BK1tmp=H1+HXI*(H2+HXI*(H3+HXI*(H4+HXI*(H5+HXI*(H6+HXI*H7)))));
		*BK1=BK1tmp/CT2;
	}
}

real erfcc(real x)
{
/*
Function to compute the complementary error function erfc(x),
with fractional error everywhere less than 1.2x10^-7.  This
algorithm uses a Chebyshev fitting. (cf: Numerical Recipes
p.164).
*/

	real z,t, erfcctmp;
	
	z = rabs(x);
	t = 1./(1.+0.5*z);
	erfcctmp = t*rexp( -z*z -1.26551223 + t*(1.00002368 + t*
				( 0.37409196 + t*(0.09678418 + t*(-0.18628806 + t*
				(0.27886807 + t*(-1.13520398 + t*(1.48851587 + t*
				(-0.82215223 + t* 0.17087277   )))))))));

	if (x < 0.0) erfcctmp = 2.0 - erfcctmp;
	return erfcctmp;
}

real expint(real x)
{
/*
Function to compute the exponential integral function E1(x),
with fractional error everywhere less than 2x10^-7.  This
algorithm is adapted from Abramowitz and Stegun, sec 5.1.
Consider using Numerical Recipes' version more general ...
*/

	real aux1, arg, aux2, aux3, expinttmp;

	arg=x;

	if (x <= 0.0) error("\narg error in expint\n");
 
	if (x <= 1.0) {
		aux1=-rlog(arg)-0.57721566+0.99999193*arg-0.24991055*rsqr(arg)+
			0.05519968*rqbe(arg)-0.00976004*rpow(arg,4)+0.00107857*rpow(arg,5);
		expinttmp=aux1;
     } else {
		aux1=rpow(arg,4)+8.5733287401*rqbe(arg)+18.0590169730*rsqr(arg)+
			8.6347608925*arg+0.2677737343;
		aux2=rpow(arg,4)+9.5733223454*rqbe(arg)+25.6329561486*rsqr(arg)+
			21.0996530827*arg+3.9584969228;
		aux3=aux1/(aux2*rexp(arg)*arg);
		expinttmp=aux3;
	}
	return expinttmp;
}

// Bessel functions from Numerical Recipes ...
// BESSI0: modified Bessel function I0(x).

double bessi0(real x)
{
    double t, tt, ti, u;
    
    t = ABS(x) / 3.75;
    tt = t * t;
    if (tt < 1.0) {
        u = 1.0 + tt * (3.5156229 + tt * (3.0899424 +
                                          tt * (1.2067492 + tt * (0.2659732 + tt * (0.0360768 +
                                                                                    tt * 0.0045813)))));
        return (u);
    } else {
        ti = 1.0 / t;
        u = 0.39894228 + ti * (0.01328592 + ti * (0.00225319 +
                                                  ti * (-0.00157565 + ti * (0.00916281 + ti * (-0.02057706 +
                                                                                               ti * (0.02635537 + ti * (-0.01647633 +
                                                                                                                        ti * 0.00392377)))))));
        return (u * exp(ABS(x)) / sqrt(ABS(x)));
    }
}

// BESSI1: modified Bessel function I1(x). 

real bessi1(real x)
{
    real t, tt, ti, u;
    
    t = x / 3.75;
    tt = t * t;
    if (tt < 1.0) {
        u = 0.5 + tt * (0.87890594 + tt * (0.51498869 +
                                           tt * (0.15084934 + tt * (0.02658733 + tt * (0.00301532 +
                                                                                       tt * 0.00032411)))));
        return (u * x);
    } else {
        if (t < 0.0)
            error("bessi1: invalid for x < -3.75\n");
        ti = 1.0 / t;
        u = 0.39894228 + ti * (-0.03988024 + ti * (-0.00362018 + 
                                                   ti * (0.00163801 + ti * (-0.01031555 + ti * (0.02282967 +
                                                                                                ti * (-0.02895312 + ti * (0.01787654 +
                                                                                                                          ti * -0.00420059)))))));
        return (u * exp(ABS(x)) / sqrt(ABS(x)));
    }
}

// BESSK0: modified Bessel function K0(x). 

real bessk0(real x)
{
    real t, tt, ti, u;
    
    if (x < 0.0)
        error("bessk0: negative argument\n");
    t = x / 2.0;
    if (t < 1.0) {
        tt = t * t;
        u = -0.57721566 + tt * (0.42278420 + tt * (0.23069756 +
                                                   tt * (0.03488590 + tt * (0.00262698 + tt * (0.00010750 +
                                                                                               tt * 0.00000740)))));
        return (u - log(t) * bessi0(x));
    } else {
        ti = 1.0 / t;
        u = 1.25331414 + ti * (-0.07832358 + ti * (0.02189568 +
                                                   ti * (-0.01062446 + ti * (0.00587872 + ti * (-0.00251540 +
                                                                                                ti * 0.00053208)))));
        return (u * exp(- x) / sqrt(x));
    }
}

// BESSJ0
double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
    
	if ((ax=rabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)));
		ans=rsqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

// BESSJ1
double bessj1(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
    
	if ((ax=rabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
		ans=rsqrt(0.636619772/ax)*(rcos(xx)*ans1-z*rsin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

// BESSJ
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double bessj(int n, double x)
{
//	double bessj0(double x);
//	double bessj1(double x);
	void nrerror(char error_text[]);
	int j,jsum,m;
	double ax,bj,bjm,bjp,sum,tox,ans;

//    if (n < 2) nrerror("Index n less than 2 in bessj");
    if (n < 2) error("Index n less than 2 in bessj");

	ax=rabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > (double) n) {
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+(int) rsqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans : ans;
}
#undef ACC
#undef BIGNO
#undef BIGNI


// BESSK1: modified Bessel function K1(x). 

real bessk1(real x)
{
    real t, tt, ti, u;
    
    if (x < 0.0)
        error("bessk1: negative argument\n");
    t = x / 2.0;
    if (t < 1.0) {
        tt = t * t;
        u = 1.0 +
	    tt * (0.15443144 + tt * (-0.67278579 + tt * (-0.18156897 +
                                                     tt * (-0.01919402 + tt * (-0.00110404 +
                                                                               tt * -0.00004686)))));
        return (u / x + log(t) * bessi1(x));
    } else {
        ti = 1.0 / t;
        u = 1.25331414 + ti * (0.23498619 + ti * (-0.03655620 +
                                                  ti * (0.01504268 + ti * (-0.00780353 + ti * (0.00325614 +
                                                                                               ti * -0.00068245)))));
        return (u * exp(- x) / sqrt(x));
    }
}


// end Bessel functions from Numerical Recipes

// Chebyshev polinomials: T(n,x) = cos(n arccos x)
// This polynomial version is quite slow. It is better to use above formula
// Recurrence relationship: T(n+1,x) = 2 x T(n,x) - T(n-1,x)
double Chebyshev(int n, double x)
{
// n must be an integer >= 0
// x must be a real number in (-1, 1)

    real T0=1.0;
    real T1;
    real Tnm1;
    real Tn;
    real Tnp1;
    int i;

    if (n==0) return T0;
    T1 = x;
    if (n==1) return T1;

    Tn = T1;
    Tnm1 = T0;

    Tnp1 = 2.0*x*Tn - Tnm1;                             // Just to avoid warning

    for (i=2; i<=n; n++) {
        Tnp1 = 2.0*x*Tn - Tnm1;
        Tn = Tnp1;
        Tnm1 = Tn;
    }
    
    return Tnp1;
}

// Use qsort from stdlib.h instead of int_piksrt
/*
 #include <stdio.h>
 #include <stdlib.h>

 // Comparison function
 int compare(const void* a, const void* b) {
    return (*(int*)a - *(int*)b);
 }

 int main() {
    int arr[] = {10, 5, 4, 6, 9};
    int n = sizeof(arr) / sizeof(arr[0]);

    qsort(arr, n, sizeof(int), compare);

    printf("Following is the sorted array: ");
    
    int i;
    for (i = 0; i < n; ++i) {
       printf("%d ", arr[i]);
    }
    printf("\n");
    return 0;
 }
 */
// Ordena un arreglo de enteros de menor a mayor
void int_piksrt(int n, int arr[])
{
    int i, j;
    int a;

    for (j=1; j<n; j++) {
        a = arr[j];
        i=j-1;
        while (i>=0 && arr[i] > a) {
            arr[i+1]=arr[i];
            i--;
        }
        arr[i+1]=a;
    }
}
