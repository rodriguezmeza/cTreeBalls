#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stdinc.h"
//#include "numrec.h"
//#include "mathfns.h"

#ifdef USEGSL
#error `USEGSL` and `COMPLEXLIB` can not used at the same time. Switched off one of them
#endif

typedef struct FCOMPLEX {double r,i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


fcomplex Cmul(fcomplex a, fcomplex b)
{
	fcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

fcomplex Complex(double re, double im)
{
	fcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

fcomplex Conjg(fcomplex z)
{
	fcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b)
{
	fcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabs(fcomplex z)
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

fcomplex Csqrt(fcomplex z)
{
	fcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

fcomplex RCmul(double x, fcomplex a)
{
	fcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}


#define NOFFSET_END 1
#define FREE_ARG char*

fcomplex *cmplx_vector(long nl, long nh)
{
    fcomplex *v;

    v=(fcomplex *)malloc((size_t) ((nh-nl+1+NOFFSET_END)*sizeof(fcomplex)));
    if (!v) error("allocation failure in cmplx_vector()");

    return v-nl+NOFFSET_END;
}

fcomplex **cmplx_matrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,
            ncol=nch-ncl+1;
    fcomplex **m;

    m=(fcomplex **) malloc((size_t)((nrow+NOFFSET_END)*sizeof(fcomplex*)));
    if (!m) error("allocation failure 1 in matrix()");
    m += NOFFSET_END;
    m -= nrl;

    m[nrl]=(fcomplex *) malloc((size_t)((nrow*ncol+NOFFSET_END)*sizeof(fcomplex)));
    if (!m[nrl]) error("allocation failure 2 in matrix()");
    m[nrl] += NOFFSET_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}


fcomplex ***cmplx_matrix3D(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    fcomplex ***t;

    /* allocate pointers to pointers to rows */
    t=(fcomplex ***) malloc((size_t)((nrow+NOFFSET_END)*sizeof(fcomplex**)));
    if (!t) error("allocation failure 1 in f3tensor()");
    t += NOFFSET_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(fcomplex **) malloc((size_t)((nrow*ncol+NOFFSET_END)*sizeof(fcomplex*)));
    if (!t[nrl]) error("allocation failure 2 in f3tensor()");
    t[nrl] += NOFFSET_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(fcomplex *) malloc((size_t)((nrow*ncol*ndep+NOFFSET_END)*sizeof(fcomplex)));
    if (!t[nrl][ncl]) error("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NOFFSET_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

void free_cmplx_matrix3D(fcomplex ***t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NOFFSET_END));
    free((FREE_ARG) (t[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (t+nrl-NOFFSET_END));
}


void free_cmplx_vector(fcomplex *v, long nl, long nh)
{
    if (nl>0)
        free((FREE_ARG) (v+nl-NOFFSET_END));
    else
        free((FREE_ARG) (v));
}

void free_cmplx_matrix(fcomplex **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (m+nrl-NOFFSET_END));
}


#undef NOFFSET_END
#undef FREE_ARG


