//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _stdinc_h
#define _stdinc_h

#include <stdio.h>
#include <stdlib.h>

#define copyright	"Copyright (c) 2000-2023 M.A. Rodriguez-Meza, MEXICO."

#if !defined(NULL)
#define NULL 0L
#endif

#define local     static

//B Here subtitution of bool definition by the ISO C standard
// ORIGINAL is defined in Makefile_fcfc_kdballtree in addons folder
#ifndef ORIGINAL
#include <stdbool.h>
#else
typedef short int bool;
#endif
//E

#if !defined(TRUE)
#define TRUE  ((bool) 1)
#define FALSE ((bool) 0)
#endif

typedef unsigned char byte;

typedef char *string;

typedef FILE *stream;

//B PRECISION SECTION

#define DOUBLEPREC
// #undef DOUBLEPREC

#define SINGLEPREC
#undef SINGLEPREC

#define MIXEDPREC
#undef MIXEDPREC

#if defined(MIXEDPREC)
#undef MIXEDPREC
#undef SINGLEPREC
#endif

#if defined(SINGLEPREC)
#undef SINGLEPREC
#undef MIXEDPREC
#endif

//B Machine parameters
#define EPSILONMACHSINGLE = 1.0E-7
#define EPSILONMACHDOUBLE = 1.0E-16
//E

//E PRECISION SECTION

#if !defined(MIXEDPREC) && !defined(SINGLEPREC) && !defined(DOUBLEPREC)
#define SINGLEPREC
#endif

#if defined(DOUBLEPREC)
#undef SINGLEPREC
#undef MIXEDPREC
typedef double real, *realptr;
#define Precision "DOUBLEPREC"
#endif

#if defined(MIXEDPREC)
#undef DOUBLEPREC
#undef SINGLEPREC
typedef float *realptr, real;
#define Precision "MIXEDPREC"
#endif

#if defined(SINGLEPREC)
#undef DOUBLEPREC
#undef MIXEDPREC
typedef float real, *realptr;
#define Precision "SINGLEPREC"
#endif

#ifdef SINGLEP
#define REAL     float
#else
#define REAL     real
#endif

//B Definition of integers:
#ifdef LONGINT
typedef long integer, *integerptr;
#define INTEGER     integer
#else
typedef int integer, *integerptr;
#define INTEGER     integer
#endif
//E

#if defined(SPEEDUPANDMEM)
typedef float FLOAT, *FLOATPTR;
#else
typedef real FLOAT, *FLOATPTR;
#endif

// Complex definitions


//B MAR: conflicts with ctelib: src/pm.c, it includes <complex.h>
typedef struct {
  real R, II;
} Cmplx;


#define CSet(a, x, y)                                       \
   a.R = x,                                                 \
   a.II = y
#define CAdd(a, b, c)                                       \
   a.R = b.R + c.R,                                         \
   a.II = b.II + c.II
#define CSub(a, b, c)                                       \
   a.R = b.R - c.R,                                         \
   a.II = b.II - c.II
#define CMul(a, b, c)                                       \
  a.R = b.R * c.R - b.II * c.II,                              \
  a.II = b.R * c.II + b.II * c.R
//E


#ifndef PI
#define PI		   3.141592653589793238462643383279502884197
#endif

#define TWO_PI     6.28318530717958647693
#define FOUR_PI   12.56637061435917295385
#define HALF_PI    1.57079632679489661923
#define FRTHRD_PI  4.18879020478639098462

#if !defined(M_LN2)
#define M_LN2		0.69314718055994530942
#endif
#if !defined(M_LN10)
#define M_LN10		2.30258509299404568402
#endif


#if !defined (HZ)
#define HZ        100
#endif


#define streq(x,y) (strcmp((x), (y)) == 0)
#define strnull(x) (strcmp((x), "") == 0)


#ifndef ABS
#define ABS(x)   (((x)<0)?-(x):(x))
#endif

/*
// Defined in class_lib/common.h :: activate if class_lib not used
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
*/

#define SignR(x,y)  (((y) >= 0) ? (x) : (- (x)))


//B Memory section

void *allocate_array(int);			// Original definition ...
void *allocate(long int);		// Correction to work with more than 16x10^6 particles...
real *AllocVecR(int);			// Define an array of reals like in fortran
int *AllocVecI(int);			// Define an array of integers like in fortran
int *AllocVecINormal(int);		// Define an array of integers like in C (starting in zero)
//void FreeVecR(real *);
//void FreeVecI(int *);
void FreeVecINormal(int *);

int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double ***dmatrix3D(long nrl, long nrh, long ncl, long nch, long nc2l, long nc2h);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix3D(double ***m, long nrl, long nrh, long ncl, long nch, long nc2l, long nc2h);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);


#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))
#define AllocMem2(a, n1, n2, t)                             \
   AllocMem (a, n1, t *);                                   \
   AllocMem (a[0], (n1) * (n2), t);                         \
   for (k = 1; k < n1; k ++) a[k] = a[k - 1] + n2;

#define Cube(x)    ((x) * (x) * (x))

#define Nint(x)                                             \
   (((x) < 0.) ? (- (int) (0.5 - (x))): ((int) (0.5 + (x))))

void error_mem_out_kd(void);

//E Memory section


//B Time section
double cputime(void);
double cpu_time(void);  // same as above but a little bit different
long rcpu_time(void);
double second(void);
double timediff(double t0,double t1);
void timer_kd(int i);
//E Time section

void error(string, ...);
void verb_print_q(int iq, int verbose, string fmt, ...);
void verb_print(int verbose, string fmt, ...);
void verb_print_warning(int verbose, string fmt, ...);
void verb_print_debug(int verbose, string fmt, ...);
void verb_log_print(int verbose, stream sout,  string fmt, ...);
void verb_print_min_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void verb_print_normal_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void verb_print_debug_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void   endrun(int);
												
void eprintf(string, ...);						

bool scanopt(string, string);					

stream stropen(string, string);					

//double second(void);
//double timediff(double t0,double t1);

//B From G
double dmax(double,double);
double dmin(double,double);
int    imax(int,int);
int    imin(int,int);
//E

#endif  /* ! _stdinc_h	*/

