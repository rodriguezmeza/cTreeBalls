//=============================================================================
//        1          2          3          4        ^ 5          6          7

#ifndef _stdinc_h
#define _stdinc_h

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define copyright	"Copyright (c) 2000-2026 M.A. Rodriguez-Meza, MEXICO."

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

/*
 * SINGLEP is a mixed-precision search profile.  It reduces the storage used
 * by positions and tree geometry.  Final body acceptance and accumulations
 * stay in the normal `real` (double-precision) type; selected tree broad-phase
 * tests use conservative storage-precision bounds.
 */
#ifdef SINGLEP
typedef float cballs_storage_real;
#define CBALLS_STORAGE_EPSILON FLT_EPSILON
#define CballsStoragePrecision "float-storage/double-compute"
#else
typedef real cballs_storage_real;
#define CBALLS_STORAGE_EPSILON DBL_EPSILON
#define CballsStoragePrecision "double-storage/double-compute"
#endif

typedef real cballs_compute_real;
typedef real cballs_accum_real;

/* Legacy code uses REAL for arithmetic and physical fields. */
#define REAL cballs_compute_real

static inline cballs_storage_real cballs_store_upper_bound(real value)
{
#ifdef SINGLEP
    float stored = (float)value;

    if (isfinite(value) && (real)stored < value)
        stored = nextafterf(stored, INFINITY);
    return stored;
#else
    return value;
#endif
}

static inline cballs_storage_real cballs_store_search_bound(real value)
{
#ifdef SINGLEP
    value += 8.0 * FLT_EPSILON * (1.0 + fabs(value));
    return cballs_store_upper_bound(value);
#else
    return value;
#endif
}

//B Definition of integers:
//
// use INTEGER_FMT as:
//
//
// use INTEGER_SCAN_FMT as:
// if (fscanf(str, "%" INTEGER_SCAN_FMT, iptr) != 1)
//    error("in_int_long: input conversion error\n");
//
#ifdef LONGINT
typedef long integer, *integerptr;
#define INTEGER     integer
#define INTEGER_FMT "ld"
#define INTEGER_SCAN_FMT INTEGER_FMT
#else
typedef int integer, *integerptr;
#define INTEGER     integer
#define INTEGER_FMT "d"
#define INTEGER_SCAN_FMT INTEGER_FMT
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
// should use param="\0" to avoid some calling problems in Cython interface...
#define strnull(x) (strcmp((x), "") == 0)


#ifndef ABS
#define ABS(x)   (((x)<0)?-(x):(x))
#endif

//B
#ifndef MIN
#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) //< the usual "min" function
#endif
#ifndef MAX
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) //< the usual "max" function
#endif
#ifndef SIGN
#define SIGN(a) (((a)>0) ? 1. : -1. )
#endif
#ifndef NRSIGN
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
//E

#define SignR(x,y)  (((y) >= 0) ? (x) : (- (x)))


//B Memory section

typedef int (*cballs_allocation_callback)(void *);

int cballs_allocation_guard(cballs_allocation_callback callback,
                            void *argument,
                            char *errmsg, size_t errmsg_size);
void cballs_allocation_failure(size_t bytes, const char *label);
int cballs_malloc_checked(void **pointer, size_t count, size_t item_size,
                          const char *label,
                          char *errmsg, size_t errmsg_size);
int cballs_calloc_checked(void **pointer, size_t count, size_t item_size,
                          const char *label,
                          char *errmsg, size_t errmsg_size);
void cballs_test_fail_allocation_after(long successful_allocations);
void cballs_test_reset_allocation_failure(void);

void *allocate_array(int);			// Original definition ...
void *allocate(long int);		// Correction to work with more than 16x10^6 particles...
real *AllocVecR(int);			// Define an array of reals like in fortran
int *AllocVecI(int);			// Define an array of integers like in fortran
int *AllocVecINormal(int);		// Define an array of integers like in C (starting in zero)
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
void verb_print_zero(int verbose, string fmt, ...);
void verb_print_warning(int verbose, string fmt, ...);
void verb_print_debug(int verbose, string fmt, ...);
void verb_log_print(int verbose, stream sout,  string fmt, ...);
void verb_print_no_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void verb_print_min_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void verb_print_normal_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void verb_print_debug_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...);
void   endrun(int);
												
void eprintf(string, ...);						

bool scanopt(string, string);					

//B from cBalls
int stropen_checked(string name, string mode, stream *out,
                    char *errmsg, size_t errmsg_size);
//E

stream stropen(string, string);

//B From G
double dmax(double,double);
double dmin(double,double);
int    imax(int,int);
int    imin(int,int);
//E

//B
// to stop using this system()...
// sprintf(outputDir,cmd->rootDir);
// sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi", outputDir, outputDir);
// system(buf);
int mkdir_p(const char *path, mode_t mode);

// to solve unsafe strcpy / sprintf usage
int copy_checked(char *dst, size_t dst_size,
                const char *src, const char *label);

int parse_int_checked(const char *text, int *value,
                      char *errmsg, size_t errmsg_size,
                      const char *label);
int parse_long_checked(const char *text, long *value,
                       char *errmsg, size_t errmsg_size,
                       const char *label);
int parse_double_checked(const char *text, double *value,
                         char *errmsg, size_t errmsg_size,
                         const char *label);
int parse_bool_checked(const char *text, bool *value,
                       char *errmsg, size_t errmsg_size,
                       const char *label);

int read_parameter_line_checked(FILE *input, char *line, size_t line_size,
                                int *has_line, unsigned long *line_number,
                                const char *filename,
                                char *errmsg, size_t errmsg_size);
int parse_parameter_line_checked(const char *line,
                                 char *name, size_t name_size,
                                 char *value, size_t value_size,
                                 int *is_data,
                                 const char *filename,
                                 unsigned long line_number,
                                 char *errmsg, size_t errmsg_size);
//E

//B CLASS type definitions
void class_protect_sprintf(char *dest, size_t dest_size, const char *tpl, ...);
//E

//B added by cBalls
#if defined(__GNUC__) || defined(__clang__)
#define CBALLS_PRINTF_FORMAT(a,b) __attribute__((format(printf,a,b)))
#else
#define CBALLS_PRINTF_FORMAT(a,b)
#endif

int format_checked(char *dst, size_t dst_size,
                   const char *label, const char *fmt, ...)
                   CBALLS_PRINTF_FORMAT(4, 5);

//E

#endif  /* ! _stdinc_h	*/

