//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include <math.h>
#include <time.h>
#include <sys/time.h>                               // To get time of the day

#include "stdinc.h"
#include "getparam.h"

                                                    // To remove:
#include <unistd.h>                                 // warning: implicit
                                                    //  declaration of function
                                                    //  ‘dup’
                                                    //  [-Wimplicit-function
                                                    //  -declaration]
                                                    //  fds = dup(fileno(inflag ?
                                                    //  stdin : stdout));
#include <sys/timeb.h>
#include <string.h>	

#include <sys/stat.h>
#include <stdarg.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>


//B Memory section

void *allocate_array(int nb)
{
    void *mem;

    mem = calloc(nb, 1);
    if (mem == NULL)
        error("allocate_array in %s: not enuf memory (%d bytes)\n",
              getargv0(), nb);
    return (mem);
}

void *allocate(long int nb)
{
    void *mem;

    mem = calloc(nb, 1);
    if (mem == NULL)
        error("allocate in %s: not enuf memory (%d bytes)\n",
              getargv0(), nb);
    return (mem);
}


void FreeVecINormal(int *v)
{
	free(v);
}

//B Following NR
#define NOFFSET_END 1
#define FREE_ARG char*

int *ivector(long nl, long nh)
{
    int *v;

    v=(int *)malloc((size_t) ((nh-nl+1+NOFFSET_END)*sizeof(int)));
    if (!v) error("allocation failure in ivector()");
    return v-nl+NOFFSET_END;
}

void free_ivector(int *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-NOFFSET_END));
}

double *dvector(long nl, long nh)
{
    double *v;

    v=(double *)malloc((size_t) ((nh-nl+1+NOFFSET_END)*sizeof(double)));
    if (!v) error("allocation failure in dvector()");

    return v-nl+NOFFSET_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,
            ncol=nch-ncl+1;
    double **m;

    m=(double **) malloc((size_t)((nrow+NOFFSET_END)*sizeof(double*)));
    if (!m) error("allocation failure 1 in matrix()");
    m += NOFFSET_END;
    m -= nrl;

    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NOFFSET_END)*sizeof(double)));
    if (!m[nrl]) error("allocation failure 2 in matrix()");
    m[nrl] += NOFFSET_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}


double ***dmatrix3D(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    double ***t;

    /* allocate pointers to pointers to rows */
    t=(double ***) malloc((size_t)((nrow+NOFFSET_END)*sizeof(double**)));
    if (!t) error("allocation failure 1 in f3tensor()");
    t += NOFFSET_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(double **) malloc((size_t)((nrow*ncol+NOFFSET_END)*sizeof(double*)));
    if (!t[nrl]) error("allocation failure 2 in f3tensor()");
    t[nrl] += NOFFSET_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NOFFSET_END)*sizeof(double)));
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

void free_dmatrix3D(double ***t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
    free((FREE_ARG) (t[nrl][ncl]+ndl-NOFFSET_END));
    free((FREE_ARG) (t[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (t+nrl-NOFFSET_END));
}


void free_dvector(double *v, long nl, long nh)
{
    if (nl>0)
        free((FREE_ARG) (v+nl-NOFFSET_END));
    else
        free((FREE_ARG) (v));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (m+nrl-NOFFSET_END));
}


#undef NOFFSET_END
#undef FREE_ARG
//E

void error_mem_out_kd(void)
{
// Memory shortage handler
  fprintf(stderr,"error_mem_out_kd: Out of memory!!\n");
  exit(1);
}


//E Memory section


//B Time section

double cputime(void)
{
    time_t ltime;
    time(&ltime);
    return(((double)ltime)/((double) 60.0));
}

double cpu_time(void)
{
    double value;
      value = (double) clock () / (double) CLOCKS_PER_SEC;
      return value;
}

// Gives time of the day
#define REALTIMECLIB   (gettimeofday(&current_time, NULL))
#define GETREALTIME(time) \
    {REALTIMECLIB; \
    (time) = current_time.tv_sec;}

long rcpu_time(void)
{
    struct timeval current_time;
    long cpurealfinal;
    GETREALTIME(cpurealfinal);
    return cpurealfinal;
}

#undef REALTIMECLIB
#undef GETREALTIME

double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

}

double timediff(double t0, double t1)
{
  double dt;
  
  dt=t1-t0;

    if(dt<0) {                                      // overflow has occured (for
                                                    //  systems with 32bit tick
                                                    //  counter)
#ifdef WALLCLOCK
        dt = 0;
#else
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}

//  Timing variables
#ifdef OPENMPCODE
#include <omp.h>
static double relbeg_kd,relend_kd,absbeg_kd,absend_kd;
#else //OPENMPCODE
#include <time.h>
static time_t relbeg_kd,relend_kd,absbeg_kd,absend_kd;
#endif //OPENMPCODE

void timer_kd(int i)
{
// Timing routine
// timer(0) -> initialize relative clock
// timer(1) -> read relative clock
// timer(2) -> read relative clock and initialize it afterwards
// timer(4) -> initialize absolute clock
// timer(5) -> read absolute clock
#ifdef OPENMPCODE
  if(i==0)
    relbeg_kd=omp_get_wtime();
  else if(i==1) {
    relend_kd=omp_get_wtime();
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend_kd-relbeg_kd));
  }
  else if(i==2) {
    relend_kd=omp_get_wtime();
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend_kd-relbeg_kd));
    relbeg_kd=omp_get_wtime();
  }
  else if(i==4)
    absbeg_kd=omp_get_wtime();
  else if(i==5) {
    absend_kd=omp_get_wtime();
    printf("    Total time ellapsed %.1lf ms\n",1000*(absend_kd-absbeg_kd));
  }
#else //OPENMPCODE
  int diff;
  
  if(i==0)
    relbeg_kd=time(NULL);
  else if(i==1) {
    relend_kd=time(NULL);
    diff=(int)(difftime(relend_kd,relbeg_kd));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
  }
  else if(i==2) {
    relend_kd=time(NULL);
    diff=(int)(difftime(relend_kd,relbeg_kd));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
    relbeg_kd=time(NULL);
  }
  else if(i==4)
    absbeg_kd=time(NULL);
  else if(i==5) {
    absend_kd=time(NULL);
    diff=(int)(difftime(absend_kd,absbeg_kd));
    printf("    Total time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
  }
#endif //OPENMPCODE
}

//E Time section


void error(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);
    exit(1);
}

void verb_print_q(int iq, int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose == iq) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_print(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose > 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_print_debug(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose > 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_log_print(int verbose, stream sout,  string fmt, ...)
{
    va_list ap;

    if (verbose > 0) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }
}

void endrun(int ierr)
{
  if(ierr)
    {
      fprintf(stdout,
		  "endrun called with an error level of %d\n\n\n", ierr);
      exit(1);
    }
  exit(0);
}


void eprintf(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);
}

bool scanopt(string opt, string key)
{
    char *op, *kp;
    
    op = (char *) opt;
    while (*op != '\0') {
        kp = key;
        while ((*op != ',' ? *op : '\0') == *kp) {
            if (*kp++ == '\0')
                return (TRUE);
            op++;
        }
        while (*op != '\0' && *op++ != ',')
            
            continue;
    }
    return (FALSE);
}


stream stropen(string name, string mode)
{
    bool inflag;
    int fds;
    stream res;
    struct stat buf;

    inflag = streq(mode, "r");
    if (name[0] == '-') {
        if (streq(name, "-")) {
            fds = dup(fileno(inflag ? stdin : stdout));
            if (fds == -1)
                error("stropen: cannot dup %s\n", inflag ? "stdin" : "stdout");
        } else
            fds = atoi(&name[1]);
        res = fdopen(fds, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen: cannot open f.d. %d for %s\n", fds, inflag ? "input" : "output");
    } else {
        if (streq(mode, "w") && stat(name, &buf) == 0)
            error("stropen: file \"%s\" already exists\n", name);
        res = fopen(name, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen: cannot open file \"%s\" for %s\n", name, inflag ? "input" : "output");
    }
    return (res);
}


//B From G

double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}

double dmin(double x,double y)
{
  if(x<y)
    return x;
  else
    return y;
}

int imax(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}

int imin(int x,int y)
{
  if(x<y)
    return x;
  else
    return y;
}

//E
