/* Own runtime globals so the archive does not pull in main.o. */
#define global
#include "globaldefs.h"
#include "stdinc.h"
#include "common_defs.h"
#include "resource_contracts.h"
#include "numrec.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
static int bad_matrix(void *unused) { (void)unused; dmatrix(1,LONG_MAX,1,LONG_MAX); return 0; }
static int bad_tensor(void *unused) { (void)unused; dmatrix3D(1,INT_MAX,1,INT_MAX,1,INT_MAX); return 0; }
static int bad_vector(void *unused) { (void)unused; dvector(LONG_MIN,LONG_MAX); return 0; }
static int over_budget(void *unused) { (void)unused; dmatrix(1,512,1,512); return 0; }
int main(void)
{
    size_t n; char error[512]; void *p=(void *)1;
    assert(!cballs_size_add(PTRDIFF_MAX,1,&n));
    assert(!cballs_size_mul(SIZE_MAX,2,&n));
    assert(!cballs_extent(LONG_MIN,LONG_MAX,&n));
    assert(!cballs_extent(2,4,&n));
    assert(cballs_shape_bytes(3,2,3,4,sizeof(double),&n));
    assert(n==25*sizeof(double)+10*sizeof(void *));
    assert(cballs_allocation_guard(bad_matrix,NULL,error,sizeof(error))==FAILURE);
    assert(strstr(error,"overflow"));
    assert(cballs_allocation_guard(bad_tensor,NULL,error,sizeof(error))==FAILURE);
    assert(cballs_allocation_guard(bad_vector,NULL,error,sizeof(error))==FAILURE);
    setenv("CBALLS_MEMORY_BUDGET_MB","1",1);
    struct cmdline_data cmd={0};struct global_data gd={0};gdhist_sincos_omp hist;
    cmd.options="only-2pcf";cmd.sizeHistN=512;cmd.mChebyshev=31;
    assert(search_init_sincos_omp(&cmd,&gd,&hist)==SUCCESS);
#ifdef TPCF
    assert(hist.histZetaMthreadcos==NULL && hist.histXithreadcos==NULL);
#endif
    assert(search_free_sincos_omp(&cmd,&gd,&hist)==SUCCESS);
    assert(cballs_allocation_guard(over_budget,NULL,error,sizeof(error))==FAILURE);
    assert(strstr(error,"budget"));
    assert(cballs_malloc_checked(&p,SIZE_MAX,8,"test",error,sizeof(error))==FAILURE && p==NULL);
    assert(cballs_calloc_checked(&p,1024*1024,8,"test",error,sizeof(error))==FAILURE && p==NULL);
    double **m=dmatrix(0,2,0,3); m[2][3]=7; assert(m[2][3]==7); free_dmatrix(m,0,2,0,3);
    double ***t=dmatrix3D(1,2,1,3,1,4);t[2][3][4]=9; assert(t[2][3][4]==9);free_dmatrix3D(t,1,2,1,3,1,4);
    double *v=dvector(0,3); v[0]=1;v[3]=2;free_dvector(v,0,3);
    v=dvector(-2,3);v[-2]=1;v[3]=2;free_dvector(v,-2,3);
    float **f=nr_matrix(0,2,0,3);f[2][3]=3;free_matrix(f,0,2,0,3);
    setenv("CBALLS_MEMORY_BUDGET_MB","0",1);
    assert(cballs_memory_preflight(1,"invalid",error,sizeof(error))==FAILURE);
    setenv("CBALLS_MEMORY_BUDGET_MB","18446744073709551615",1);
    assert(cballs_memory_preflight(1,"invalid",error,sizeof(error))==FAILURE);
    unsetenv("CBALLS_MEMORY_BUDGET_MB");
    puts("PASS: checked dimensions, compound preflight, budget errors, zero/one-based storage");
    return 0;
}
