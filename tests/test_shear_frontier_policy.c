/* Standalone policy tests: compile once normally and once with
 * -DSHEAR_MPI_ENABLED. No MPI launcher or large histogram allocation needed. */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define INTEGER long
#define real double
#define SUCCESS 0
#define FAILURE 1
#define _ERRORMSGSIZE_ 1024
#define SHEAR_ENGINE_NAME "policy-test"
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
struct cmdline_data {char error_message[_ERRORMSGSIZE_];};
#include "../support/shear_kernel/shear_frontier_policy.h"

int main(void) {
    struct cmdline_data cmd = {0};
    INTEGER target;
    size_t limit;
    unsetenv("CBALLS_SHEAR_TASKS");
    unsetenv("CBALLS_SHEAR_TASK_MEMORY_MB");
    assert(shear_frontier_policy(&cmd,64,2,4096,256,&target,&limit)==SUCCESS);
#ifdef SHEAR_MPI_ENABLED
    assert(target==1024); /* Eight tasks per worker, beyond the old cap. */
#else
    assert(target==256);
#endif
    setenv("CBALLS_SHEAR_TASK_MEMORY_MB","1",1);
    setenv("CBALLS_SHEAR_TASKS","4096",1);
    assert(shear_frontier_policy(&cmd,64,4,4096,256,&target,&limit)==SUCCESS);
    assert(target==128 && limit==32);
    shear_frontier_priority order[128];
    for(int i=0;i<128;i++) { order[i].task=i; order[i].cost=i==0?1.e12L:0; }
    assert(shear_frontier_assign(&cmd,order,128,4,limit)==SUCCESS);
    int counts[4]={0},slots[4][32]={{0}};
    for(int i=0;i<128;i++) {
        assert(order[i].owner>=0 && order[i].owner<4);
        assert(order[i].local_index>=0 && order[i].local_index<32);
        assert(!slots[order[i].owner][order[i].local_index]++);
        counts[order[i].owner]++;
    }
    for(int i=0;i<4;i++)assert(counts[i]==32);
    assert(shear_frontier_assign(&cmd,order,128,4,31)==FAILURE);
    /* A single large histogram is an explicit unavoidable minimum. */
    assert(shear_frontier_policy(&cmd,64,4,1048576,256,&target,&limit)==SUCCESS);
    assert(limit==1 && target==4);
    assert(shear_frontier_policy(&cmd,64,4,SIZE_MAX,256,&target,&limit)==FAILURE);
    for(int i=0;i<6;i++) {
        const char *bad[]={"-1","nan","3.4","","99999999999999999999999999","8x"};
        setenv("CBALLS_SHEAR_TASKS",bad[i],1);
        assert(shear_frontier_policy(&cmd,2,2,4096,64,&target,&limit)==FAILURE);
    }
    puts("SHEAR_FRONTIER_POLICY_OK");
    return 0;
}
