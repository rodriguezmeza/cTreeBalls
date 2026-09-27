#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct {
    int64_t forest_id;
    size_t leg_bin,ordinal;
    double radius,payload[8];
} lya_neighbor;
#include "../addons/lya_forest_omp/lya_neighbor_sort.h"

static int oracle(const void *aa,const void *bb)
{
    const lya_neighbor *a=aa,*b=bb;
    if(a->forest_id!=b->forest_id)return a->forest_id<b->forest_id?-1:1;
    if(a->leg_bin!=b->leg_bin)return a->leg_bin<b->leg_bin?-1:1;
    if(a->radius!=b->radius)return a->radius<b->radius?-1:1;
    return (a->ordinal>b->ordinal)-(a->ordinal<b->ordinal);
}

int main(void)
{
    uint64_t state=190926;
    for(size_t count=0;count<=8192;count=count<32?count+1:count*2) {
        lya_neighbor *input=calloc(count+1,sizeof(*input));
        lya_neighbor *want=malloc((count+1)*sizeof(*want));
        lya_neighbor *actual=malloc((count+1)*sizeof(*actual));
        assert(input && want && actual);
        for(int pattern=0;pattern<6;pattern++) {
            for(size_t i=0;i<count;i++) {
                state=state*UINT64_C(6364136223846793005)+1;
                input[i].forest_id=pattern<2?(int64_t)(state%47)-23:INT64_MIN;
                input[i].leg_bin=pattern<3?state%29:0;
                input[i].radius=pattern<4?(double)(state%71):1.;
                input[i].ordinal=pattern%2?count-i:i;
                for(int k=0;k<8;k++)input[i].payload[k]=(double)i+k*.25;
            }
            memcpy(want,input,count*sizeof(*input));qsort(want,count,sizeof(*want),oracle);
            for(int force_heap=0;force_heap<2;force_heap++) {
                memcpy(actual,input,count*sizeof(*input));
                if(force_heap) lya_neighbor_intro(actual,count,0);else lya_sort_neighbors(actual,count);
                assert(memcmp(actual,want,count*sizeof(*actual))==0);
            }
        }
        free(input);free(want);free(actual);
    }
    return 0;
}
