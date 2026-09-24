#ifndef CBALLS_RESOURCE_CONTRACTS_H
#define CBALLS_RESOURCE_CONTRACTS_H
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
/* All shapes are checked before multiplication or signed conversion. */
static inline int cballs_size_add(size_t a, size_t b, size_t *out)
{
    if (a > (size_t)PTRDIFF_MAX || b > (size_t)PTRDIFF_MAX-a) return 0;
    *out = a+b; return 1;
}
static inline int cballs_size_mul(size_t a, size_t b, size_t *out)
{
    if (a && b > (size_t)PTRDIFF_MAX/a) return 0;
    *out = a*b; return 1;
}
static inline int cballs_extent(long lo, long hi, size_t *out)
{
    /* NR shifted pointers support the zero/one/negative origins in this code.
     * Reject origins >1: they would form a pointer before the allocation. */
    if (hi < lo || hi < 0 || lo > 1 || lo == LONG_MIN) return 0;
    uintmax_t distance = (uintmax_t)hi-(uintmax_t)lo;
    if (distance >= (uintmax_t)PTRDIFF_MAX) return 0;
    *out = (size_t)distance+1; return 1;
}
/* Includes NR padding and all pointer tables. rank is 1, 2 or 3. */
static inline int cballs_shape_bytes(int rank, size_t a, size_t b, size_t c,
                                     size_t item, size_t *total)
{
    size_t n, padded, data, pointers=0, planes=0;
    if (rank < 1 || rank > 3 || !a || (rank>1 && !b) || (rank>2 && !c)) return 0;
    n=a;
    if (rank>1 && (!cballs_size_add(a,1,&padded)
                  || !cballs_size_mul(padded,sizeof(void *),&planes)
                  || !cballs_size_mul(a,b,&n))) return 0;
    if (rank>2 && (!cballs_size_add(n,1,&padded)
                  || !cballs_size_mul(padded,sizeof(void *),&pointers)
                  || !cballs_size_mul(n,c,&n))) return 0;
    return cballs_size_add(n,1,&padded) && cballs_size_mul(padded,item,&data)
        && cballs_size_add(data,pointers,total) && cballs_size_add(*total,planes,total);
}
#endif
