#ifndef _mask_utils_h
#define _mask_utils_h

#include <stdint.h>
#include <string.h>

#define MASK_NODE_EMPTY  (-1)
#define MASK_NODE_MASKED 0
#define MASK_NODE_VALID  1
#define MASK_NODE_MIXED  2

static inline short mask_node_combine(short aggregate, short child)
{
    if (aggregate == MASK_NODE_EMPTY)
        return child;
    if (aggregate == child)
        return aggregate;
    return MASK_NODE_MIXED;
}

static inline int mask_node_can_approximate(short state)
{
    return state == MASK_NODE_VALID;
}

/* Bit tests remain effective when the main build uses -ffast-math. */
static inline int mask_float_is_finite(float value)
{
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

static inline int mask_real_is_finite(double value)
{
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT64_C(0x7ff0000000000000))
            != UINT64_C(0x7ff0000000000000);
}

#endif
