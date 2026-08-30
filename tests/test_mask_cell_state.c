#include <assert.h>
#include <math.h>
#include "mask_utils.h"

int main(void)
{
    short state = MASK_NODE_EMPTY;

    state = mask_node_combine(state, MASK_NODE_VALID);
    assert(state == MASK_NODE_VALID);
    assert(mask_node_can_approximate(state));

    state = mask_node_combine(state, MASK_NODE_MASKED);
    assert(state == MASK_NODE_MIXED);
    assert(!mask_node_can_approximate(state));

    state = mask_node_combine(MASK_NODE_MASKED, MASK_NODE_MASKED);
    assert(state == MASK_NODE_MASKED);
    assert(!mask_node_can_approximate(state));

    assert(mask_float_is_finite(1.0f));
    assert(!mask_float_is_finite(INFINITY));
    assert(mask_real_is_finite(1.0));
    assert(!mask_real_is_finite(NAN));

    return 0;
}
