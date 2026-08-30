#include <assert.h>
#include <math.h>
#include <string.h>

#include "globaldefs.h"
#include "tree_contracts.h"

int main(void)
{
    cell parent;
    cell child_cell;
    body child_body;
    compute_vector center_of_mass_sum;
    short mask = MASK_NODE_EMPTY;

    memset(&parent, 0, sizeof(parent));
    memset(&child_cell, 0, sizeof(child_cell));
    memset(&child_body, 0, sizeof(child_body));
    CLRV(center_of_mass_sum);

    Type(&parent) = CELL;
    Type(&child_cell) = CELL;
    Type(&child_body) = BODY;
    Mask(&child_cell) = MASK_NODE_VALID;
    Mask(&child_body) = MASK_NODE_VALID;

    Nb(&child_cell) = 2;
    Weight(&child_cell) = 5.0;
#ifndef NOWKAvg
    Kappa(&child_cell) = 16.0;
    Weight(&child_body) = 1.0;
#else
    Kappa(&child_cell) = 15.0;
    Weight(&child_body) = 1.0;
#endif
    Kappa(&child_body) = 30.0;

    cballs_cell_accumulate_child((nodeptr)&parent, (nodeptr)&child_cell,
                                 FALSE, &mask, center_of_mass_sum);
    cballs_cell_accumulate_child((nodeptr)&parent, (nodeptr)&child_body,
                                 FALSE, &mask, center_of_mass_sum);
    cballs_cell_finalize_averages((nodeptr)&parent);

    assert(Nb(&parent) == 3);
#ifndef NOWKAvg
    assert(fabs(Kappa(&parent) - (110.0/6.0)) < 1.0e-12);
    assert(fabs(Weight(&parent) - 6.0) < 1.0e-12);
#else
    assert(fabs(Kappa(&parent) - 20.0) < 1.0e-12);
    assert(fabs(Weight(&parent) - 2.0) < 1.0e-12);
#endif

    return 0;
}
