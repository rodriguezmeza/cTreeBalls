#include <assert.h>
#include <math.h>
#include <string.h>

#include "globaldefs.h"
#include "tree_contracts.h"

static void set_position(nodeptr p, real x, real y, real z)
{
    Pos(p)[0] = x;
    Pos(p)[1] = y;
#if NDIM == 3
    Pos(p)[2] = z;
#else
    (void)z;
#endif
}

int main(void)
{
    struct cmdline_data cmd;
    struct global_data gd;
    cell parent;
    body first;
    body second;
    body pivot;
    compute_vector center_of_mass_sum;
    compute_vector dr;
    short mask = MASK_NODE_EMPTY;
    real distance;

    memset(&cmd, 0, sizeof(cmd));
    memset(&gd, 0, sizeof(gd));
    memset(&parent, 0, sizeof(parent));
    memset(&first, 0, sizeof(first));
    memset(&second, 0, sizeof(second));
    memset(&pivot, 0, sizeof(pivot));
    CLRV(center_of_mass_sum);

    Type(&parent) = CELL;
    Type(&first) = BODY;
    Type(&second) = BODY;
    Type(&pivot) = BODY;
    Mask(&first) = MASK_NODE_VALID;
    Mask(&second) = MASK_NODE_VALID;
    Mass(&first) = 1.0;
    Mass(&second) = 3.0;
    Weight(&first) = 2.0;
    Weight(&second) = 3.0;
    Kappa(&first) = 10.0;
    Kappa(&second) = 20.0;
    Selected(&second) = TRUE;
    Update(&first) = TRUE;
    set_position((nodeptr)&first, 1.0, 0.0, 0.0);
    set_position((nodeptr)&second, 3.0, 0.0, 0.0);

    assert(!cballs_cell_accumulate_child((nodeptr)&parent,
                                         (nodeptr)&first, FALSE,
                                         &mask, center_of_mass_sum));
    assert(!cballs_cell_accumulate_child((nodeptr)&parent,
                                         (nodeptr)&second, FALSE,
                                         &mask, center_of_mass_sum));
    cballs_cell_finalize_averages((nodeptr)&parent);

    assert(Nb(&parent) == 2);
    assert(fabs(Weight(&parent) - 5.0) < 1.0e-12);
    assert(fabs(Kappa(&parent) - 16.0) < 1.0e-12);
    assert(fabs(Mass(&parent) - 4.0) < 1.0e-12);
    assert(fabs(center_of_mass_sum[0] - 10.0) < 1.0e-12);
    assert(Selected(&parent));
    assert(Update(&parent));

    gd.Rcut = 1.0;
    set_position((nodeptr)&pivot, 0.0, 0.0, 0.0);
    set_position((nodeptr)&first, 0.5, 0.0, 0.0);
    assert(cballs_accept_body_contract(&cmd, &gd, &pivot,
                                       (nodeptr)&first, &distance, dr));
    assert(fabs(distance - 0.5) < 1.0e-12);

    set_position((nodeptr)&first, 1.0, 0.0, 0.0);
    assert(!cballs_accept_body_contract(&cmd, &gd, &pivot,
                                        (nodeptr)&first, &distance, dr));

    cmd.usePeriodic = TRUE;
    gd.Rcut = 0.5;
    gd.Box[0] = 10.0;
    gd.Box[1] = 10.0;
#if NDIM == 3
    gd.Box[2] = 10.0;
#endif
    set_position((nodeptr)&pivot, 4.8, 0.0, 0.0);
    set_position((nodeptr)&first, -4.8, 0.0, 0.0);
    assert(cballs_accept_body_contract(&cmd, &gd, &pivot,
                                       (nodeptr)&first, &distance, dr));
#ifdef SINGLEP
    assert(fabs(distance - 0.4) < 8.0 * FLT_EPSILON);
#else
    assert(fabs(distance - 0.4) < 1.0e-12);
#endif

    return 0;
}
