#ifndef CTREEBALLS_KDTREE_SHEAR_SPHERE_TREE_H
#define CTREEBALLS_KDTREE_SHEAR_SPHERE_TREE_H

#include "fcfc_balltree.h"

int kdtree_shear_sphere_build(
    struct cmdline_data *, struct global_data *, bodyptr, INTEGER, int,
    bool, fcfc_balltreeptr *);
int kdtree_shear_sphere_frontier(
    struct cmdline_data *, fcfc_balltreeptr, INTEGER, INTEGER **, INTEGER *);
void kdtree_shear_sphere_free(fcfc_balltreeptr);

#endif
