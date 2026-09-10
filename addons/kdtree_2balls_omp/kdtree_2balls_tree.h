#ifndef _kdtree_2balls_tree_h
#define _kdtree_2balls_tree_h

#include "fcfc_balltree.h"

int kdtree_2balls_tree_build_pivot(
    struct cmdline_data *, struct global_data *, bodyptr, INTEGER, int,
    fcfc_balltreeptr *);
int kdtree_2balls_tree_build_neighbor(
    struct cmdline_data *, struct global_data *, bodyptr, INTEGER, int,
    fcfc_balltreeptr *);
int kdtree_2balls_tree_frontier(
    struct cmdline_data *, fcfc_balltreeptr, INTEGER, INTEGER **, INTEGER *);
void kdtree_2balls_tree_free(fcfc_balltreeptr);

#endif /* !_kdtree_2balls_tree_h */
