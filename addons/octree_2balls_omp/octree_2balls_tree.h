#ifndef _octree_2balls_tree_h
#define _octree_2balls_tree_h

#include "fcfc_balltree.h"

int octree_2balls_tree_build(struct cmdline_data *, struct global_data *,
                             bodyptr, INTEGER, int, fcfc_balltreeptr *);
int octree_2balls_tree_frontier(struct cmdline_data *, fcfc_balltreeptr,
                                INTEGER, INTEGER **, INTEGER *);
void octree_2balls_tree_free(fcfc_balltreeptr);

#endif /* !_octree_2balls_tree_h */
