/* ==============================================================================
 MODULE: kdtree.h        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Based on: zeno lib
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: include "kdtree.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _kdtree_h
#define _kdtree_h

// Structure representing the bounds of a kd node
typedef struct {
    vector minb;                                    // min pos of bounding box
    vector maxb;                                    // max pos of bounding box
    vector width;                                   // width of bounding box
    vector center;                                  // center of bounding box
    real radius;                                    // radius of bounding box
} bound;

//  Structure representing a node of a kd-tree
typedef struct {
    bound bnd;                                      // min, max bnd of ball-cell
    int dim;                                        // splitting dimension
    real split;                                     // median coordinate value
    int first;                                      // index of frst body in node
    int last;                                       // index of last body in node
    real kappa;                                     // scalar field at cmpos
    real weight;                                    // point weight
    vector cmpos;                                   // center of mass pos
    matrix Ixy;                                     // inertia tensor
    real etaxy;                                     // deformation factor xy
#if NDIM == 3
    real etaxz;                                     // deformation factor xz
    real etayz;                                     // deformation factor yz
#endif
} kdnode;

// Macros for indicies into kd tree
#define KDROOT		1
#define Lower(i)	((i)<<1)
#define Upper(i)	(((i)<<1)+1)
#define Parent(i)	((i)>>1)
#define Sibling(i) 	(((i)&1)?(i)-1:(i)+1)

// Q: i is an odd number? If no divide i by 2 and go back...
#define SetNext(i) {				\
  while (i&1)					    \
    i=i>>1;					        \
  ++i;						        \
}

// Structure to represent context of KD tree
typedef struct {
    INTEGER npoint;
    bodyptr *bptr;
    bound bnd;
    int nnode;
    int nsplit;
    kdnode *ntab;
} kdcontext, *kdxptr;

kdxptr init_kdtree(struct cmdline_data* cmd,
                   struct  global_data* gd,
                   bodyptr, INTEGER);
void build_kdtree(struct cmdline_data* cmd,
                  struct  global_data* gd,
                  kdxptr, int);
void finish_kdtree(kdxptr);

#endif  /* ! _kdtree_h */
