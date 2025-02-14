// Use:
//#include "protodefs_kdtree_omp.h"

// NMultipoles has been switched off for kdtree_omp
//  NMultipoles -> NMultipoles_kdtree

#ifndef _protodefs_kdtree_omp_h
#define _protodefs_kdtree_omp_h

global int searchcalc_kdtree_omp(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 bodyptr *btable,
                                 INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2);

#endif	// ! _protodefs_kdtree_omp_h
