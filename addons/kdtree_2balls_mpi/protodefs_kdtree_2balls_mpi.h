#ifndef _protodefs_kdtree_2balls_mpi_h
#define _protodefs_kdtree_2balls_mpi_h

#include "protodefs_kdtree_2balls_omp.h"
#include "fcfc_kdtree_2balls_mpi.h"

#define KDTREE2BALLSMPIMETHOD 198

global int searchcalc_kdtree_2balls_mpi(
    struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *,
    INTEGER, INTEGER *, int, int);

#endif /* !_protodefs_kdtree_2balls_mpi_h */
