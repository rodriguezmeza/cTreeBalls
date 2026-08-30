#ifndef _protodefs_octree_2balls_mpi_h
#define _protodefs_octree_2balls_mpi_h

#define OCTREE2BALLSMPIMETHOD 177

#include "fcfc_octree_2balls_mpi.h"

global int searchcalc_octree_2balls_mpi(struct cmdline_data *,
                                        struct global_data *, bodyptr *,
                                        INTEGER *, INTEGER, INTEGER *,
                                        int, int);

#endif /* !_protodefs_octree_2balls_mpi_h */
