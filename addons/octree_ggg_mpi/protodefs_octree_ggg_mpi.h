#ifndef _protodefs_octree_ggg_mpi_h
#define _protodefs_octree_ggg_mpi_h

#include "fcfc_octree_ggg_mpi.h"

global int searchcalc_octree_ggg_mpi(struct cmdline_data *,
                                     struct global_data *, bodyptr *,
                                     INTEGER *, INTEGER, INTEGER *, int, int);
global int fcfc_octree_ggg_mpi_prepare(struct cmdline_data *,
                                       struct global_data *);
global int fcfc_octree_ggg_mpi_finalize(struct cmdline_data *);
global int fcfc_octree_ggg_mpi_output_enabled(struct cmdline_data *);
global int fcfc_octree_ggg_mpi_consensus(struct cmdline_data *, int,
                                         const char *);

#endif
