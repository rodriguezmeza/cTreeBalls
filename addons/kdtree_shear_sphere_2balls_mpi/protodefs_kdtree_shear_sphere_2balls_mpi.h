#ifndef _protodefs_kdtree_shear_sphere_2balls_mpi_h
#define _protodefs_kdtree_shear_sphere_2balls_mpi_h

#include "protodefs_octree_shear_sphere_omp.h"
#include "fcfc_kdtree_shear_sphere_2balls_mpi.h"

#define KDTREESHEARSPHERE2BALLSMPIMETHOD 204

global int prepare_kdtree_shear_sphere_2balls_mpi_catalogs(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btable, INTEGER *nbody);

global int searchcalc_kdtree_shear_sphere_2balls_mpi(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btable, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2, int cat3);

#endif
