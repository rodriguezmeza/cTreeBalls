#ifndef _protodefs_balltree_shear_sphere_2balls_mpi_h
#define _protodefs_balltree_shear_sphere_2balls_mpi_h

#include "protodefs_octree_shear_sphere_omp.h"
#include "fcfc_balltree_shear_sphere_2balls_mpi.h"

#define BALLTREESHEARSPHERE2BALLSMPIMETHOD 205

global int prepare_balltree_shear_sphere_2balls_mpi_catalogs(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btable, INTEGER *nbody);

global int searchcalc_balltree_shear_sphere_2balls_mpi(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btable, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2, int cat3);

#endif
