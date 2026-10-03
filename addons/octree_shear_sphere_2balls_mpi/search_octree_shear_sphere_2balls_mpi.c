/* MPI specialization of the full-sky shear two-ball estimator. */
#define OCTREE_SHEAR_SPHERICAL 1
#define SHEAR_SPHERE_FAST_KERNEL 1
#define OCTREE_SHEAR_SPHERICAL_TWO_BALLS 1
#define SHEAR_SPHERE_BODY_POSITIONS_UNIT 1
#ifdef BALLS4SCANLEV
#define SHEAR_SPHERE_NATIVE_FRONTIER_SCHEDULER 1
#endif
#define SHEAR_MPI_ENABLED 1
#define SHEAR_ENGINE_NAME "octree-shear-sphere-2balls-mpi"
#define prepare_octree_shear_catalogs \
    prepare_octree_shear_sphere_2balls_mpi_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_octree_shear_sphere_2balls_mpi
#include "../../support/shear_kernel/search_octree_shear_omp.c"
