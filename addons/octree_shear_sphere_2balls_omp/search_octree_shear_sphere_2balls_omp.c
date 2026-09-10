/* dual-node-style two-cell specialization of the full-sky shear estimator. */
#define OCTREE_SHEAR_SPHERICAL 1
#define OCTREE_SHEAR_SPHERICAL_TWO_BALLS 1
#define SHEAR_ENGINE_NAME "octree-shear-sphere-2balls-omp"
#define prepare_octree_shear_catalogs \
    prepare_octree_shear_sphere_2balls_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_octree_shear_sphere_2balls_omp
#include "../octree_shear_omp/search_octree_shear_omp.c"
