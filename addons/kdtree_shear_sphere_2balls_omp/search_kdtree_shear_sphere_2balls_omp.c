/* Full-sky spin-2 specialization over a median KD tree. */
#define OCTREE_SHEAR_SPHERICAL 1
#define SHEAR_SPHERE_BINARY_TWO_BALLS 1
#define SHEAR_SPHERE_BINARY_TREE_BUILD kdtree_shear_sphere_build
#define SHEAR_ENGINE_NAME "kdtree-shear-sphere-2balls-omp"
#define prepare_octree_shear_catalogs \
    prepare_kdtree_shear_sphere_2balls_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_kdtree_shear_sphere_2balls_omp
#include "../octree_shear_omp/search_octree_shear_omp.c"
