/* Spherical specialization of the shared spin-2 shear estimator. */
#define OCTREE_SHEAR_SPHERICAL 1
#define prepare_octree_shear_catalogs prepare_octree_shear_sphere_catalogs
#define searchcalc_octree_shear_omp searchcalc_octree_shear_sphere_omp
#include "../octree_shear_omp/search_octree_shear_omp.c"
