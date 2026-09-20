/* Spherical specialization of the shared spin-2 shear estimator. */
#define OCTREE_SHEAR_SPHERICAL 1
#define prepare_octree_shear_catalogs prepare_octree_shear_sphere_catalogs
#define searchcalc_octree_shear_omp searchcalc_octree_shear_sphere_omp
#include "../shear_sphere_shared/shear_sphere_engine.c"
