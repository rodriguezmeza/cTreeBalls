/* dual-node-style two-cell specialization of the full-sky shear estimator. */
#define OCTREE_SHEAR_SPHERICAL 1
#define OCTREE_SHEAR_SPHERICAL_TWO_BALLS 1
#define SHEAR_SPHERE_BODY_POSITIONS_UNIT 1
#ifdef BALLS4SCANLEV
#define SHEAR_SPHERE_NATIVE_FRONTIER_SCHEDULER 1
#endif
#define SHEAR_ENGINE_NAME "octree-shear-sphere-2balls-omp"
#define prepare_octree_shear_catalogs \
    prepare_octree_shear_sphere_2balls_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_octree_shear_sphere_2balls_omp
#include "../shear_sphere_shared/shear_sphere_engine.c"
