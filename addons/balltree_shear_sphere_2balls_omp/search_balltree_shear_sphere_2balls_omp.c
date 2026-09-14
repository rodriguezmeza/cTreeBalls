/* Full-sky spin-2 specialization over an FCFC-style PCA ball tree. */
#define OCTREE_SHEAR_SPHERICAL 1
#define SHEAR_SPHERE_BINARY_TWO_BALLS 1
#define SHEAR_SPHERE_BINARY_POSITIONS_UNIT 1
#ifdef BALLS4SCANLEV
#define SHEAR_SPHERE_BINARY_FRONTIER_SCHEDULER 1
#endif
#define SHEAR_SPHERE_BINARY_TREE_HEADER "fcfc_balltree.h"
#define SHEAR_SPHERE_BINARY_TREE_BUILD fcfc_balltree_build_shear_sphere
#define SHEAR_SPHERE_BINARY_TREE_FRONTIER fcfc_balltree_frontier
#define SHEAR_SPHERE_BINARY_TREE_FREE fcfc_balltree_free
#define SHEAR_ENGINE_NAME "balltree-shear-sphere-2balls-omp"
#define prepare_octree_shear_catalogs \
    prepare_balltree_shear_sphere_2balls_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_balltree_shear_sphere_2balls_omp
#include "../shear_sphere_shared/shear_sphere_engine.c"
