/* Full-sky spin-2 specialization over a median KD tree. */
#define OCTREE_SHEAR_SPHERICAL 1
#define SHEAR_SPHERE_BINARY_TWO_BALLS 1
#define SHEAR_SPHERE_BINARY_POSITIONS_UNIT 1
/* Keep these arithmetic/traversal specializations private to this engine. */
#define SHEAR_SPHERE_FAST_KERNEL 1
#define SHEAR_SPHERE_BINARY_GEODESIC_PAIRS 1
#ifdef BALLS4SCANLEV
#define SHEAR_SPHERE_BINARY_FRONTIER_SCHEDULER 1
#define SHEAR_SPHERE_PIVOT_REUSE 1
#define SHEAR_SPHERE_BINARY_PIVOT_REUSE 1
#endif
#define SHEAR_SPHERE_BINARY_TREE_HEADER "kdtree_shear_sphere_tree.h"
#define SHEAR_SPHERE_BINARY_TREE_BUILD kdtree_shear_sphere_build
#define SHEAR_SPHERE_BINARY_TREE_FRONTIER kdtree_shear_sphere_frontier
#define SHEAR_SPHERE_BINARY_TREE_FREE kdtree_shear_sphere_free
#define SHEAR_ENGINE_NAME "kdtree-shear-sphere-2balls-omp"
#define prepare_octree_shear_catalogs \
    prepare_kdtree_shear_sphere_2balls_catalogs
#define searchcalc_octree_shear_omp \
    searchcalc_kdtree_shear_sphere_2balls_omp
#include "../../support/shear_kernel/search_octree_shear_omp.c"
