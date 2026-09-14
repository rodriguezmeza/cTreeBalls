/* Runtime IDs for the maintained active add-on profile. */
#ifndef _startrun_include_11_h
#define _startrun_include_11_h

#ifdef KDTREE2BALLSOMP
#include "startrun_kdtree_2balls_omp_11.h"
#endif
#ifdef KDTREE2BALLSMPI
#include "startrun_kdtree_2balls_mpi_11.h"
#endif
#ifdef BALLTREE2BALLSOMP
#include "startrun_balltree_2balls_omp_11.h"
#endif
#ifdef BALLTREE2BALLSMPI
#include "startrun_balltree_2balls_mpi_11.h"
#endif
#ifdef OCTREE2BALLSOMP
#include "startrun_octree_2balls_omp_11.h"
#endif
#ifdef OCTREE2BALLSMPI
#include "startrun_octree_2balls_mpi_11.h"
#endif

#ifdef OCTREESHEARSPHERE2BALLSOMP
#include "startrun_octree_shear_sphere_2balls_omp_11.h"
#endif
#ifdef KDTREESHEARSPHERE2BALLSOMP
#include "startrun_kdtree_shear_sphere_2balls_omp_11.h"
#endif
#ifdef BALLTREESHEARSPHERE2BALLSOMP
#include "startrun_balltree_shear_sphere_2balls_omp_11.h"
#endif

#ifdef KDTREEBOXOMP
#include "startrun_kdtree_box_omp_11.h"
#endif
#ifdef NEIGHBORBOXESOMP
#include "startrun_neighbor_boxes_omp_11.h"
#endif

#ifdef OCTREE3PCF3DOMP
#include "startrun_octree_3pcf_3d_omp_11.h"
#endif
#ifdef OCTREE3PCF3DMPI
#include "startrun_octree_3pcf_3d_mpi_11.h"
#endif

#ifdef LYAFORESTOMP
#include "startrun_lya_forest_omp_11.h"
#endif
#ifdef LYAFORESTMPI
#include "startrun_lya_forest_mpi_11.h"
#endif

/* The core octree-sincos method remains part of every standard build. */
#ifdef OCTREESINCOSOMP
#include "startrun_octree_sincos_omp_11.h"
#endif

#endif
