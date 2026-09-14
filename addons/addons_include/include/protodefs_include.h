/* Prototypes for the maintained active add-on profile. */
#ifndef _protodefs_include_h
#define _protodefs_include_h

#ifdef KDTREE2BALLSOMP
#include "protodefs_kdtree_2balls_omp.h"
#endif
#ifdef KDTREE2BALLSMPI
#include "protodefs_kdtree_2balls_mpi.h"
#endif
#ifdef BALLTREE2BALLSOMP
#include "protodefs_balltree_2balls_omp.h"
#endif
#ifdef BALLTREE2BALLSMPI
#include "protodefs_balltree_2balls_mpi.h"
#endif
#ifdef OCTREE2BALLSOMP
#include "protodefs_octree_2balls_omp.h"
#endif
#ifdef OCTREE2BALLSMPI
#include "protodefs_octree_2balls_mpi.h"
#endif

#ifdef OCTREESHEARSPHERE2BALLSOMP
#include "protodefs_octree_shear_sphere_2balls_omp.h"
#endif
#ifdef KDTREESHEARSPHERE2BALLSOMP
#include "protodefs_kdtree_shear_sphere_2balls_omp.h"
#endif
#ifdef BALLTREESHEARSPHERE2BALLSOMP
#include "protodefs_balltree_shear_sphere_2balls_omp.h"
#endif

#ifdef KDTREEBOXOMP
#include "protodefs_kdtree_box_omp.h"
#endif
#ifdef NEIGHBORBOXESOMP
#include "protodefs_neighbor_boxes_omp.h"
#endif

#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
#include "protodefs_octree_3pcf_3d_omp.h"
#endif
#ifdef OCTREE3PCF3DMPI
#include "cb3d_mpi.h"
#endif
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
#include "cb3d_parallel.h"
#endif

#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
#include "protodefs_lya_forest_omp.h"
#endif
#ifdef LYAFORESTMPI
#include "lya_forest_mpi.h"
#endif

#ifdef COSMOLIB
#include "protodefs_cosmolib.h"
#endif
#ifdef PXD
#include "protodefs_pxd.h"
#endif
#ifdef OCTREESINCOSOMP
#include "protodefs_octree_sincos_omp.h"
#endif

#endif
