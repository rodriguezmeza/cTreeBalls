#ifndef CTREEBALLS_PROTODEFS_KDTREE_SHEAR_SPHERE_2BALLS_OMP_H
#define CTREEBALLS_PROTODEFS_KDTREE_SHEAR_SPHERE_2BALLS_OMP_H

#define KDTREESHEARSPHERE2BALLSOMPMETHOD 201

global int prepare_kdtree_shear_sphere_2balls_catalogs(
        struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *);
global int searchcalc_kdtree_shear_sphere_2balls_omp(
        struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *,
        INTEGER, INTEGER *, int, int, int);

#endif
