#ifndef CTREEBALLS_PROTODEFS_BALLTREE_SHEAR_SPHERE_2BALLS_OMP_H
#define CTREEBALLS_PROTODEFS_BALLTREE_SHEAR_SPHERE_2BALLS_OMP_H

#define BALLTREESHEARSPHERE2BALLSOMPMETHOD 202

global int prepare_balltree_shear_sphere_2balls_catalogs(
        struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *);
global int searchcalc_balltree_shear_sphere_2balls_omp(
        struct cmdline_data *, struct global_data *, bodyptr *, INTEGER *,
        INTEGER, INTEGER *, int, int, int);

#endif
