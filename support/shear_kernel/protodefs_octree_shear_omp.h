#ifndef _protodefs_octree_shear_omp_h
#define _protodefs_octree_shear_omp_h

#define OCTREESHEARMETHOD 173

global int prepare_octree_shear_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btable, INTEGER *nbody);

global int searchcalc_octree_shear_omp(struct cmdline_data *cmd,
                                       struct global_data *gd,
                                       bodyptr *btable, INTEGER *nbody,
                                       INTEGER ipmin, INTEGER *ipmax,
                                       int cat1, int cat2, int cat3);

#endif
