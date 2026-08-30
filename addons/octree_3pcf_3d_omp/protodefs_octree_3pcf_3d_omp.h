// Use:
//#include "protodefs_octree_3pcf_3d_omp.h"

#ifndef _protodefs_octree_3pcf_3d_omp_h
#define _protodefs_octree_3pcf_3d_omp_h

global int searchcalc_octree_3pcf_3d_omp(struct cmdline_data* cmd,
                                         struct global_data* gd,
                                         bodyptr *btable, INTEGER *nbody,
                                         INTEGER ipmin, INTEGER *ipmax,
                                         int cat1, int cat2);

global int cb3d_prepare_survey_catalogs(struct cmdline_data* cmd,
                                        struct global_data* gd,
                                        bodyptr *btable, INTEGER *nbody,
                                        int data_cat, int random_cat,
                                        int *dmr_cat, int *normalized_random_cat,
                                        REAL *random_scale);

global int cb3d_prepare_common_frame(struct cmdline_data* cmd,
                                     struct global_data* gd,
                                     bodyptr *btable, INTEGER *nbody);

global int searchcalc_octree_3pcf_3d_omp_survey(
    struct cmdline_data* cmd, struct global_data* gd,
    bodyptr *btable, INTEGER *nbody,
    int dmr_cat, int normalized_random_cat, REAL random_scale);

#endif // !_protodefs_octree_3pcf_3d_omp_h
