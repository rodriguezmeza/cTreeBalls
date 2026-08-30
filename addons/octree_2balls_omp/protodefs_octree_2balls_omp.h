#ifndef _protodefs_octree_2balls_omp_h
#define _protodefs_octree_2balls_omp_h

#define OCTREE2BALLSMETHOD 176

global int searchcalc_octree_2balls_omp(struct cmdline_data *,
                                        struct global_data *, bodyptr *,
                                        INTEGER *, INTEGER, INTEGER *,
                                        int, int);

#endif /* !_protodefs_octree_2balls_omp_h */
