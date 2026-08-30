#ifndef _protodefs_balltree_omp_h
#define _protodefs_balltree_omp_h

global int searchcalc_balltree_omp(struct cmdline_data *,
                                   struct global_data *, bodyptr *,
                                   INTEGER *, INTEGER, INTEGER *, int, int);
#ifdef BALLTREEMPI
global int searchcalc_balltree_mpi(struct cmdline_data *,
                                   struct global_data *, bodyptr *,
                                   INTEGER *, INTEGER, INTEGER *, int, int);
global int fcfc_balltree_mpi_prepare(struct cmdline_data *,
                                     struct global_data *);
global int fcfc_balltree_mpi_finalize(struct cmdline_data *);
global int fcfc_balltree_mpi_output_enabled(struct cmdline_data *);
global int fcfc_balltree_mpi_consensus(struct cmdline_data *, int,
                                       const char *);
#endif

#endif /* !_protodefs_balltree_omp_h */
