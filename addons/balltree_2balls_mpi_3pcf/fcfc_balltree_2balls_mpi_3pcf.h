#ifndef _fcfc_balltree_2balls_mpi_3pcf_h
#define _fcfc_balltree_2balls_mpi_3pcf_h

#ifdef BALLTREE2BALLSMPI3PCF

#include <mpi.h>

int fcfc_balltree_2balls_mpi_3pcf_prepare(struct cmdline_data *,
                                           struct global_data *);
int fcfc_balltree_2balls_mpi_3pcf_finalize(struct cmdline_data *);
int fcfc_balltree_2balls_mpi_3pcf_active(void);
int fcfc_balltree_2balls_mpi_3pcf_is_root(void);
int fcfc_balltree_2balls_mpi_3pcf_rank(void);
int fcfc_balltree_2balls_mpi_3pcf_size(void);
int fcfc_balltree_2balls_mpi_3pcf_output_enabled(struct cmdline_data *);
int fcfc_balltree_2balls_mpi_3pcf_consensus(struct cmdline_data *, int,
                                            const char *);
int fcfc_balltree_2balls_mpi_3pcf_task_owned(INTEGER);
int fcfc_balltree_2balls_mpi_3pcf_reduce_reals(struct cmdline_data *,
                                               real *, size_t);
int fcfc_balltree_2balls_mpi_3pcf_reduce_integers(struct cmdline_data *,
                                                  INTEGER *, size_t);

#endif /* BALLTREE2BALLSMPI3PCF */

#endif /* !_fcfc_balltree_2balls_mpi_3pcf_h */
