#ifndef _fcfc_kdtree_2balls_mpi_h
#define _fcfc_kdtree_2balls_mpi_h

#ifdef KDTREE2BALLSMPI
#include <mpi.h>

int fcfc_kdtree_2balls_mpi_prepare(
    struct cmdline_data *, struct global_data *);
int fcfc_kdtree_2balls_mpi_finalize(struct cmdline_data *);
int fcfc_kdtree_2balls_mpi_is_root(void);
int fcfc_kdtree_2balls_mpi_size(void);
int fcfc_kdtree_2balls_mpi_output_enabled(struct cmdline_data *);
int fcfc_kdtree_2balls_mpi_consensus(
    struct cmdline_data *, int, const char *);
int fcfc_kdtree_2balls_mpi_task_owned(INTEGER);
int fcfc_kdtree_2balls_mpi_reduce_reals(
    struct cmdline_data *, real *, size_t);
int fcfc_kdtree_2balls_mpi_reduce_integers(
    struct cmdline_data *, INTEGER *, size_t);

#endif
#endif /* !_fcfc_kdtree_2balls_mpi_h */
