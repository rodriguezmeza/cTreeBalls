/* FCFC-style MPI support for the cTreeBalls ball-tree addon. */

#ifndef _fcfc_balltree_mpi_h
#define _fcfc_balltree_mpi_h

#ifdef BALLTREEMPI

#include <stdint.h>
#include <mpi.h>

typedef struct {
    MPI_Win window;
    uint64_t *counter;
    uint64_t local_next;
    uint64_t task_count;
    uint64_t step;
    int ready;
} fcfc_balltree_mpi_scheduler;

int fcfc_balltree_mpi_prepare(struct cmdline_data *, struct global_data *);
int fcfc_balltree_mpi_finalize(struct cmdline_data *);
int fcfc_balltree_mpi_active(void);
int fcfc_balltree_mpi_is_root(void);
int fcfc_balltree_mpi_rank(void);
int fcfc_balltree_mpi_size(void);
int fcfc_balltree_mpi_output_enabled(struct cmdline_data *);
int fcfc_balltree_mpi_consensus(struct cmdline_data *, int, const char *);

int fcfc_balltree_mpi_broadcast_catalogs(struct cmdline_data *,
                                         struct global_data *, bodyptr *,
                                         INTEGER *, int, int);
int fcfc_balltree_mpi_build(struct cmdline_data *, struct global_data *,
                            bodyptr, INTEGER, int, fcfc_balltreeptr *);
int fcfc_balltree_mpi_frontier(struct cmdline_data *, fcfc_balltreeptr,
                               INTEGER, INTEGER **, INTEGER *);

int fcfc_balltree_mpi_scheduler_init(struct cmdline_data *,
                                     fcfc_balltree_mpi_scheduler *,
                                     INTEGER, INTEGER);
int fcfc_balltree_mpi_scheduler_claim(struct cmdline_data *,
                                      fcfc_balltree_mpi_scheduler *,
                                      INTEGER *, INTEGER *);
int fcfc_balltree_mpi_scheduler_destroy(struct cmdline_data *,
                                        fcfc_balltree_mpi_scheduler *);

int fcfc_balltree_mpi_reduce_reals(struct cmdline_data *, real *, size_t);
int fcfc_balltree_mpi_reduce_integers(struct cmdline_data *, INTEGER *, size_t);

#endif /* BALLTREEMPI */

#endif /* !_fcfc_balltree_mpi_h */
