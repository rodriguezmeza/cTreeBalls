#ifndef _fcfc_octree_ggg_mpi_h
#define _fcfc_octree_ggg_mpi_h

#ifdef OCTREEGGGMPI

#include <mpi.h>
#include <stdint.h>

typedef struct {
    MPI_Win window;
    uint64_t *counter;
    uint64_t task_count;
    uint64_t step;
    uint64_t local_next;
    int ready;
} fcfc_octree_ggg_mpi_scheduler;

int fcfc_octree_ggg_mpi_prepare(struct cmdline_data *, struct global_data *);
int fcfc_octree_ggg_mpi_finalize(struct cmdline_data *);
int fcfc_octree_ggg_mpi_active(void);
int fcfc_octree_ggg_mpi_is_root(void);
int fcfc_octree_ggg_mpi_rank(void);
int fcfc_octree_ggg_mpi_size(void);
int fcfc_octree_ggg_mpi_output_enabled(struct cmdline_data *);
int fcfc_octree_ggg_mpi_consensus(struct cmdline_data *, int, const char *);

int fcfc_octree_ggg_mpi_scheduler_init(
    struct cmdline_data *, fcfc_octree_ggg_mpi_scheduler *, INTEGER, INTEGER);
int fcfc_octree_ggg_mpi_scheduler_claim(
    struct cmdline_data *, fcfc_octree_ggg_mpi_scheduler *, INTEGER *, INTEGER *);
int fcfc_octree_ggg_mpi_scheduler_destroy(
    struct cmdline_data *, fcfc_octree_ggg_mpi_scheduler *);

int fcfc_octree_ggg_mpi_reduce_reals(struct cmdline_data *, real *, size_t);
int fcfc_octree_ggg_mpi_reduce_integers(struct cmdline_data *, INTEGER *, size_t);

#endif /* OCTREEGGGMPI */

#endif /* !_fcfc_octree_ggg_mpi_h */
