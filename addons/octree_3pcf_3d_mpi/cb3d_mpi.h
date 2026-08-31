#ifndef CBALLS_CB3D_MPI_H
#define CBALLS_CB3D_MPI_H
#ifdef OCTREE3PCF3DMPI
#include <mpi.h>
int cb3d_mpi_prepare(struct cmdline_data *, struct global_data *);
int cb3d_mpi_finalize(struct cmdline_data *);
int cb3d_mpi_active(void);
int cb3d_mpi_is_root(void);
int cb3d_mpi_rank(void);
int cb3d_mpi_size(void);
int cb3d_mpi_output_enabled(struct cmdline_data *);
int cb3d_mpi_consensus(struct cmdline_data *, int, const char *);
int cb3d_mpi_reduce_reals(struct cmdline_data *, real *, size_t);
int cb3d_mpi_reduce_integers(struct cmdline_data *, INTEGER *, size_t);
#endif
#endif
