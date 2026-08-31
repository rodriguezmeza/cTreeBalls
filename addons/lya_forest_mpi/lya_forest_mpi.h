#ifndef CBALLS_LYA_FOREST_MPI_H
#define CBALLS_LYA_FOREST_MPI_H
#ifdef LYAFORESTMPI
#include <mpi.h>
#include <stdint.h>
int lya_forest_mpi_prepare(struct cmdline_data *, struct global_data *);
int lya_forest_mpi_finalize(struct cmdline_data *);
int lya_forest_mpi_active(void);
int lya_forest_mpi_is_root(void);
int lya_forest_mpi_rank(void);
int lya_forest_mpi_size(void);
int lya_forest_mpi_output_enabled(struct cmdline_data *);
int lya_forest_mpi_consensus(struct cmdline_data *, int, const char *);
int lya_forest_mpi_reduce_reals(struct cmdline_data *, real *, size_t);
int lya_forest_mpi_reduce_integers(struct cmdline_data *, INTEGER *, size_t);
int lya_forest_mpi_reduce_long_doubles(struct cmdline_data *, long double *, size_t);
int lya_forest_mpi_reduce_uint64(struct cmdline_data *, uint64_t *, size_t);
#endif
#endif
