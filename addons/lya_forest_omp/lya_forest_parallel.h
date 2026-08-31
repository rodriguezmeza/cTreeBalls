#ifndef CBALLS_LYA_FOREST_PARALLEL_H
#define CBALLS_LYA_FOREST_PARALLEL_H
#include <stdint.h>

static inline int lya_parallel_consensus(struct cmdline_data *cmd, int status,
                                         const char *operation)
{
#ifdef LYAFORESTMPI
    return lya_forest_mpi_consensus(cmd, status, operation);
#else
    (void)cmd; (void)operation;
    return status;
#endif
}

static inline int lya_parallel_publish(struct cmdline_data *cmd)
{
#ifdef LYAFORESTMPI
    return lya_forest_mpi_output_enabled(cmd);
#else
    (void)cmd;
    return TRUE;
#endif
}

static inline size_t lya_parallel_first(struct cmdline_data *cmd)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return (size_t)lya_forest_mpi_rank();
#endif
    (void)cmd;
    return 0;
}

static inline size_t lya_parallel_stride(struct cmdline_data *cmd)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return (size_t)lya_forest_mpi_size();
#endif
    (void)cmd;
    return 1;
}

static inline int lya_parallel_reduce_reals(struct cmdline_data *cmd,
                                                REAL *values, size_t count)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return lya_forest_mpi_reduce_reals(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int lya_parallel_reduce_integers(struct cmdline_data *cmd,
                                                INTEGER *values, size_t count)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return lya_forest_mpi_reduce_integers(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int lya_parallel_reduce_long_doubles(struct cmdline_data *cmd,
                                                long double *values, size_t count)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return lya_forest_mpi_reduce_long_doubles(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int lya_parallel_reduce_uint64(struct cmdline_data *cmd,
                                                uint64_t *values, size_t count)
{
#ifdef LYAFORESTMPI
    if (lya_forest_is_mpi_method(cmd->searchMethod))
        return lya_forest_mpi_reduce_uint64(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}
#endif
