#ifndef CBALLS_CB3D_PARALLEL_H
#define CBALLS_CB3D_PARALLEL_H

static inline int cb3d_is_mpi_method(const char *name)
{
    return name != NULL && (strcmp(name, "octree-3pcf-3d-mpi") == 0
                            || strcmp(name, "octree-ggg-3d-mpi") == 0);
}

static inline int cb3d_is_method(const char *name)
{
    return cb3d_is_mpi_method(name)
        || (name != NULL && (strcmp(name, "octree-3pcf-3d-omp") == 0
                            || strcmp(name, "octree-ggg-3d-omp") == 0));
}

static inline int cb3d_parallel_consensus(struct cmdline_data *cmd, int status,
                                         const char *operation)
{
#ifdef OCTREE3PCF3DMPI
    return cb3d_mpi_consensus(cmd, status, operation);
#else
    (void)cmd; (void)operation;
    return status;
#endif
}

static inline int cb3d_parallel_publish(struct cmdline_data *cmd)
{
#ifdef OCTREE3PCF3DMPI
    return cb3d_mpi_output_enabled(cmd);
#else
    (void)cmd;
    return TRUE;
#endif
}

static inline INTEGER cb3d_parallel_first(struct cmdline_data *cmd)
{
#ifdef OCTREE3PCF3DMPI
    if (cb3d_is_mpi_method(cmd->searchMethod)) return cb3d_mpi_rank();
#endif
    (void)cmd;
    return 0;
}

static inline INTEGER cb3d_parallel_stride(struct cmdline_data *cmd)
{
#ifdef OCTREE3PCF3DMPI
    if (cb3d_is_mpi_method(cmd->searchMethod)) return cb3d_mpi_size();
#endif
    (void)cmd;
    return 1;
}

static inline int cb3d_parallel_reduce_reals(struct cmdline_data *cmd,
                                             REAL *values, size_t count)
{
#ifdef OCTREE3PCF3DMPI
    if (cb3d_is_mpi_method(cmd->searchMethod))
        return cb3d_mpi_reduce_reals(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int cb3d_parallel_reduce_integers(struct cmdline_data *cmd,
                                                INTEGER *values, size_t count)
{
#ifdef OCTREE3PCF3DMPI
    if (cb3d_is_mpi_method(cmd->searchMethod))
        return cb3d_mpi_reduce_integers(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}
#endif
