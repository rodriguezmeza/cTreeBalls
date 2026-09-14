#ifndef _cballs_mpi_dispatch_h
#define _cballs_mpi_dispatch_h

#ifdef CBALLS_MPI_ENABLED

static inline int cballs_mpi_prepare(struct cmdline_data *cmd,
                                     struct global_data *gd)
{
    int status = SUCCESS;
#ifdef OCTREE3PCF3DMPI
    status = cb3d_mpi_prepare(cmd, gd);
    if (status == FAILURE) return FAILURE;
#endif
#ifdef LYAFORESTMPI
    status = lya_forest_mpi_prepare(cmd, gd);
    if (status == FAILURE) return FAILURE;
#endif
#ifdef KDTREE2BALLSMPI
    status = fcfc_kdtree_2balls_mpi_prepare(cmd, gd);
    if (status == FAILURE) return FAILURE;
#endif
#ifdef BALLTREE2BALLSMPI
    if (status == SUCCESS)
        status = fcfc_balltree_2balls_mpi_prepare(cmd, gd);
#endif
#ifdef OCTREE2BALLSMPI
    if (status == SUCCESS)
        status = fcfc_octree_2balls_mpi_prepare(cmd, gd);
#endif
    return status;
}

static inline int cballs_mpi_finalize(struct cmdline_data *cmd)
{
    int status = SUCCESS;
#ifdef OCTREE3PCF3DMPI
    if (cb3d_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef LYAFORESTMPI
    if (lya_forest_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef KDTREE2BALLSMPI
    if (fcfc_kdtree_2balls_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef BALLTREE2BALLSMPI
    if (fcfc_balltree_2balls_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef OCTREE2BALLSMPI
    if (fcfc_octree_2balls_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
    return status;
}

static inline int cballs_mpi_output_enabled(struct cmdline_data *cmd)
{
    int enabled = TRUE;
#ifdef OCTREE3PCF3DMPI
    enabled = enabled && cb3d_mpi_output_enabled(cmd);
#endif
#ifdef LYAFORESTMPI
    enabled = enabled && lya_forest_mpi_output_enabled(cmd);
#endif
#ifdef KDTREE2BALLSMPI
    enabled = enabled && fcfc_kdtree_2balls_mpi_output_enabled(cmd);
#endif
#ifdef BALLTREE2BALLSMPI
    enabled = enabled && fcfc_balltree_2balls_mpi_output_enabled(cmd);
#endif
#ifdef OCTREE2BALLSMPI
    enabled = enabled && fcfc_octree_2balls_mpi_output_enabled(cmd);
#endif
    return enabled;
}

static inline int cballs_mpi_consensus(struct cmdline_data *cmd,
                                       int local_status,
                                       const char *operation)
{
    int status = local_status;
#ifdef OCTREE3PCF3DMPI
    status = cb3d_mpi_consensus(cmd, status, operation);
#endif
#ifdef LYAFORESTMPI
    status = lya_forest_mpi_consensus(cmd, status, operation);
#endif
#ifdef KDTREE2BALLSMPI
    status = fcfc_kdtree_2balls_mpi_consensus(cmd, status, operation);
#endif
#ifdef BALLTREE2BALLSMPI
    status = fcfc_balltree_2balls_mpi_consensus(cmd, status, operation);
#endif
#ifdef OCTREE2BALLSMPI
    status = fcfc_octree_2balls_mpi_consensus(cmd, status, operation);
#endif
    return status;
}

#endif /* CBALLS_MPI_ENABLED */

#endif /* !_cballs_mpi_dispatch_h */
