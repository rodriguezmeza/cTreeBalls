#ifndef _cballs_mpi_dispatch_h
#define _cballs_mpi_dispatch_h

#ifdef CBALLS_MPI_ENABLED

static inline int cballs_mpi_prepare(struct cmdline_data *cmd,
                                     struct global_data *gd)
{
    int status = SUCCESS;
#ifdef BALLTREEMPI
    status = fcfc_balltree_mpi_prepare(cmd, gd);
#endif
#ifdef BALLTREE2BALLSMPI
    if (status == SUCCESS)
        status = fcfc_balltree_2balls_mpi_prepare(cmd, gd);
#endif
#ifdef OCTREE2BALLSMPI
    if (status == SUCCESS)
        status = fcfc_octree_2balls_mpi_prepare(cmd, gd);
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    if (status == SUCCESS)
        status = fcfc_balltree_2balls_mpi_3pcf_prepare(cmd, gd);
#endif
#ifdef OCTREEGGGMPI
    if (status == SUCCESS) status = fcfc_octree_ggg_mpi_prepare(cmd, gd);
#endif
    return status;
}

static inline int cballs_mpi_finalize(struct cmdline_data *cmd)
{
    int status = SUCCESS;
#ifdef BALLTREEMPI
    if (fcfc_balltree_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef BALLTREE2BALLSMPI
    if (fcfc_balltree_2balls_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef OCTREE2BALLSMPI
    if (fcfc_octree_2balls_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    if (fcfc_balltree_2balls_mpi_3pcf_finalize(cmd) == FAILURE)
        status = FAILURE;
#endif
#ifdef OCTREEGGGMPI
    if (fcfc_octree_ggg_mpi_finalize(cmd) == FAILURE) status = FAILURE;
#endif
    return status;
}

static inline int cballs_mpi_output_enabled(struct cmdline_data *cmd)
{
    int enabled = TRUE;
#ifdef BALLTREEMPI
    enabled = enabled && fcfc_balltree_mpi_output_enabled(cmd);
#endif
#ifdef BALLTREE2BALLSMPI
    enabled = enabled && fcfc_balltree_2balls_mpi_output_enabled(cmd);
#endif
#ifdef OCTREE2BALLSMPI
    enabled = enabled && fcfc_octree_2balls_mpi_output_enabled(cmd);
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    enabled = enabled
        && fcfc_balltree_2balls_mpi_3pcf_output_enabled(cmd);
#endif
#ifdef OCTREEGGGMPI
    enabled = enabled && fcfc_octree_ggg_mpi_output_enabled(cmd);
#endif
    return enabled;
}

static inline int cballs_mpi_consensus(struct cmdline_data *cmd,
                                       int local_status,
                                       const char *operation)
{
    int status = local_status;
#ifdef BALLTREEMPI
    status = fcfc_balltree_mpi_consensus(cmd, status, operation);
#endif
#ifdef BALLTREE2BALLSMPI
    status = fcfc_balltree_2balls_mpi_consensus(cmd, status, operation);
#endif
#ifdef OCTREE2BALLSMPI
    status = fcfc_octree_2balls_mpi_consensus(cmd, status, operation);
#endif
#ifdef BALLTREE2BALLSMPI3PCF
    status = fcfc_balltree_2balls_mpi_3pcf_consensus(
        cmd, status, operation);
#endif
#ifdef OCTREEGGGMPI
    status = fcfc_octree_ggg_mpi_consensus(cmd, status, operation);
#endif
    return status;
}

#endif /* CBALLS_MPI_ENABLED */

#endif /* !_cballs_mpi_dispatch_h */
