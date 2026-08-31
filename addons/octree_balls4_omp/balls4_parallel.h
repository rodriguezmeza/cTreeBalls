#ifndef CBALLS_BALLS4_PARALLEL_H
#define CBALLS_BALLS4_PARALLEL_H

static inline bool balls4_distributed(const struct cmdline_data *cmd)
{
#ifdef OCTREEBALLS4MPI
    return cmd->searchMethod && strcmp(cmd->searchMethod, "octree-balls4-mpi") == 0;
#else
    (void)cmd;
    return FALSE;
#endif
}

static inline int balls4_consensus(struct cmdline_data *cmd, int status,
                                   const char *operation)
{
#ifdef OCTREEBALLS4MPI
    if (balls4_distributed(cmd))
        return fcfc_octree_balls4_mpi_consensus(cmd, status, operation);
#endif
    (void)cmd; (void)operation;
    return status;
}

static inline bool balls4_publish(struct cmdline_data *cmd)
{
#ifdef OCTREEBALLS4MPI
    if (balls4_distributed(cmd)) return fcfc_octree_balls4_mpi_is_root();
#endif
    (void)cmd;
    return TRUE;
}

static inline bool balls4_task_owned(struct cmdline_data *cmd, INTEGER task)
{
#ifdef OCTREEBALLS4MPI
    if (balls4_distributed(cmd)) return fcfc_octree_balls4_mpi_task_owned(task);
#endif
    (void)cmd; (void)task;
    return TRUE;
}

static inline int balls4_reduce(struct cmdline_data *cmd, real *values, size_t count)
{
#ifdef OCTREEBALLS4MPI
    if (balls4_distributed(cmd))
        return fcfc_octree_balls4_mpi_reduce_reals(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int balls4_reduce_counts(struct cmdline_data *cmd, INTEGER *values,
                                       size_t count)
{
#ifdef OCTREEBALLS4MPI
    if (balls4_distributed(cmd))
        return fcfc_octree_balls4_mpi_reduce_integers(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}
#endif
