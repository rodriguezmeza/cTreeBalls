#ifndef CBALLS_KDTREE_PARALLEL_H
#define CBALLS_KDTREE_PARALLEL_H

static inline bool kdtree_distributed(const struct cmdline_data *cmd)
{
    if (cmd == NULL || cmd->searchMethod == NULL) return FALSE;
#ifdef KDTREEMPI
    if (strcmp(cmd->searchMethod, "kdtree-mpi") == 0) return TRUE;
#endif
#ifdef KDTREE2BALLSMPI
    if (strcmp(cmd->searchMethod, "kdtree-2balls-mpi") == 0
        && cballs_opt_legacy_one_ball(cmd))
        return TRUE;
#endif
    return FALSE;
}

static inline bool kdtree_two_ball_legacy_distributed(
        const struct cmdline_data *cmd)
{
#ifdef KDTREE2BALLSMPI
    return cmd != NULL && cmd->searchMethod != NULL
        && strcmp(cmd->searchMethod, "kdtree-2balls-mpi") == 0
        && cballs_opt_legacy_one_ball(cmd);
#else
    (void)cmd;
    return FALSE;
#endif
}

static inline int kdtree_consensus(struct cmdline_data *cmd, int status,
                                   const char *operation)
{
#ifdef KDTREE2BALLSMPI
    if (kdtree_two_ball_legacy_distributed(cmd))
        return fcfc_kdtree_2balls_mpi_consensus(cmd, status, operation);
#endif
#ifdef KDTREEMPI
    if (kdtree_distributed(cmd))
        return fcfc_kdtree_mpi_consensus(cmd, status, operation);
#endif
    (void)cmd; (void)operation;
    return status;
}

static inline bool kdtree_publish(const struct cmdline_data *cmd)
{
#ifdef KDTREE2BALLSMPI
    if (kdtree_two_ball_legacy_distributed(cmd))
        return fcfc_kdtree_2balls_mpi_is_root();
#endif
#ifdef KDTREEMPI
    if (kdtree_distributed(cmd)) return fcfc_kdtree_mpi_is_root();
#endif
    (void)cmd;
    return TRUE;
}

static inline bool kdtree_task_owned(const struct cmdline_data *cmd,
                                     INTEGER task)
{
#ifdef KDTREE2BALLSMPI
    if (kdtree_two_ball_legacy_distributed(cmd))
        return fcfc_kdtree_2balls_mpi_task_owned(task);
#endif
#ifdef KDTREEMPI
    if (kdtree_distributed(cmd)) return fcfc_kdtree_mpi_task_owned(task);
#endif
    (void)cmd; (void)task;
    return TRUE;
}

static inline int kdtree_reduce(struct cmdline_data *cmd, real *values,
                                size_t count)
{
#ifdef KDTREE2BALLSMPI
    if (kdtree_two_ball_legacy_distributed(cmd))
        return fcfc_kdtree_2balls_mpi_reduce_reals(cmd, values, count);
#endif
#ifdef KDTREEMPI
    if (kdtree_distributed(cmd))
        return fcfc_kdtree_mpi_reduce_reals(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

static inline int kdtree_reduce_counts(struct cmdline_data *cmd,
                                       INTEGER *values, size_t count)
{
#ifdef KDTREE2BALLSMPI
    if (kdtree_two_ball_legacy_distributed(cmd))
        return fcfc_kdtree_2balls_mpi_reduce_integers(cmd, values, count);
#endif
#ifdef KDTREEMPI
    if (kdtree_distributed(cmd))
        return fcfc_kdtree_mpi_reduce_integers(cmd, values, count);
#endif
    (void)cmd; (void)values; (void)count;
    return SUCCESS;
}

#endif
