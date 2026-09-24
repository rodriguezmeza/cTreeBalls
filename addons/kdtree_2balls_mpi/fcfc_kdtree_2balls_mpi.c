/* MPI lifecycle and deterministic reductions for KD two-ball scans. */

#include <limits.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_kdtree_2balls_mpi.h"

#ifdef KDTREE2BALLSMPI

#define KDTREE_2BALLS_MPI_ROOT 0

#if defined(SINGLEPREC)
#define KDTREE_2BALLS_MPI_REAL MPI_FLOAT
#else
#define KDTREE_2BALLS_MPI_REAL MPI_DOUBLE
#endif

#ifdef LONGINT
#define KDTREE_2BALLS_MPI_INTEGER MPI_LONG
#else
#define KDTREE_2BALLS_MPI_INTEGER MPI_INT
#endif

/* Per-context rank view; process MPI ownership lives in mpi_runtime.c. */
#define CBALLS_MPI_STATE (cballs_runtime_current()->mpi[CBALLS_MPI_FCFC_KDTREE_2BALLS_MPI])
#define mpi_active (CBALLS_MPI_STATE.active)
#define mpi_rank (CBALLS_MPI_STATE.rank)
#define mpi_size (CBALLS_MPI_STATE.size ? CBALLS_MPI_STATE.size : 1)



static int kdtree_2balls_mpi_error(
        struct cmdline_data *cmd, const char *operation, int status)
{
    return cballs_mpi_shared_error(cmd, operation, status);
}

static bool kdtree_2balls_mpi_method(const struct cmdline_data *cmd)
{
    return cmd != NULL && cmd->searchMethod != NULL
        && strcmp(cmd->searchMethod, "kdtree-2balls-mpi") == 0;
}

int fcfc_kdtree_2balls_mpi_prepare(
        struct cmdline_data *cmd, struct global_data *gd)
{
    return cballs_mpi_shared_prepare(&CBALLS_MPI_STATE, cmd, gd, !(!kdtree_2balls_mpi_method(cmd)));
}

int fcfc_kdtree_2balls_mpi_finalize(struct cmdline_data *cmd)
{
    return cballs_mpi_shared_finalize(&CBALLS_MPI_STATE, cmd);
}

int fcfc_kdtree_2balls_mpi_is_root(void)
{
    return !mpi_active || mpi_rank == KDTREE_2BALLS_MPI_ROOT;
}

int fcfc_kdtree_2balls_mpi_size(void)
{
    return mpi_size;
}

int fcfc_kdtree_2balls_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !kdtree_2balls_mpi_method(cmd) || !mpi_active
        || mpi_rank == KDTREE_2BALLS_MPI_ROOT;
}

int fcfc_kdtree_2balls_mpi_consensus(
        struct cmdline_data *cmd, int local_status, const char *operation)
{
    return cballs_mpi_shared_consensus(&CBALLS_MPI_STATE, cmd, !(!kdtree_2balls_mpi_method(cmd)), local_status, operation);
}

int fcfc_kdtree_2balls_mpi_task_owned(INTEGER task)
{
    return !mpi_active || mpi_size <= 1
        || (int)(task % (INTEGER)mpi_size) == mpi_rank;
}

int fcfc_kdtree_2balls_mpi_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == KDTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(
                MPI_IN_PLACE, values + offset, chunk,
                KDTREE_2BALLS_MPI_REAL, MPI_SUM,
                KDTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(
                values + offset, NULL, chunk,
                KDTREE_2BALLS_MPI_REAL, MPI_SUM,
                KDTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return kdtree_2balls_mpi_error(
                cmd, "MPI KD two-ball real reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int fcfc_kdtree_2balls_mpi_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == KDTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(
                MPI_IN_PLACE, values + offset, chunk,
                KDTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                KDTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(
                values + offset, NULL, chunk,
                KDTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                KDTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return kdtree_2balls_mpi_error(
                cmd, "MPI KD two-ball integer reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* KDTREE2BALLSMPI */
