/* Deterministic MPI runtime for the native-octree two-ball estimator. */

#include <limits.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_octree_2balls_mpi.h"

#ifdef OCTREE2BALLSMPI

#define OCTREE_2BALLS_MPI_ROOT 0

#if defined(SINGLEPREC)
#define OCTREE_2BALLS_MPI_REAL MPI_FLOAT
#else
#define OCTREE_2BALLS_MPI_REAL MPI_DOUBLE
#endif

#ifdef LONGINT
#define OCTREE_2BALLS_MPI_INTEGER MPI_LONG
#else
#define OCTREE_2BALLS_MPI_INTEGER MPI_INT
#endif

/* Per-context rank view; process MPI ownership lives in mpi_runtime.c. */
#define CBALLS_MPI_STATE (cballs_runtime_current()->mpi[CBALLS_MPI_FCFC_OCTREE_2BALLS_MPI])
#define mpi_active (CBALLS_MPI_STATE.active)
#define mpi_rank (CBALLS_MPI_STATE.rank)
#define mpi_size (CBALLS_MPI_STATE.size ? CBALLS_MPI_STATE.size : 1)

static bool fcfc_octree_2balls_native_mpi_selected(
        const struct cmdline_data *cmd)
{
    return cmd != NULL && cmd->searchMethod != NULL
        && strcmp(cmd->searchMethod, "octree-2balls-mpi") == 0
        && !cballs_opt_legacy_one_ball(cmd);
}



static int mpi_error(struct cmdline_data *cmd, const char *operation,
                     int status)
{
    return cballs_mpi_shared_error(cmd, operation, status);
}

int fcfc_octree_2balls_mpi_prepare(struct cmdline_data *cmd,
                                   struct global_data *gd)
{
    return cballs_mpi_shared_prepare(&CBALLS_MPI_STATE, cmd, gd, !(!fcfc_octree_2balls_native_mpi_selected(cmd)));
}

int fcfc_octree_2balls_mpi_finalize(struct cmdline_data *cmd)
{
    return cballs_mpi_shared_finalize(&CBALLS_MPI_STATE, cmd);
}

int fcfc_octree_2balls_mpi_active(void) { return mpi_active; }
int fcfc_octree_2balls_mpi_is_root(void)
{
    return !mpi_active || mpi_rank == OCTREE_2BALLS_MPI_ROOT;
}
int fcfc_octree_2balls_mpi_rank(void) { return mpi_rank; }
int fcfc_octree_2balls_mpi_size(void) { return mpi_size; }

int fcfc_octree_2balls_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !fcfc_octree_2balls_native_mpi_selected(cmd)
        || !mpi_active || mpi_rank == OCTREE_2BALLS_MPI_ROOT;
}

int fcfc_octree_2balls_mpi_consensus(struct cmdline_data *cmd,
                                     int local_status,
                                     const char *operation)
{
    return cballs_mpi_shared_consensus(&CBALLS_MPI_STATE, cmd, !(!fcfc_octree_2balls_native_mpi_selected(cmd)), local_status, operation);
}

int fcfc_octree_2balls_mpi_task_owned(INTEGER task)
{
    if (!mpi_active || mpi_size <= 1) return TRUE;
    return (int)(task % (INTEGER)mpi_size) == mpi_rank;
}

int fcfc_octree_2balls_mpi_reduce_reals(struct cmdline_data *cmd,
                                        real *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == OCTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                OCTREE_2BALLS_MPI_REAL, MPI_SUM,
                                OCTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                OCTREE_2BALLS_MPI_REAL, MPI_SUM,
                                OCTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI two-ball real reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int fcfc_octree_2balls_mpi_reduce_integers(struct cmdline_data *cmd,
                                           INTEGER *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == OCTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                OCTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                                OCTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                OCTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                                OCTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI two-ball integer reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* OCTREE2BALLSMPI */
