/* Replicated-catalog MPI runtime for the 3D scalar estimators. */

#include <limits.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "cb3d_mpi.h"

#ifdef OCTREE3PCF3DMPI

#define CB3D_MPI_ROOT 0

#if defined(SINGLEPREC)
#define CB3D_MPI_REAL MPI_FLOAT
#else
#define CB3D_MPI_REAL MPI_DOUBLE
#endif

#ifdef LONGINT
#define CB3D_MPI_INTEGER MPI_LONG
#else
#define CB3D_MPI_INTEGER MPI_INT
#endif

/* Per-context rank view; process MPI ownership lives in mpi_runtime.c. */
#define CBALLS_MPI_STATE (cballs_runtime_current()->mpi[CBALLS_MPI_CB3D_MPI])
#define mpi_active (CBALLS_MPI_STATE.active)
#define mpi_rank (CBALLS_MPI_STATE.rank)
#define mpi_size (CBALLS_MPI_STATE.size ? CBALLS_MPI_STATE.size : 1)



static int mpi_error(struct cmdline_data *cmd, const char *operation,
                     int status)
{
    return cballs_mpi_shared_error(cmd, operation, status);
}

int cb3d_mpi_prepare(struct cmdline_data *cmd,
                                   struct global_data *gd)
{
    return cballs_mpi_shared_prepare(&CBALLS_MPI_STATE, cmd, gd, !(!cb3d_is_mpi_method(cmd->searchMethod)));
}

int cb3d_mpi_finalize(struct cmdline_data *cmd)
{
    return cballs_mpi_shared_finalize(&CBALLS_MPI_STATE, cmd);
}

int cb3d_mpi_active(void) { return mpi_active; }
int cb3d_mpi_is_root(void)
{
    return !mpi_active || mpi_rank == CB3D_MPI_ROOT;
}
int cb3d_mpi_rank(void) { return mpi_rank; }
int cb3d_mpi_size(void) { return mpi_size; }

int cb3d_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !cb3d_is_mpi_method(cmd->searchMethod)
        || !mpi_active || mpi_rank == CB3D_MPI_ROOT;
}

int cb3d_mpi_consensus(struct cmdline_data *cmd,
                                     int local_status,
                                     const char *operation)
{
    return cballs_mpi_shared_consensus(&CBALLS_MPI_STATE, cmd, !(!cb3d_is_mpi_method(cmd->searchMethod)), local_status, operation);
}

int cb3d_mpi_reduce_reals(struct cmdline_data *cmd,
                                        real *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == CB3D_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                CB3D_MPI_REAL, MPI_SUM,
                                CB3D_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                CB3D_MPI_REAL, MPI_SUM,
                                CB3D_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI 3D scalar real reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int cb3d_mpi_reduce_integers(struct cmdline_data *cmd,
                                           INTEGER *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == CB3D_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                CB3D_MPI_INTEGER, MPI_SUM,
                                CB3D_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                CB3D_MPI_INTEGER, MPI_SUM,
                                CB3D_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI 3D scalar integer reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif
