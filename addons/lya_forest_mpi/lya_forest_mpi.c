/* Replicated-catalog MPI runtime for the Ly-alpha estimators. */

#include <limits.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "lya_forest_mpi.h"

#ifdef LYAFORESTMPI

#define LYA_FOREST_MPI_ROOT 0

#if defined(SINGLEPREC)
#define LYA_FOREST_MPI_REAL MPI_FLOAT
#else
#define LYA_FOREST_MPI_REAL MPI_DOUBLE
#endif

#ifdef LONGINT
#define LYA_FOREST_MPI_INTEGER MPI_LONG
#else
#define LYA_FOREST_MPI_INTEGER MPI_INT
#endif

/* Per-context rank view; process MPI ownership lives in mpi_runtime.c. */
#define CBALLS_MPI_STATE (cballs_runtime_current()->mpi[CBALLS_MPI_LYA_FOREST_MPI])
#define mpi_active (CBALLS_MPI_STATE.active)
#define mpi_rank (CBALLS_MPI_STATE.rank)
#define mpi_size (CBALLS_MPI_STATE.size ? CBALLS_MPI_STATE.size : 1)



static int mpi_error(struct cmdline_data *cmd, const char *operation,
                     int status)
{
    return cballs_mpi_shared_error(cmd, operation, status);
}

int lya_forest_mpi_prepare(struct cmdline_data *cmd,
                                   struct global_data *gd)
{
    return cballs_mpi_shared_prepare(&CBALLS_MPI_STATE, cmd, gd, !(!lya_forest_is_mpi_method(cmd->searchMethod)));
}

int lya_forest_mpi_finalize(struct cmdline_data *cmd)
{
    return cballs_mpi_shared_finalize(&CBALLS_MPI_STATE, cmd);
}

int lya_forest_mpi_active(void) { return mpi_active; }
int lya_forest_mpi_is_root(void)
{
    return !mpi_active || mpi_rank == LYA_FOREST_MPI_ROOT;
}
int lya_forest_mpi_rank(void) { return mpi_rank; }
int lya_forest_mpi_size(void) { return mpi_size; }

int lya_forest_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !lya_forest_is_mpi_method(cmd->searchMethod)
        || !mpi_active || mpi_rank == LYA_FOREST_MPI_ROOT;
}

int lya_forest_mpi_consensus(struct cmdline_data *cmd,
                                     int local_status,
                                     const char *operation)
{
    return cballs_mpi_shared_consensus(&CBALLS_MPI_STATE, cmd, !(!lya_forest_is_mpi_method(cmd->searchMethod)), local_status, operation);
}

int lya_forest_mpi_reduce_reals(struct cmdline_data *cmd,
                                        real *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == LYA_FOREST_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                LYA_FOREST_MPI_REAL, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                LYA_FOREST_MPI_REAL, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI Ly-alpha real reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int lya_forest_mpi_reduce_integers(struct cmdline_data *cmd,
                                           INTEGER *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == LYA_FOREST_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                LYA_FOREST_MPI_INTEGER, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                LYA_FOREST_MPI_INTEGER, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI Ly-alpha integer reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int lya_forest_mpi_reduce_long_doubles(struct cmdline_data *cmd,
                                        long double *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == LYA_FOREST_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                MPI_LONG_DOUBLE, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                MPI_LONG_DOUBLE, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI Ly-alpha long-double reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int lya_forest_mpi_reduce_uint64(struct cmdline_data *cmd,
                                        uint64_t *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == LYA_FOREST_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                MPI_UINT64_T, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                MPI_UINT64_T, MPI_SUM,
                                LYA_FOREST_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI Ly-alpha uint64 reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* LYAFORESTMPI */
