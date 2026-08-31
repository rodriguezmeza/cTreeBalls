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

static int mpi_active = FALSE;
static int mpi_rank = 0;
static int mpi_size = 1;
static int mpi_owned = FALSE;
static int mpi_finalized = FALSE;

static void finalize_at_exit(void)
{
    int finalized = FALSE;

    if (!mpi_owned || mpi_finalized) return;
    if (MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized)
        MPI_Finalize();
    mpi_finalized = TRUE;
    mpi_active = FALSE;
}

static int mpi_error(struct cmdline_data *cmd, const char *operation,
                     int status)
{
    char detail[MPI_MAX_ERROR_STRING];
    int length = 0;

    detail[0] = '\0';
    MPI_Error_string(status, detail, &length);
    snprintf(cmd->error_message, _ERRORMSGSIZE_, "%s failed%s%s",
             operation, length ? ": " : "", length ? detail : "");
    return FAILURE;
}

int lya_forest_mpi_prepare(struct cmdline_data *cmd,
                                   struct global_data *gd)
{
    int initialized = FALSE;
    int finalized = FALSE;
    int provided = MPI_THREAD_SINGLE;
    int status;

    if (!lya_forest_is_mpi_method(cmd->searchMethod))
        return SUCCESS;
    if ((status = MPI_Finalized(&finalized)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Finalized", status);
    if (finalized) {
        mpi_active = FALSE;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Ly-alpha MPI cannot start after MPI_Finalize");
        return FAILURE;
    }
    if (mpi_active) {
        if (mpi_rank != LYA_FOREST_MPI_ROOT) {
            cmd->verbose = 0;
            cmd->verbose_log = 0;
            gd->flagPrint = FALSE;
        }
        return SUCCESS;
    }

    if ((status = MPI_Initialized(&initialized)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Initialized", status);
    if (!initialized) {
        status = MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI_Init_thread", status);
        mpi_owned = TRUE;
        if (atexit(finalize_at_exit) != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "Ly-alpha MPI could not register MPI cleanup");
            finalize_at_exit();
            return FAILURE;
        }
    } else if ((status = MPI_Query_thread(&provided)) != MPI_SUCCESS) {
        return mpi_error(cmd, "MPI_Query_thread", status);
    }
    if (provided < MPI_THREAD_FUNNELED) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Ly-alpha MPI requires MPI_THREAD_FUNNELED support");
        return FAILURE;
    }
    if ((status = MPI_Comm_set_errhandler(MPI_COMM_WORLD,
                                          MPI_ERRORS_RETURN)) != MPI_SUCCESS
        || (status = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) != MPI_SUCCESS
        || (status = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI communicator setup", status);

    mpi_active = TRUE;
    if (mpi_rank != LYA_FOREST_MPI_ROOT) {
        cmd->verbose = 0;
        cmd->verbose_log = 0;
        gd->flagPrint = FALSE;
    }
    return SUCCESS;
}

int lya_forest_mpi_finalize(struct cmdline_data *cmd)
{
    int finalized = FALSE;
    int status;

    if (!mpi_active || !mpi_owned || mpi_finalized) return SUCCESS;
    if ((status = MPI_Finalized(&finalized)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Finalized", status);
    if (!finalized && (status = MPI_Finalize()) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Finalize", status);
    mpi_finalized = TRUE;
    mpi_active = FALSE;
    return SUCCESS;
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
    int local_success = local_status == SUCCESS;
    int all_success = FALSE;
    int status;

    if (!mpi_active || !lya_forest_is_mpi_method(cmd->searchMethod))
        return local_status;
    status = MPI_Allreduce(&local_success, &all_success, 1, MPI_INT, MPI_MIN,
                           MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) return mpi_error(cmd, operation, status);
    if (!all_success) {
        if (local_success || cmd->error_message[0] == '\0')
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s failed on at least one MPI rank", operation);
        return FAILURE;
    }
    return SUCCESS;
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
