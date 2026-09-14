/* Deterministic MPI runtime for the FCFC dual/triple-node estimator. */

#include <limits.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_balltree_2balls_mpi.h"

#ifdef BALLTREE2BALLSMPI

#define BALLTREE_2BALLS_MPI_ROOT 0
#define BALLTREE_2BALLS_MPI_METHOD "balltree-2balls-mpi"

#if defined(SINGLEPREC)
#define BALLTREE_2BALLS_MPI_REAL MPI_FLOAT
#else
#define BALLTREE_2BALLS_MPI_REAL MPI_DOUBLE
#endif

#ifdef LONGINT
#define BALLTREE_2BALLS_MPI_INTEGER MPI_LONG
#else
#define BALLTREE_2BALLS_MPI_INTEGER MPI_INT
#endif

static int mpi_active = FALSE;
static int mpi_rank = 0;
static int mpi_size = 1;
static int mpi_owned = FALSE;
static int mpi_finalized = FALSE;

static int method_selected(const struct cmdline_data *cmd)
{
    return cmd->searchMethod != NULL
        && strcmp(cmd->searchMethod, BALLTREE_2BALLS_MPI_METHOD) == 0;
}

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

int fcfc_balltree_2balls_mpi_prepare(struct cmdline_data *cmd,
                                     struct global_data *gd)
{
    int initialized = FALSE;
    int finalized = FALSE;
    int provided = MPI_THREAD_SINGLE;
    int status;

    if (!method_selected(cmd)) return SUCCESS;
    if (mpi_active) {
        if (mpi_rank != BALLTREE_2BALLS_MPI_ROOT) {
            cmd->verbose = 0;
            cmd->verbose_log = 0;
            gd->flagPrint = FALSE;
        }
        return SUCCESS;
    }

    if ((status = MPI_Finalized(&finalized)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Finalized", status);
    if (finalized) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s cannot start after MPI_Finalize",
                 BALLTREE_2BALLS_MPI_METHOD);
        return FAILURE;
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
                     "%s could not register MPI cleanup",
                     BALLTREE_2BALLS_MPI_METHOD);
            finalize_at_exit();
            return FAILURE;
        }
    } else if ((status = MPI_Query_thread(&provided)) != MPI_SUCCESS) {
        return mpi_error(cmd, "MPI_Query_thread", status);
    }
    if (provided < MPI_THREAD_FUNNELED) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires MPI_THREAD_FUNNELED support",
                 BALLTREE_2BALLS_MPI_METHOD);
        return FAILURE;
    }
    if ((status = MPI_Comm_set_errhandler(MPI_COMM_WORLD,
                                          MPI_ERRORS_RETURN)) != MPI_SUCCESS
        || (status = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) != MPI_SUCCESS
        || (status = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI communicator setup", status);

    mpi_active = TRUE;
    if (mpi_rank != BALLTREE_2BALLS_MPI_ROOT) {
        cmd->verbose = 0;
        cmd->verbose_log = 0;
        gd->flagPrint = FALSE;
    }
    return SUCCESS;
}

int fcfc_balltree_2balls_mpi_finalize(struct cmdline_data *cmd)
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

int fcfc_balltree_2balls_mpi_active(void) { return mpi_active; }
int fcfc_balltree_2balls_mpi_is_root(void)
{
    return !mpi_active || mpi_rank == BALLTREE_2BALLS_MPI_ROOT;
}
int fcfc_balltree_2balls_mpi_rank(void) { return mpi_rank; }
int fcfc_balltree_2balls_mpi_size(void) { return mpi_size; }

int fcfc_balltree_2balls_mpi_output_enabled(struct cmdline_data *cmd)
{
    return !method_selected(cmd) || !mpi_active
        || mpi_rank == BALLTREE_2BALLS_MPI_ROOT;
}

int fcfc_balltree_2balls_mpi_consensus(
        struct cmdline_data *cmd, int local_status, const char *operation)
{
    int local_success = local_status == SUCCESS;
    int all_success = FALSE;
    int status;

    if (!mpi_active || !method_selected(cmd)) return local_status;
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

int fcfc_balltree_2balls_mpi_task_owned(INTEGER task)
{
    if (!mpi_active || mpi_size <= 1) return TRUE;
    return (int)(task % (INTEGER)mpi_size) == mpi_rank;
}

int fcfc_balltree_2balls_mpi_reduce_reals(
        struct cmdline_data *cmd, real *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == BALLTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                BALLTREE_2BALLS_MPI_REAL, MPI_SUM,
                                BALLTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                BALLTREE_2BALLS_MPI_REAL, MPI_SUM,
                                BALLTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI two-ball real reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int fcfc_balltree_2balls_mpi_reduce_integers(
        struct cmdline_data *cmd, INTEGER *values, size_t count)
{
    size_t offset = 0;

    if (!mpi_active || mpi_size <= 1) return SUCCESS;
    while (offset < count) {
        const int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;

        if (mpi_rank == BALLTREE_2BALLS_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                BALLTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                                BALLTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                BALLTREE_2BALLS_MPI_INTEGER, MPI_SUM,
                                BALLTREE_2BALLS_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI two-ball integer reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* BALLTREE2BALLSMPI */
