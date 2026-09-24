/* One process-level MPI owner; each runtime context has its own rank views.
 * Engine reduction buffers, task ownership and summation order stay local to
 * the backends. All entry points run on the MPI main thread (FUNNELED).
 */
#include "globaldefs.h"
#include "mpi_runtime.h"
#ifdef CBALLS_MPI_ENABLED
#include <mpi.h>
#include <limits.h>

typedef struct { int owned, finalized, cleanup_registered; } cballs_mpi_process;
static cballs_mpi_process mpi_process;

static void cballs_mpi_at_exit(void)
{
    int finalized = FALSE;
    if (!mpi_process.owned || mpi_process.finalized) return;
    if (MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized)
        MPI_Finalize();
    mpi_process.finalized = TRUE;
}

int cballs_mpi_shared_error(struct cmdline_data *cmd, const char *operation,
                             int status)
{
    char detail[MPI_MAX_ERROR_STRING] = {0};
    int length = 0;
    MPI_Error_string(status, detail, &length);
    snprintf(cmd->error_message, _ERRORMSGSIZE_, "%s failed%s%s",
             operation, length ? ": " : "", length ? detail : "");
    return FAILURE;
}

int cballs_mpi_shared_prepare(cballs_mpi_engine_state *view,
                               struct cmdline_data *cmd,
                               struct global_data *gd, int selected)
{
    int initialized = FALSE, finalized = FALSE, provided = MPI_THREAD_SINGLE;
    int status, main_thread = FALSE;
    if (!selected) return SUCCESS;
    if ((status = MPI_Finalized(&finalized)) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI_Finalized", status);
    if (finalized) {
        view->active = FALSE;
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s cannot start after MPI_Finalize", cmd->searchMethod);
        return FAILURE;
    }
    if ((status = MPI_Initialized(&initialized)) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI_Initialized", status);
    if (!initialized) {
        status = MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
        if (status != MPI_SUCCESS)
            return cballs_mpi_shared_error(cmd, "MPI_Init_thread", status);
        mpi_process.owned = TRUE;
        if (!mpi_process.cleanup_registered) {
            if (atexit(cballs_mpi_at_exit) != 0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_, "could not register MPI cleanup");
                cballs_mpi_at_exit();
                return FAILURE;
            }
            mpi_process.cleanup_registered = TRUE;
        }
    } else if ((status = MPI_Query_thread(&provided)) != MPI_SUCCESS) {
        return cballs_mpi_shared_error(cmd, "MPI_Query_thread", status);
    }
    if (provided < MPI_THREAD_FUNNELED) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires MPI_THREAD_FUNNELED support", cmd->searchMethod);
        return FAILURE;
    }
    if ((status = MPI_Is_thread_main(&main_thread)) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI_Is_thread_main", status);
    if (!main_thread) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "MPI estimators must run on the MPI main thread");
        return FAILURE;
    }
    if ((status = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN)) != MPI_SUCCESS
        || (status = MPI_Comm_rank(MPI_COMM_WORLD, &view->rank)) != MPI_SUCCESS
        || (status = MPI_Comm_size(MPI_COMM_WORLD, &view->size)) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI communicator setup", status);
    view->active = TRUE;
    if (view->rank != 0) {
        cmd->verbose = cmd->verbose_log = 0;
        gd->flagPrint = FALSE;
    }
    return SUCCESS;
}

int cballs_mpi_shared_finalize(cballs_mpi_engine_state *view,
                                struct cmdline_data *cmd)
{
    int finalized = FALSE, status;
    if (!view->active) return SUCCESS;
    view->active = FALSE;
    /* Borrowed MPI (e.g. mpi4py) is never finalized by this library. */
    if (!mpi_process.owned || mpi_process.finalized) return SUCCESS;
    if ((status = MPI_Finalized(&finalized)) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI_Finalized", status);
    if (!finalized && (status = MPI_Finalize()) != MPI_SUCCESS)
        return cballs_mpi_shared_error(cmd, "MPI_Finalize", status);
    mpi_process.finalized = TRUE;
    return SUCCESS;
}

int cballs_mpi_shared_consensus(cballs_mpi_engine_state *view,
                                 struct cmdline_data *cmd, int selected,
                                 int local_status, const char *operation)
{
    int first_failure, candidate, status;
    if (!selected || !view->active) return local_status;
    candidate = local_status == SUCCESS ? INT_MAX : view->rank;
    status = MPI_Allreduce(&candidate, &first_failure, 1, MPI_INT, MPI_MIN,
                           MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) return cballs_mpi_shared_error(cmd, operation, status);
    if (first_failure == INT_MAX) return SUCCESS;
    /* All ranks retain the same actionable message from the first failed rank. */
    if (view->rank == first_failure && cmd->error_message[0] == '\0')
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s failed on MPI rank %d", operation, first_failure);
    status = MPI_Bcast(cmd->error_message, _ERRORMSGSIZE_, MPI_CHAR,
                       first_failure, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) return cballs_mpi_shared_error(cmd, operation, status);
    return FAILURE;
}
#endif
