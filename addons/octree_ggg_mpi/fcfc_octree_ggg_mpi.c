/* FCFC-style dynamic MPI scheduling for the octree-GGG estimator. */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_octree_ggg_mpi.h"

#ifdef OCTREEGGGMPI

#define FCFC_MPI_ROOT 0

#define FCFC_MPI_REAL MPI_DOUBLE

#ifdef LONGINT
#define FCFC_MPI_INTEGER MPI_LONG
#else
#define FCFC_MPI_INTEGER MPI_INT
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

int fcfc_octree_ggg_mpi_prepare(struct cmdline_data *cmd,
                                struct global_data *gd)
{
    int initialized = FALSE;
    int finalized = FALSE;
    int provided = MPI_THREAD_SINGLE;
    int status;

    if (cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "octree-ggg-mpi") != 0)
        return SUCCESS;
    if (mpi_active) {
        if (mpi_rank != FCFC_MPI_ROOT) {
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
                 "octree-ggg-mpi cannot start after MPI_Finalize");
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
                     "octree-ggg-mpi could not register MPI cleanup");
            finalize_at_exit();
            return FAILURE;
        }
    } else if ((status = MPI_Query_thread(&provided)) != MPI_SUCCESS) {
        return mpi_error(cmd, "MPI_Query_thread", status);
    }
    if (provided < MPI_THREAD_FUNNELED) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-mpi requires MPI_THREAD_FUNNELED support");
        return FAILURE;
    }
    if ((status = MPI_Comm_set_errhandler(MPI_COMM_WORLD,
                                          MPI_ERRORS_RETURN)) != MPI_SUCCESS
        || (status = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank)) != MPI_SUCCESS
        || (status = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI communicator setup", status);

    mpi_active = TRUE;
    if (mpi_rank != FCFC_MPI_ROOT) {
        cmd->verbose = 0;
        cmd->verbose_log = 0;
        gd->flagPrint = FALSE;
    }
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_finalize(struct cmdline_data *cmd)
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

int fcfc_octree_ggg_mpi_active(void) { return mpi_active; }
int fcfc_octree_ggg_mpi_is_root(void) {
    return !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}
int fcfc_octree_ggg_mpi_rank(void) { return mpi_rank; }
int fcfc_octree_ggg_mpi_size(void) { return mpi_size; }

int fcfc_octree_ggg_mpi_output_enabled(struct cmdline_data *cmd)
{
    return cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "octree-ggg-mpi") != 0
        || !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}

int fcfc_octree_ggg_mpi_consensus(struct cmdline_data *cmd,
                                  int local_status,
                                  const char *operation)
{
    int local_success = local_status == SUCCESS;
    int all_success = FALSE;
    int status;

    if (!mpi_active || cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "octree-ggg-mpi") != 0)
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

int fcfc_octree_ggg_mpi_scheduler_init(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler,
    INTEGER task_count, INTEGER step)
{
    MPI_Aint bytes = 0;
    int root_status = SUCCESS;
    int status;

    memset(scheduler, 0, sizeof(*scheduler));
    scheduler->window = MPI_WIN_NULL;
    scheduler->task_count = (uint64_t)task_count;
    scheduler->step = step > 0 ? (uint64_t)step : 1;
    if (mpi_size == 1) {
        scheduler->ready = TRUE;
        return SUCCESS;
    }
    if (mpi_rank == FCFC_MPI_ROOT) {
        bytes = sizeof(*scheduler->counter);
        if ((status = MPI_Alloc_mem(bytes, MPI_INFO_NULL,
                                    &scheduler->counter)) != MPI_SUCCESS)
            root_status = mpi_error(cmd, "MPI_Alloc_mem", status);
        else
            *scheduler->counter = 0;
    }
    if ((status = MPI_Bcast(&root_status, 1, MPI_INT, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI scheduler-status broadcast", status);
    if (root_status == FAILURE) return FAILURE;
    status = MPI_Win_create(scheduler->counter, bytes, sizeof(uint64_t),
                            MPI_INFO_NULL, MPI_COMM_WORLD, &scheduler->window);
    if (status != MPI_SUCCESS) {
        if (mpi_rank == FCFC_MPI_ROOT) MPI_Free_mem(scheduler->counter);
        scheduler->counter = NULL;
        return mpi_error(cmd, "MPI_Win_create", status);
    }
    scheduler->ready = TRUE;
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_scheduler_claim(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler,
    INTEGER *first, INTEGER *last)
{
    uint64_t start;
    int status;

    if (!scheduler->ready) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-mpi scheduler is not initialized");
        return FAILURE;
    }
    if (mpi_size == 1) {
        start = scheduler->local_next;
        scheduler->local_next += scheduler->step;
    } else {
        int locked = FALSE;
        status = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, FCFC_MPI_ROOT, 0,
                              scheduler->window);
        if (status == MPI_SUCCESS) {
            locked = TRUE;
            status = MPI_Fetch_and_op(&scheduler->step, &start, MPI_UINT64_T,
                                      FCFC_MPI_ROOT, 0, MPI_SUM,
                                      scheduler->window);
        }
        if (locked) {
            const int unlock_status = MPI_Win_unlock(FCFC_MPI_ROOT,
                                                      scheduler->window);
            if (status == MPI_SUCCESS) status = unlock_status;
        }
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI RMA task claim", status);
    }
    if (start >= scheduler->task_count) {
        *first = *last = (INTEGER)scheduler->task_count;
        return SUCCESS;
    }
    uint64_t end = start + scheduler->step;
    if (end > scheduler->task_count) end = scheduler->task_count;
    *first = (INTEGER)start;
    *last = (INTEGER)end;
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_scheduler_destroy(
    struct cmdline_data *cmd, fcfc_octree_ggg_mpi_scheduler *scheduler)
{
    int status = MPI_SUCCESS;
    int free_status = MPI_SUCCESS;

    if (!scheduler->ready) return SUCCESS;
    if (mpi_size > 1) {
        status = MPI_Win_free(&scheduler->window);
        if (mpi_rank == FCFC_MPI_ROOT && scheduler->counter != NULL)
            free_status = MPI_Free_mem(scheduler->counter);
    }
    memset(scheduler, 0, sizeof(*scheduler));
    if (status != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Win_free", status);
    if (free_status != MPI_SUCCESS)
        return mpi_error(cmd, "MPI_Free_mem", free_status);
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_reduce_reals(struct cmdline_data *cmd,
                                     real *values, size_t count)
{
    size_t offset = 0;
    while (offset < count) {
        int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;
        if (mpi_rank == FCFC_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                FCFC_MPI_REAL, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk, FCFC_MPI_REAL,
                                MPI_SUM, FCFC_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI histogram reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

int fcfc_octree_ggg_mpi_reduce_integers(struct cmdline_data *cmd,
                                        INTEGER *values, size_t count)
{
    size_t offset = 0;
    while (offset < count) {
        int chunk = count - offset > (size_t)INT_MAX
            ? INT_MAX : (int)(count - offset);
        int status;
        if (mpi_rank == FCFC_MPI_ROOT)
            status = MPI_Reduce(MPI_IN_PLACE, values + offset, chunk,
                                FCFC_MPI_INTEGER, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        else
            status = MPI_Reduce(values + offset, NULL, chunk,
                                FCFC_MPI_INTEGER, MPI_SUM, FCFC_MPI_ROOT,
                                MPI_COMM_WORLD);
        if (status != MPI_SUCCESS)
            return mpi_error(cmd, "MPI counter reduction", status);
        offset += (size_t)chunk;
    }
    return SUCCESS;
}

#endif /* OCTREEGGGMPI */
