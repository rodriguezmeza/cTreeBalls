/*
 * MPI scheduling for the cTreeBalls FCFC-style ball tree.
 *
 * The root-built tree broadcast, breadth-first task frontier, dynamic RMA
 * scheduler, and final reduction follow FCFC's MPI design:
 * https://github.com/cheng-zhao/FCFC
 * Copyright (c) 2020--2022 Cheng Zhao, used under the MIT license included in
 * fcfc_balltree.c.
 */

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#include "globaldefs.h"
#include "fcfc_balltree.h"
#include "fcfc_balltree_mpi.h"

#ifdef BALLTREEMPI

#define FCFC_MPI_ROOT 0
#define FCFC_MPI_BODY_CHUNK 65536

#define FCFC_MPI_REAL MPI_DOUBLE

#ifdef LONGINT
#define FCFC_MPI_INTEGER MPI_LONG
#else
#define FCFC_MPI_INTEGER MPI_INT
#endif

typedef struct {
    real pos[NDIM];
    real mass;
    real kappa;
    real weight;
    real radius;
    short type;
} fcfc_mpi_body_wire;

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
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s failed%s%s", operation, length ? ": " : "",
             length ? detail : "");
    return FAILURE;
}

int fcfc_balltree_mpi_prepare(struct cmdline_data *cmd,
                              struct global_data *gd)
{
    int initialized = FALSE;
    int finalized = FALSE;
    int provided = MPI_THREAD_SINGLE;
    int status;

    if (cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "balltree-mpi") != 0)
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
                 "balltree-mpi cannot start after MPI_Finalize");
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
                     "balltree-mpi could not register MPI cleanup");
            finalize_at_exit();
            return FAILURE;
        }
    } else {
        if ((status = MPI_Query_thread(&provided)) != MPI_SUCCESS)
            return mpi_error(cmd, "MPI_Query_thread", status);
    }
    if (provided < MPI_THREAD_FUNNELED) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi requires MPI_THREAD_FUNNELED support");
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

int fcfc_balltree_mpi_finalize(struct cmdline_data *cmd)
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

int fcfc_balltree_mpi_active(void) { return mpi_active; }
int fcfc_balltree_mpi_is_root(void) {
    return !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}
int fcfc_balltree_mpi_rank(void) { return mpi_rank; }
int fcfc_balltree_mpi_size(void) { return mpi_size; }

int fcfc_balltree_mpi_output_enabled(struct cmdline_data *cmd)
{
    return cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "balltree-mpi") != 0
        || !mpi_active || mpi_rank == FCFC_MPI_ROOT;
}

int fcfc_balltree_mpi_consensus(struct cmdline_data *cmd, int local_status,
                                const char *operation)
{
    int local_success = local_status == SUCCESS;
    int all_success = FALSE;
    int status;

    if (!mpi_active || cmd->searchMethod == NULL
        || strcmp(cmd->searchMethod, "balltree-mpi") != 0)
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

static void pack_body(fcfc_mpi_body_wire *wire, bodyptr p)
{
    int k;
    DO_COORD(k) wire->pos[k] = Pos(p)[k];
    wire->mass = Mass(p);
    wire->kappa = Kappa(p);
    wire->weight = Weight(p);
    wire->radius = Radius(p);
    wire->type = Type(p);
}

static void unpack_body(const fcfc_mpi_body_wire *wire, bodyptr p)
{
    int k;
    DO_COORD(k) Pos(p)[k] = wire->pos[k];
    Mass(p) = wire->mass;
    Kappa(p) = wire->kappa;
    Weight(p) = wire->weight;
    Radius(p) = wire->radius;
    Type(p) = wire->type;
}

static int broadcast_catalog(struct cmdline_data *cmd, bodyptr btab,
                             INTEGER nbody)
{
    INTEGER root_count = nbody;
    fcfc_mpi_body_wire *buffer = NULL;
    int status;

    if ((status = MPI_Bcast(&root_count, 1, FCFC_MPI_INTEGER, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI catalog-size broadcast", status);
    int catalog_valid = root_count == nbody && nbody > 0 && btab != NULL;
    if (!catalog_valid)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi catalog sizes differ across ranks");
    if (fcfc_balltree_mpi_consensus(cmd,
            catalog_valid ? SUCCESS : FAILURE,
            "MPI catalog validation") == FAILURE)
        return FAILURE;

    buffer = calloc(FCFC_MPI_BODY_CHUNK, sizeof(*buffer));
    if (fcfc_balltree_mpi_consensus(cmd,
            buffer == NULL ? FAILURE : SUCCESS,
            "MPI catalog buffer allocation") == FAILURE) {
        free(buffer);
        return FAILURE;
    }

    for (INTEGER first = 0; first < nbody; first += FCFC_MPI_BODY_CHUNK) {
        INTEGER count = nbody - first;
        if (count > FCFC_MPI_BODY_CHUNK) count = FCFC_MPI_BODY_CHUNK;
        if (mpi_rank == FCFC_MPI_ROOT)
            for (INTEGER i = 0; i < count; i++)
                pack_body(&buffer[i], nthBody(btab, first + i));
        status = MPI_Bcast(buffer, (int)((size_t)count * sizeof(*buffer)),
                           MPI_BYTE, FCFC_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
            free(buffer);
            return mpi_error(cmd, "MPI catalog broadcast", status);
        }
        if (mpi_rank != FCFC_MPI_ROOT)
            for (INTEGER i = 0; i < count; i++)
                unpack_body(&buffer[i], nthBody(btab, first + i));
    }
    free(buffer);
    return SUCCESS;
}

int fcfc_balltree_mpi_broadcast_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btab, INTEGER *nbody,
                                         int cat1, int cat2)
{
    (void)gd;
    if (!mpi_active) return SUCCESS;
    if (broadcast_catalog(cmd, btab[cat1], nbody[cat1]) == FAILURE)
        return FAILURE;
    if (cat2 != cat1
        && broadcast_catalog(cmd, btab[cat2], nbody[cat2]) == FAILURE)
        return FAILURE;
    return SUCCESS;
}

static int broadcast_bytes(struct cmdline_data *cmd, void *data, size_t bytes,
                           const char *operation)
{
    unsigned char *cursor = data;
    while (bytes != 0) {
        int count = bytes > (size_t)INT_MAX ? INT_MAX : (int)bytes;
        int status = MPI_Bcast(cursor, count, MPI_BYTE, FCFC_MPI_ROOT,
                               MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) return mpi_error(cmd, operation, status);
        cursor += count;
        bytes -= (size_t)count;
    }
    return SUCCESS;
}

int fcfc_balltree_mpi_build(struct cmdline_data *cmd, struct global_data *gd,
                            bodyptr btab, INTEGER nbody, int nleaf,
                            fcfc_balltreeptr *result)
{
    fcfc_balltreeptr tree = NULL;
    uint64_t metadata[4] = {0, 0, 0, 0};
    INTEGER *indices = NULL;
    int build_status = SUCCESS;
    int status;

    if (!mpi_active)
        return fcfc_balltree_build(cmd, gd, btab, nbody, nleaf, result);
    *result = NULL;
    if (mpi_rank == FCFC_MPI_ROOT) {
        build_status = fcfc_balltree_build(cmd, gd, btab, nbody, nleaf,
                                           &tree);
        if (build_status == SUCCESS) {
            metadata[0] = (uint64_t)tree->npoint;
            metadata[1] = (uint64_t)tree->nnode;
            metadata[2] = (uint64_t)tree->capacity;
            metadata[3] = (uint64_t)tree->max_depth;
        }
    }
    if ((status = MPI_Bcast(&build_status, 1, MPI_INT, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS) {
        fcfc_balltree_free(tree);
        return mpi_error(cmd, "MPI tree-status broadcast", status);
    }
    if (build_status == FAILURE) {
        if (mpi_rank != FCFC_MPI_ROOT)
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "balltree-mpi tree construction failed on rank 0");
        return FAILURE;
    }
    if ((status = MPI_Bcast(metadata, 4, MPI_UINT64_T, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS) {
        fcfc_balltree_free(tree);
        return mpi_error(cmd, "MPI tree-metadata broadcast", status);
    }
    if (metadata[0] != (uint64_t)nbody || metadata[1] < 1
        || metadata[1] > metadata[2]
        || metadata[0] > (uint64_t)SIZE_MAX / sizeof(bodyptr)
        || metadata[1] > (uint64_t)SIZE_MAX / sizeof(fcfc_ballnode)) {
        fcfc_balltree_free(tree);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi received invalid tree metadata");
        return FAILURE;
    }

    if (mpi_rank != FCFC_MPI_ROOT) {
        tree = calloc(1, sizeof(*tree));
        if (tree != NULL) {
            tree->npoint = (INTEGER)metadata[0];
            tree->nnode = (INTEGER)metadata[1];
            tree->capacity = tree->nnode;
            tree->max_depth = (int)metadata[3];
            tree->bptr = malloc((size_t)tree->npoint * sizeof(*tree->bptr));
            tree->nodes = calloc((size_t)tree->nnode, sizeof(*tree->nodes));
        }
    }
    if (fcfc_balltree_mpi_consensus(cmd,
            tree != NULL && tree->bptr != NULL && tree->nodes != NULL
                ? SUCCESS : FAILURE,
            "MPI tree allocation") == FAILURE) {
        fcfc_balltree_free(tree);
        return FAILURE;
    }

    if (broadcast_bytes(cmd, tree->nodes,
                        (size_t)tree->nnode * sizeof(*tree->nodes),
                        "MPI tree-node broadcast") == FAILURE) {
        fcfc_balltree_free(tree);
        return FAILURE;
    }

    indices = malloc(FCFC_MPI_BODY_CHUNK * sizeof(*indices));
    if (fcfc_balltree_mpi_consensus(cmd,
            indices == NULL ? FAILURE : SUCCESS,
            "MPI tree-index allocation") == FAILURE) {
        free(indices);
        fcfc_balltree_free(tree);
        return FAILURE;
    }
    for (INTEGER first = 0; first < nbody; first += FCFC_MPI_BODY_CHUNK) {
        INTEGER count = nbody - first;
        int indices_valid = TRUE;
        if (count > FCFC_MPI_BODY_CHUNK) count = FCFC_MPI_BODY_CHUNK;
        if (mpi_rank == FCFC_MPI_ROOT)
            for (INTEGER i = 0; i < count; i++)
                indices[i] = tree->bptr[first + i] - btab;
        status = MPI_Bcast(indices, (int)count, FCFC_MPI_INTEGER,
                           FCFC_MPI_ROOT, MPI_COMM_WORLD);
        if (status != MPI_SUCCESS) {
            free(indices);
            fcfc_balltree_free(tree);
            return mpi_error(cmd, "MPI tree-index broadcast", status);
        }
        for (INTEGER i = 0; i < count; i++)
            if (indices[i] < 0 || indices[i] >= nbody) {
                indices_valid = FALSE;
                break;
            }
        if (!indices_valid)
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "balltree-mpi received an invalid body index");
        if (fcfc_balltree_mpi_consensus(cmd,
                indices_valid ? SUCCESS : FAILURE,
                "MPI tree-index validation") == FAILURE) {
            free(indices);
            fcfc_balltree_free(tree);
            return FAILURE;
        }
        if (mpi_rank != FCFC_MPI_ROOT)
            for (INTEGER i = 0; i < count; i++)
                tree->bptr[first + i] = nthBody(btab, indices[i]);
    }
    free(indices);
    if (mpi_rank != FCFC_MPI_ROOT)
        gd->bytes_tot += sizeof(*tree)
            + (size_t)tree->npoint * sizeof(*tree->bptr)
            + (size_t)tree->nnode * sizeof(*tree->nodes);
    *result = tree;
    return SUCCESS;
}

int fcfc_balltree_mpi_frontier(struct cmdline_data *cmd,
                               fcfc_balltreeptr tree, INTEGER minimum_count,
                               INTEGER **result, INTEGER *result_count)
{
    INTEGER *frontier = NULL;
    INTEGER count = 0;
    int local_status = SUCCESS;
    int status;

    if (mpi_rank == FCFC_MPI_ROOT)
        local_status = fcfc_balltree_frontier(cmd, tree, minimum_count,
                                              &frontier, &count);
    if ((status = MPI_Bcast(&local_status, 1, MPI_INT, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
        return mpi_error(cmd, "MPI frontier-status broadcast", status);
    if (local_status == FAILURE) {
        if (mpi_rank != FCFC_MPI_ROOT)
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "balltree-mpi frontier construction failed on rank 0");
        return FAILURE;
    }
    if ((status = MPI_Bcast(&count, 1, FCFC_MPI_INTEGER, FCFC_MPI_ROOT,
                            MPI_COMM_WORLD)) != MPI_SUCCESS) {
        free(frontier);
        return mpi_error(cmd, "MPI frontier-size broadcast", status);
    }
    if (count < 1 || count > INT_MAX
        || (size_t)count > SIZE_MAX / sizeof(*frontier)) {
        free(frontier);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi received an invalid frontier size");
        return FAILURE;
    }
    if (mpi_rank != FCFC_MPI_ROOT)
        frontier = malloc((size_t)count * sizeof(*frontier));
    if (fcfc_balltree_mpi_consensus(cmd,
            frontier == NULL ? FAILURE : SUCCESS,
            "MPI frontier allocation") == FAILURE) {
        free(frontier);
        return FAILURE;
    }
    if ((status = MPI_Bcast(frontier, (int)count, FCFC_MPI_INTEGER,
                            FCFC_MPI_ROOT, MPI_COMM_WORLD)) != MPI_SUCCESS) {
        free(frontier);
        return mpi_error(cmd, "MPI frontier broadcast", status);
    }
    *result = frontier;
    *result_count = count;
    return SUCCESS;
}

int fcfc_balltree_mpi_scheduler_init(struct cmdline_data *cmd,
                                     fcfc_balltree_mpi_scheduler *scheduler,
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

int fcfc_balltree_mpi_scheduler_claim(struct cmdline_data *cmd,
                                      fcfc_balltree_mpi_scheduler *scheduler,
                                      INTEGER *first, INTEGER *last)
{
    uint64_t start;
    int status;

    if (!scheduler->ready) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "balltree-mpi scheduler is not initialized");
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

int fcfc_balltree_mpi_scheduler_destroy(struct cmdline_data *cmd,
                                        fcfc_balltree_mpi_scheduler *scheduler)
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

int fcfc_balltree_mpi_reduce_reals(struct cmdline_data *cmd, real *values,
                                   size_t count)
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

int fcfc_balltree_mpi_reduce_integers(struct cmdline_data *cmd,
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

#endif /* BALLTREEMPI */
