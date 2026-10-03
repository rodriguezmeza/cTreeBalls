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
    if (!selected || !view->active || cballs_runtime_current()->mpi_local_depth)
        return local_status;
    cballs_runtime_state *runtime = cballs_runtime_current();
    size_t trace_len = strlen(runtime->mpi_boundary_trace);
    if (trace_len + strlen(operation) + 2 < sizeof(runtime->mpi_boundary_trace)) {
        strcat(runtime->mpi_boundary_trace, operation);
        strcat(runtime->mpi_boundary_trace, "\n");
    }
    if (runtime->mpi_test_failure[0]
        && !strcmp(runtime->mpi_test_failure, operation)) {
        runtime->mpi_test_failure[0] = '\0';
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "injected application failure at %s on MPI rank %d", operation, view->rank);
        local_status = FAILURE;
    }
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

/* Startup consists of guarded local phases. Nested reader/helper consensus is
 * deferred to the phase boundary, so a local early return cannot skip a
 * collective that peers enter. Never use this scope around a kernel reduction. */
void cballs_mpi_local_begin(void) {
#ifdef CBALLS_MPI_ENABLED
    cballs_runtime_current()->mpi_local_depth++;
#endif
}
void cballs_mpi_local_end(void) {
#ifdef CBALLS_MPI_ENABLED
    cballs_runtime_current()->mpi_local_depth--;
#endif
}
void cballs_mpi_test_boundary(const char *operation) {
#ifdef CBALLS_MPI_ENABLED
    cballs_runtime_state *r = cballs_runtime_current();
    snprintf(r->mpi_test_failure, sizeof(r->mpi_test_failure), "%s", operation);
    r->mpi_boundary_trace[0] = '\0';
#endif
}
const char *cballs_mpi_boundary_trace(void) {
#ifdef CBALLS_MPI_ENABLED
    return cballs_runtime_current()->mpi_boundary_trace;
#else
    return "";
#endif
}
int cballs_mpi_bootstrap(struct cmdline_data *cmd, struct global_data *gd,
                          const char *method) {
#ifdef CBALLS_MPI_ENABLED
    const int launched = getenv("OMPI_COMM_WORLD_SIZE") || getenv("PMI_SIZE") || getenv("PMIX_RANK");
    const int known_serial = method && !strstr(method,"-mpi") && cballs_search_method_id(method)>=0;
    const int selected = !known_serial && ((method && strstr(method, "-mpi")) || launched);
    if (!selected) {
        /* MPI drivers may execute an OpenMP engine on rank zero alone. Context
         * views end here; process ownership remains with its original owner. */
        for (int i=0;i<CBALLS_MPI_ENGINE_COUNT;i++)
            cballs_runtime_current()->mpi[i].active=FALSE;
    }
    char *saved = cmd->searchMethod;
    cmd->searchMethod = (char *)(method ? method : "MPI startup");
    int status = cballs_mpi_shared_prepare(
        &cballs_runtime_current()->mpi[CBALLS_MPI_STARTUP], cmd, gd, selected);
    cmd->searchMethod = saved;
    return status;
#else
    return SUCCESS;
#endif
}
int cballs_mpi_context_consensus(struct cmdline_data *cmd, int status,
                                  const char *operation) {
#ifdef CBALLS_MPI_ENABLED
    cballs_runtime_state *r = cballs_runtime_current();
    if (r->mpi[CBALLS_MPI_STARTUP].active)
        return cballs_mpi_shared_consensus(&r->mpi[CBALLS_MPI_STARTUP], cmd, TRUE, status, operation);
    for (int i = 0; i < CBALLS_MPI_ENGINE_COUNT; i++)
        if (r->mpi[i].active)
            return cballs_mpi_shared_consensus(&r->mpi[i], cmd, TRUE, status, operation);
#endif
    return status;
}
int cballs_mpi_agree_text(struct cmdline_data *cmd, const char *text,
                          const char *operation) {
#ifdef CBALLS_MPI_ENABLED
    cballs_runtime_state *r = cballs_runtime_current();
    int active = 0;
    for (int i = 0; i < CBALLS_MPI_ENGINE_COUNT; i++) active |= r->mpi[i].active;
    if (!active) return SUCCESS;
    unsigned long long hash = 1469598103934665603ULL, lo, hi;
    for (const unsigned char *p = (const unsigned char *)text; *p; p++)
        hash = (hash ^ *p)*1099511628211ULL;
    int status = MPI_Allreduce(&hash, &lo, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    if (status == MPI_SUCCESS)
        status = MPI_Allreduce(&hash, &hi, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) return cballs_mpi_shared_error(cmd, operation, status);
    if (lo != hi) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "%s: MPI ranks supplied inconsistent run controls", operation);
        return FAILURE;
    }
#endif
    return SUCCESS;
}

int cballs_mpi_agree_catalogs(struct cmdline_data *cmd, struct global_data *gd)
{
#ifdef CBALLS_MPI_ENABLED
    int active = FALSE;
    for (int i=0;i<CBALLS_MPI_ENGINE_COUNT;i++) active |= cballs_runtime_current()->mpi[i].active;
    if (!active) return SUCCESS;
    unsigned long long hash = 1469598103934665603ULL;
#define HASH_VALUE(value) do { const unsigned char *bytes=(const unsigned char *)&(value); \
    for(size_t b=0;b<sizeof(value);b++) hash=(hash^bytes[b])*1099511628211ULL; } while(0)
    HASH_VALUE(gd->ninfiles);
    for (int cat=0;cat<gd->ninfiles;cat++) {
        /* Companion masks modify catalog zero and do not own catalog one. */
        if (!bodytable[cat]) continue;
        HASH_VALUE(gd->nbodyTable[cat]);
        for (INTEGER row=0;row<gd->nbodyTable[cat];row++) {
            bodyptr p=bodytable[cat]+row;
            for(int axis=0;axis<NDIM;axis++) HASH_VALUE(Pos(p)[axis]);
            HASH_VALUE(Kappa(p)); HASH_VALUE(Weight(p)); HASH_VALUE(Mask(p));
#ifdef THREEPCFSHEAR
            if (strstr(cmd->searchMethod,"shear")) { HASH_VALUE(Gamma1(p)); HASH_VALUE(Gamma2(p)); }
#endif
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
            if (strstr(cmd->searchMethod,"3pcf-3d")) HASH_VALUE(Octree3pcf3dLosId(p));
#endif
#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
            if (!strncmp(cmd->searchMethod,"lya-",4)) HASH_VALUE(LyaForestId(p));
#endif
        }
    }
#undef HASH_VALUE
    char text[32]; snprintf(text,sizeof(text),"%016llx",hash);
    return cballs_mpi_agree_text(cmd,text,"MPI replicated catalog agreement");
#else
    return SUCCESS;
#endif
}
