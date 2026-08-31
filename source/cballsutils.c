/* ==============================================================================
 MODULE: cballsutils.c			    [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "tree_contracts.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

#include <errno.h>
#include <sys/wait.h>

local cballs_runtime_state cballs_default_runtime;
local cballs_runtime_state *cballs_active_runtime;

#ifdef USEGSL
local _Thread_local gsl_error_handler_t *cballs_previous_gsl_handler;

local void cballs_gsl_allocation_handler(const char *reason,
                                          const char *file,
                                          int line,
                                          int gsl_errno)
{
    if (gsl_errno == GSL_ENOMEM)
        cballs_allocation_failure(0, reason);

    if (cballs_previous_gsl_handler != NULL) {
        cballs_previous_gsl_handler(reason, file, line, gsl_errno);
        return;
    }

    gsl_stream_printf("ERROR", file, line, reason);
    fflush(stdout);
    fprintf(stderr, "Default GSL error handler invoked.\n");
    fflush(stderr);
    abort();
}
#endif

local void cballs_runtime_store(cballs_runtime_state *state)
{
#ifdef USEGSL
    state->r_gsl = r_gsl;
    state->histZetaMatrix = histZetaMatrix;
#else
    state->idum = idum;
#endif
    memcpy(state->errmsg, errmsg, sizeof(errmsg));
    memcpy(state->bodytable, bodytable, sizeof(bodytable));
    memcpy(state->nodetablescanlev, nodetablescanlev,
           sizeof(nodetablescanlev));
    memcpy(state->nodetablescanlev_root, nodetablescanlev_root,
           sizeof(nodetablescanlev_root));
    memcpy(state->roottable, roottable, sizeof(roottable));
#ifdef CBALLS_NEEDS_BALLS4_SCAN
    memcpy(state->nodetablescanlevB4, nodetablescanlevB4,
           sizeof(nodetablescanlevB4));
#endif
#ifndef MACONLY
    memcpy(state->celltable, celltable, sizeof(celltable));
#endif
    memcpy(state->tree_is_threaded, tree_is_threaded,
           sizeof(tree_is_threaded));
    state->histXi2pcf_omp = histXi2pcf_omp;
    state->rootnode = rootnode;
    state->inout_xval = inout_xval;
    state->inout_yval = inout_yval;
    state->inout_zval = inout_zval;
    state->inout_uval = inout_uval;
    state->inout_vval = inout_vval;
    state->inout_wval = inout_wval;
#ifdef CBALLS_RUNTIME_BALLS_GLOBALS
    state->bodytabbf = bodytabbf;
    state->bodytabsm = bodytabsm;
    state->bodytabSel = bodytabSel;
    state->nodetab = nodetab;
    state->nodetabscanlev = nodetabscanlev;
    state->nodetabscanlev_root = nodetabscanlev_root;
    state->nodetable = nodetable;
    state->nodetable_root = nodetable_root;
#endif
}

local void cballs_runtime_load(const cballs_runtime_state *state)
{
#ifdef USEGSL
    r_gsl = state->r_gsl;
    histZetaMatrix = state->histZetaMatrix;
#else
    idum = state->idum;
#endif
    memcpy(errmsg, state->errmsg, sizeof(errmsg));
    memcpy(bodytable, state->bodytable, sizeof(bodytable));
    memcpy(nodetablescanlev, state->nodetablescanlev,
           sizeof(nodetablescanlev));
    memcpy(nodetablescanlev_root, state->nodetablescanlev_root,
           sizeof(nodetablescanlev_root));
    memcpy(roottable, state->roottable, sizeof(roottable));
#ifdef CBALLS_NEEDS_BALLS4_SCAN
    memcpy(nodetablescanlevB4, state->nodetablescanlevB4,
           sizeof(nodetablescanlevB4));
#endif
#ifndef MACONLY
    memcpy(celltable, state->celltable, sizeof(celltable));
#endif
    memcpy(tree_is_threaded, state->tree_is_threaded,
           sizeof(tree_is_threaded));
    histXi2pcf_omp = state->histXi2pcf_omp;
    rootnode = state->rootnode;
    inout_xval = state->inout_xval;
    inout_yval = state->inout_yval;
    inout_zval = state->inout_zval;
    inout_uval = state->inout_uval;
    inout_vval = state->inout_vval;
    inout_wval = state->inout_wval;
#ifdef CBALLS_RUNTIME_BALLS_GLOBALS
    bodytabbf = state->bodytabbf;
    bodytabsm = state->bodytabsm;
    bodytabSel = state->bodytabSel;
    nodetab = state->nodetab;
    nodetabscanlev = state->nodetabscanlev;
    nodetabscanlev_root = state->nodetabscanlev_root;
    nodetable = state->nodetable;
    nodetable_root = state->nodetable_root;
#endif
}

global cballs_runtime_state *cballs_runtime_create(void)
{
    return (cballs_runtime_state *) calloc(1, sizeof(cballs_runtime_state));
}

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    char *filename;
} cballs_guarded_stage_context;

local int cballs_call_start_run_common(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return StartRun_Common(context->cmd, context->gd);
}

local int cballs_call_print_parameter_file(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return PrintParameterFile(context->cmd, context->gd, context->filename);
}

local int cballs_call_set_number_threads(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return SetNumberThreads(context->cmd);
}

local int cballs_call_main_loop(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return MainLoop(context->cmd, context->gd);
}

local int cballs_call_end_run(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return EndRun(context->cmd, context->gd);
}

local int cballs_call_end_run_free_memory(void *argument)
{
    cballs_guarded_stage_context *context = argument;
    return EndRun_FreeMemory(context->cmd, context->gd);
}

local int cballs_guard_stage(struct cmdline_data *cmd,
                             struct global_data *gd,
                             char *filename,
                             cballs_allocation_callback callback)
{
    cballs_guarded_stage_context context;
    int status;
#ifdef USEGSL
    gsl_error_handler_t *previous_gsl_handler;
#endif

    if (cmd == NULL)
        return FAILURE;

    context.cmd = cmd;
    context.gd = gd;
    context.filename = filename;
#ifdef USEGSL
    previous_gsl_handler =
        gsl_set_error_handler(cballs_gsl_allocation_handler);
    cballs_previous_gsl_handler = previous_gsl_handler;
#endif
    status = cballs_allocation_guard(callback, &context,
                                     cmd->error_message,
                                     sizeof(cmd->error_message));
#ifdef USEGSL
    gsl_set_error_handler(previous_gsl_handler);
    cballs_previous_gsl_handler = NULL;
#endif
    return status;
}

global int cballs_start_run_common_guarded(struct cmdline_data *cmd,
                                           struct global_data *gd)
{
    return cballs_guard_stage(cmd, gd, NULL,
                              cballs_call_start_run_common);
}

global int cballs_print_parameter_file_guarded(struct cmdline_data *cmd,
                                               struct global_data *gd,
                                               char *filename)
{
    return cballs_guard_stage(cmd, gd, filename,
                              cballs_call_print_parameter_file);
}

global int cballs_set_number_threads_guarded(struct cmdline_data *cmd)
{
    return cballs_guard_stage(cmd, NULL, NULL,
                              cballs_call_set_number_threads);
}

global int cballs_main_loop_guarded(struct cmdline_data *cmd,
                                    struct global_data *gd)
{
    return cballs_guard_stage(cmd, gd, NULL, cballs_call_main_loop);
}

global int cballs_end_run_guarded(struct cmdline_data *cmd,
                                  struct global_data *gd)
{
    return cballs_guard_stage(cmd, gd, NULL, cballs_call_end_run);
}

global int cballs_end_run_free_memory_guarded(struct cmdline_data *cmd,
                                              struct global_data *gd)
{
    return cballs_guard_stage(cmd, gd, NULL,
                              cballs_call_end_run_free_memory);
}

global int cballs_runtime_activate(cballs_runtime_state *state)
{
    if (state == NULL)
        return FAILURE;
    if (state == cballs_active_runtime)
        return SUCCESS;

    if (cballs_active_runtime != NULL)
        cballs_runtime_store(cballs_active_runtime);
    else
        cballs_runtime_store(&cballs_default_runtime);

    cballs_runtime_load(state);
    cballs_active_runtime = state;
    return SUCCESS;
}

global const void *cballs_runtime_bodytable_at(
    const cballs_runtime_state *state, int ifile)
{
    if (state == NULL || ifile < 0 || ifile >= MAXITEMS)
        return NULL;
    if (state == cballs_active_runtime)
        return bodytable[ifile];
    return state->bodytable[ifile];
}

global int cballs_compiled_ndim(void)
{
    return NDIM;
}

global int cballs_max_memory_catalogs(void)
{
    return MAXITEMS;
}

global int cballs_load_memory_catalog(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      int ifile,
                                      const double *positions,
                                      size_t nbody,
                                      int ndim,
                                      const double *kappa,
                                      const double *weights,
                                      const unsigned char *mask,
                                      const double *gamma1,
                                      const double *gamma2)
{
    bodyptr catalog;
    size_t i;
    int k;

    if (cmd == NULL || gd == NULL) return FAILURE;
    if (ifile < 0 || ifile >= MAXITEMS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog index %d is outside [0, %d)",
                 ifile, MAXITEMS);
        return FAILURE;
    }
    if (ifile == 0) {
        for (k = 0; k < MAXITEMS; k++) {
            if (bodytable[k] != NULL) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "cannot replace a live C catalog; clean the run first");
                return FAILURE;
            }
        }
        gd->ninfiles = 0;
    }
    if (ifile != gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalogs must be loaded contiguously from index 0");
        return FAILURE;
    }
    if (positions == NULL || ndim != NDIM) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory positions must have compiled dimension %d", NDIM);
        return FAILURE;
    }
    if (nbody < 3) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d needs at least 3 bodies", ifile);
        return FAILURE;
    }
#ifdef LONGINT
    if (nbody > (size_t)LONG_MAX) {
#else
    if (nbody > (size_t)INT_MAX) {
#endif
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d is too large for INTEGER", ifile);
        return FAILURE;
    }
    if (nbody > SIZE_MAX / sizeof(body)
        || nbody > (size_t)(LONG_MAX / (long)sizeof(body))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d byte size is not representable", ifile);
        return FAILURE;
    }
    if ((gamma1 == NULL) != (gamma2 == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "gamma1 and gamma2 must either both be supplied or both omitted");
        return FAILURE;
    }
#ifndef THREEPCFSHEAR
    if (gamma1 != NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "gamma arrays require a THREEPCFSHEAR build");
        return FAILURE;
    }
#endif

    for (i = 0; i < nbody; i++) {
        for (k = 0; k < NDIM; k++) {
            if (!isfinite(positions[i * (size_t)NDIM + (size_t)k])) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "in-memory catalog %d contains a non-finite position",
                         ifile);
                return FAILURE;
            }
        }
        if (kappa != NULL && !isfinite(kappa[i])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite kappa", ifile);
            return FAILURE;
        }
        if (weights != NULL && !isfinite(weights[i])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite weights", ifile);
            return FAILURE;
        }
        if (mask != NULL && mask[i] > 1) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d mask values must be 0 or 1", ifile);
            return FAILURE;
        }
#ifdef THREEPCFSHEAR
        if (gamma1 != NULL
            && (!isfinite(gamma1[i]) || !isfinite(gamma2[i]))) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite shear", ifile);
            return FAILURE;
        }
#endif
    }

    catalog = (bodyptr)calloc(nbody, sizeof(body));
    if (catalog == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "not enough memory for in-memory catalog %d", ifile);
        return FAILURE;
    }

    for (i = 0; i < nbody; i++) {
        bodyptr p = catalog + i;
        Type(p) = BODY;
        Update(p) = TRUE;
        Update2(p) = TRUE;
        Mask(p) = mask == NULL || mask[i] ? MASK_NODE_VALID
                                          : MASK_NODE_MASKED;
        Mass(p) = 1.0;
        Kappa(p) = kappa == NULL ? 1.0 : (REAL)kappa[i];
        Weight(p) = weights == NULL ? 1.0 : (REAL)weights[i];
        for (k = 0; k < NDIM; k++)
            Pos(p)[k] = (REAL)positions[i * (size_t)NDIM + (size_t)k];
#ifdef THREEPCFSHEAR
        Gamma1(p) = gamma1 == NULL ? 1.0 : (REAL)gamma1[i];
        Gamma2(p) = gamma2 == NULL ? 1.0 : (REAL)gamma2[i];
#endif
        Id(p) = (INTEGER)i + 1;
    }

    bodytable[ifile] = catalog;
    gd->nbodyTable[ifile] = (INTEGER)nbody;
    gd->ninfiles = ifile + 1;
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += (INTEGER)(nbody * sizeof(body));
    cmd->nbody = (INTEGER)nbody;
    for (k = 0; k < NDIM; k++)
        gd->Box[k] = cmd->lengthBox;
    gd->input_comment = "Python in-memory catalog";
    return SUCCESS;
}

global int cballs_set_memory_forest_ids(struct cmdline_data *cmd,
                                        struct global_data *gd, int ifile,
                                        const int64_t *forest_ids, size_t nbody)
{
#if (defined(LYAFORESTOMP) || defined(LYAFORESTMPI)) && NDIM == 3
    size_t i;
    int k;
    REAL minimum[NDIM], maximum[NDIM];
    if (cmd == NULL || gd == NULL) return FAILURE;
    if (ifile < 0 || ifile >= gd->ninfiles || forest_ids == NULL
        || bodytable[ifile] == NULL || nbody != (size_t)gd->nbodyTable[ifile]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "invalid in-memory forest catalog");
        return FAILURE;
    }
    for (k = 0; k < NDIM; k++)
        minimum[k] = maximum[k] = Pos(bodytable[ifile])[k];
    for (i = 0; i < nbody; i++) {
        bodyptr p = bodytable[ifile] + i;
        REAL distance = hypot(hypot(Pos(p)[0], Pos(p)[1]), Pos(p)[2]);
#ifdef LONGINT
        const int64_t id_min = LONG_MIN, id_max = LONG_MAX;
#else
        const int64_t id_min = INT_MIN, id_max = INT_MAX;
#endif
        if (forest_ids[i] < id_min || forest_ids[i] > id_max
            || !isfinite(distance) || distance <= 0.0
            || !isfinite(Weight(p)) || Weight(p) < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "forest pixels require representable INTEGER IDs, "
                     "positive observer distance and non-negative weights");
            return FAILURE;
        }
        LyaForestId(p) = (INTEGER)forest_ids[i];
        LyaDistance(p) = distance;
        for (k = 0; k < NDIM; k++) {
            LyaLOS(p)[k] = Pos(p)[k] / distance;
            minimum[k] = MIN(minimum[k], Pos(p)[k]);
            maximum[k] = MAX(maximum[k], Pos(p)[k]);
        }
    }
    for (k = 0; k < NDIM; k++) gd->Box[k] = maximum[k] - minimum[k];
    gd->input_comment = "Python in-memory Lyman-alpha forest catalog";
    return SUCCESS;
#else
    (void)gd; (void)ifile; (void)forest_ids; (void)nbody;
    if (cmd != NULL)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "forest catalogs require LYAFORESTOMP or LYAFORESTMPI and NDIM=3");
    return FAILURE;
#endif
}

global void cballs_runtime_destroy(cballs_runtime_state *state)
{
    if (state == NULL)
        return;

    if (state == cballs_active_runtime) {
        cballs_runtime_store(state);
        cballs_runtime_load(&cballs_default_runtime);
        cballs_active_runtime = NULL;
    }
    free(state);
}

global real cballs_normalize_or_zero(real numerator, real denominator)
{
    return denominator == 0.0 ? 0.0 : numerator / denominator;
}

global int cballs_system_checked(struct cmdline_data *cmd,
                                 string routineName,
                                 string command)
{
    int status = system(command);

    if (status == -1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: system command failed: %s: %s",
                 routineName, command, strerror(errno));
        return FAILURE;
    }

    if (WIFEXITED(status)) {
        int exit_status = WEXITSTATUS(status);
        if (exit_status == 0)
            return SUCCESS;

        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: system command exited with status %d: %s",
                 routineName, exit_status, command);
        return FAILURE;
    }

    if (WIFSIGNALED(status)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: system command killed by signal %d: %s",
                 routineName, WTERMSIG(status), command);
        return FAILURE;
    }

    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "%s: system command returned unexpected status %d: %s",
             routineName, status, command);
    return FAILURE;
}

global int cballs_stream_close_checked(struct cmdline_data *cmd,
                                       string routineName,
                                       stream *stream_ptr,
                                       string filename)
{
    int failed = 0;
    int saved_errno = 0;

    if (stream_ptr == NULL || *stream_ptr == NULL)
        return SUCCESS;

    if (ferror(*stream_ptr)) {
        failed = 1;
        saved_errno = errno;
    }

    if (fclose(*stream_ptr) != 0) {
        failed = 1;
        if (saved_errno == 0)
            saved_errno = errno;
    }

    *stream_ptr = NULL;

    if (failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: error closing '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}


#define CLOSE_STREAM(s)        \
    do {                       \
        if ((s) != NULL) {     \
            fclose(s);         \
            (s) = NULL;        \
        }                      \
    } while (0)


typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    gdhistptr_sincos_omp hist;
} sincos_hist_init_context;

local int search_init_sincos_omp_unguarded(void *argument)
{
    sincos_hist_init_context *context = argument;
    struct cmdline_data *cmd = context->cmd;
    gdhistptr_sincos_omp hist = context->hist;
    int n;
    int m;

#ifdef TPCF
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
#endif
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
//B 2pcf
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
//E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
#ifdef TPCF
        hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        
        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,
                      1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

        hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif

   for (n = 1; n <= cmd->sizeHistN; n++) {
       hist->histNthread[n] = 0.0;
       hist->histNNSubthread[n] = 0.0;
       hist->histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
       hist->histNNSubXi2pcfthreadp[n] = 0.0;
       hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
       hist->histXi2pcfthread[n] = 0.0;
       hist->histXi2pcfthreadsub[n] = 0.0;
   }

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        }
#endif

    return SUCCESS;
}

global int search_init_sincos_omp(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  gdhistptr_sincos_omp hist)
{
    sincos_hist_init_context context;
    ErrorMsg allocation_error;

    memset(hist, 0, sizeof(*hist));
    context.cmd = cmd;
    context.gd = gd;
    context.hist = hist;
    if (cballs_allocation_guard(search_init_sincos_omp_unguarded,
                                &context, allocation_error,
                                sizeof(allocation_error)) == FAILURE) {
        search_free_sincos_omp(cmd, gd, hist);
        return FAILURE;
    }
    return SUCCESS;
}

global int search_free_sincos_omp(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
#define FREE_DVECTOR_IF_SET(p,nl,nh) \
    do { if ((p) != NULL) { free_dvector((p),(nl),(nh)); (p) = NULL; } } while (0)
#define FREE_DMATRIX_IF_SET(p,nrl,nrh,ncl,nch) \
    do { if ((p) != NULL) { free_dmatrix((p),(nrl),(nrh),(ncl),(nch)); (p) = NULL; } } while (0)
#define FREE_DMATRIX3D_IF_SET(p,nrl,nrh,ncl,nch,ndl,ndh) \
    do { if ((p) != NULL) { free_dmatrix3D((p),(nrl),(nrh),(ncl),(nch),(ndl),(ndh)); (p) = NULL; } } while (0)
#ifdef TPCF
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX_IF_SET(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadcossin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadsincos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX3D_IF_SET(hist->histZetaMthreadcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        FREE_DMATRIX_IF_SET(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif
    FREE_DVECTOR_IF_SET(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histXi2pcfthread,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
    FREE_DVECTOR_IF_SET(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNNSubthread,1,cmd->sizeHistN);
    FREE_DVECTOR_IF_SET(hist->histNthread,1,cmd->sizeHistN);
#ifdef TPCF
        FREE_DVECTOR_IF_SET(hist->ChebsU,1,cmd->mChebyshev+1);
        FREE_DVECTOR_IF_SET(hist->ChebsT,1,cmd->mChebyshev+1);
#endif

#undef FREE_DVECTOR_IF_SET
#undef FREE_DMATRIX_IF_SET
#undef FREE_DMATRIX3D_IF_SET

    return SUCCESS;
}

global int computeBodyProperties_sincos(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp hist)
{
    int n;
    int m;
    real xi = 0.0;
    real xi_2p = 0.0;

    //B Normalization of histograms
    if (Type(p) == BODY) {
#ifdef NOSTANDARNORMHIST
        xi = Kappa(p);
        xi_2p = Kappa(p);
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
            xi_2p = KappaRmin(p);
            xi = NbRmin(p)*xi_2p;
    }
#endif
#else // ! NOSTANDARNORMHIST
        xi = Kappa(p)/nbody;
#ifdef BALLS4SCANLEV
        xi_2p = (Weight(p)/Nb(p))*Kappa(p);
#else
        xi_2p = Weight(p)*Kappa(p);
#endif
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
#ifdef BALLS4SCANLEV
            xi_2p = KappaRmin(p);
#endif
            xi = NbRmin(p)*xi_2p/nbody;
    }
#endif
#endif // ! NOSTANDARNORMHIST
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
    //E Normalization of histograms
    const bool raw_multipoles = cballs_raw_legacy_multipoles(cmd);
    if (raw_multipoles) {
        xi = Weight(p)*Kappa(p);
        xi_2p = xi;
    }

#ifdef TPCF
        for (m=1; m<=cmd->mChebyshev+1; m++)
            //B Normalization of histograms
#ifdef NOSTANDARNORMHIST
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= 1.0;
                hist->histXithreadsin[m][n] /= 1.0;
            }
#else
            for (n=1; n<=cmd->sizeHistN; n++) {
                if (!raw_multipoles) {
                    hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                    hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                }
            }
#endif
            //E
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos,
                      hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin,
                      hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos,
                      hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            OUTVP_ext(hist->xiOUTVPcossin,
                      hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMtmpcossin,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,
                      hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            MULMS_ext(hist->histZetaMtmpcossin,
                      hist->xiOUTVPcossin,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],
                     hist->histZetaMthreadcos[m],
                     hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],
                     hist->histZetaMthreadsin[m],
                     hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],
                     hist->histZetaMthreadsincos[m],
                     hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            ADDM_ext(hist->histZetaMthreadcossin[m],
                     hist->histZetaMthreadcossin[m],
                     hist->histZetaMtmpcossin,cmd->sizeHistN);
        }
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
    }

    return SUCCESS;
}

global int search_init_gd_hist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;
    int n1, n2, l;

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaM[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
    }
    
#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaGmRe[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaGmIm[m], cmd->sizeHistN);
        }
#endif

    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

global int search_init_gd_hist_sincos(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
#ifdef TPCF
            for (m = 1; m <= cmd->mChebyshev+1; m++) {
                // HERE MUST BE gd->histXicos and gd->histXisin
            }
#endif
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

//B Computation of histogram of all B-B encounters
// The correlation function is estimated as:
//    xi=(V/v(r))*(DD(r)/N^2)
// where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
// total # particles.
//
// Note: only rminHistN = 0 works and agree with CUTE_BOX
//  you may try options=cute-box-rmin to correct a bit the results...
//      but for biger values of rminHistN differences grow...
//
local int search_compute_Xi(struct  cmdline_data* cmd,
                            struct  global_data* gd, int nbody)
{
    int k;
    int n;
    real normFac;
    real Vol;
    //B correct cute-box-rmin
    real deltaR = gd->deltaR;
    if ((cballs_opt_cute_box_rmin(cmd)))
        deltaR = cmd->rangeN/cmd->sizeHistN;
    //E

    Vol = 1.0;
    DO_COORD(k)
        Vol = Vol*gd->Box[k];

if (!cmd->useLogHist) {
    if ((cballs_opt_cute_box(cmd))) {
        gd->histNN[1]-=nbody;
    }
}
    real *edd;
    real *corr;
    real *ercorr;
    edd = dvector(1,cmd->sizeHistN);
    corr = dvector(1,cmd->sizeHistN);
    ercorr = dvector(1,cmd->sizeHistN);
    real rho_av=(real)nbody/Vol;

    for (n = 1; n <= cmd->sizeHistN; n++)
        edd[n] = 1./rsqrt(gd->histNN[n]);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        if(gd->histNN[n]==0) {
            corr[n]=0;
            ercorr[n]=0;
        } else {
            double r0,r1,vr,rho_r;
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    r0 = rpow(10.0, ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                    r1 = rpow(10.0, ((real)(n+1-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                } else {
                    r0 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n))*gd->deltaR );
                    r1 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n+1))*gd->deltaR );
                }
            } else {
                //B correct cute-box-rmin
                if ((cballs_opt_cute_box_rmin(cmd))) {
                    r0=(real)n*deltaR;
                    r1=(real)(n+1)*deltaR;
                } else {
                    r0=(real)n*gd->deltaR;
                    r1=(real)(n+1)*gd->deltaR;
                }
            }

#if (NDIM==3)
            if (cballs_opt_cute_box(cmd)) {
                //B this version does not give same results as CB
                //      although the programming is the same...
                vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
                rho_r=gd->histNN[n]/((real)nbody*vr);
                corr[n]=rho_r/rho_av-1;             // Correlation function
                ercorr[n]=(1+corr[n])*edd[n];       // Poisson errors
                gd->histCF[n] = corr[n];            // Original line
                //E
            } else {
                if (cmd->useLogHist) {
                    vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
// rho_r/rho_av = ( histNN[n]/(nbody*vr) ) / (nbody/Vol)
                    normFac = Vol/(vr*((real)(nbody*nbody)));
                    gd->histCF[n] = gd->histNN[n] * normFac - 1.0;
                } else {
                    //B correct cute-box-rmin
                    if ((cballs_opt_cute_box_rmin(cmd))) {
                        normFac = Vol/(2.0*PI*rpow(deltaR,3.0)*nbody*nbody);
                    } else {
                        normFac = Vol/(2.0*PI*rpow(gd->deltaR,3.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5);
// This line gives results in agreement with CB:
                    }
                    gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5) -1.0;
                    //E
                }
            }
#else
            if (cballs_opt_cute_box(cmd)) {
                // This should be CB version...
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            } else {
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5);
// This line gives results in agreement with CB:
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            }
#endif // ! NDIM
        }
    }

    free_dvector(ercorr,1,cmd->sizeHistN);
    free_dvector(corr,1,cmd->sizeHistN);
    free_dvector(edd,1,cmd->sizeHistN);

    return SUCCESS;
}


global int search_compute_HistN(struct  cmdline_data* cmd, 
                                struct  global_data* gd, int nbody)
{
    int n;
    real normFac;

//B Check this factor is correct...
// to agree with cute_box normalization commented out these lines
//B these does not work!!
//    normFac = 1.0;
//E
    normFac = 0.5;
    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNN[n] *= normFac;
//E
    if (cballs_opt_and_cf(cmd))
        search_compute_Xi(cmd, gd, nbody);

    return SUCCESS;
}


// No se usa qsize: update
global bool reject_cell(struct  cmdline_data* cmd, struct
                        global_data* gd, nodeptr p, nodeptr q, real qsize)
{
    real drpq, drpq2;
    compute_vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}


//B 2023.11.22
global bool reject_balls(struct  cmdline_data* cmd,
                         struct  global_data* gd, nodeptr p, nodeptr q,
                         real *drpq, compute_vector dr)
{
    real drpq2;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    *drpq = rsqrt(drpq2);

    if ( *drpq >= gd->Rcut + Radius(p) + Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}
//E

global bool reject_cell_balls(struct  cmdline_data* cmd,
                              struct  global_data* gd, nodeptr p, nodeptr q,
                              real *drpq, compute_vector dr)
{
    real drpq2;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    *drpq = rsqrt(drpq2);

    if ( *drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

global bool reject_bodycell(struct  cmdline_data* cmd,
                            struct  global_data* gd, nodeptr p, nodeptr q)
{
    real drpq, drpq2;
    compute_vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

global bool reject_cellcell(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q)
{
    real drpq, drpq2;
    compute_vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(p)+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

global bool accept_body(struct  cmdline_data* cmd, struct  global_data* gd,
                        bodyptr p, nodeptr q, real *drpq, compute_vector dr)
{
    return cballs_accept_body_contract(cmd, gd, p, q, drpq, dr);
}

#ifdef SMOOTHPIVOT
local size_t smooth_hash_cell(const INTEGER cell[NDIM], size_t mask)
{
    unsigned long long hash = 1469598103934665603ULL;
    int k;

    DO_COORD(k) {
        hash ^= (unsigned long long)cell[k];
        hash *= 1099511628211ULL;
    }

    return (size_t)hash & mask;
}

local int smooth_body_cell(struct cmdline_data* cmd,
                           struct global_data* gd,
                           bodyptr p,
                           const real cell_width[NDIM],
                           const INTEGER periodic_cells[NDIM],
                           INTEGER cell[NDIM])
{
    int k;

    DO_COORD(k) {
        real coordinate;

        if (!isfinite(Pos(p)[k])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "prepare_smooth_pivots: non-finite body coordinate[%d]=%g",
                     k, Pos(p)[k]);
            return FAILURE;
        }

        if (cmd->usePeriodic) {
            real box = gd->Box[k];

            if (!isfinite(box) || box <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "prepare_smooth_pivots: invalid periodic box[%d]=%g",
                         k, box);
                return FAILURE;
            }

            coordinate = Pos(p)[k] + 0.5*box;
            coordinate -= floor(coordinate/box)*box;
            cell[k] = (INTEGER)floor(coordinate/cell_width[k]);
            if (cell[k] >= periodic_cells[k])
                cell[k] = periodic_cells[k] - 1;
        } else {
            cell[k] = (INTEGER)floor(Pos(p)[k]/cell_width[k]);
        }
    }

    return SUCCESS;
}

local int smooth_neighbor_cells(INTEGER center, INTEGER periodic_cells,
                                bool periodic, INTEGER neighbors[3])
{
    int count = 0;
    int offset;

    for (offset = -1; offset <= 1; offset++) {
        INTEGER candidate = center + offset;
        int i;

        if (periodic) {
            candidate %= periodic_cells;
            if (candidate < 0)
                candidate += periodic_cells;
        }

        for (i = 0; i < count; i++)
            if (neighbors[i] == candidate)
                break;
        if (i == count)
            neighbors[count++] = candidate;
    }

    return count;
}

local void smooth_claim_cell(struct cmdline_data* cmd,
                             struct global_data* gd,
                             bodyptr p, bodyptr scan_table,
                             const INTEGER target[NDIM],
                             const INTEGER *body_cells,
                             const INTEGER *bucket_heads,
                             const INTEGER *bucket_next,
                             size_t bucket_mask)
{
    INTEGER iq;
    size_t bucket = smooth_hash_cell(target, bucket_mask);

    for (iq = bucket_heads[bucket]; iq >= 0; iq = bucket_next[iq]) {
        bodyptr q;
        real dr1;
        compute_vector dr;
        int k;

        DO_COORD(k)
            if (body_cells[(size_t)iq*NDIM + k] != target[k])
                break;
        if (k != NDIM)
            continue;

        q = scan_table + iq;
        if (p == q)
            continue;
        if (cballs_opt_read_mask(cmd)
            && Mask(q) == MASK_NODE_MASKED)
            continue;
        if (!accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)
            || dr1 > gd->rsmooth[0])
            continue;

        if (Update(q) == TRUE) {
            Update(q) = FALSE;
            NbRmin(p) += 1;
#ifndef NOWKAvg
            KappaRmin(p) += Weight(q)*Kappa(q);
#else
            KappaRmin(p) += Kappa(q);
#endif
            WeightRmin(p) += Weight(q);
        } else {
            NbRminOverlap(p) += 1;
        }
    }
}

/*
 * Build smoothing groups in stable pivot order before an OpenMP search starts.
 * Search workers may read the resulting body fields but must not mutate them.
 */
global int prepare_smooth_pivots(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 bodyptr *btable, INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2)
{
    INTEGER first, last, npivot, nscan;
    INTEGER *body_cells = NULL;
    INTEGER *bucket_heads = NULL;
    INTEGER *bucket_next = NULL;
    INTEGER periodic_cells[NDIM];
    real cell_width[NDIM];
    size_t bucket_count = 16;
    size_t bucket_mask;
    bodyptr p, q;
    int k;
    int status = FAILURE;

    if (btable == NULL || nbody == NULL || ipmax == NULL
        || cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid catalog arguments");
        return FAILURE;
    }

    first = ipmin - 1;
    last = ipmax[cat1];
    npivot = nbody[cat1];
    nscan = nbody[cat2];
    if (first < 0 || last < first || last > npivot || nscan < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid pivot or scan bounds");
        return FAILURE;
    }
    if (!isfinite(gd->rsmooth[0]) || gd->rsmooth[0] < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid smoothing radius %g",
                 gd->rsmooth[0]);
        return FAILURE;
    }
    if (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0) {
        DO_BODY(q, btable[cat2], btable[cat2] + nscan)
            Update(q) = TRUE;
        if (cat1 != cat2)
            DO_BODY(p, btable[cat1], btable[cat1] + npivot)
                Update(p) = TRUE;
        DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            KappaRmin(p) = Kappa(p);
            WeightRmin(p) = Weight(p);
        }
        return SUCCESS;
    }
    if ((size_t)nscan > ((size_t)-1)/(NDIM*sizeof(INTEGER))
        || (size_t)nscan > ((size_t)-1)/2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: catalog is too large");
        return FAILURE;
    }

    while (bucket_count < (size_t)nscan*2) {
        if (bucket_count > ((size_t)-1)/2) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "prepare_smooth_pivots: hash table size overflow");
            return FAILURE;
        }
        bucket_count <<= 1;
    }
    if (bucket_count > ((size_t)-1)/sizeof(*bucket_heads)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: hash table allocation overflow");
        return FAILURE;
    }
    bucket_mask = bucket_count - 1;

    body_cells = malloc((size_t)nscan*NDIM*sizeof(*body_cells));
    bucket_next = malloc((size_t)nscan*sizeof(*bucket_next));
    bucket_heads = malloc(bucket_count*sizeof(*bucket_heads));
    if (body_cells == NULL || bucket_next == NULL || bucket_heads == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: unable to allocate temporary ownership map");
        goto cleanup;
    }

    DO_COORD(k) {
        if (cmd->usePeriodic) {
            real box = gd->Box[k];
            if (!isfinite(box) || box <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "prepare_smooth_pivots: invalid periodic box[%d]=%g",
                         k, box);
                goto cleanup;
            }
            periodic_cells[k] = (INTEGER)floor(box/gd->rsmooth[0]);
            if (periodic_cells[k] < 1)
                periodic_cells[k] = 1;
            cell_width[k] = box/(real)periodic_cells[k];
        } else {
            periodic_cells[k] = 0;
            cell_width[k] = gd->rsmooth[0];
        }
    }

    for (size_t ibucket = 0; ibucket < bucket_count; ibucket++)
        bucket_heads[ibucket] = -1;

    DO_BODY(q, btable[cat2], btable[cat2] + nscan) {
        INTEGER iq = q - btable[cat2];
        INTEGER *cell = body_cells + (size_t)iq*NDIM;
        size_t bucket;

        if (smooth_body_cell(cmd, gd, q, cell_width, periodic_cells, cell)
            == FAILURE)
            goto cleanup;
        bucket = smooth_hash_cell(cell, bucket_mask);
        bucket_next[iq] = bucket_heads[bucket];
        bucket_heads[bucket] = iq;
        Update(q) = TRUE;
    }

    if (cat1 != cat2)
        DO_BODY(p, btable[cat1], btable[cat1] + npivot)
            Update(p) = TRUE;

    DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
        NbRmin(p) = 1;
        NbRminOverlap(p) = 0;
        KappaRmin(p) = Kappa(p);
        WeightRmin(p) = Weight(p);
    }

    DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
        INTEGER center[NDIM];
        INTEGER neighbors[NDIM][3];
        int neighbor_count[NDIM];

        if (Update(p) == FALSE)
            continue;
        if (smooth_body_cell(cmd, gd, p, cell_width, periodic_cells, center)
            == FAILURE)
            goto cleanup;
        DO_COORD(k)
            neighbor_count[k] = smooth_neighbor_cells(center[k],
                                                       periodic_cells[k],
                                                       cmd->usePeriodic,
                                                       neighbors[k]);

#if NDIM == 3
        for (int i0 = 0; i0 < neighbor_count[0]; i0++)
            for (int i1 = 0; i1 < neighbor_count[1]; i1++)
                for (int i2 = 0; i2 < neighbor_count[2]; i2++) {
                    INTEGER target[NDIM] = {
                        neighbors[0][i0], neighbors[1][i1], neighbors[2][i2]
                    };
                    smooth_claim_cell(cmd, gd, p, btable[cat2], target,
                                      body_cells, bucket_heads, bucket_next,
                                      bucket_mask);
                }
#elif NDIM == 2
        for (int i0 = 0; i0 < neighbor_count[0]; i0++)
            for (int i1 = 0; i1 < neighbor_count[1]; i1++) {
                INTEGER target[NDIM] = {
                    neighbors[0][i0], neighbors[1][i1]
                };
                smooth_claim_cell(cmd, gd, p, btable[cat2], target,
                                  body_cells, bucket_heads, bucket_next,
                                  bucket_mask);
            }
#else
#error prepare_smooth_pivots supports only NDIM=2 or NDIM=3
#endif
    }

    status = SUCCESS;

cleanup:
    free(bucket_heads);
    free(bucket_next);
    free(body_cells);
    return status;
}
#endif


//B Calculate number of threads
// Use in the terminal:
// export OMP_NUM_THREADS=8
// to set the maximum number of threads that can be used in a run
global int ThreadCount(struct  cmdline_data* cmd, struct  global_data* gd,
                       INTEGER nbody, int cat1) {
    int nthreads=0;
#pragma omp parallel
    {
#pragma omp atomic
        nthreads++;
    }
    verb_print(cmd->verbose, "\nUsing %d threads \n",nthreads);

//B socket:
#ifdef ADDONS
#include "cballsutils_include_01.h"
#endif
//E

    return SUCCESS;
}
//E


//B coordinate transformations routines

//
// All angles and input distances are in radians.
// Will be useful to convert to degrees, hours, arcmin and arcsec.
// And write this info in log file.
//
// deg = rads * 180/PI // rads to degrees.
// deg = rads * 180/(PI*15.0) // rads to hours.
// deg = (rads * 180/PI) * 60 // rads to arcmin.
// deg = (rads * 180/PI) * 60 * 60 // rads to arcsec.
//
// Also consider wrapping angles to the range (-pi, pi]:
// if ( angle <= -PI ) angle += 2.0*PI
// if ( angle >   PI ) angle -= 2.0*PI
//
// Transform to x,y,z coordinates:
// ra = rhs.getX();
// dec = rhs.getY();
//    x = cosdec * cosra;
//    y = cosdec * sinra;
//    z = sindec;
//

#define ARFKEN              0
#define NOARFKEN            1
#define GALACTIC            2
#define ECLIPTIC            3
#define CELESTIAL           4                   // also known as equatorial

local int spherical_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                  real theta, real phi, vector xyz);
local int galactic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                 real latitud, real longitud, vector xyz);
local int ecliptic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                 real DEC, real RA, vector xyz);
local int celestial_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                  real DEC, real RA, vector xyz);

global int coordinate_string_to_int(struct cmdline_data* cmd,
                                    struct  global_data* gd)
{
    gd->coordTag = -1;
    gd->coordTag = ARFKEN;
    if (cballs_opt_arfken(cmd)) gd->coordTag = ARFKEN;
    if (cballs_opt_no_arfken(cmd)) gd->coordTag = NOARFKEN;
    if (cballs_opt_galactic(cmd)) gd->coordTag = GALACTIC;
    if (cballs_opt_ecliptic(cmd)) gd->coordTag = ECLIPTIC;
    if (cballs_opt_celestial(cmd)) gd->coordTag = CELESTIAL;

    return SUCCESS;
}

global int coordinate_transformation(struct cmdline_data* cmd, struct  global_data* gd,
                                    real theta, real phi, vector xyz)
{
    switch(gd->coordTag) {
        case ARFKEN:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case NOARFKEN:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case GALACTIC:
            galactic_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case ECLIPTIC:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case CELESTIAL:
            celestial_to_cartesians(cmd, gd, theta, phi, xyz); break;
        default:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
    }

    return SUCCESS;
}

local int spherical_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real theta, real phi, vector xyz)
{
    real ra, dec;

    if (cballs_opt_no_arfken(cmd)) {
        //B this would be a galactic->cartesian tranformation...
        //      if RA is theta and DEC is phi
        //      more preciselly: RA is the longitud (l)
        //                      DEC is the latitud (b)
        ra = theta;
        dec = PIO2 - phi;
        xyz[0] = rcos(dec)*rcos(ra);
        xyz[1] = rcos(dec)*rsin(ra);
        xyz[2] = rsin(dec);
        //E
    } else {
        ra = phi;
        dec = theta;
// Standard transformation. See Arfken
        xyz[0] = rsin(dec)*rcos(ra);
        xyz[1] = rsin(dec)*rsin(ra);
        xyz[2] = rcos(dec);
    }

    return SUCCESS;
}

// latitud (b): the angle of an object northward of the galactic equator
// longitud (l): the angular distance of an object eastward along the galactic equator
local int galactic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real latitud, real longitud, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\ngalactic_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    phi = longitud;
    theta = latitud;
    //B Standard transformation. See Arfken
    xyz[0] = rsin(theta)*rcos(phi);
    xyz[1] = rsin(theta)*rsin(phi);
    xyz[2] = rcos(theta);
    //E

    return SUCCESS;
}

local int ecliptic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real DEC, real RA, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\necliptic_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    if (cballs_opt_ra_reversed(cmd)) {
        phi = TWOPI - RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    } else {
        phi = RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    }
    
    return SUCCESS;
}

local int celestial_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real DEC, real RA, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\ncelestial_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    if (cballs_opt_ra_reversed(cmd)) {
        phi = TWOPI - RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    } else {
        phi = RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    }
    
    return SUCCESS;
}

global int spherical_periodic_condition(real *thetaL, real *thetaR,
                                        real *phiL, real *phiR)
{
    while(*thetaL<0) {
        *thetaL += 2.0*PI;
    }
    while(*thetaL>2.0*PI) {
        *thetaL -= 2.0*PI;
    }
    while(*thetaR<0) {
        *thetaR += 2.0*PI;
    }
    while(*thetaR>2.0*PI) {
        *thetaR -= 2.0*PI;
    }

    while(*phiL<0) {
        *phiL += PI;
    }
    while(*phiL>PI) {
        *phiL -= PI;
    }
    while(*phiR<0) {
        *phiR += PI;
    }
    while(*phiR>PI) {
        *phiR -= PI;
    }

    return SUCCESS;
}

#undef ARFKEN
#undef NOARFKEN
#undef GALACTIC
#undef ECLIPTIC
#undef CELESTIAL                                // also known as equatorial
#undef NULLTRANSFORM

//E coordinate transformations routines


//B section of several routines to do pre/post processing

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETASTDDEV \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETAHEADER \
"# [1] rBins; [2] diagonal; [3] theta2=Nbins/4.0; [4] theta2=2.0*Nbins/4.0; \
[5] theta2=3.0*Nbins/4.0; [6] theta2=4.0*Nbins/4.0 - 1.0\n"

#define MHISTZETAHEADERSTDDEV \
"# [1] rBins; [2] diagonal &SD; [3] theta2=Nbins/4.0 &SD; [4] theta2=2.0*Nbins/4.0 &SD; \
[5] theta2=3.0*Nbins/4.0 &SD; [6] theta2=4.0*Nbins/4.0 - 1.0 &SD\n"

// routine to compute covariance matriz of correlations of Takahashi simulations
//  test were done on:
//  run/Cosma/Takahasi_nres12_balls-omp/zs9_balls-omp/rxxx/Output files
global int statHistogram(struct cmdline_data* cmd, struct  global_data* gd)
{
    string routine_name = "statHistogram";
    string routineName = routine_name;
    char namebuf1[256];
    char namebuf2[256];
    char namebuf3[256];
    char namebuf4[256];
    char namebuf5[256];
    char namebuf6[256];
    struct stat buf;
    char rootDirPath[MAXLENGTHOFFILES];
    int nrealization = 0;
    int status1 = 1;
    int status2 = 1;
    int status3 = 1;
    int ifile = 0;
    int i;
    int j;
    int m;
    int n1;
    int n2;
    int npts1;
    int npts2;
    int npts3;
    real matElement;
    real matElement2;
    real matElementIm;
    int sizeHistN;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    real ZetaStdDev;
    real Zeta2StdDev;
    real Zeta3StdDev;
    real Zeta4StdDev;
    real Zeta5StdDev;
    int Nbins;
    int NR;

    //B
    stream instr1 = NULL, instr2 = NULL, instr3 = NULL;
    stream outstr1 = NULL, outstr2 = NULL;

    real **mat1 = NULL, **mat2 = NULL, **mattmp = NULL;
    real **Zv = NULL, **ZvAvg = NULL;

    real ***mat3 = NULL;
    real ***matAvg = NULL;
    real ***matStdDev = NULL;
    real ***matCovMat = NULL;

    real *rBin = NULL;

    int Nv = 0;
    int N = 0;
    int rc = FAILURE;
    //E
    
    //B
#define STAT_FAIL(...)                                      \
    do {                                                    \
        snprintf(cmd->error_message, _ERRORMSGSIZE_,        \
                 __VA_ARGS__);                              \
        goto fail;                                          \
    } while (0)

#define STAT_FORMAT_OR_FAIL(dst, size, label, ...)          \
    do {                                                    \
        if (format_checked((dst), (size), (label),          \
                           __VA_ARGS__) != 0) {             \
            STAT_FAIL("%s: could not format %s",            \
                      routine_name, (label));               \
        }                                                   \
    } while (0)
    //E
    
    sizeHistN = cmd->sizeHistN;
    mat1 = dmatrix(1,sizeHistN,1,sizeHistN);
    mat2 = dmatrix(1,sizeHistN,1,sizeHistN);

//B read one file to got needed info to go...
    verb_print(cmd->verbose,
               "\n%s: first read two file to get some needed info:",
               routine_name);
    STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                        "rootDirPath", "%s/%s%03d/%s/",
                        cmd->rootDir, "r", nrealization,
                        cmd->suffixOutFiles);
    m=1;
    //B file 1
    STAT_FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                        "namebuf1", "%s%s_%d%s",
                        rootDirPath, gd->infilenames[ifile], m, EXTFILES);

    verb_print(cmd->verbose,
                "\n%s: opening file %s...",
               routine_name,namebuf1);
    if (stat(namebuf1, &buf) == 0)               // no input file exists?
        OPEN_OUTPUT_OR_FAIL(instr1, namebuf1, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name,namebuf1);
        status1 = 0;
    }
    //E file 1
    //B file 2
    STAT_FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                        "namebuf2", "%s%s_%d%s",
                        rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);

    verb_print(cmd->verbose,
               "\n%s: opening file %s...",
               routine_name,namebuf2);
    if (stat(namebuf2, &buf) == 0)               // no input file exists?
        OPEN_OUTPUT_OR_FAIL(instr2, namebuf2, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name,namebuf2);
        status2 = 0;
    }
    //E file 2
        
    if (status1 == 0 || status2 == 0) {
        if (status1 != 0)
            CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);

        if (status2 != 0)
            CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);

        STAT_FAIL("\n%s: one of the input files does not exist: %s %s",
                  routine_name, namebuf1, namebuf2);
    } else {
        //B processing input files
        inout_InputDataMatrix(namebuf1, mat1, &npts1,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        inout_InputDataMatrix(namebuf2, mat2, &npts2,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts1 != npts2) {
            STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                      routine_name, nrealization, "npts are different:", npts1,npts2);
        }
        cmd->sizeHistN = npts1;
        CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);

        CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
        //E
    } // ! status
    verb_print(cmd->verbose,
               "\n%s: mChebyshev and sizeHistN: %d %d\n",
               routine_name, cmd->mChebyshev, cmd->sizeHistN);
    
    free_dmatrix(mat2, 1, sizeHistN, 1, sizeHistN);
    mat2 = NULL;
    free_dmatrix(mat1, 1, sizeHistN, 1, sizeHistN);
    mat1 = NULL;
    
//E read one file


//B read rBin file
    rBin = dvector(1,sizeHistN);
    //B file 3
    verb_print(cmd->verbose,
               "\n%s: reading rbins...",
               routine_name, cmd->mChebyshev, cmd->sizeHistN);
    STAT_FORMAT_OR_FAIL(namebuf5, sizeof(namebuf5),
                        "namebuf5", "%s%s%s",
                        rootDirPath, "rbins", EXTFILES);

    verb_print(cmd->verbose,
               "\n%s: opening file %s...",routine_name, namebuf5);
    if (stat(namebuf5, &buf) == 0)               // no input file exists?
        OPEN_OUTPUT_OR_FAIL(instr3, namebuf5, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name, namebuf5);
        status3 = 0;
    }
    //E file 3
    if (status3 == 0) {
        STAT_FAIL("\n%s: the input file does not exist: %s",
                  routine_name, namebuf5);
    } else {
        //B processing input file
        inout_InputDataVector(namebuf5, rBin, &npts3,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts3 != cmd->sizeHistN) {
            STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                      routine_name, nrealization, "npts is not equal to sizeHistN:",
                      npts3,cmd->sizeHistN);
        }
        CLOSE_OUTPUT_OR_FAIL(instr3, namebuf5);
        //E
    } // ! status3
    verb_print(cmd->verbose, "\ndone.\n");
//B read rBin file

    Nbins = cmd->sizeHistN;

    mat1 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat2 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat3 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    matAvg =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    matStdDev =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);


    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(matAvg[m], cmd->sizeHistN);
        CLRM_ext(matStdDev[m], cmd->sizeHistN);
    }

//B read realization files and compute mean values
    verb_print(cmd->verbose,
               "\n%s: mean value computation:", routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/%s%03d/%s/",
                            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);

        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            STAT_FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                                "namebuf1", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) == 0)               // no input file exists?
                OPEN_OUTPUT_OR_FAIL(instr1, namebuf1, "r");
            else {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            STAT_FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                            "namebuf2", "%s%s_%d%s",
                            rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",
                       routine_name, namebuf2);
            if (stat(namebuf2, &buf) == 0)               // no input file exists?
                OPEN_OUTPUT_OR_FAIL(instr2, namebuf2, "r");
            else {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
                if (status1 != 0) //fclose(instr1);
                CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);

                if (status2 != 0) //fclose(instr2);
                CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
                break;
            } else {
                //B output file 1
                STAT_FORMAT_OR_FAIL(namebuf3, sizeof(namebuf3),
                                    "namebuf3", "%s%s_%d%s",
                                    rootDirPath, cmd->outfile, m, EXTFILES);

                verb_print(cmd->verbose,
                           "\n%s: opening file %s... to save statistics",
                           routine_name, namebuf3);
                OPEN_OUTPUT_OR_FAIL(outstr1, namebuf3, "w!");
                //E output file 1
                //B output file 2
                STAT_FORMAT_OR_FAIL(namebuf4, sizeof(namebuf4),
                                    "namebuf4", "%s%s%s_%d%s",
                                    rootDirPath, "m", cmd->outfile, m, EXTFILES);

                verb_print(cmd->verbose,
                           "\n%s: opening file %s... to save statistics",
                           routine_name, namebuf4);
                OPEN_OUTPUT_OR_FAIL(outstr2, namebuf4, "w!");
                //E output file 2

                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                              routine_name, nrealization,
                              "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        matElement = mat1[n1][n2]+mat2[n1][n2];
                        if (cballs_opt_same_infiles(cmd))
                            matElement *= 0.5;
                        mat3[m][n1][n2] = matElement;
                        matAvg[m][n1][n2] += matElement;
                        WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3,
                                             "%16.8e ",matElement);
                    }
                    WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3, "\n");
                }

                WRITE_OUTPUT_OR_FAIL(outstr2, namebuf4, MHISTZETAHEADER);
                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    Zeta = mat3[m][n1][n1];
                    Zeta2 = mat3[m][n1][(int)(Nbins/4.0)];
                    Zeta3 = mat3[m][n1][(int)(2.0*Nbins/4.0)];
                    Zeta4 = mat3[m][n1][(int)(3.0*Nbins/4.0)];
                    Zeta5 = mat3[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                    WRITE_OUTPUT_OR_FAIL(outstr2, namebuf4, MHISTZETA,
                                         rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
                }

                CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
                CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
                CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf3);
                CLOSE_OUTPUT_OR_FAIL(outstr2, namebuf4);
                //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/", cmd->rootDir);
        //E

        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            STAT_FORMAT_OR_FAIL(namebuf3, sizeof(namebuf3),
                                "namebuf3", "%s%s_%s_%d%s",
                                rootDirPath, cmd->outfile, "Avg", m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            OPEN_OUTPUT_OR_FAIL(outstr1, namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    matAvg[m][n1][n2] /= (nrealization-1);
                    WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3,
                                         "%16.8e ",matAvg[m][n1][n2]);
                }
                WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3, "\n");
            }
            CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf3);

            //B output file 2
            STAT_FORMAT_OR_FAIL(namebuf6, sizeof(namebuf6),
                            "namebuf6", "%s%s%s_%s_%d%s",
                            rootDirPath, "m", cmd->outfile, "Avg", m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf6);
            OPEN_OUTPUT_OR_FAIL(outstr2, namebuf6, "w!");
            //E output file 2
            WRITE_OUTPUT_OR_FAIL(outstr2, namebuf6, MHISTZETAHEADER);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                Zeta = matAvg[m][n1][n1];
                Zeta2 = matAvg[m][n1][(int)(Nbins/4.0)];
                Zeta3 = matAvg[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = matAvg[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = matAvg[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                WRITE_OUTPUT_OR_FAIL(outstr2, namebuf6,
                                MHISTZETA,rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            CLOSE_OUTPUT_OR_FAIL(outstr2, namebuf6);

        } // ! end m loop
    } // ! nrealization > 3
//E read realization files and compute mean values


//
//B read realization files and compute std values
//
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
               "\n\n%s: standard deviation computation:", routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/%s%03d/%s/",
                            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E

        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            STAT_FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                                "namebuf1", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            STAT_FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                                "namebuf2", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                              routine_name, nrealization,
                              "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (cballs_opt_same_infiles(cmd))
                            matElement = rabs(
                                           0.5*(mat1[n1][n2]+mat2[n1][n2])
                                           -matAvg[m][n1][n2]
                                           );
                        else
                            matElement = rabs(
                                           mat1[n1][n2]+mat2[n1][n2]
                                           -matAvg[m][n1][n2]
                                           );

                        matStdDev[m][n1][n2] += matElement;
                    }
                }
                //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            STAT_FORMAT_OR_FAIL(namebuf3, sizeof(namebuf3),
                                "namebuf3", "%s%s_%s_%d%s",
                            rootDirPath, cmd->outfile, "StdDev", m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            OPEN_OUTPUT_OR_FAIL(outstr1, namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    matStdDev[m][n1][n2] /= (nrealization-1);
                    WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3,
                                         "%16.8e ",matStdDev[m][n1][n2]);
                }
                WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3, "\n");
            }
            CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf3);

            //B output file 2
            STAT_FORMAT_OR_FAIL(namebuf6, sizeof(namebuf6),
                        "namebuf6", "%s%s%s_%s_%d%s",
                        rootDirPath, "m", cmd->outfile, "StdDev", m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf6);
            OPEN_OUTPUT_OR_FAIL(outstr2, namebuf6, "w!");
            //E output file 2
            WRITE_OUTPUT_OR_FAIL(outstr2, namebuf6, MHISTZETAHEADERSTDDEV);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                Zeta = matAvg[m][n1][n1];
                ZetaStdDev = matStdDev[m][n1][n1];
                Zeta2 = matAvg[m][n1][(int)(Nbins/4.0)];
                Zeta2StdDev = matStdDev[m][n1][(int)(Nbins/4.0)];
                Zeta3 = matAvg[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta3StdDev = matStdDev[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = matAvg[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta4StdDev = matStdDev[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = matAvg[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                Zeta5StdDev = matStdDev[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                WRITE_OUTPUT_OR_FAIL(outstr2, namebuf6, MHISTZETASTDDEV,
                                     rBin[n1],Zeta,ZetaStdDev,
                                     Zeta2,Zeta2StdDev,
                                     Zeta3,Zeta3StdDev,
                                     Zeta4,Zeta4StdDev,
                                     Zeta5,Zeta5StdDev
                                     );
            }
            CLOSE_OUTPUT_OR_FAIL(outstr2, namebuf6);

        } // ! end m loop
    } // ! nrealization > 3

//
//E read realization files and compute std values
//

//
//B read realization files and compute covariance matrices values
//
    Nbins = cmd->sizeHistN*cmd->sizeHistN;

    N = cmd->sizeHistN;
    Nv = N * N;
    
    Zv = dmatrix(1,cmd->mChebyshev+1,1,Nv);
    ZvAvg = dmatrix(1,cmd->mChebyshev+1,1,Nv);

    mattmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    matCovMat =
            dmatrix3D(1,cmd->mChebyshev+1,1,Nv,1,Nv);

    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    CLRM_ext(mattmp, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(matCovMat[m], Nv);
    }
    CLRM_ext_ext(Zv, cmd->mChebyshev+1, Nv);
    CLRM_ext_ext(ZvAvg, cmd->mChebyshev+1, Nv);

    int n3, n4;
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
               "\n\n%s: covariance matrices computation:", routine_name);
    verb_print(cmd->verbose,
               "\n%s: reading files and computing covariance matrices...",
               routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/%s%03d/%s/",
                            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            STAT_FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                                "namebuf1", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            STAT_FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                                "namebuf2", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                            "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2

            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                              routine_name, nrealization,
                              "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (cballs_opt_same_infiles(cmd))
                            Zv[m][N*(n1-1)+n2] += 0.5*(mat1[n1][n2]+mat2[n1][n2]);
                        else
                            Zv[m][N*(n1-1)+n2] += mat1[n1][n2]+mat2[n1][n2];
                    }
                }
                    //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    NR = nrealization-1;

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, NR);
    
    verb_print(cmd->verbose,
               "\n%s: computing mean vectors...\n",
               routine_name);
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (i=1; i<=Nv; i++) {
            ZvAvg[m][i] = Zv[m][i]/((real)(NR-1));          // sample mean
        }
    }


    //B computing covariance matrices
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
    "\n%s: reading files again and final computation of covariance matrices...",
               routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/%s%03d/%s/",
                            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            STAT_FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                                "namebuf1", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            STAT_FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                                "namebuf2", "%s%s_%d%s",
                                rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                            "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2

            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    STAT_FAIL("\n%s: in realization r%03d: %s %d %d",
                              routine_name, nrealization,
                              "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (cballs_opt_same_infiles(cmd)) {
                            matElement=
                            0.5*(mat1[n1][n2]+mat2[n1][n2])-ZvAvg[m][N*(n1-1)+n2];
                            for (n3=1; n3<=cmd->sizeHistN; n3++) {
                                for (n4=1; n4<=cmd->sizeHistN; n4++) {
                                    matElement2=
                                    0.5*(mat1[n3][n4]+mat2[n3][n4])-ZvAvg[m][N*(n3-1)+n4];
                                    matCovMat[m][N*(n1-1)+n2][N*(n3-1)+n4] +=
                                    matElement*matElement2;
                                }
                            }
                        } else {
                            matElement=
                            mat1[n1][n2]+mat2[n1][n2]-ZvAvg[m][N*(n1-1)+n2];
                            for (n3=1; n3<=cmd->sizeHistN; n3++) {
                                for (n4=1; n4<=cmd->sizeHistN; n4++) {
                                    matElement2=
                                    mat1[n3][n4]+mat2[n3][n4]-ZvAvg[m][N*(n3-1)+n4];
                                    matCovMat[m][N*(n1-1)+n2][N*(n3-1)+n4] +=
                                    matElement*matElement2;
                                }
                            }

                        }
                    }
                }
                    //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    NR = nrealization-1;

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed (again): %d\n",
               routine_name, NR);

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (i=1; i<=Nv; i++) {
            for (int j=1; j<=Nv; j++) {
                matCovMat[m][i][j] /= ((real)(NR-1));          // sample mean
            }
        }
    }
    //E computing covariance matrices

    verb_print(cmd->verbose,
               "\n%s: saving covariance matrices...",
               routine_name);

    if (NR > 2) {
        //B we use rootDir string to construct input/output file names
        STAT_FORMAT_OR_FAIL(rootDirPath, sizeof(rootDirPath),
                            "rootDirPath", "%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            STAT_FORMAT_OR_FAIL(namebuf3, sizeof(namebuf3),
                            "namebuf3", "%s%s_%s_%d%s",
                            rootDirPath, cmd->outfile, "CovMat", m, EXTFILES);

            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            OPEN_OUTPUT_OR_FAIL(outstr1, namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=Nv; n1++) {
                for (n2=1; n2<=Nv; n2++) {
                    WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3,
                                         "%16.8e ",matCovMat[m][n1][n2]);
                }
                WRITE_OUTPUT_OR_FAIL(outstr1, namebuf3, "\n");
            }
            CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf3);
        } // ! end m loop
    } // ! nrealization > 3
    verb_print(cmd->verbose," done.\n");

//
//E read realization files and compute covariance matrices values
//
    
    //B
    rc = SUCCESS;

    fail:
        CLOSE_STREAM(instr1);
        CLOSE_STREAM(instr2);
        CLOSE_STREAM(instr3);
        CLOSE_STREAM(outstr1);
        CLOSE_STREAM(outstr2);

        if (ZvAvg != NULL && Nv > 0)
            free_dmatrix(ZvAvg, 1, cmd->mChebyshev+1, 1, Nv);
        if (Zv != NULL && Nv > 0)
            free_dmatrix(Zv, 1, cmd->mChebyshev+1, 1, Nv);

        if (matCovMat != NULL && Nv > 0)
            free_dmatrix3D(matCovMat, 1, cmd->mChebyshev+1, 1, Nv, 1, Nv);
        if (mattmp != NULL)
            free_dmatrix(mattmp, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (matStdDev != NULL)
            free_dmatrix3D(matStdDev, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (matAvg != NULL)
            free_dmatrix3D(matAvg, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat3 != NULL)
            free_dmatrix3D(mat3, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat2 != NULL)
            free_dmatrix(mat2, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat1 != NULL)
            free_dmatrix(mat1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (rBin != NULL)
            free_dvector(rBin, 1, cmd->sizeHistN);

    #undef STAT_FORMAT_OR_FAIL
    #undef STAT_FAIL

        return rc;
    //E

}

// routine to compute edge corrections using two saved histZetaM histograms
global int computeEdgeCorrections(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    string routineName = "computeEdgeCorrections";
    char namebuf1[256];
    char namebuf2[256];
    char namebuf3[256];
    char namebuf4[256];
    char namebuf5[256];
    struct stat buf;
    int status1 = 1;
    int status2 = 1;
    int status3 = 1;
    int ifile = 0;
    int m;
    int n1;
    int n2;
    int npts1;
    int npts2;
    int npts3;
    int mChebyshev;
    real rBin, rbinlog;
    
    //B initialize everything:
    stream instr1 = NULL, instr2 = NULL, instr3 = NULL;
    stream outstr1 = NULL, outstr2 = NULL;
    real **mat1 = NULL, **mat2 = NULL;
    real ***mat3 = NULL, ***mat4 = NULL, ***mat5 = NULL;
    real *rBins = NULL;
    //E

    char rootDirPath1[MAXLENGTHOFFILES];
    char preFileName1[MAXLENGTHOFFILES];
    char rootDirPath2[MAXLENGTHOFFILES];
    char preFileName2[MAXLENGTHOFFILES];
    
#define FORMAT_OR_FAIL(dst, size, label, ...)                         \
    do {                                                              \
        if (format_checked((dst), (size), (label), __VA_ARGS__) != 0) {\
            snprintf(cmd->error_message, _ERRORMSGSIZE_,              \
                     "\ncomputeEdgeCorrections: could not format %s: %s", \
                     (label), strerror(errno));                                        \
            goto fail;                                                \
        }                                                             \
    } while (0)
    
    if (extractInputRootDir(gd->infilenames[ifile],
                            rootDirPath1, sizeof(rootDirPath1),
                            preFileName1, sizeof(preFileName1), ifile,
                            cmd->verbose, cmd->verbose_log, gd->outlog) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "computeEdgeCorrections: invalid input path '%s'",
                 gd->infilenames[ifile]);
        goto fail;
    }
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: rootDir input file(%d) %s\n",
                ifile, rootDirPath1);
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: preFileName input file (%d) %s\n",
                 ifile, preFileName1);
    if (extractInputRootDir(gd->infilenames[ifile],
                            rootDirPath2, sizeof(rootDirPath2),
                            preFileName2, sizeof(preFileName2), ifile + 1,
                            cmd->verbose, cmd->verbose_log, gd->outlog) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "computeEdgeCorrections: invalid input path '%s'",
                 gd->infilenames[ifile]);
        goto fail;
    }
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: rootDir input file(%d) %s\n",
                ifile+1, rootDirPath2);
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: preFileName input file (%d) %s\n",
                 ifile+1, preFileName2);
    
    //B read one file to get needed info...
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: first read one file to get some needed info:");
    m=1;
    //B file 1
    FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                   "namebuf1", "%s/%s_%s_%d%s",
                   rootDirPath1, preFileName1, "cos", m, EXTFILES);


    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
    if (stat(namebuf1, &buf) == 0)               // no input file exists?
        OPEN_OUTPUT_OR_FAIL(instr1, namebuf1, "r");
    else {
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf1);
        status1 = 0;
    }
    //E file 1
    int nrow, ncol;
        
    if (status1 == 0) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: input file does not exist: %s",
                       namebuf1);
    } else {
        //B processing input files
        inout_InputDataMatrix_info(namebuf1, &nrow, &ncol,
                                   cmd->verbose, cmd->verbose_log, gd->outlog
                                   );
        if (nrow != ncol) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: in input files: %s %d %d",
                             "nrow and ncol are different:", nrow, ncol);
        }
        CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
        //E
    } // ! status
    cmd->sizeHistN = nrow;
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: sizeHistN: %d\n",
               cmd->sizeHistN);
    //E read one file to get needed info...

    // should know the size of the matrix
    mat1 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat2 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

//B read rBin file
    rBins = dvector(1,cmd->sizeHistN);
    //B file 3
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: reading rbins...");
    
    FORMAT_OR_FAIL(namebuf3, sizeof(namebuf3),
                   "namebuf3", "rbins");

    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: opening file %s...",namebuf3);
    if (stat(namebuf3, &buf) == 0)               // no input file exists?
        OPEN_OUTPUT_OR_FAIL(instr3, namebuf3, "r");
    else {
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf3);
        status3 = 0;
    }
    //E file 3
    if (status3 == 0) {
        COSMO_FAIL_GOTO(cmd, fail,
                   "\ncomputeEdgeCorrections: the input file does not exist: %s",
                         namebuf3);
    } else {
        //B processing input file
        inout_InputDataVector(namebuf3, rBins, &npts3,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts3 != cmd->sizeHistN) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: in rBin file: %s %d %d",
                             "npts is not equal to sizeHistN:",
                             npts3,cmd->sizeHistN);
        }
        CLOSE_OUTPUT_OR_FAIL(instr3, namebuf3);
        //E
    } // ! status3
    verb_print(cmd->verbose, "\ndone.\n");
//E read rBin file

    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: touching multipole files:");
    mChebyshev = 1;
    status1 = 1;
    do {
        FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                       "namebuf1", "%s/%s_%s_%d%s",
                       rootDirPath1, preFileName1, "cos", mChebyshev, EXTFILES);

        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0) {            // no input file exists?
            CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
        } else {
            verb_print(cmd->verbose,
                       "\nstatHistograms: Input file does not exist: %s\n",
                       namebuf1);
            status1 = 0;
            mChebyshev--;
        }
        mChebyshev++;
    } while (status1);

    cmd->mChebyshev = --mChebyshev - 1;
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: found %d multipole files\n",
               cmd->mChebyshev+1);

    //B reading histograms files set 1
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: reading multipole files set 1:");
    mat3 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    CLRM_ext(mat1, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
    }
    //B input monopoles
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        //B input files cos type
        status1 = 1;
        FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                       "namebuf1", "%s/%s_%s_%d%s",
                       rootDirPath1, preFileName1, "cos", m, EXTFILES);

        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0)               // no input file exists?
            OPEN_OUTPUT_OR_FAIL(instr1, namebuf1, "r");
        else {
            verb_print(cmd->verbose,
                    "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf1);
            status1 = 0;
        }
        if (status1 == 0) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: input file does not exist: %s",
                             namebuf1);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts1 != cmd->sizeHistN) {
                COSMO_FAIL_GOTO(cmd, fail,
                           "\ncomputeEdgeCorrections: in file %s: %s %d %d",
                                 namebuf1, "npts, sizeHistN are different:",
                                 npts1, cmd->sizeHistN);
            }
            CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
            //E
        } // ! status
        //E input files cos type
        
        //B input files sin type
        status2 = 1;
        FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                       "namebuf2", "%s/%s_%s_%d%s",
                       rootDirPath1, preFileName1, "sin", m, EXTFILES);

        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf2);
        if (stat(namebuf2, &buf) == 0)               // no input file exists?
            OPEN_OUTPUT_OR_FAIL(instr2, namebuf2, "r");
        else {
            verb_print(cmd->verbose,
                    "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf2);
            status2 = 0;
        }
        if (status2 == 0) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: input file does not exist: %s",
                             namebuf2);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts2 != cmd->sizeHistN) {
                COSMO_FAIL_GOTO(cmd, fail,
                           "\ncomputeEdgeCorrections: in file %s: %s %d %d",
                                 namebuf1, "npts, sizeHistN are different:",
                                 npts2, cmd->sizeHistN);
            }
            CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
            //E
        } // ! status
        //E input files sin type

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                mat3[m][n1][n2] = mat1[n1][n2] + mat2[n1][n2];
            }
        }
    } // ! loop monopoles
    //E input monopoles
    //E reading histograms files set 1

    //B reading histograms files set 2 multipoles of N
    verb_print(cmd->verbose,
               "\n\ncomputeEdgeCorrections: reading multipole files set 2:");
    mat4 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat4[m], cmd->sizeHistN);
    }
    //B input monopoles
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        //B input files cos type
        status1 = 1;
        FORMAT_OR_FAIL(namebuf1, sizeof(namebuf1),
                       "namebuf1", "%s/%s_%s_%d%s",
                       rootDirPath2, preFileName2, "cos_N", m, EXTFILES);

        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0)               // no input file exists?
            OPEN_OUTPUT_OR_FAIL(instr1, namebuf1, "r");
        else {
            verb_print(cmd->verbose,
                    "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf1);
            status1 = 0;
        }
        if (status1 == 0) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: input file does not exist: %s",
                             namebuf1);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts1 != cmd->sizeHistN) {
                COSMO_FAIL_GOTO(cmd, fail,
                           "\ncomputeEdgeCorrections: in file %s: %s %d %d",
                                 namebuf1, "npts, sizeHistN are different:",
                                 npts1, cmd->sizeHistN);
            }
            CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
            //E
        } // ! status
        //E input files cos type
        
        //B input files sin type
        status2 = 1;
        FORMAT_OR_FAIL(namebuf2, sizeof(namebuf2),
                       "namebuf2", "%s/%s_%s_%d%s",
                       rootDirPath2, preFileName2, "sin_N", m, EXTFILES);

        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf2);
        if (stat(namebuf2, &buf) == 0)               // no input file exists?
            OPEN_OUTPUT_OR_FAIL(instr2, namebuf2, "r");
        else {
            verb_print(cmd->verbose,
                    "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf2);
            status2 = 0;
        }
        if (status2 == 0) {
            COSMO_FAIL_GOTO(cmd, fail,
                       "\ncomputeEdgeCorrections: input file does not exist: %s",
                             namebuf2);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts2 != cmd->sizeHistN) {
                COSMO_FAIL_GOTO(cmd, fail,
                           "\ncomputeEdgeCorrections: in file %s: %s %d %d",
                                 namebuf2, "npts, sizeHistN are different:",
                                 npts2, cmd->sizeHistN);
            }
            CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
            //E
        } // ! status
        //E input files sin type

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                mat4[m][n1][n2] = mat1[n1][n2] + mat2[n1][n2];
            }
        }
    } // ! loop monopoles
    //E input monopoles
    //E reading histograms files set 2

    mat5 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat5[m], cmd->sizeHistN);
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            matrixClm(cmd, gd, mat3, mat4, n1, n2, mat5);
            if (cmd->verbose_log>=3)  {
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "\n\nhistZetaM elements again (%d, %d):\n\n",
                               n1, n2);
                for (m=1; m<=cmd->mChebyshev+1; m++) {
                        verb_log_print(cmd->verbose_log, gd->outlog,
                                       "%g\n", mat5[m][n1][n2]);
                }
            }
        }
    }

    verb_print_q(2, cmd->verbose, "\n\n");
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        FORMAT_OR_FAIL(namebuf4, sizeof(namebuf4),
                       "namebuf4", "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                       "_EE", m, EXTFILES);

        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf4);
        OPEN_OUTPUT_OR_FAIL(outstr1, namebuf4, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                WRITE_OUTPUT_OR_FAIL(outstr1, namebuf4, "%16.8e ",
                                     mat5[m][n1][n2]);
            }
            WRITE_OUTPUT_OR_FAIL(outstr1, namebuf4, "\n");
        }
        CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf4);
    }

    //B  and saves matrix ZetaM for each m multipole at a set of theta2 angles
    if (cballs_opt_out_m_histzeta(cmd)) {
        real Zeta;
        real Zeta2;
        real Zeta3;
        real Zeta4;
        real Zeta5;
        int Nbins;
        
        verb_print_q(2, cmd->verbose, "Printing : logscale ... %d %g %g\n",
                     cmd->useLogHist, cmd->rminHist, cmd->rangeN);
        Nbins = cmd->sizeHistN;
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            FORMAT_OR_FAIL(namebuf5, sizeof(namebuf5),
                        "namebuf5", "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                        "_EE", m, EXTFILES);

            verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",
                         namebuf5);
            OPEN_OUTPUT_OR_FAIL(outstr2, namebuf5, "w!");
            WRITE_OUTPUT_OR_FAIL(outstr2, namebuf5, MHISTZETAHEADER);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                if (cmd->useLogHist) {
                gd->deltaR = rlog10(cmd->rangeN/cmd->rminHist)/cmd->sizeHistN;
                    if (cmd->rminHist==0) {
                        rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                        + rlog10(cmd->rangeN);
                    } else {
                        rbinlog = rlog10(cmd->rminHist)
                                    + ((real)(n1)-0.5)*gd->deltaR;
                    }
                    rBin=rpow(10.0,rbinlog);
                } else {
                    gd->deltaR = (cmd->rangeN-cmd->rminHist)/cmd->sizeHistN;
                    rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
                }
                Zeta = mat5[m][n1][n1];
                Zeta2 = mat5[m][n1][(int)(Nbins/4.0)];
                Zeta3 = mat5[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = mat5[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = mat5[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                WRITE_OUTPUT_OR_FAIL(outstr2, namebuf5,
                                    MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            CLOSE_OUTPUT_OR_FAIL(outstr2, namebuf5);
        }
    }
    //E

    free_dmatrix3D(mat5,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat4,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat3,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat2,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(rBins,1,cmd->sizeHistN);

    return SUCCESS;

    fail:
        CLOSE_OUTPUT_OR_FAIL(instr1, namebuf1);
        CLOSE_OUTPUT_OR_FAIL(instr2, namebuf2);
        CLOSE_OUTPUT_OR_FAIL(instr3, namebuf3);
        CLOSE_OUTPUT_OR_FAIL(outstr1, namebuf4);
        CLOSE_OUTPUT_OR_FAIL(outstr2, namebuf5);

        if (mat5 != NULL) free_dmatrix3D(mat5, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat4 != NULL) free_dmatrix3D(mat4, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat3 != NULL) free_dmatrix3D(mat3, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat2 != NULL) free_dmatrix(mat2, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (mat1 != NULL) free_dmatrix(mat1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        if (rBins != NULL) free_dvector(rBins, 1, cmd->sizeHistN);

        return FAILURE;

#undef FORMAT_OR_FAIL

}

#ifdef USEGSL
global int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    double ***mat3, double ***mat4,
                    int n1, int n2, double ***mat5)
{
    // 0 <= l, m < 2*mChebyshev + 1
    // 1 <= n1, n2 <= sizeHistN
    int l, m;
    int lmx;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    real C1;
    int s;

    gsl_vector * bl = gsl_vector_alloc (neqs);
    gsl_matrix * Clm = gsl_matrix_alloc (neqs, neqs);
    gsl_matrix * ClmChk = gsl_matrix_alloc (neqs, neqs);
    gsl_vector *x = gsl_vector_alloc (neqs);
    gsl_permutation * p = gsl_permutation_alloc (neqs);

    gsl_vector *t = gsl_vector_alloc (neqs);
    gsl_matrix * u = gsl_matrix_alloc (neqs, neqs);
    real v;

    for (l=0; l<neqs; l++) {
        gsl_vector_set(bl, l, 0.0);
        for (m=0; m<neqs; m++) {
            gsl_matrix_set(Clm, l, m, 0.0);
            gsl_matrix_set(ClmChk, l, m, 0.0);
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nMatrix and b elements (%d, %d):\n\n", n1, n2);
    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;
//B PIVOTLOOP
        gsl_vector_set(bl, l, mat3[lmx][n1][n2]/mat4[1][n1][n2]);
//E
        if (l<=mx) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), gsl_vector_get(bl, l), l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, gsl_vector_get(bl, l), l, lmx);
        }
        for (m=0; m<neqs; m++) {
            if (l-m>=-mx && l-m<0) {
//B PIVOTLOOP
                C1 = mat4[m-l+1][n1][n2]/mat4[1][n1][n2];
//E
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                //B
                if (cballs_opt_full_sky(cmd)) {
                    if (l!=m) {
                        gsl_matrix_set(Clm, l, m, 0.0);
                        gsl_matrix_set( ClmChk, l, m, 0.0);
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            } // ! l-m >= -mx && l-m < 0
            if (l-m>=0 && l-m<=mx) {
                C1 = mat4[l-m+1][n1][n2]/mat4[1][n1][n2];
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                //B
                if (cballs_opt_full_sky(cmd)) {
                    if (l!=m) {
                        gsl_matrix_set(Clm, l, m, 0.0);
                        gsl_matrix_set( ClmChk, l, m, 0.0);
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            } // ! l-m >= 0 && l-m <= mx
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", gsl_matrix_get(Clm, l, m));
        } // ! loop m
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    } // ! loop l

    gsl_linalg_LU_decomp (Clm, p, &s);
    gsl_linalg_LU_solve (Clm, p, bl, x);
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"x = \n");
        gsl_vector_fprintf (gd->outlog, x, "%g");
    }

    // check A x = b
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=0; l<neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ( gsl_matrix_get(ClmChk,l,m)*gsl_vector_get(x,m) );
            }
            gsl_vector_set(t, l, v);
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",
                           gsl_vector_get(bl,l),gsl_vector_get(t,l));
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nhistZetaM elements:\n\n");
    for (l=0; l<neqs; l++) {
        if (l>=mx) {
            lmx = l+1-mx;
            mat5[lmx][n1][n2] = gsl_vector_get(x, l);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%g\n", mat5[lmx][n1][n2]);
        }
    }

    gsl_matrix_free (u);
    gsl_vector_free (t);
    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_matrix_free (ClmChk);
    gsl_matrix_free (Clm);
    gsl_vector_free (bl);

    return SUCCESS;
}
#else // ! USEGSL
// routine to compute matriz elements Clm and its inverse
//  Be careful with the use of NONORMHISTON and NMultipoles switches:
//  NMultipolesON = 1 and NONORMHISTON = 1
global int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    double ***mat3, double ***mat4,
                    int n1, int n2, double ***mat5)
{
// 1 <= l, m <= 2*mChebyshev + 1
// 1 <= n1, n2 <= sizeHistN
    int l, m;
    int j;
    int *indx;
    real p;
    real **Clm;
    real **ClmChk;
    real **u;
    real *bl;
    real *blChk;
    real *t;
    real C1;

    int lm, lp;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    int lmx;

    Clm = dmatrix(1,neqs,1,neqs);
    ClmChk = dmatrix(1,neqs,1,neqs);
    u = dmatrix(1,neqs,1,neqs);
    bl = dvector(1,neqs);
    blChk = dvector(1,neqs);
    t = dvector(1,neqs);
    indx = ivector(1,neqs);

    CLRM_ext(Clm,neqs);
    CLRM_ext(ClmChk,neqs);
    CLRV_ext(bl,neqs);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nMatrix and b elements (%d, %d):\n\n", n1, n2);
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        bl[l] = mat3[lmx][n1][n2]/mat4[1][n1][n2];
        blChk[l] = bl[l];

        if (l<=mx+1) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), bl[l], l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, bl[l], l, lmx);
        }
        for (m=1; m<=neqs; m++) {
            if (l-m>=-mx && l-m<0) {

                //B check indexs m-l+1 ! m start at 1 and GSL version starts at 0
                C1 = mat4[m-l+1][n1][n2]/mat4[1][n1][n2];

                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];

                //B
                if (cballs_opt_full_sky(cmd)) {
                    if (l!=m) {
                        Clm[l][m] = 0.0;
                        ClmChk[l][m] = 0.0;
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = mat4[l-m+1][n1][n2]/mat4[1][n1][n2];

                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];

                //B
                if (cballs_opt_full_sky(cmd)) {
                    if (l!=m) {
                        Clm[l][m] = 0.0;
                        ClmChk[l][m] = 0.0;
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", Clm[l][m]);
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    ludcmp(Clm,neqs,indx,&p);
    lubksb(Clm,neqs,indx,bl);

    // vector solutions
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nVector solution:\n");
        for (l=1;l<=neqs;l++) {
            verb_log_print(cmd->verbose_log, gd->outlog,"%8s %g\n"," ", bl[l]);
        }
    }

    // check A x = b
    real v;
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=1; l<=neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ClmChk[l][m]*bl[m];
            }
            t[l] = v;
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",blChk[l],t[l]);
        }
    }

//B correction 2025-04-06
    for (l=1; l<=neqs; l++) {
        if (l>=mx+2) {
            lmx = (l-2)+2-mx;
            mat5[lmx][n1][n2] = bl[l];
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%g\n", mat5[lmx][n1][n2]);
        }

    }
//E correction 2025-04-06

    free_ivector(indx,1,cmd->mChebyshev+1);
    free_dvector(t,1,cmd->mChebyshev+1);
    free_dvector(blChk,1,cmd->mChebyshev+1);
    free_dvector(bl,1,cmd->mChebyshev+1);
    free_dmatrix(u,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(ClmChk,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(Clm,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);

    return SUCCESS;
}

#endif // ! USEGSL

#undef MHISTZETAHEADER
#undef MHISTZETA
#undef MHISTZETAHEADERSTDDEV
#undef MHISTZETASTDDEV

//E

//B socket:
#ifdef ADDONS
#include "cballsutils_include_02.h"
#endif
//E
