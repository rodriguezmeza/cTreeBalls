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
#include "resource_contracts.h"
#include "../support/pca_tree/fcfc_balltree.h"

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
#ifndef USEGSL
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
#ifndef USEGSL
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

global cballs_runtime_state *cballs_runtime_current(void)
{
    return cballs_active_runtime != NULL ? cballs_active_runtime : &cballs_default_runtime;
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
    if (callback != cballs_call_end_run_free_memory)
        status = cballs_mpi_context_consensus(cmd, status,
            callback == cballs_call_print_parameter_file ? "MPI parameter output" :
            callback == cballs_call_set_number_threads ? "MPI thread configuration" :
            callback == cballs_call_start_run_common ? "MPI startup completion" :
            callback == cballs_call_main_loop ? "MPI computation completion" : "MPI run finalization");
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
    int status = cballs_guard_stage(cmd, gd, NULL, cballs_call_main_loop);
    if (status == FAILURE && gd != NULL) {
        gd->histogram_results_ready = FALSE;
        gd->histogram_products = 0;
    }
    return status;
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

global void cballs_runtime_destroy(cballs_runtime_state *state)
{
    if (state == NULL)
        return;

    if (state == cballs_active_runtime) {
        cballs_runtime_store(state);
        cballs_runtime_load(&cballs_default_runtime);
        cballs_active_runtime = NULL;
    }
    cballs_runtime_clear_results(state);
    cballs_runtime_clear_caches(state);
    free(state);
}


/* Results and backend attachments are object-owned. Legacy activation remains
 * serialized; this does not make concurrent direct native entry safe. */
void *cballs_runtime_attach(int slot, size_t size, void (*destroy)(void *), size_t (*bytes)(const void *))
{
    if(slot<0 || slot>=2 || !size || !destroy || !bytes)
        cballs_resource_failure("invalid runtime attachment");
    cballs_runtime_attachment *a=&cballs_runtime_current()->attachments[slot];
    if (!a->data) {
        a->data=calloc(1,size);
        if (!a->data) cballs_allocation_failure(size,"runtime attachment");
        a->destroy=destroy; a->bytes=bytes;
    }
    return a->data;
}
size_t cballs_runtime_cache_bytes(const cballs_runtime_state *state)
{
    size_t total=0;
    for (int i=0;i<2;i++) if (state->attachments[i].data)
        total+=state->attachments[i].bytes(state->attachments[i].data);
    return total;
}
void cballs_runtime_clear_caches(cballs_runtime_state *state)
{
    for (int i=0;i<2;i++) {
        if (state->attachments[i].data) state->attachments[i].destroy(state->attachments[i].data);
        memset(&state->attachments[i],0,sizeof(state->attachments[i]));
    }
}
void cballs_runtime_clear_results(cballs_runtime_state *state)
{
    for (int i=0;i<8;i++) { free(state->results[i].data); memset(&state->results[i],0,sizeof(state->results[i])); }
}
size_t cballs_runtime_result_bytes(const cballs_runtime_state *state)
{
    size_t total=0; for(int i=0;i<8;i++) total+=state->results[i].bytes; return total;
}
void cballs_runtime_retain_results(cballs_runtime_state *state,int enabled) {state->retain_results=enabled;}
void cballs_result_adopt(const char *name,void **data,int kind,int rank,const size_t *shape)
{
    cballs_runtime_state *state=cballs_runtime_current();
    if (!state->retain_results || !*data) return;
    cballs_result_array *a=NULL;
    for(int i=0;i<8;i++) if(!state->results[i].data) {a=&state->results[i];break;}
    if(!a || rank<1 || rank>5 || kind<0 || kind>2)
        cballs_resource_failure("invalid native result ownership request");
    size_t count=1,bytes;
    for(int i=0;i<rank;i++) if(!cballs_size_mul(count,shape[i],&count))
        cballs_resource_failure("native result dimensions overflow");
    if(!cballs_size_mul(count,kind==1?sizeof(long double):kind==2?sizeof(uint64_t):sizeof(REAL),&bytes))
        cballs_resource_failure("native result size overflow");
    snprintf(a->name,sizeof(a->name),"%s",name);a->kind=kind;a->rank=rank;a->count=count;a->bytes=bytes;
    memcpy(a->shape,shape,rank*sizeof(size_t));a->data=*data;*data=NULL;
}
const cballs_result_array *cballs_result_at(const cballs_runtime_state *state,int index)
{return index>=0 && index<8 && state->results[index].data ? &state->results[index] : NULL;}
double cballs_result_value(const cballs_result_array *a,size_t i)
{return a->kind==1?(double)((long double*)a->data)[i]:a->kind==2?(double)((uint64_t*)a->data)[i]:(double)((REAL*)a->data)[i];}
uint64_t cballs_result_count_value(const cballs_result_array *a,size_t i) {return ((uint64_t*)a->data)[i];}
size_t cballs_resource_type_size(int kind)
{
    switch(kind) {case 0:return sizeof(body);case 1:return sizeof(cell);case 2:return sizeof(fcfc_ballnode);
    case 3:return sizeof(fcfc_ballpoint);case 4:return sizeof(REAL);case 5:return sizeof(size_t);case 6:return sizeof(long double);default:return 0;}
}
size_t cballs_resource_default_mib(void) {return CBALLS_DEFAULT_MEMORY_BUDGET_MIB;}
int cballs_resource_base(struct cmdline_data *cmd,struct global_data *gd,size_t *out)
{
    size_t total=gd->common_histogram_bytes,bytes;
    for(int i=0;i<gd->ninfiles;i++) {
        if(gd->nbodyTable[i]<0 || !cballs_size_mul((size_t)gd->nbodyTable[i],sizeof(body),&bytes)
            || !cballs_size_add(total,bytes,&total)) goto invalid;
    }
    if(!cballs_size_add(total,cballs_runtime_cache_bytes(cballs_runtime_current()),&total)) goto invalid;
    *out=total;
    return cballs_memory_preflight(total,"catalogs, common histograms and retained caches per rank",cmd->error_message,sizeof(cmd->error_message));
invalid:
    snprintf(cmd->error_message,sizeof(cmd->error_message),"aggregate resource dimension overflow");return FAILURE;
}
