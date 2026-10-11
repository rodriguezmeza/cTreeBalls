/* Coarse task telemetry; histogram arithmetic and scheduling are unchanged. */
#ifndef CBALLS_LYA_PROGRESS_H
#define CBALLS_LYA_PROGRESS_H

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    const char *phase;
    int enabled;
    size_t rank, total, completed, reported, interval;
    double started, last_report;
} lya_progress;

static double lya_progress_wall_time(void)
{
#ifdef OPENMPCODE
    return omp_get_wtime();
#else
    return CPUTIME;
#endif
}

static size_t lya_progress_local_tasks(size_t total, size_t first, size_t stride)
{
    return first < total ? 1 + (total - 1 - first) / stride : 0;
}

static void lya_progress_print(lya_progress *p, double now)
{
    verb_print_min_info(p->cmd->verbose, p->cmd->verbose_log, p->gd->outlog,
        "%s: %s progress: rank %zu completed tasks %zu / %zu (%.1f%%); elapsed %.2f s\n",
        p->cmd->searchMethod, p->phase, p->rank, p->completed, p->total,
        p->total ? 100.0*(double)p->completed/(double)p->total : 100.0,
        now - p->started);
    p->reported = p->completed;
    p->last_report = now;
}

/* Initialize on the calling thread, before OpenMP; MPI state is read only here.
 * The normal MPI logger selects rank zero. Its denominator is rank-local. */
static void lya_progress_init(lya_progress *p, struct cmdline_data *cmd,
                              struct global_data *gd, const char *phase,
                              size_t local_tasks)
{
    *p = (lya_progress){0};
    p->cmd = cmd; p->gd = gd; p->phase = phase;
    p->rank = lya_parallel_first(cmd);
    p->total = local_tasks;
    p->interval = local_tasks / 100 + (local_tasks % 100 != 0);
    p->enabled = cmd->verbose >= VERBOSEMININFO
        || (cmd->verbose_log >= VERBOSEMININFO && gd->outlog != NULL);
    if (p->enabled) {
        p->started = lya_progress_wall_time();
        lya_progress_print(p, p->started);
    }
}

/* Called only after a successful commit in the existing omp ordered region.
 * No added locks, MPI collectives, or per-pair/per-triplet instrumentation.
 * Print about every 1%, or after five seconds at the next task boundary. */
static void lya_progress_step(lya_progress *p)
{
    if (!p->enabled) return;
    p->completed++;
    double now = lya_progress_wall_time();
    if (p->completed == p->total || p->completed - p->reported >= p->interval
        || now - p->last_report >= 5.0)
        lya_progress_print(p, now);
}

static void lya_progress_finish(const lya_progress *p, int success)
{
    if (p->enabled && success && p->completed == p->total)
        verb_print_min_info(p->cmd->verbose, p->cmd->verbose_log, p->gd->outlog,
            "%s: %s counting complete on rank %zu; continuing reductions/output\n",
            p->cmd->searchMethod, p->phase, p->rank);
}

#endif
