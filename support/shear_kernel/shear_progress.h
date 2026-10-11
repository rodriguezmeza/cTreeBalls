/* Task completion telemetry, separate from scientific accumulation. */
#ifndef CTREEBALLS_SHEAR_PROGRESS_H
#define CTREEBALLS_SHEAR_PROGRESS_H

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    bool enabled;
    int rank;
    INTEGER total;
    INTEGER completed;
    double started;
} shear_task_progress;

static void shear_print_task_progress(const shear_task_progress *progress)
{
    verb_print_min_info(
        progress->cmd->verbose, progress->cmd->verbose_log,
        progress->gd->outlog,
        SHEAR_ENGINE_NAME ": 3PCF progress: rank %d completed tasks %"
        INTEGER_FMT " / %" INTEGER_FMT " (%.1f%%); elapsed %.2f s\n",
        progress->rank, progress->completed, progress->total,
        100.0*(double)progress->completed/(double)progress->total,
        shear_profile_wall_time() - progress->started);
}

static void shear_complete_progress_task(shear_task_progress *progress)
{
    if (!progress->enabled) return;
    /* One update per coarse frontier task, never per neighbor or coefficient.
     * Keep count and flushed output together so workers cannot reorder lines.
     * This counts completed tasks, not elapsed-time fraction or pivot bodies. */
#pragma omp critical(cballs_shear_task_progress)
    {
        progress->completed++;
        shear_print_task_progress(progress);
    }
}

#endif
