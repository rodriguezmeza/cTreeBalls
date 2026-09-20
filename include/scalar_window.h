#ifndef CBALLS_SCALAR_WINDOW_H
#define CBALLS_SCALAR_WINDOW_H

#include <errno.h>
#include <stdint.h>

/* One status per (radial bin 1, radial bin 2), shared by all solved modes. */
enum cballs_scalar_window_status {
    CBALLS_WINDOW_NOT_COMPUTED = 0,
    CBALLS_WINDOW_VALID = 1,
    CBALLS_WINDOW_EMPTY = 2,
    CBALLS_WINDOW_SINGULAR = 3,
    CBALLS_WINDOW_NONFINITE = 4
};

static inline void cballs_scalar_window_free(struct global_data *gd)
{
    free(gd->scalar_window_status);
    free(gd->scalar_window_w0);
    free(gd->scalar_window_pivot_ratio);
    gd->scalar_window_status = NULL;
    gd->scalar_window_w0 = gd->scalar_window_pivot_ratio = NULL;
    gd->scalar_window_bins = 0;
    gd->scalar_window_ready = FALSE;
}

static inline int cballs_scalar_window_begin(struct cmdline_data *cmd,
                                              struct global_data *gd)
{
    size_t bins = (size_t)cmd->sizeHistN, count;
    cballs_scalar_window_free(gd);
    if (cmd->sizeHistN <= 0 || bins > SIZE_MAX / bins) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "scalar window: invalid radial dimensions");
        return FAILURE;
    }
    count = bins * bins;
    if (cballs_calloc_checked((void **)&gd->scalar_window_status, count,
            sizeof(*gd->scalar_window_status), "scalar window status",
            cmd->error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_malloc_checked((void **)&gd->scalar_window_w0, count,
            sizeof(real), "scalar window monopole",
            cmd->error_message, _ERRORMSGSIZE_) == FAILURE
        || cballs_malloc_checked((void **)&gd->scalar_window_pivot_ratio, count,
            sizeof(real), "scalar window pivot ratio",
            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        cballs_scalar_window_free(gd);
        return FAILURE;
    }
    for (size_t i = 0; i < count; i++) {
        gd->scalar_window_w0[i] = NAN;
        gd->scalar_window_pivot_ratio[i] = NAN;
    }
    gd->scalar_window_bins = cmd->sizeHistN;
    return SUCCESS;
}

/* The ratio of smallest to largest accepted elimination pivots is a proxy,
   not a matrix condition number. An unavailable ratio remains NaN. */
static inline int cballs_scalar_window_write(struct cmdline_data *cmd,
                                              struct global_data *gd)
{
    string routineName = "scalar window diagnostics";
    char path[MAXLENGTHOFFILES + 80];
    stream output = NULL;
    if (cballs_opt_no_out_hist(cmd)) return SUCCESS;
    int length = snprintf(path, sizeof(path), "%s_window_diagnostics%s",
                          gd->fpfnamehistZetaMFileName, EXTFILES);
    if (length < 0 || (size_t)length >= sizeof(path)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "scalar window: diagnostics path too long");
        return FAILURE;
    }
    OPEN_OUTPUT_OR_FAIL(output, path, "w!");
    WRITE_OUTPUT_OR_FAIL(output, path,
        "# bin1 bin2 status window_monopole pivot_ratio\n"
        "# status: 1=valid 2=empty/nonpositive 3=singular 4=nonfinite\n");
    for (int i = 0; i < gd->scalar_window_bins; i++) {
        for (int j = 0; j < gd->scalar_window_bins; j++) {
            size_t bin = (size_t)i * gd->scalar_window_bins + j;
            WRITE_OUTPUT_OR_FAIL(output, path, "%d %d %u %.17g %.17g\n",
                i+1, j+1, (unsigned)gd->scalar_window_status[bin],
                (double)gd->scalar_window_w0[bin],
                (double)gd->scalar_window_pivot_ratio[bin]);
        }
    }
    CLOSE_OUTPUT_OR_FAIL(output, path);
    return SUCCESS;
}

#endif
