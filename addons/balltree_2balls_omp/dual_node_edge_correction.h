/* Scalar angular window deconvolution shared by all two-ball traversals.
 * W_l = sum w_i w_j w_k exp(i*l*(phi_ij-phi_ik)), with distinct j,k.
 * Solve sum_n W_(l-n)/W_0 * zeta_n = S_l/W_0 for -M <= l,n <= M.
 * Window orders through 2*M are required, including their imaginary parts.
 */
#ifndef CBALLS_DUAL_NODE_EDGE_CORRECTION_H
#define CBALLS_DUAL_NODE_EDGE_CORRECTION_H

#include "scalar_window_solver.h"

static inline double dual_node_edge_timer_now(void)
{
#ifdef OPENMPCODE
    return omp_get_wtime();
#else
    return CPUTIME;
#endif
}

static bool dual_node_triple_values(size_t stride, int orders, int window_orders,
                                   size_t *values)
{
    size_t planes;
    if (!stride || orders <= 0 || window_orders < 0
        || stride > SIZE_MAX / stride)
        return FALSE;
    if ((size_t)orders > (SIZE_MAX - 1) / DUAL_NODE_ZETA_COMPONENTS)
        return FALSE;
    planes = DUAL_NODE_ZETA_COMPONENTS * (size_t)orders + 1;
    if ((size_t)window_orders > (SIZE_MAX - planes) / 2) return FALSE;
    planes += 2 * (size_t)window_orders;
    if (planes > SIZE_MAX / (stride * stride)) return FALSE;
    *values = planes * stride * stride;
    return TRUE;
}

static int dual_node_write_edge_matrix(
        struct cmdline_data *cmd, struct global_data *gd,
        const char *suffix, int order, const real *flat, real **matrix,
        size_t stride)
{
    string routineName = "two-ball edge output";
    char path[MAXLENGTHOFFILES + 80];
    stream output = NULL;
    const int length = snprintf(path, sizeof(path), "%s_%s_%d%s",
        gd->fpfnamehistZetaMFileName, suffix, order + 1, EXTFILES);
    if (length < 0 || (size_t)length >= sizeof(path)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge output path too long", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    OPEN_OUTPUT_OR_FAIL(output, path, "w!");
    for (int i = 1; i <= cmd->sizeHistN; i++) {
        for (int j = 1; j <= cmd->sizeHistN; j++)
            WRITE_OUTPUT_OR_FAIL(output, path, "%.17g ",
                flat ? flat[(size_t)i * stride + (size_t)j] : matrix[i][j]);
        WRITE_OUTPUT_OR_FAIL(output, path, "\n");
    }
    CLOSE_OUTPUT_OR_FAIL(output, path);
    return SUCCESS;
}

static int dual_node_publish_edge(
        struct cmdline_data *cmd, struct global_data *gd,
        const real *tasks, INTEGER task_count, size_t values_per_task,
        size_t stride, int orders)
{
    const int window_orders = dual_node_window_orders(cmd);
    const int mmax = orders - 1;
    const size_t plane = stride * stride;
    size_t window_values;
    real *window = NULL;
    double complex *matrix = NULL, *rhs = NULL, *modes = NULL;
    double solve_cpu_started = 0.0, solve_wall_started = 0.0;
    int singular = 0, empty = 0;
    int solve_timing_active = FALSE;
    int status = FAILURE;

    if (!cballs_opt_edge_corrections(cmd)) return SUCCESS;
    if (window_orders <= 0 || gd->histZetaM_EE == NULL
        || gd->histZetaM_EE_Im == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge correction requires allocated 3PCF histograms",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if ((size_t)window_orders > SIZE_MAX / 2 / plane / sizeof(real)
        || (size_t)window_orders > SIZE_MAX / (size_t)window_orders
                                                   / sizeof(*matrix)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge workspace size overflow", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (cballs_scalar_window_begin(cmd, gd) == FAILURE) return FAILURE;
    window_values = 2 * (size_t)window_orders * plane;
    window = calloc(window_values, sizeof(*window));
    matrix = calloc((size_t)window_orders * (size_t)window_orders, sizeof(*matrix));
    rhs = calloc((size_t)window_orders, sizeof(*rhs));
    modes = calloc((size_t)window_orders + orders, sizeof(*modes));
    if (!window || !matrix || !rhs || !modes) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge workspace allocation failed", DUAL_NODE_METHOD_NAME);
        goto cleanup;
    }

    /* Keep the same ascending task order as the signal reduction. */
    for (INTEGER task = 0; task < task_count; task++) {
        const real *source = tasks + (size_t)task * values_per_task
            + (DUAL_NODE_ZETA_COMPONENTS * (size_t)orders + 1) * plane;
        for (size_t i = 0; i < window_values; i++) window[i] += source[i];
    }
    solve_cpu_started = CPUTIME;
    solve_wall_started = dual_node_edge_timer_now();
    solve_timing_active = TRUE;
    for (int i = 1; i <= cmd->sizeHistN; i++) {
        for (int j = 1; j <= cmd->sizeHistN; j++) {
            const size_t bin = (size_t)i * stride + (size_t)j;
            const double wzero = window[bin];
            const size_t diagnostic_bin = (size_t)(i-1)*cmd->sizeHistN + j-1;
            gd->scalar_window_w0[diagnostic_bin] = wzero;
            gd->scalar_window_status[diagnostic_bin] = CBALLS_WINDOW_EMPTY;
            for (int m = 1; m <= orders; m++) {
                gd->histZetaM_EE[m][i][j] = NAN;
                gd->histZetaM_EE_Im[m][i][j] = NAN;
            }
            for (int m=0; m<orders; m++) {
                const double re = gd->histZetaMcos[m+1][i][j]
                                + gd->histZetaMsin[m+1][i][j];
                const double im = gd->histZetaMsincos[m+1][i][j]
                                - gd->histZetaMcossin[m+1][i][j];
                modes[m] = cballs_scalar_complex(re, im);
            }
            for (int m=0; m<window_orders; m++)
                modes[orders+m] = cballs_scalar_complex(window[(size_t)m*plane+bin],
                    window[((size_t)window_orders+m)*plane+bin]);
            double pivot_ratio;
            const int solve_status = cballs_scalar_edge_deconvolve(
                modes, modes+orders, mmax, matrix, rhs, &pivot_ratio);
            gd->scalar_window_status[diagnostic_bin] = solve_status;
            gd->scalar_window_pivot_ratio[diagnostic_bin] = pivot_ratio;
            if (solve_status != CBALLS_WINDOW_VALID) {
                if (solve_status == CBALLS_WINDOW_EMPTY) empty++;
                else singular++;
                continue;
            }
            for (int m = 0; m < orders; m++) {
                gd->histZetaM_EE[m + 1][i][j] = creal(rhs[mmax + m]);
                gd->histZetaM_EE_Im[m + 1][i][j] = cimag(rhs[mmax + m]);
            }
        }
    }
    gd->cpu_edge_correction += CPUTIME - solve_cpu_started;
    gd->wall_edge_correction += dual_node_edge_timer_now() - solve_wall_started;
    solve_timing_active = FALSE;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: edge correction uses window modes 0..%d; "
        "%d empty and %d rejected radial-bin pairs set to NaN\n",
        DUAL_NODE_METHOD_NAME, window_orders - 1, empty, singular);
    if (!scanopt(cmd->options, "no-out-Hist")) {
        if (cballs_scalar_edge_manifest(cmd, gd, FALSE) == FAILURE) goto cleanup;
        for (int order = 0; order < window_orders; order++) {
            if (dual_node_write_edge_matrix(cmd, gd, "window_Re", order,
                    window + (size_t)order * plane, NULL, stride) == FAILURE
                || dual_node_write_edge_matrix(cmd, gd, "window_Im", order,
                    window + ((size_t)window_orders + order) * plane,
                    NULL, stride) == FAILURE)
                goto cleanup;
        }
        for (int order = 0; order < orders; order++) {
            if (dual_node_write_edge_matrix(cmd, gd, "edge_cos", order, NULL,
                    gd->histZetaMcos[order+1], stride) == FAILURE
                || dual_node_write_edge_matrix(cmd, gd, "edge_sin", order, NULL,
                    gd->histZetaMsin[order+1], stride) == FAILURE
                || dual_node_write_edge_matrix(cmd, gd, "edge_sincos", order, NULL,
                    gd->histZetaMsincos[order+1], stride) == FAILURE
                || dual_node_write_edge_matrix(cmd, gd, "edge_cossin", order, NULL,
                    gd->histZetaMcossin[order+1], stride) == FAILURE) goto cleanup;
            if (dual_node_write_edge_matrix(cmd, gd, "EE", order, NULL,
                    gd->histZetaM_EE[order + 1], stride) == FAILURE
                || dual_node_write_edge_matrix(cmd, gd, "EE_Im", order, NULL,
                    gd->histZetaM_EE_Im[order + 1], stride) == FAILURE)
                goto cleanup;
        }
    }
    if (cballs_scalar_window_write(cmd, gd) == FAILURE) goto cleanup;
    if (!cballs_opt_no_out_hist(cmd)
        && cballs_scalar_edge_manifest(cmd, gd, TRUE) == FAILURE) goto cleanup;
    gd->scalar_window_ready = TRUE;
    status = SUCCESS;
cleanup:
    if (solve_timing_active) {
        gd->cpu_edge_correction += CPUTIME - solve_cpu_started;
        gd->wall_edge_correction += dual_node_edge_timer_now() - solve_wall_started;
    }
    free(modes);
    free(rhs);
    free(matrix);
    free(window);
    return status;
}

#endif
