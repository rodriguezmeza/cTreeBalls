/* Scalar angular window deconvolution shared by all two-ball traversals.
 * W_l = sum w_i w_j w_k exp(i*l*(phi_ij-phi_ik)), with distinct j,k.
 * Solve sum_n W_(l-n)/W_0 * zeta_n = S_l/W_0 for -M <= l,n <= M.
 * Window orders through 2*M are required, including their imaginary parts.
 */
#ifndef CBALLS_TREECORR_EDGE_CORRECTION_H
#define CBALLS_TREECORR_EDGE_CORRECTION_H

#include <complex.h>

static bool treecorr_triple_values(size_t stride, int orders, int window_orders,
                                   size_t *values)
{
    size_t planes;
    if (!stride || orders <= 0 || window_orders < 0
        || stride > SIZE_MAX / stride)
        return FALSE;
    if ((size_t)orders > (SIZE_MAX - 1) / TREECORR_ZETA_COMPONENTS)
        return FALSE;
    planes = TREECORR_ZETA_COMPONENTS * (size_t)orders + 1;
    if ((size_t)window_orders > (SIZE_MAX - planes) / 2) return FALSE;
    planes += 2 * (size_t)window_orders;
    if (planes > SIZE_MAX / (stride * stride)) return FALSE;
    *values = planes * stride * stride;
    return TRUE;
}

static int treecorr_write_edge_matrix(
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
                 "%s: edge output path too long", TREECORR_METHOD_NAME);
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

/* Scaled complex LU with partial pivoting. Singular windows have no unique
 * truncated estimate; match the GGG policy by returning zero for those bins.
 */
static bool treecorr_edge_solve(double complex *a, double complex *rhs, int n)
{
    const double tolerance = 128.0 * DBL_EPSILON * n;
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double largest = cabs(a[(size_t)col * n + col]);
        for (int row = col + 1; row < n; row++) {
            const double value = cabs(a[(size_t)row * n + col]);
            if (value > largest) {
                pivot = row;
                largest = value;
            }
        }
        if (!isfinite(largest) || largest <= tolerance) return FALSE;
        if (pivot != col) {
            for (int j = col; j < n; j++) {
                const double complex swap = a[(size_t)col * n + j];
                a[(size_t)col * n + j] = a[(size_t)pivot * n + j];
                a[(size_t)pivot * n + j] = swap;
            }
            const double complex swap = rhs[col];
            rhs[col] = rhs[pivot];
            rhs[pivot] = swap;
        }
        for (int row = col + 1; row < n; row++) {
            const double complex factor = a[(size_t)row * n + col]
                                        / a[(size_t)col * n + col];
            for (int j = col + 1; j < n; j++)
                a[(size_t)row * n + j] -= factor * a[(size_t)col * n + j];
            rhs[row] -= factor * rhs[col];
        }
    }
    for (int row = n - 1; row >= 0; row--) {
        for (int j = row + 1; j < n; j++)
            rhs[row] -= a[(size_t)row * n + j] * rhs[j];
        rhs[row] /= a[(size_t)row * n + row];
        if (!isfinite(creal(rhs[row])) || !isfinite(cimag(rhs[row])))
            return FALSE;
    }
    return TRUE;
}

static int treecorr_publish_edge(
        struct cmdline_data *cmd, struct global_data *gd,
        const real *tasks, INTEGER task_count, size_t values_per_task,
        size_t stride, int orders)
{
    const int window_orders = treecorr_window_orders(cmd);
    const int mmax = orders - 1;
    const size_t plane = stride * stride;
    size_t window_values;
    real *window = NULL;
    double complex *matrix = NULL, *rhs = NULL;
    int singular = 0, empty = 0;
    int status = FAILURE;

    if (!cballs_opt_edge_corrections(cmd)) return SUCCESS;
    if (window_orders <= 0 || gd->histZetaM_EE == NULL
        || gd->histZetaM_EE_Im == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge correction requires allocated 3PCF histograms",
                 TREECORR_METHOD_NAME);
        return FAILURE;
    }
    if ((size_t)window_orders > SIZE_MAX / 2 / plane / sizeof(real)
        || (size_t)window_orders > SIZE_MAX / (size_t)window_orders
                                                   / sizeof(*matrix)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge workspace size overflow", TREECORR_METHOD_NAME);
        return FAILURE;
    }
    window_values = 2 * (size_t)window_orders * plane;
    window = calloc(window_values, sizeof(*window));
    matrix = calloc((size_t)window_orders * (size_t)window_orders, sizeof(*matrix));
    rhs = calloc((size_t)window_orders, sizeof(*rhs));
    if (!window || !matrix || !rhs) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: edge workspace allocation failed", TREECORR_METHOD_NAME);
        goto cleanup;
    }

    /* Keep the same ascending task order as the signal reduction. */
    for (INTEGER task = 0; task < task_count; task++) {
        const real *source = tasks + (size_t)task * values_per_task
            + (TREECORR_ZETA_COMPONENTS * (size_t)orders + 1) * plane;
        for (size_t i = 0; i < window_values; i++) window[i] += source[i];
    }
    for (size_t i = 0; i < window_values; i++) {
        if (!isfinite(window[i])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: nonfinite angular window", TREECORR_METHOD_NAME);
            goto cleanup;
        }
    }

    for (int i = 1; i <= cmd->sizeHistN; i++) {
        for (int j = 1; j <= cmd->sizeHistN; j++) {
            const size_t bin = (size_t)i * stride + (size_t)j;
            const double wzero = window[bin];
            for (int m = 1; m <= orders; m++) {
                gd->histZetaM_EE[m][i][j] = 0.0;
                gd->histZetaM_EE_Im[m][i][j] = 0.0;
            }
            if (!(wzero > 0.0)) {
                empty++;
                continue;
            }
            for (int row = 0; row < window_orders; row++) {
                const int ell = row - mmax;
                const int m = abs(ell) + 1;
                const double re = gd->histZetaMcos[m][i][j]
                                + gd->histZetaMsin[m][i][j];
                const double im = gd->histZetaMsincos[m][i][j]
                                - gd->histZetaMcossin[m][i][j];
                if (!isfinite(re) || !isfinite(im)) {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "%s: nonfinite signal multipole", TREECORR_METHOD_NAME);
                    goto cleanup;
                }
                rhs[row] = (re + I * (ell < 0 ? -im : im)) / wzero;
                for (int col = 0; col < window_orders; col++) {
                    const int difference = row - col;
                    const size_t wi = (size_t)abs(difference) * plane + bin;
                    const double wr = window[wi];
                    const double wm = window[(size_t)window_orders * plane + wi];
                    matrix[(size_t)row * window_orders + col] =
                        (wr + I * (difference < 0 ? -wm : wm)) / wzero;
                }
            }
            if (!treecorr_edge_solve(matrix, rhs, window_orders)) {
                singular++;
                continue;
            }
            for (int m = 0; m < orders; m++) {
                gd->histZetaM_EE[m + 1][i][j] = creal(rhs[mmax + m]);
                gd->histZetaM_EE_Im[m + 1][i][j] = cimag(rhs[mmax + m]);
            }
        }
    }
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "%s: edge correction uses window modes 0..%d; "
        "%d empty and %d singular radial-bin pairs set to zero\n",
        TREECORR_METHOD_NAME, window_orders - 1, empty, singular);
    if (!scanopt(cmd->options, "no-out-Hist")) {
        for (int order = 0; order < window_orders; order++) {
            if (treecorr_write_edge_matrix(cmd, gd, "window_Re", order,
                    window + (size_t)order * plane, NULL, stride) == FAILURE
                || treecorr_write_edge_matrix(cmd, gd, "window_Im", order,
                    window + ((size_t)window_orders + order) * plane,
                    NULL, stride) == FAILURE)
                goto cleanup;
        }
        for (int order = 0; order < orders; order++) {
            if (treecorr_write_edge_matrix(cmd, gd, "EE", order, NULL,
                    gd->histZetaM_EE[order + 1], stride) == FAILURE
                || treecorr_write_edge_matrix(cmd, gd, "EE_Im", order, NULL,
                    gd->histZetaM_EE_Im[order + 1], stride) == FAILURE)
                goto cleanup;
        }
    }
    status = SUCCESS;
cleanup:
    free(rhs);
    free(matrix);
    free(window);
    return status;
}

#endif
