#ifndef CBALLS_SCALAR_EDGE_CORRECTION_H
#define CBALLS_SCALAR_EDGE_CORRECTION_H
/* Private numerical implementation: include after the engine types.
 * Arithmetic and iteration order are preserved; see engineering_contracts.rst. */
#if defined(NMultipoles) && defined(NONORMHIST) \
    && defined(THREEPCFCONVERGENCE)
typedef struct {
    real re;
    real im;
} ggg_complex;

local ggg_complex ggg_complex_make(real re, real im)
{
    ggg_complex value = {re, im};
    return value;
}

local ggg_complex ggg_complex_sub(ggg_complex a, ggg_complex b)
{
    return ggg_complex_make(a.re-b.re, a.im-b.im);
}

local ggg_complex ggg_complex_mul(ggg_complex a, ggg_complex b)
{
    return ggg_complex_make(a.re*b.re-a.im*b.im,
                            a.re*b.im+a.im*b.re);
}

local real ggg_complex_abs2(ggg_complex value)
{
    return value.re*value.re + value.im*value.im;
}

local ggg_complex ggg_complex_div(ggg_complex numerator,
                                  ggg_complex denominator)
{
    real norm = ggg_complex_abs2(denominator);
    return ggg_complex_make(
        (numerator.re*denominator.re + numerator.im*denominator.im)/norm,
        (numerator.im*denominator.re - numerator.re*denominator.im)/norm);
}

local ggg_complex ggg_histogram_mode(real ***coscos, real ***sinsin,
                                     real ***sincos, real ***cossin,
                                     int order, int n1, int n2)
{
    int index = abs(order) + 1;
    ggg_complex value = ggg_complex_make(
        coscos[index][n1][n2] + sinsin[index][n1][n2],
        sincos[index][n1][n2] - cossin[index][n1][n2]);
    if (order < 0)
        value.im = -value.im;
    return value;
}

local int compute_edge_corrections_ggg(struct cmdline_data *cmd,
                                       struct global_data *gd,
                                       gdlptr_sincos_omp_ggg gdl,
                                       gdlptr_sincos_omp_ggg_N gdlN)
{
    int nmax = cmd->mChebyshev;
    int multipoles = 2*nmax + 1;
    size_t matrix_count;
    ggg_complex *matrix = NULL;
    ggg_complex *rhs = NULL;
    int n1, n2;
    int singular_bins = 0;

    if ((size_t)multipoles > SIZE_MAX/(size_t)multipoles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-omp: edge-correction workspace overflow");
        return FAILURE;
    }
    if (cballs_scalar_window_begin(cmd, gd) == FAILURE) return FAILURE;
    matrix_count = (size_t)multipoles*(size_t)multipoles;
    matrix = calloc(matrix_count, sizeof(*matrix));
    rhs = calloc((size_t)multipoles, sizeof(*rhs));
    if (matrix == NULL || rhs == NULL) {
        free(matrix);
        free(rhs);
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-ggg-omp: cannot allocate edge-correction workspace");
        return FAILURE;
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            ggg_complex nzero = ggg_histogram_mode(
                gdlN->histZetaMcos, gdlN->histZetaMsin,
                gdlN->histZetaMsincos, gdlN->histZetaMcossin,
                0, n1, n2);
            const size_t bin = (size_t)(n1-1)*cmd->sizeHistN + n2-1;
            real smallest_pivot = DBL_MAX, largest_pivot = 0.0;
            gd->scalar_window_w0[bin] = nzero.re;
            gd->scalar_window_status[bin] = CBALLS_WINDOW_NONFINITE;
            real tolerance;
            real matrix_scale = 0.0;
            bool singular = FALSE;
            int row, column;

            for (row=1; row<=nmax+1; row++) {
                gdl->histZetaM_EE[row][n1][n2] = NAN;
                gdl->histZetaM_EE_Im[row][n1][n2] = NAN;
            }

            memset(matrix, 0, matrix_count*sizeof(*matrix));
            memset(rhs, 0, (size_t)multipoles*sizeof(*rhs));
            if (!isfinite((double)nzero.re) || !isfinite((double)nzero.im)) {
                singular_bins++;
                continue;
            }
            if (nzero.re <= 0.0 || ggg_complex_abs2(nzero) <= DBL_MIN) {
                gd->scalar_window_status[bin] = CBALLS_WINDOW_EMPTY;
                singular_bins++;
                continue;
            }

            for (row=0; row<multipoles; row++) {
                int ell = row-nmax;
                for (column=0; column<multipoles; column++) {
                    int order = ell-(column-nmax);
                    ggg_complex value = ggg_histogram_mode(
                        gdlN->histZetaMcos, gdlN->histZetaMsin,
                        gdlN->histZetaMsincos, gdlN->histZetaMcossin,
                        order, n1, n2);
                    matrix[(size_t)row*multipoles+column] =
                        ggg_complex_div(value, nzero);
                    if (!isfinite((double)matrix[
                            (size_t)row*multipoles+column].re)
                        || !isfinite((double)matrix[
                            (size_t)row*multipoles+column].im)) {
                        singular = TRUE;
                        break;
                    }
                    if (ggg_complex_abs2(
                            matrix[(size_t)row*multipoles+column])
                        > matrix_scale)
                        matrix_scale = ggg_complex_abs2(
                            matrix[(size_t)row*multipoles+column]);
                }
                rhs[row] = ggg_complex_div(ggg_histogram_mode(
                    gdl->histZetaMcos, gdl->histZetaMsin,
                    gdl->histZetaMsincos, gdl->histZetaMcossin,
                    ell, n1, n2), nzero);
                if (singular || !isfinite((double)rhs[row].re)
                    || !isfinite((double)rhs[row].im)) {
                    gd->scalar_window_status[bin] = CBALLS_WINDOW_NONFINITE;
                    singular = TRUE;
                    break;
                }
            }
            if (singular || !isfinite((double)matrix_scale)) {
                singular_bins++;
                continue;
            }
            tolerance = 128.0*DBL_EPSILON*(1.0+sqrt(matrix_scale));
            gd->scalar_window_status[bin] = CBALLS_WINDOW_SINGULAR;

            for (column=0; column<multipoles; column++) {
                int pivot = column;
                real pivot_norm = ggg_complex_abs2(
                    matrix[(size_t)column*multipoles+column]);
                int candidate;
                if (!isfinite((double)pivot_norm)) {
                    gd->scalar_window_status[bin] = CBALLS_WINDOW_NONFINITE;
                    singular = TRUE;
                    break;
                }
                for (candidate=column+1; candidate<multipoles; candidate++) {
                    real candidate_norm = ggg_complex_abs2(
                        matrix[(size_t)candidate*multipoles+column]);
                    if (!isfinite((double)candidate_norm)) {
                        gd->scalar_window_status[bin] = CBALLS_WINDOW_NONFINITE;
                        singular = TRUE;
                        break;
                    }
                    if (candidate_norm > pivot_norm) {
                        pivot = candidate;
                        pivot_norm = candidate_norm;
                    }
                }
                if (singular)
                    break;
                if (pivot_norm <= tolerance*tolerance) {
                    gd->scalar_window_pivot_ratio[bin] = 0.0;
                    singular = TRUE;
                    break;
                }
                smallest_pivot = fmin(smallest_pivot, sqrt(pivot_norm));
                largest_pivot = fmax(largest_pivot, sqrt(pivot_norm));
                if (pivot != column) {
                    int k;
                    for (k=0; k<multipoles; k++) {
                        ggg_complex swap = matrix[
                            (size_t)column*multipoles+k];
                        matrix[(size_t)column*multipoles+k] = matrix[
                            (size_t)pivot*multipoles+k];
                        matrix[(size_t)pivot*multipoles+k] = swap;
                    }
                    {
                        ggg_complex swap = rhs[column];
                        rhs[column] = rhs[pivot];
                        rhs[pivot] = swap;
                    }
                }
                {
                    ggg_complex divisor = matrix[
                        (size_t)column*multipoles+column];
                    int k;
                    for (k=0; k<multipoles; k++)
                        matrix[(size_t)column*multipoles+k] =
                            ggg_complex_div(
                                matrix[(size_t)column*multipoles+k], divisor);
                    rhs[column] = ggg_complex_div(rhs[column], divisor);
                }
                for (row=0; row<multipoles; row++) {
                    ggg_complex factor;
                    int k;
                    if (row == column)
                        continue;
                    factor = matrix[(size_t)row*multipoles+column];
                    if (ggg_complex_abs2(factor) == 0.0)
                        continue;
                    for (k=0; k<multipoles; k++)
                        matrix[(size_t)row*multipoles+k] = ggg_complex_sub(
                            matrix[(size_t)row*multipoles+k],
                            ggg_complex_mul(factor,
                                matrix[(size_t)column*multipoles+k]));
                    rhs[row] = ggg_complex_sub(
                        rhs[row], ggg_complex_mul(factor, rhs[column]));
                }
            }
            if (singular) {
                singular_bins++;
                continue;
            }
            for (row=nmax; row<multipoles; row++) {
                if (!isfinite((double)rhs[row].re)
                    || !isfinite((double)rhs[row].im)) {
                    gd->scalar_window_status[bin] = CBALLS_WINDOW_NONFINITE;
                    singular = TRUE;
                    break;
                }
            }
            if (singular) {
                singular_bins++;
                continue;
            }
            gd->scalar_window_status[bin] = CBALLS_WINDOW_VALID;
            gd->scalar_window_pivot_ratio[bin] = smallest_pivot/largest_pivot;
            for (row=nmax; row<multipoles; row++) {
                int output_index = row-nmax+1;
                gdl->histZetaM_EE[output_index][n1][n2] = rhs[row].re;
                gdl->histZetaM_EE_Im[output_index][n1][n2] = rhs[row].im;
            }
        }
    }

    if (singular_bins > 0)
        verb_print_normal_info(
            cmd->verbose, cmd->verbose_log, gd->outlog,
            "octree-ggg-omp: %d edge-correction radial-bin pairs "
            "were unsupported and set to NaN\n", singular_bins);
    free(matrix);
    free(rhs);
    if (gd->rootDirFlag && cballs_scalar_window_write(cmd, gd) == FAILURE)
        return FAILURE;
    gd->scalar_window_ready = TRUE;
    return SUCCESS;
}
#endif


#endif
