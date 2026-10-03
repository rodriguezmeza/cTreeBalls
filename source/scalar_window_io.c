/* Versioned raw scalar-window exports and checked saved-file correction.
 * File indices are mode+1; the manifest is the last-written completion marker.
 * See docs/saved_scalar_edge.md for the scientific and file-format contract. */
#include "globaldefs.h"
#include "scalar_window_solver.h"
#include <ctype.h>
#include <limits.h>
#include <unistd.h>

#define EDGE_PATH_SIZE (MAXLENGTHOFFILES + 96)
#define EDGE_MAGIC "CBALLS_SCALAR_EDGE"
#define EDGE_NORM "raw-ordered-distinct"
#define EDGE_PHASE "cc+ss+i(sc-cs)"

typedef struct {
    int bins, mmax, wmax, dimensions, periodic, logarithmic;
    double lower_cutoff, box[3], *edges;
} edge_metadata;

static int edge_error(struct cmdline_data *cmd, const char *what, const char *path)
{
    snprintf(cmd->error_message, _ERRORMSGSIZE_, "saved scalar edge: %s: %s",
             what, path);
    return FAILURE;
}

static int edge_path(struct cmdline_data *cmd, char *path,
                     const char *prefix, const char *suffix, int mode)
{
    const int length = mode < 0
        ? snprintf(path, EDGE_PATH_SIZE, "%s_%s%s", prefix, suffix, EXTFILES)
        : snprintf(path, EDGE_PATH_SIZE, "%s_%s_%d%s", prefix, suffix, mode+1, EXTFILES);
    return length < 0 || length >= EDGE_PATH_SIZE
        ? edge_error(cmd, "path too long", prefix) : SUCCESS;
}

static int edge_remove(struct cmdline_data *cmd, const char *path)
{
    if (unlink(path) != 0 && errno != ENOENT)
        return edge_error(cmd, "cannot invalidate completion marker", path);
    return SUCCESS;
}

/* Parse tokens with complete conversion and reject surplus input. No atoi,
 * inferred mode counts, unbounded line buffers, or current-directory rbins. */
static int edge_token(FILE *file, const char *expected)
{
    char token[256];
    return fscanf(file, "%255s", token) == 1 && !strcmp(token, expected);
}

static int edge_integer(FILE *file, int *value)
{
    char token[256], *end;
    if (fscanf(file, "%255s", token) != 1) return FALSE;
    errno = 0;
    long number = strtol(token, &end, 10);
    if (errno || end == token || *end || number < 0 || number > INT_MAX-1)
        return FALSE;
    *value = (int)number;
    return TRUE;
}

static int edge_number(FILE *file, double *value)
{
    char token[256], *end;
    if (fscanf(file, "%255s", token) != 1) return FALSE;
    errno = 0;
    *value = strtod(token, &end);
    /* Explicit nan/inf payloads are permitted, and get per-bin status 4.
     * Representable subnormals are retained; overflow or underflow to zero
     * in text is a malformed export. */
    return (errno == 0 || (errno == ERANGE && isfinite(*value) && *value != 0))
        && end != token && !*end;
}

static int edge_eof(FILE *file)
{
    int ch;
    do { ch = fgetc(file); } while (ch != EOF && isspace((unsigned char)ch));
    return ch == EOF && !ferror(file);
}

static int edge_read_metadata(struct cmdline_data *cmd, const char *prefix,
                              edge_metadata *metadata)
{
    char path[EDGE_PATH_SIZE];
    FILE *file;
    int version, status = FAILURE;
    if (edge_path(cmd, path, prefix, "edge_manifest", -1) == FAILURE) return FAILURE;
    file = fopen(path, "r");
    if (!file) return edge_error(cmd,
        "cannot open v1 manifest; regenerate legacy/incomplete exports", path);
    if (!edge_token(file, EDGE_MAGIC) || !edge_integer(file, &version) || version != 1
        || !edge_token(file, "bins") || !edge_integer(file, &metadata->bins)
        || !edge_token(file, "signal_mmax") || !edge_integer(file, &metadata->mmax)
        || !edge_token(file, "window_mmax") || !edge_integer(file, &metadata->wmax)
        || !edge_token(file, "dimensions") || !edge_integer(file, &metadata->dimensions)
        || !edge_token(file, "periodic") || !edge_integer(file, &metadata->periodic)
        || metadata->bins < 1 || metadata->mmax > (INT_MAX-1)/2
        || metadata->wmax < 2*metadata->mmax || metadata->wmax > (INT_MAX-1)/2
        || (metadata->dimensions != 2 && metadata->dimensions != 3)
        || metadata->periodic > 1 || !edge_token(file, "box")) goto malformed;
    for (int k=0; k<metadata->dimensions; k++)
        if (!edge_number(file, &metadata->box[k]) || !isfinite(metadata->box[k])
            || (metadata->periodic && metadata->box[k] <= 0)) goto malformed;
    if (!edge_token(file, "geometry") || !edge_token(file,
            metadata->dimensions == 2 ? "planar-fourier" : "observer-tangent-fourier")
        || !edge_token(file, "normalization") || !edge_token(file, EDGE_NORM)
        || !edge_token(file, "phase") || !edge_token(file, EDGE_PHASE)
        || !edge_token(file, "logarithmic") || !edge_integer(file, &metadata->logarithmic)
        || metadata->logarithmic > 1
        || !edge_token(file, "lower_cutoff") || !edge_number(file, &metadata->lower_cutoff)
        || !isfinite(metadata->lower_cutoff) || metadata->lower_cutoff < 0
        || !edge_token(file, "edges")) goto malformed;
    size_t cells, minimum_bytes;
    if (!cballs_size_mul((size_t)metadata->bins, (size_t)metadata->bins, &cells)
        || !cballs_size_mul(cells, sizeof(double complex), &minimum_bytes)) {
        edge_error(cmd, "resource workspace size overflow", path); goto cleanup;
    }
    if (cballs_memory_preflight(minimum_bytes, "saved scalar minimum workspace",
            cmd->error_message, _ERRORMSGSIZE_) == FAILURE) goto cleanup;
    if (cballs_malloc_checked((void **)&metadata->edges, (size_t)metadata->bins+1,
            sizeof(double), "saved scalar radial edges", cmd->error_message,
            _ERRORMSGSIZE_) == FAILURE) goto cleanup;
    for (int i=0; i<=metadata->bins; i++)
        if (!edge_number(file, &metadata->edges[i]) || !isfinite(metadata->edges[i])
            || metadata->edges[i] < 0
            || (i && metadata->edges[i] <= metadata->edges[i-1])) goto malformed;
    if (metadata->lower_cutoff >= metadata->edges[metadata->bins]
        || !edge_token(file, "complete") || !edge_eof(file)) goto malformed;
    status = SUCCESS;
    goto cleanup;
malformed:
    edge_error(cmd, "malformed or unsupported manifest", path);
cleanup:
    if (fclose(file) != 0 && status == SUCCESS)
        status = edge_error(cmd, "manifest close failed", path);
    return status;
}

/* Independent windows may have different objects/weights and more modes.
 * Their measurement convention and radial grid must nevertheless agree. */
static int edge_compatible(struct cmdline_data *cmd, const edge_metadata *s,
                           const edge_metadata *w)
{
    if (s->bins != w->bins || s->dimensions != w->dimensions
        || s->periodic != w->periodic || s->logarithmic != w->logarithmic
        || s->lower_cutoff != w->lower_cutoff || w->wmax < 2*s->mmax)
        return edge_error(cmd, "incompatible dimensions, modes or grid", "input prefixes");
    for (int k=0; k<s->dimensions; k++)
        if (s->periodic && s->box[k] != w->box[k])
            return edge_error(cmd, "incompatible periodic boxes", "input prefixes");
    for (int i=0; i<=s->bins; i++)
        if (s->edges[i] != w->edges[i])
            return edge_error(cmd, "incompatible radial edges", "input prefixes");
    return SUCCESS;
}

/* A plane has exactly bins*bins whitespace-separated values, row-major. */
static int edge_read_plane(struct cmdline_data *cmd, const char *prefix,
        const char *suffix, int mode, size_t plane, double complex *target,
        bool imaginary, double sign)
{
    char path[EDGE_PATH_SIZE];
    FILE *file;
    int status = FAILURE;
    if (edge_path(cmd, path, prefix, suffix, mode) == FAILURE) return FAILURE;
    file = fopen(path, "r");
    if (!file) return edge_error(cmd, "cannot open required mode", path);
    for (size_t i=0; i<plane; i++) {
        double value;
        if (!edge_number(file, &value)) goto malformed;
        target[i] = imaginary ? cballs_scalar_complex(creal(target[i]), cimag(target[i])+sign*value)
                              : cballs_scalar_complex(creal(target[i])+sign*value, cimag(target[i]));
    }
    if (!edge_eof(file)) goto malformed;
    status = SUCCESS;
    goto cleanup;
malformed:
    edge_error(cmd, "malformed matrix or wrong value count", path);
cleanup:
    if (fclose(file) != 0 && status == SUCCESS)
        status = edge_error(cmd, "matrix close failed", path);
    return status;
}

static int edge_write_plane(struct cmdline_data *cmd, struct global_data *gd,
        const char *suffix, int mode, const double complex *values, bool imaginary)
{
    char path[EDGE_PATH_SIZE];
    int status = SUCCESS;
    if (edge_path(cmd, path, gd->fpfnamehistZetaMFileName, suffix, mode) == FAILURE)
        return FAILURE;
    FILE *file = fopen(path, "w");
    if (!file) return edge_error(cmd, "cannot open output", path);
    for (int i=0; i<cmd->sizeHistN && status == SUCCESS; i++) {
        for (int j=0; j<cmd->sizeHistN; j++) {
            double complex value = values[(size_t)i*cmd->sizeHistN+j];
            if (fprintf(file, "%.17g ", imaginary ? cimag(value) : creal(value)) < 0)
                status = FAILURE;
        }
        if (fputc('\n', file) == EOF) status = FAILURE;
    }
    if (fclose(file) != 0) status = FAILURE;
    return status == SUCCESS ? SUCCESS : edge_error(cmd, "output write failed", path);
}

/* Export metadata contains the actual bin edges, including the legacy
 * zero-cutoff logarithmic first-bin extension used by scalar traversals. */
int cballs_scalar_edge_manifest(struct cmdline_data *cmd,
                                struct global_data *gd, bool publish)
{
    char path[EDGE_PATH_SIZE];
    if (edge_path(cmd, path, gd->fpfnamehistZetaMFileName, "edge_manifest", -1)
        == FAILURE) return FAILURE;
    if (!publish) return edge_remove(cmd, path);
    FILE *file = fopen(path, "w");
    if (!file) return edge_error(cmd, "cannot publish manifest", path);
    int status = SUCCESS;
#define EM_WRITE(...) do { if (fprintf(file, __VA_ARGS__) < 0) status = FAILURE; } while (0)
    EM_WRITE(EDGE_MAGIC " 1\nbins %d\nsignal_mmax %d\nwindow_mmax %d\n"
        "dimensions %d\nperiodic %d\nbox", cmd->sizeHistN, cmd->mChebyshev,
        2*cmd->mChebyshev, NDIM, !!cmd->usePeriodic);
    for (int k=0; k<NDIM; k++) EM_WRITE(" %.17g", (double)gd->Box[k]);
    EM_WRITE("\ngeometry %s\nnormalization " EDGE_NORM "\nphase " EDGE_PHASE
        "\nlogarithmic %d\nlower_cutoff %.17g\nedges",
        NDIM == 2 ? "planar-fourier" : "observer-tangent-fourier",
        !!cmd->useLogHist, (double)cmd->rminHist);
    for (int i=0; i<=cmd->sizeHistN; i++) {
        double edge;
        if (!cmd->useLogHist)
            edge = cmd->rminHist+(cmd->rangeN-cmd->rminHist)*(double)i/cmd->sizeHistN;
        else if (cmd->rminHist > 0)
            edge = cmd->rminHist*pow(cmd->rangeN/cmd->rminHist, (double)i/cmd->sizeHistN);
        else edge = cmd->rangeN*pow(10., ((i == 0 ? -1 : i)-cmd->sizeHistN)/(double)cmd->logHistBinsPD);
        if (!isfinite(edge)) status = FAILURE;
        EM_WRITE(" %.17g", edge);
    }
    EM_WRITE("\ncomplete\n");
#undef EM_WRITE
    if (fclose(file) != 0) status = FAILURE;
    if (status == FAILURE) {
        unlink(path);
        return edge_error(cmd, "manifest write failed", path);
    }
    return SUCCESS;
}

int computeEdgeCorrections(struct cmdline_data *cmd, struct global_data *gd)
{
    edge_metadata signal_meta = {0}, window_meta = {0};
    const int saved_bins = cmd->sizeHistN, saved_mmax = cmd->mChebyshev;
    double complex *signal = NULL, *window = NULL, *result = NULL;
    double complex *matrix = NULL, *rhs = NULL, *modes = NULL;
    size_t plane, scount, wcount, matrix_count, total, count;
    int status = FAILURE, bins, mmax, orders, n;
    char marker[EDGE_PATH_SIZE] = {0};
    bool writing = FALSE;
    if (gd->ninfiles != 2)
        return edge_error(cmd, "exactly two prefixes required", "in=signal_prefix,window_prefix");
    if (cballs_opt_no_out_hist(cmd) || cballs_opt_full_sky(cmd))
        return edge_error(cmd, "no-out-Hist and full-sky are unsupported for this solve", "options");
    if (edge_read_metadata(cmd, gd->infilenames[0], &signal_meta) == FAILURE
        || edge_read_metadata(cmd, gd->infilenames[1], &window_meta) == FAILURE
        || edge_compatible(cmd, &signal_meta, &window_meta) == FAILURE) goto cleanup;
    bins = signal_meta.bins; mmax = signal_meta.mmax; orders = mmax+1; n = 2*mmax+1;
    /* Bound the complete live workspace, including diagnostics, before allocation. */
    if (!cballs_size_mul((size_t)bins, (size_t)bins, &plane)
        || !cballs_size_mul(plane, (size_t)orders, &scount)
        || !cballs_size_mul(plane, (size_t)n, &wcount)
        || !cballs_size_mul((size_t)n, (size_t)n, &matrix_count)
        || !cballs_size_add(scount, scount, &total)
        || !cballs_size_add(total, wcount, &total)
        || !cballs_size_add(total, matrix_count, &total)
        || !cballs_size_add(total, (size_t)2*n+orders, &total)
        || !cballs_size_mul(total, sizeof(double complex), &total)
        || !cballs_size_mul(plane, 2*sizeof(real)+sizeof(unsigned char), &count)
        || !cballs_size_add(total, count, &total)) {
        edge_error(cmd, "workspace size overflow", "manifest dimensions"); goto cleanup;
    }
    if (cballs_memory_preflight(total, "saved scalar edge workspace", cmd->error_message,
                                _ERRORMSGSIZE_) == FAILURE) goto cleanup;
#define EM_ALLOC(pointer, values) do { \
    if (cballs_calloc_checked((void **)&pointer, values, sizeof(*pointer), \
        "saved scalar edge workspace", cmd->error_message, _ERRORMSGSIZE_) == FAILURE) \
        goto cleanup; \
} while (0)
    EM_ALLOC(signal, scount); EM_ALLOC(window, wcount); EM_ALLOC(result, scount);
    EM_ALLOC(matrix, matrix_count); EM_ALLOC(rhs, n); EM_ALLOC(modes, (size_t)n+orders);
#undef EM_ALLOC
    for (int m=0; m<orders; m++) {
        double complex *target = signal+(size_t)m*plane;
        if (edge_read_plane(cmd, gd->infilenames[0], "edge_cos", m, plane, target, FALSE, 1) == FAILURE
            || edge_read_plane(cmd, gd->infilenames[0], "edge_sin", m, plane, target, FALSE, 1) == FAILURE
            || edge_read_plane(cmd, gd->infilenames[0], "edge_sincos", m, plane, target, TRUE, 1) == FAILURE
            || edge_read_plane(cmd, gd->infilenames[0], "edge_cossin", m, plane, target, TRUE, -1) == FAILURE)
            goto cleanup;
    }
    for (int m=0; m<n; m++) {
        double complex *target = window+(size_t)m*plane;
        if (edge_read_plane(cmd, gd->infilenames[1], "window_Re", m, plane, target, FALSE, 1) == FAILURE
            || edge_read_plane(cmd, gd->infilenames[1], "window_Im", m, plane, target, TRUE, 1) == FAILURE)
            goto cleanup;
    }
    for (size_t b=0; b<plane; b++)
        if ((isfinite(cimag(signal[b])) && cimag(signal[b]) != 0)
            || (isfinite(cimag(window[b])) && cimag(window[b]) != 0)) {
            edge_error(cmd, "monopoles must be real", "input matrices"); goto cleanup;
        }
    cmd->sizeHistN = bins; cmd->mChebyshev = mmax;
    if (cballs_scalar_window_begin(cmd, gd) == FAILURE) goto cleanup;
    for (size_t b=0; b<plane; b++) {
        for (int m=0; m<orders; m++) modes[m] = signal[(size_t)m*plane+b];
        for (int m=0; m<n; m++) modes[orders+m] = window[(size_t)m*plane+b];
        double ratio;
        const int solve_status = cballs_scalar_edge_deconvolve(
            modes, modes+orders, mmax, matrix, rhs, &ratio);
        gd->scalar_window_w0[b] = creal(window[b]);
        gd->scalar_window_status[b] = solve_status;
        gd->scalar_window_pivot_ratio[b] = ratio;
        for (int m=0; m<orders; m++) result[(size_t)m*plane+b] =
            solve_status == CBALLS_WINDOW_VALID ? rhs[mmax+m] : cballs_scalar_complex(NAN, NAN);
    }
    if (edge_path(cmd, marker, gd->fpfnamehistZetaMFileName, "edge_result", -1) == FAILURE
        || edge_remove(cmd, marker) == FAILURE) goto cleanup;
    writing = TRUE;
    for (int m=0; m<orders; m++)
        if (edge_write_plane(cmd, gd, "EE", m, result+(size_t)m*plane, FALSE) == FAILURE
            || edge_write_plane(cmd, gd, "EE_Im", m, result+(size_t)m*plane, TRUE) == FAILURE)
            goto cleanup;
    if (cballs_scalar_window_write(cmd, gd) == FAILURE) goto cleanup;
    FILE *file = fopen(marker, "w");
    if (!file) { edge_error(cmd, "cannot open completion marker", marker); goto cleanup; }
    int written = fprintf(file, "CBALLS_SCALAR_EDGE_RESULT 1\n"
        "signal_prefix %s\nwindow_prefix %s\nbins %d\nmmax %d\n"
        "equation sum_n W_(ell-n)*zeta_n=S_ell; -M<=ell,n<=M\n"
        "invalid_bins NaN; see window_diagnostics\nedges",
        gd->infilenames[0], gd->infilenames[1], bins, mmax);
    for (int i=0; i<=bins; i++)
        if (fprintf(file, " %.17g", signal_meta.edges[i]) < 0) written = -1;
    if (fprintf(file, "\ncomplete\n") < 0) written = -1;
    if (fclose(file) != 0) written = -1;
    if (written < 0) { edge_error(cmd, "completion marker write failed", marker); goto cleanup; }
    /* This is a file-producing preprocessing task, not a native histogram run. */
    status = SUCCESS;
cleanup:
    if (status == FAILURE && writing) unlink(marker);
    cballs_scalar_window_free(gd);
    /* Startup controls are compared across MPI ranks after root preprocessing.
     * File dimensions are local workspace settings, not a rank-local command change. */
    cmd->sizeHistN = saved_bins; cmd->mChebyshev = saved_mmax;
    free(modes); free(rhs); free(matrix); free(result); free(window); free(signal);
    free(window_meta.edges); free(signal_meta.edges);
    return status;
}
