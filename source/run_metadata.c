/* Shared, checked JSON provenance for native products and Python snapshots. */
#include "globaldefs.h"
#include "cballs_build_fingerprint.h"
#include <stdarg.h>
#include <errno.h>
#ifdef CBALLS_MPI_ENABLED
#include <mpi.h>
#endif
#ifdef OPENMPCODE
#include <omp.h>
#endif

const char *cballs_build_id(void) { return CBALLS_BUILD_ID; }
const char *cballs_build_json(void) { return CBALLS_BUILD_JSON; }

typedef struct { char *data; size_t size, used; int failed; } metadata_buffer;

static void append(metadata_buffer *b, const char *format, ...)
{
    va_list args, copy;
    if (b->failed) return;
    va_start(args, format); va_copy(copy, args);
    int length = vsnprintf(NULL, 0, format, copy);
    va_end(copy);
    if (length < 0 || (size_t)length > SIZE_MAX-b->used-1) b->failed = TRUE;
    else if (b->used+(size_t)length+1 > b->size) {
        size_t size = b->used+(size_t)length+1;
        if (size <= SIZE_MAX/2) size *= 2;
        /* Metadata is also a bounded allocation; do not grow unchecked when
         * a caller supplies extreme axis sizes. */
        if (cballs_memory_preflight(size,"run metadata",NULL,0)==FAILURE) {
            b->failed=TRUE; va_end(args); return;
        }
        char *data = realloc(b->data, size);
        if (!data) b->failed = TRUE;
        else { b->data = data; b->size = size; }
    }
    if (!b->failed) {
        vsnprintf(b->data+b->used, b->size-b->used, format, args);
        b->used += (size_t)length;
    }
    va_end(args);
}

static void string_value(metadata_buffer *b, const char *value)
{
    if (!value) { append(b, "null"); return; }
    append(b, "\"");
    for (const unsigned char *p=(const unsigned char *)value; *p; p++) {
        if (*p == '"' || *p == '\\') append(b, "\\%c", *p);
        else if (*p < 32) append(b, "\\u%04x", *p);
        else append(b, "%c", *p);
    }
    append(b, "\"");
}

static void number(metadata_buffer *b, double value)
{ if (isfinite(value)) append(b, "%.17g", value); else append(b, "null"); }

static void named_string(metadata_buffer *b, const char *key, const char *value)
{ append(b, "\"%s\":", key); string_value(b, value); }

static void axis(metadata_buffer *b, const char *key, double lo, double hi, size_t bins)
{
    size_t bytes;
    if (!cballs_size_add(bins,1,&bytes) || !cballs_size_mul(bytes,32,&bytes)
        || cballs_memory_preflight(bytes,"metadata axis",NULL,0)==FAILURE) {
        b->failed=TRUE; return;
    }
    append(b, "\"%s\":[", key);
    if (bins > 0) for (size_t i=0; i<=bins && !b->failed; i++) {
        if (i) append(b, ",");
        number(b, lo+(hi-lo)*(double)i/bins);
    }
    append(b, "]");
}

int cballs_run_metadata(struct cmdline_data *cmd, struct global_data *gd, char **result)
{
    metadata_buffer b = {0};
    const char *method = cmd->searchMethod ? cmd->searchMethod : "";
    const int forest = !strncmp(method, "lya-", 4);
    const int radial_forest = forest && strstr(method, "-1d-") != NULL;
    const int shear = strstr(method, "shear") != NULL;
    const int spherical_shear = shear && strstr(method,"-sphere-") != NULL;
    const int physical = strstr(method, "3pcf-3d") != NULL || strstr(method, "ggg-3d") != NULL;
    const int boxes = strstr(method, "box") != NULL;
    const int neighbor = !strcmp(method, "neighbor-boxes-omp");
    const int survey = scanopt(cmd->options, "survey-estimator-3d");
    int rank=0, ranks=1, mpi_active=FALSE;
    *result = NULL;
#ifdef CBALLS_MPI_ENABLED
    int initialized=0, finalized=0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Finalized(&finalized);
    if (initialized && !finalized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        mpi_active = TRUE;
    }
#endif
    append(&b, "{\"schema_version\":1,\"build\":%s,", cballs_build_json());
    named_string(&b, "engine", method);
    append(&b, ",\"engine_id\":%d,", gd->searchMethod_int);
    named_string(&b, "product_kind", gd->stopflag ? "catalog-or-startup" : "correlation");
    append(&b, ",");
    named_string(&b, "estimator", forest ? (radial_forest ? (strstr(method, "same-los") ? "equal-forest mean of within-forest weighted radial 2PCF" : "weighted forest radial correlations") : "weighted anisotropic forest correlations")
        : shear ? (spherical_shear ? "full-sky spin-2 pair correlations and window-corrected natural 3PCF components" : "flat-sky spin-2 pair correlations and natural 3PCF components")
        : physical ? (survey ? "D-alpha*R survey Legendre 2PCF/3PCF with random-window correction" : "weighted Euclidean Legendre 2PCF/3PCF")
        : boxes ? (neighbor ? "periodic ordered pair counts" : "periodic unordered pair counts and optional N-squared shell density contrast")
        : cballs_opt_edge_corrections(cmd) ? "scalar Fourier angular window deconvolution"
        : cballs_opt_no_normalize_histzeta(cmd) ? "raw scalar distinct-neighbor tangent-Fourier moments"
        : "weight-normalized scalar Fourier moments");
    append(&b, ",\"options\":"); string_value(&b, cmd->options);
    append(&b, ",\"bin_edges\":{");
    if (!forest) {
        int logarithmic = cmd->useLogHist;
        double minimum = cmd->rminHist;
        if (neighbor) {
            minimum = 0;
#ifdef LOGBINCBON
#ifdef _LOGBIN_
            logarithmic = TRUE;
#else
            logarithmic = FALSE;
#endif
#endif
        }
        append(&b, "\"radial\":[");
        for (int i=0; i<=cmd->sizeHistN; i++) {
            double edge;
            if (i) append(&b, ",");
            if (!logarithmic) edge = minimum+(cmd->rangeN-minimum)*(double)i/cmd->sizeHistN;
            else if (minimum > 0) edge = minimum*pow(cmd->rangeN/minimum, (double)i/cmd->sizeHistN);
            else edge = cmd->rangeN*pow(10., ((i == 0 ? -1 : i)-cmd->sizeHistN)/(double)cmd->logHistBinsPD);
            number(&b, edge);
        }
        append(&b, "],\"logarithmic\":%s,\"zero_cutoff_truncation_extension\":%s",
               logarithmic ? "true":"false", logarithmic && minimum == 0 ? "true":"false");
        if (shear) { append(&b, ","); axis(&b, "phi_radians", -PI, PI, gd->shearAngularBins); }
    }
#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
    if (forest) {
        const int kind=lya_forest_method_kind(method);
        const int pair=kind==0 || kind==2 || kind==3 || kind==5 || kind==6 || kind==8;
        const int triple=kind==1 || kind==2 || kind==4 || kind==5 || kind==7;
        if (pair) {
            axis(&b,"two_point_abs_parallel",0,cmd->lya2RpMax,(size_t)cmd->lya2RpBins);
            if (!radial_forest) {
                append(&b,","); axis(&b,"two_point_transverse",0,cmd->lya2RtMax,(size_t)cmd->lya2RtBins);
            }
        }
        if (triple) {
            if (pair) append(&b,",");
            if (!radial_forest) {
                axis(&b,"three_point_side",0,cmd->lya3RMax,(size_t)cmd->lya3RBins);
                append(&b,","); axis(&b,"three_point_los_angle",0,PI,(size_t)cmd->lya3ThetaBins);
                append(&b,","); axis(&b,"three_point_opening_cosine",-1,1,(size_t)cmd->lya3MuBins);
            } else {
                size_t signed_bins;
                if (!cballs_size_mul(2,(size_t)cmd->lya3RBins,&signed_bins)) b.failed=TRUE;
                else axis(&b,"three_point_signed_parallel",-cmd->lya3RMax,cmd->lya3RMax,signed_bins);
            }
        }
    }
#endif
    append(&b, "},\"multipole_max\":%d,\"coordinate_convention\":", cmd->mChebyshev);
    string_value(&b, forest ? "observer-centered Cartesian comoving positions; radial norm/LOS separation; no periodic wrapping"
        : shear ? (spherical_shear ? "unit-sphere Cartesian directions; chord-distance bins; parallel-transported spin-2 tangent components" : "Cartesian tangent-plane positions; planar distances and spin-2 components; no spherical transport")
        : physical ? "Cartesian physical separations; interior 3D opening angle; Legendre multipoles"
        : boxes || cmd->usePeriodic ? "Cartesian periodic box; minimum-image separations"
        : NDIM == 2 ? "Cartesian plane; Euclidean separations and planar Fourier bearings"
        : "original observer frame; Cartesian chords and tangent bearings; undefined radial/antipodal bearings excluded from angular moments");
    append(&b, ",\"geometry\":{\"dimensions\":%d,\"spherical_shear\":%s,\"scalar_observer_frame\":%s},",
        NDIM,spherical_shear?"true":"false",!shear && cballs_observer_frame(cmd)?"true":"false");
    append(&b, "\"resources\":{\"common_histogram_bytes\":%zu,\"common_scalar_3pcf\":%s,\"memory_budget_bytes\":%zu,\"budget_scope\":\"per-rank common histogram plan and individual common allocations; not total RSS\"}",
        gd->common_histogram_bytes,gd->common_scalar_3pcf?"true":"false",gd->memory_budget_bytes);
    append(&b, ",\"box\":[");
    for (int k=0; k<NDIM; k++) { if (k) append(&b, ","); number(&b, gd->Box[k]); }
    append(&b, "],\"weights\":{");
    named_string(&b, "contract", forest || shear || physical ? "supplied catalog weights; unit weights when loader defaults them"
                  : boxes ? "count estimator; catalog scalar field is not the density contrast"
                  : "weights-norm selects supplied catalog weights; otherwise unit weights");
    append(&b, ",\"weights_norm_option\":%s},\"masks\":{\"read_mask\":%s,",
        cballs_opt_weights_norm(cmd)?"true":"false", cballs_opt_read_mask(cmd)?"true":"false");
    named_string(&b, "selection", forest ? (strstr(method, "same-los") ? "same-forest pairs" : "exclude same-forest pairs; three distinct forests for triplets")
        : "reader/embedded mask selection and engine options recorded with input descriptors");
    append(&b, "},\"effective_smoothing\":{\"supported\":%s,\"enabled\":%s,\"radius\":",
        cballs_run_supports_smooth_pivot(cmd)?"true":"false", cballs_opt_smooth_pivot(cmd)?"true":"false");
    number(&b, cballs_opt_smooth_pivot(cmd) ? gd->rsmooth[0] : 0);
    append(&b, ",\"requested_radius_units\":"); string_value(&b,spherical_shear?"arcmin":"Cartesian catalog units");
    append(&b, ",\"effective_radius_units\":"); string_value(&b,spherical_shear?"unit-sphere chord":"Cartesian catalog units");
    append(&b, ",\"requested\":"); string_value(&b, cmd->rsmooth);
    append(&b, ",\"nsmooth\":%d},\"opening_tolerance\":{\"theta\":", cmd->nsmooth);
    number(&b, cmd->theta);
    append(&b, ",\"no_one_ball\":%s,\"no_two_balls\":%s,\"legacy_one_ball\":%s},",
        cballs_opt_no_one_ball(cmd)?"true":"false", cballs_opt_no_two_balls(cmd)?"true":"false",
        cballs_opt_legacy_one_ball(cmd)?"true":"false");
    append(&b, "\"precision\":{\"storage_bits\":%u,\"compute_bits\":%u,\"accumulator_bits\":%u,\"integer_bits\":%u,\"long_double_bits\":%u},",
        (unsigned)(CHAR_BIT*sizeof(cballs_storage_real)), (unsigned)(CHAR_BIT*sizeof(cballs_compute_real)),
        (unsigned)(CHAR_BIT*sizeof(cballs_accum_real)), (unsigned)(CHAR_BIT*sizeof(INTEGER)), (unsigned)(CHAR_BIT*sizeof(long double)));
    append(&b, "\"parallel\":{\"mpi_initialized\":%s,\"rank\":%d,\"world_ranks\":%d,\"estimator_ranks\":%d,\"threads_requested\":%d",
        mpi_active?"true":"false", rank, ranks, strstr(method,"-mpi")?ranks:1, cmd->numthreads);
#ifdef OPENMPCODE
    int observed_threads = 1;
#pragma omp parallel
    {
#pragma omp single
        observed_threads = omp_get_num_threads();
    }
    append(&b, ",\"openmp_max_threads\":%d,\"openmp_dynamic\":%s,\"openmp_probe_threads\":%d",
           omp_get_max_threads(), omp_get_dynamic()?"true":"false", observed_threads);
#else
    append(&b, ",\"openmp_max_threads\":1,\"openmp_dynamic\":false,\"openmp_probe_threads\":1");
#endif
    append(&b, "},\"inputs\":{");
    named_string(&b, "files", cmd->infile); append(&b, ",");
    named_string(&b, "formats", cmd->infilefmt); append(&b, ",");
    named_string(&b, "catalog_selection", cmd->iCatalogs); append(&b, ",");
#ifdef IOLIB
    named_string(&b, "columns", cmd->columns); append(&b, ",");
#endif
    named_string(&b, "test_model", cmd->testmodel);
    append(&b, ",\"seed\":%d,\"catalog_sizes\":[", cmd->seed);
    for (int i=0; i<gd->ninfiles; i++) { if (i) append(&b, ","); append(&b, "%lld", (long long)gd->nbodyTable[i]); }
    append(&b, "]},\"output\":{"); named_string(&b, "root", cmd->rootDir);
    append(&b, ","); named_string(&b, "catalog", cmd->outfile);
    append(&b, ","); named_string(&b, "catalog_format", cmd->outfilefmt);
    append(&b, "},\"native_timing_semantics\":\"CPUTIME and getCPUTime are process CPU seconds on this rank; no division or MPI sum/max is implied\",\"scalar_window_ready\":%s}\n", gd->scalar_window_ready?"true":"false");
    if (b.failed) {
        free(b.data);
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "run metadata: allocation or formatting failed");
        return FAILURE;
    }
    *result = b.data;
    return SUCCESS;
}

int cballs_write_run_metadata(struct cmdline_data *cmd, struct global_data *gd)
{
    char path[MAXLENGTHOFFILES+64], *json = NULL;
    FILE *file = NULL;
    int status = SUCCESS;
    if (!gd->rootDirFlag || (cballs_opt_no_out_hist(cmd) && strnull(cmd->outfile))) return SUCCESS;
#ifdef CBALLS_MPI_ENABLED
    if (!cballs_mpi_output_enabled(cmd)) return SUCCESS;
#endif
    int length = snprintf(path, sizeof(path), "%s/run-metadata.json", cmd->rootDir);
    if (length < 0 || (size_t)length >= sizeof(path)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "run metadata: output path too long");
        return FAILURE;
    }
    if (cballs_run_metadata(cmd, gd, &json) == FAILURE) return FAILURE;
    file = fopen(path, "w");
    if (!file || fputs(json, file) == EOF) status = FAILURE;
    if (file && fclose(file) != 0) status = FAILURE;
    if (status == FAILURE) snprintf(cmd->error_message, _ERRORMSGSIZE_,
        "run metadata: cannot write '%s': %s", path, strerror(errno));
    free(json);
    return status;
}
