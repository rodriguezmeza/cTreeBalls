#ifndef CBALLS_LYA_DISTORTION_H
#define CBALLS_LYA_DISTORTION_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Folded observer-bisector (rp,rt), rp-major/rt-minor. Coordinates are the
 * original observer-centred Cartesian positions. Model support may differ
 * from data support. Projection order: -1 identity, 0 mean, 1 mean+slope. */
typedef struct {
    int rp_bins, rt_bins, model_rp_bins, model_rt_bins;
    int projection_order, threads;
    double rp_max, rt_max, model_rp_max, model_rt_max;
    size_t max_bytes;
} cballs_lya_distortion_config;

int cballs_lya_distortion_available(void);
/* Conservative native workspace + caller inputs/outputs, before allocation. */
int cballs_lya_distortion_plan(size_t n, const cballs_lya_distortion_config *config,
                               size_t *bytes, char *error, size_t error_size);
/* weights select/weight DATA pixels. projection_weights define each forest's
 * fit over ALL supplied pixels, including zero-data-weight model pixels.
 * NULL projection_weights uses weights. log_lambda is required for order 1.
 * Output numerator is unnormalised: D[A,B] = numerator[A,B] / row_weights[A].
 * Empty rows have zero numerator/weight; the Python API marks D rows as NaN.
 * Arrays remain caller-owned. Returns 0 on success, nonzero on failure.
 * stats: pixels, forests, forest pairs, pruned pairs, pairs with data,
 * positive-weight data pixel pairs, degenerate-slope forests, planned bytes.
 * Native calls use no cballs context or borrowed Python pointers after return.
 */
int cballs_lya_distortion_compute(const double *positions, const double *weights,
    const int64_t *forest_ids, const double *log_lambda,
    const double *projection_weights, size_t n,
    const cballs_lya_distortion_config *config,
    double *numerator, size_t numerator_length,
    double *row_weights, size_t row_weights_length,
    uint64_t stats[8], char *error, size_t error_size);

#ifdef __cplusplus
}
#endif
#endif
