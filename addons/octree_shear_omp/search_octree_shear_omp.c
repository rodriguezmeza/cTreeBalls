/*
 * Flat-sky weak-lensing shear 2PCF and 3PCF estimator.
 *
 * The octree is used for exact radial pruning.  Accepted cells are descended
 * to leaves because a cell-averaged spin-2 field does not preserve the phase
 * factors required by the natural shear three-point functions.
 */

#include "globaldefs.h"
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdarg.h>

#if NDIM < 2

global int prepare_octree_shear_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btable, INTEGER *nbody)
{
    (void)gd;
    (void)btable;
    (void)nbody;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-shear-omp requires NDIM >= 2");
    return FAILURE;
}

global int searchcalc_octree_shear_omp(struct cmdline_data *cmd,
                                       struct global_data *gd,
                                       bodyptr *btable, INTEGER *nbody,
                                       INTEGER ipmin, INTEGER *ipmax,
                                       int cat1, int cat2, int cat3)
{
    (void)gd;
    (void)btable;
    (void)nbody;
    (void)ipmin;
    (void)ipmax;
    (void)cat1;
    (void)cat2;
    (void)cat3;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-shear-omp requires NDIM >= 2");
    return FAILURE;
}

#else

typedef struct {
    real re;
    real im;
} shear_complex;

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    bodyptr pivot;
    int bins;
    int ring_max;
    shear_complex *g_ring_first;
    shear_complex *w_ring_first;
    shear_complex *g_ring_second;
    shear_complex *w_ring_second;
    shear_complex *g_ring_active;
    shear_complex *w_ring_active;
    shear_complex *diag_g6;
    shear_complex *diag_g2;
    shear_complex *diag_abs2;
    real *diag_w2;
    shear_complex *xi_plus;
    shear_complex *xi_minus;
    real *xi_weight;
    bool collect_first_leg;
    bool same_neighbor_catalog;
} shear_pivot_workspace;

static shear_complex shear_make(real re, real im)
{
    shear_complex value;
    value.re = re;
    value.im = im;
    return value;
}

static shear_complex shear_add(shear_complex a, shear_complex b)
{
    return shear_make(a.re + b.re, a.im + b.im);
}

static shear_complex shear_sub(shear_complex a, shear_complex b)
{
    return shear_make(a.re - b.re, a.im - b.im);
}

static shear_complex shear_mul(shear_complex a, shear_complex b)
{
    return shear_make(a.re*b.re - a.im*b.im,
                      a.re*b.im + a.im*b.re);
}

static shear_complex shear_conj(shear_complex value)
{
    return shear_make(value.re, -value.im);
}

static shear_complex shear_scale(shear_complex value, real scale)
{
    return shear_make(value.re*scale, value.im*scale);
}

static real shear_abs2(shear_complex value)
{
    return value.re*value.re + value.im*value.im;
}

static shear_complex shear_div(shear_complex numerator,
                               shear_complex denominator)
{
    real norm = shear_abs2(denominator);
    shear_complex product = shear_mul(numerator, shear_conj(denominator));
    return shear_scale(product, 1.0/norm);
}

static int shear_size_mul(size_t a, size_t b, size_t *result)
{
    if (a != 0 && b > SIZE_MAX/a)
        return FAILURE;
    *result = a*b;
    return SUCCESS;
}

static int shear_calloc(struct cmdline_data *cmd, void **pointer,
                        size_t count, size_t item_size, const char *label)
{
    size_t bytes;

    *pointer = NULL;
    if (shear_size_mul(count, item_size, &bytes) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: allocation size overflow for %s", label);
        return FAILURE;
    }
    if (bytes == 0)
        bytes = 1;
    *pointer = calloc(1, bytes);
    if (*pointer == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: cannot allocate %zu bytes for %s",
                 bytes, label);
        return FAILURE;
    }
    return SUCCESS;
}

static void shear_clear_results(struct global_data *gd)
{
    free(gd->histShearGammaIm);
    free(gd->histShearGammaRe);
    free(gd->histShearDenominatorIm);
    free(gd->histShearDenominatorRe);
    free(gd->histShearGammaMultipoleIm);
    free(gd->histShearGammaMultipoleRe);
    free(gd->histShearGammaNumeratorIm);
    free(gd->histShearGammaNumeratorRe);
    free(gd->histShearXiWeight);
    free(gd->histShearXiMinusIm);
    free(gd->histShearXiMinusRe);
    free(gd->histShearXiPlusIm);
    free(gd->histShearXiPlusRe);
    gd->histShearGammaIm = NULL;
    gd->histShearGammaRe = NULL;
    gd->histShearDenominatorIm = NULL;
    gd->histShearDenominatorRe = NULL;
    gd->histShearGammaMultipoleIm = NULL;
    gd->histShearGammaMultipoleRe = NULL;
    gd->histShearGammaNumeratorIm = NULL;
    gd->histShearGammaNumeratorRe = NULL;
    gd->histShearXiWeight = NULL;
    gd->histShearXiMinusIm = NULL;
    gd->histShearXiMinusRe = NULL;
    gd->histShearXiPlusIm = NULL;
    gd->histShearXiPlusRe = NULL;
    gd->shearMultipoleMax = 0;
    gd->shearAngularBins = 0;
}

static size_t shear_ring_index(int order, int bin, int ring_max, int bins)
{
    return (size_t)(order + ring_max)*(size_t)bins + (size_t)bin;
}

static size_t shear_gamma_index(int component, int order_index,
                                int bin1, int bin2,
                                int multipoles, int bins)
{
    return ((((size_t)component*(size_t)multipoles
              + (size_t)order_index)*(size_t)bins
             + (size_t)bin1)*(size_t)bins + (size_t)bin2);
}

static size_t shear_denominator_index(int order_index, int bin1, int bin2,
                                      int bins)
{
    return (((size_t)order_index*(size_t)bins + (size_t)bin1)
            *(size_t)bins + (size_t)bin2);
}

static size_t shear_angular_index(int component, int phi_bin,
                                  int bin1, int bin2,
                                  int phi_bins, int bins)
{
    return ((((size_t)component*(size_t)phi_bins + (size_t)phi_bin)
             *(size_t)bins + (size_t)bin1)*(size_t)bins + (size_t)bin2);
}

static int shear_radial_bin(struct cmdline_data *cmd,
                            struct global_data *gd, real distance)
{
    int bin;

    if (!isfinite(distance) || distance <= cmd->rminHist
        || distance >= cmd->rangeN)
        return -1;

    if (cmd->useLogHist) {
        if (cmd->rminHist == 0.0)
            bin = (int)(cmd->logHistBinsPD
                        *(rlog10(distance) - rlog10(cmd->rangeN))
                        + cmd->sizeHistN);
        else
            bin = (int)(rlog10(distance/cmd->rminHist)*gd->i_deltaR);
    } else {
        bin = (int)((distance - cmd->rminHist)*gd->i_deltaR);
    }

    return bin >= 0 && bin < cmd->sizeHistN ? bin : -1;
}

static void shear_accumulate_neighbor(shear_pivot_workspace *work, nodeptr q)
{
    struct cmdline_data *cmd = work->cmd;
    struct global_data *gd = work->gd;
    bodyptr p = work->pivot;
    real distance;
    compute_vector dr;
    real cos_phi;
    real sin_phi;
    real weight;
    int bin;
    int order;
    shear_complex phase;
    shear_complex phase_power;
    shear_complex gamma;
    shear_complex weighted_gamma;
    shear_complex pivot_gamma;
    shear_complex pivot_weighted_gamma;
    shear_complex z2;
    shear_complex z4;
    shear_complex z6;

    if ((nodeptr)p == q)
        return;
    if (cballs_opt_read_mask(cmd)
        && Mask(q) != MASK_NODE_VALID)
        return;
    if (!accept_body(cmd, gd, p, q, &distance, dr))
        return;
    bin = shear_radial_bin(cmd, gd, distance);
    if (bin < 0)
        return;

    /* accept_body gives p-q; the polar phase used here is q-p. */
    cos_phi = -dr[0]/distance;
    sin_phi = -dr[1]/distance;
    phase = shear_make(cos_phi, sin_phi);
    gamma = shear_make(Gamma1(q), Gamma2(q));
    weight = Weight(q);
    weighted_gamma = shear_scale(gamma, weight);

    phase_power = shear_make(1.0, 0.0);
    work->g_ring_active[shear_ring_index(
        0, bin, work->ring_max, work->bins)] =
        shear_add(work->g_ring_active[
                      shear_ring_index(0, bin, work->ring_max, work->bins)],
                  weighted_gamma);
    work->w_ring_active[shear_ring_index(
        0, bin, work->ring_max, work->bins)].re +=
        weight;
    for (order = 1; order <= work->ring_max; order++) {
        size_t positive;
        size_t negative;
        shear_complex conjugate_power;

        phase_power = shear_mul(phase_power, phase);
        conjugate_power = shear_conj(phase_power);
        positive = shear_ring_index(order, bin, work->ring_max, work->bins);
        negative = shear_ring_index(-order, bin, work->ring_max, work->bins);
        work->g_ring_active[positive] = shear_add(
            work->g_ring_active[positive],
            shear_mul(weighted_gamma, phase_power));
        work->g_ring_active[negative] = shear_add(
            work->g_ring_active[negative],
            shear_mul(weighted_gamma, conjugate_power));
        work->w_ring_active[positive] = shear_add(
            work->w_ring_active[positive], shear_scale(phase_power, weight));
        work->w_ring_active[negative] = shear_add(
            work->w_ring_active[negative],
            shear_scale(conjugate_power, weight));
    }

    if (!work->collect_first_leg)
        return;

    z2 = shear_mul(phase, phase);
    z4 = shear_mul(z2, z2);
    z6 = shear_mul(z4, z2);
    work->diag_g6[bin] = shear_add(
        work->diag_g6[bin],
        shear_mul(shear_mul(weighted_gamma, weighted_gamma), shear_conj(z6)));
    work->diag_g2[bin] = shear_add(
        work->diag_g2[bin],
        shear_mul(shear_mul(weighted_gamma, weighted_gamma), shear_conj(z2)));
    work->diag_abs2[bin] = shear_add(
        work->diag_abs2[bin],
        shear_scale(shear_conj(z2), weight*weight*shear_abs2(gamma)));
    work->diag_w2[bin] += weight*weight;

    pivot_gamma = shear_make(Gamma1(p), Gamma2(p));
    pivot_weighted_gamma = shear_scale(pivot_gamma, Weight(p));
    work->xi_plus[bin] = shear_add(
        work->xi_plus[bin],
        shear_mul(pivot_weighted_gamma, shear_conj(weighted_gamma)));
    work->xi_minus[bin] = shear_add(
        work->xi_minus[bin],
        shear_mul(shear_mul(pivot_weighted_gamma, weighted_gamma),
                  shear_conj(z4)));
    work->xi_weight[bin] += Weight(p)*weight;
}

static void shear_walk_tree(shear_pivot_workspace *work, nodeptr q, real qsize)
{
    nodeptr child;

    if (q == NULL || q == (nodeptr)work->pivot)
        return;
    if (Type(q) == CELL) {
        if (cballs_opt_read_mask(work->cmd)
            && Mask(q) == MASK_NODE_MASKED)
            return;
        if (reject_cell(work->cmd, work->gd, (nodeptr)work->pivot, q, qsize))
            return;
        for (child = More(q); child != Next(q); child = Next(child))
            shear_walk_tree(work, child, qsize/2.0);
        return;
    }
    shear_accumulate_neighbor(work, q);
}

static int shear_validate_catalog(struct cmdline_data *cmd, bodyptr table,
                                  INTEGER count, const char *role)
{
    INTEGER i;
    int axis;

    if (table == NULL || count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: empty %s catalog", role);
        return FAILURE;
    }
    for (i = 0; i < count; i++) {
        bodyptr body_i = table + i;
        if (!isfinite(Pos(body_i)[0]) || !isfinite(Pos(body_i)[1])
            || !isfinite(Gamma1(body_i)) || !isfinite(Gamma2(body_i))
            || !isfinite(Weight(body_i)) || Weight(body_i) < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-shear-omp: invalid %s catalog row %" INTEGER_FMT,
                     role, i);
            return FAILURE;
        }
        for (axis = 2; axis < NDIM; axis++) {
            real reference = Pos(table)[axis];
            real tolerance = 64.0*DBL_EPSILON
                *(1.0 + rabs(reference) + rabs(Pos(body_i)[axis]));
            if (!isfinite(Pos(body_i)[axis])
                || rabs(Pos(body_i)[axis] - reference) > tolerance) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "octree-shear-omp: %s catalog is not flat in axis %d",
                         role, axis);
                return FAILURE;
            }
        }
    }
    return SUCCESS;
}

static int shear_validate_common_plane(struct cmdline_data *cmd,
                                       bodyptr first, bodyptr second,
                                       const char *first_role,
                                       const char *second_role)
{
    int axis;

    for (axis = 2; axis < NDIM; axis++) {
        real first_coordinate = Pos(first)[axis];
        real second_coordinate = Pos(second)[axis];
        real tolerance = 64.0*DBL_EPSILON
            *(1.0 + rabs(first_coordinate) + rabs(second_coordinate));

        if (rabs(first_coordinate - second_coordinate) > tolerance) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-shear-omp: %s and %s catalogs do not share "
                     "a tangent plane in axis %d",
                     first_role, second_role, axis);
            return FAILURE;
        }
    }
    return SUCCESS;
}

global int prepare_octree_shear_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btable, INTEGER *nbody)
{
    compute_vector minimum;
    compute_vector maximum;
    compute_vector center;
    bodyptr p;
    int cat1;
    int cat2;
    int cat3;
    int ifile;
    int axis;

    if (cmd == NULL)
        return FAILURE;
    if (gd == NULL || btable == NULL || nbody == NULL || gd->ninfiles <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid catalog preparation state");
        return FAILURE;
    }

    cat1 = gd->iCatalogs[0];
    cat2 = gd->ninfiles >= 2 ? gd->iCatalogs[1] : cat1;
    cat3 = gd->ninfiles >= 3 ? gd->iCatalogs[2] : cat2;
    if (cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles
        || cat3 < 0 || cat3 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid catalog selection %d:%d:%d",
                 cat1, cat2, cat3);
        return FAILURE;
    }
    if (shear_validate_catalog(cmd, btable[cat1], nbody[cat1], "pivot")
        == FAILURE)
        return FAILURE;
    if (cat2 != cat1
        && shear_validate_catalog(cmd, btable[cat2], nbody[cat2],
                                  "first neighbor") == FAILURE)
        return FAILURE;
    if (cat3 != cat1 && cat3 != cat2
        && shear_validate_catalog(cmd, btable[cat3], nbody[cat3],
                                  "second neighbor") == FAILURE)
        return FAILURE;
    if (cat2 != cat1
        && shear_validate_common_plane(cmd, btable[cat1], btable[cat2],
                                       "pivot", "first neighbor") == FAILURE)
        return FAILURE;
    if (cat3 != cat1
        && shear_validate_common_plane(cmd, btable[cat1], btable[cat3],
                                       "pivot", "second neighbor") == FAILURE)
        return FAILURE;

    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        if (btable[ifile] == NULL || nbody[ifile] <= 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "octree-shear-omp: empty input catalog %d", ifile);
            return FAILURE;
        }
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile]) {
            DO_COORD(axis) {
                real coordinate = Pos(p)[axis];
                if (ifile == 0 && p == btable[0]) {
                    minimum[axis] = coordinate;
                    maximum[axis] = coordinate;
                } else {
                    if (coordinate < minimum[axis])
                        minimum[axis] = coordinate;
                    if (coordinate > maximum[axis])
                        maximum[axis] = coordinate;
                }
            }
        }
    }
    DO_COORD(axis)
        center[axis] = 0.5*(minimum[axis] + maximum[axis]);
    for (ifile = 0; ifile < gd->ninfiles; ifile++)
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile])
            DO_COORD(axis)
                Pos(p)[axis] -= center[axis];

    return SUCCESS;
}

static int shear_validate(struct cmdline_data *cmd, struct global_data *gd,
                          bodyptr *btable, INTEGER *nbody,
                          INTEGER ipmin, INTEGER *ipmax,
                          int cat1, int cat2, int cat3)
{
    if (cmd == NULL)
        return FAILURE;
    if (gd == NULL || btable == NULL || nbody == NULL || ipmax == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid search state");
        return FAILURE;
    }
    if (cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles
        || cat3 < 0 || cat3 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid catalog selection %d:%d:%d",
                 cat1, cat2, cat3);
        return FAILURE;
    }
    if (cmd->usePeriodic) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: periodic geometry is not a flat-sky survey");
        return FAILURE;
    }
    if (cmd->sizeHistN <= 0 || cmd->sizeHistPhi <= 0
        || cmd->mChebyshev < 0
        || cmd->mChebyshev > (INT_MAX - 3)/4
        || !isfinite(cmd->rangeN)
        || !isfinite(cmd->rminHist) || cmd->rangeN <= cmd->rminHist
        || cmd->rminHist < 0.0 || !isfinite(gd->i_deltaR)
        || gd->i_deltaR <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid histogram domain");
        return FAILURE;
    }
    if (cmd->useLogHist && cmd->rminHist == 0.0
        && cmd->logHistBinsPD <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: logHistBinsPD must be positive");
        return FAILURE;
    }
    if (ipmin < 1 || ipmax[cat1] < ipmin || ipmax[cat1] > nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: invalid pivot interval");
        return FAILURE;
    }
    if (roottable[cat2] == NULL || roottable[cat3] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: a neighbor tree is not available");
        return FAILURE;
    }
    if (shear_validate_catalog(cmd, btable[cat1], nbody[cat1], "pivot")
        == FAILURE)
        return FAILURE;
    if (cat2 != cat1
        && shear_validate_catalog(cmd, btable[cat2], nbody[cat2],
                                  "first neighbor") == FAILURE)
        return FAILURE;
    if (cat3 != cat1 && cat3 != cat2
        && shear_validate_catalog(cmd, btable[cat3], nbody[cat3],
                                  "second neighbor") == FAILURE)
        return FAILURE;
    if (cat2 != cat1
        && shear_validate_common_plane(cmd, btable[cat1], btable[cat2],
                                       "pivot", "first neighbor") == FAILURE)
        return FAILURE;
    if (cat3 != cat1
        && shear_validate_common_plane(cmd, btable[cat1], btable[cat3],
                                       "pivot", "second neighbor") == FAILURE)
        return FAILURE;
    return SUCCESS;
}

static int shear_allocate_results(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  size_t bins, size_t gamma_count,
                                  size_t denominator_count,
                                  size_t angular_count)
{
#define SHEAR_ALLOC_RESULT(member, count)                                \
    do {                                                                 \
        if (shear_calloc(cmd, (void **)&gd->member, (count),             \
                         sizeof(*gd->member), #member) == FAILURE)       \
            goto fail;                                                   \
    } while (0)

    shear_clear_results(gd);
    SHEAR_ALLOC_RESULT(histShearXiPlusRe, bins);
    SHEAR_ALLOC_RESULT(histShearXiPlusIm, bins);
    SHEAR_ALLOC_RESULT(histShearXiMinusRe, bins);
    SHEAR_ALLOC_RESULT(histShearXiMinusIm, bins);
    SHEAR_ALLOC_RESULT(histShearXiWeight, bins);
    SHEAR_ALLOC_RESULT(histShearGammaNumeratorRe, gamma_count);
    SHEAR_ALLOC_RESULT(histShearGammaNumeratorIm, gamma_count);
    SHEAR_ALLOC_RESULT(histShearGammaMultipoleRe, gamma_count);
    SHEAR_ALLOC_RESULT(histShearGammaMultipoleIm, gamma_count);
    SHEAR_ALLOC_RESULT(histShearDenominatorRe, denominator_count);
    SHEAR_ALLOC_RESULT(histShearDenominatorIm, denominator_count);
    SHEAR_ALLOC_RESULT(histShearGammaRe, angular_count);
    SHEAR_ALLOC_RESULT(histShearGammaIm, angular_count);
#undef SHEAR_ALLOC_RESULT
    return SUCCESS;

fail:
#undef SHEAR_ALLOC_RESULT
    shear_clear_results(gd);
    return FAILURE;
}

static void shear_reduce_pivot(struct cmdline_data *cmd,
                               struct global_data *gd,
                               shear_pivot_workspace *work,
                               int nmax)
{
    int bin;
    int bin1;
    int bin2;
    int order;
    int multipoles = 2*nmax + 1;
    int denominator_offset = 2*nmax;
    shear_complex pivot_gamma = shear_make(Gamma1(work->pivot),
                                           Gamma2(work->pivot));
    shear_complex pivot_weighted_gamma = shear_scale(pivot_gamma,
                                                     Weight(work->pivot));

    for (bin = 0; bin < work->bins; bin++) {
        gd->histShearXiPlusRe[bin] += work->xi_plus[bin].re;
        gd->histShearXiPlusIm[bin] += work->xi_plus[bin].im;
        gd->histShearXiMinusRe[bin] += work->xi_minus[bin].re;
        gd->histShearXiMinusIm[bin] += work->xi_minus[bin].im;
        gd->histShearXiWeight[bin] += work->xi_weight[bin];
    }

    for (bin1 = 0; bin1 < work->bins; bin1++) {
        for (bin2 = 0; bin2 < work->bins; bin2++) {
            for (order = -2*nmax; order <= 2*nmax; order++) {
                shear_complex value = shear_mul(
                    work->w_ring_first[shear_ring_index(
                        order, bin1, work->ring_max, work->bins)],
                    shear_conj(work->w_ring_second[shear_ring_index(
                        order, bin2, work->ring_max, work->bins)]));
                size_t index = shear_denominator_index(
                    order + denominator_offset, bin1, bin2, work->bins);
                if (work->same_neighbor_catalog && bin1 == bin2)
                    value.re -= work->diag_w2[bin1];
                value = shear_scale(value, Weight(work->pivot));
                gd->histShearDenominatorRe[index] += value.re;
                gd->histShearDenominatorIm[index] += value.im;
            }

            for (order = -nmax; order <= nmax; order++) {
                shear_complex g_n3 = work->g_ring_first[shear_ring_index(
                    order - 3, bin1, work->ring_max, work->bins)];
                shear_complex g_mn3 = work->g_ring_second[shear_ring_index(
                    -order - 3, bin2, work->ring_max, work->bins)];
                shear_complex g_n1 = work->g_ring_first[shear_ring_index(
                    order - 1, bin1, work->ring_max, work->bins)];
                shear_complex g_mn1 = work->g_ring_second[shear_ring_index(
                    -order - 1, bin2, work->ring_max, work->bins)];
                shear_complex g_mn1_b1 = work->g_ring_first[shear_ring_index(
                    -order - 1, bin1, work->ring_max, work->bins)];
                shear_complex g_n1_b2 = work->g_ring_second[shear_ring_index(
                    order - 1, bin2, work->ring_max, work->bins)];
                shear_complex products[4];
                shear_complex pivots[4];
                int component;

                products[0] = shear_mul(g_n3, g_mn3);
                products[1] = shear_mul(g_n1, g_mn1);
                products[2] = shear_mul(shear_conj(g_mn1_b1), g_mn3);
                products[3] = shear_mul(g_n3, shear_conj(g_n1_b2));
                if (work->same_neighbor_catalog && bin1 == bin2) {
                    products[0] = shear_sub(products[0], work->diag_g6[bin1]);
                    products[1] = shear_sub(products[1], work->diag_g2[bin1]);
                    products[2] = shear_sub(products[2], work->diag_abs2[bin1]);
                    products[3] = shear_sub(products[3], work->diag_abs2[bin1]);
                }
                pivots[0] = pivot_weighted_gamma;
                pivots[1] = shear_conj(pivot_weighted_gamma);
                pivots[2] = pivot_weighted_gamma;
                pivots[3] = pivot_weighted_gamma;
                for (component = 0; component < 4; component++) {
                    shear_complex value = shear_scale(
                        shear_mul(pivots[component], products[component]), -1.0);
                    size_t index = shear_gamma_index(
                        component, order + nmax, bin1, bin2,
                        multipoles, work->bins);
                    gd->histShearGammaNumeratorRe[index] += value.re;
                    gd->histShearGammaNumeratorIm[index] += value.im;
                }
            }
        }
    }
    (void)cmd;
}

static int shear_solve_mode_coupling(struct cmdline_data *cmd,
                                     struct global_data *gd,
                                     int bins, int nmax)
{
    int multipoles = 2*nmax + 1;
    int denominator_offset = 2*nmax;
    size_t matrix_count;
    size_t rhs_count;
    shear_complex *matrix = NULL;
    shear_complex *rhs = NULL;
    int bin1;
    int bin2;

    if (shear_size_mul((size_t)multipoles, (size_t)multipoles,
                       &matrix_count) == FAILURE
        || shear_size_mul((size_t)multipoles, 4, &rhs_count) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: mode-coupling workspace overflow");
        return FAILURE;
    }
    if (shear_calloc(cmd, (void **)&matrix, matrix_count,
                     sizeof(*matrix), "mode-coupling matrix") == FAILURE
        || shear_calloc(cmd, (void **)&rhs, rhs_count,
                        sizeof(*rhs), "mode-coupling right-hand side")
           == FAILURE) {
        free(matrix);
        free(rhs);
        return FAILURE;
    }

    for (bin1 = 0; bin1 < bins; bin1++) {
        for (bin2 = 0; bin2 < bins; bin2++) {
            shear_complex n_zero;
            real matrix_scale = 0.0;
            real tolerance;
            int row;
            int column;
            bool singular = FALSE;
            size_t n0_index = shear_denominator_index(
                denominator_offset, bin1, bin2, bins);

            n_zero = shear_make(gd->histShearDenominatorRe[n0_index],
                                gd->histShearDenominatorIm[n0_index]);
            memset(matrix, 0, matrix_count*sizeof(*matrix));
            memset(rhs, 0, rhs_count*sizeof(*rhs));
            if (shear_abs2(n_zero) <= DBL_MIN)
                continue;

            for (row = 0; row < multipoles; row++) {
                int ell = row - nmax;
                int component;
                for (column = 0; column < multipoles; column++) {
                    int order = ell - (column - nmax);
                    size_t index = shear_denominator_index(
                        order + denominator_offset, bin1, bin2, bins);
                    shear_complex value = shear_make(
                        gd->histShearDenominatorRe[index],
                        gd->histShearDenominatorIm[index]);
                    matrix[(size_t)row*multipoles + column] =
                        shear_div(value, n_zero);
                    if (shear_abs2(value) > matrix_scale)
                        matrix_scale = shear_abs2(value);
                }
                for (component = 0; component < 4; component++) {
                    size_t index = shear_gamma_index(
                        component, row, bin1, bin2, multipoles, bins);
                    rhs[(size_t)row*4 + component] = shear_div(
                        shear_make(gd->histShearGammaNumeratorRe[index],
                                   gd->histShearGammaNumeratorIm[index]),
                        n_zero);
                }
            }
            tolerance = 128.0*DBL_EPSILON
                *(1.0 + rsqrt(matrix_scale/shear_abs2(n_zero)));

            for (column = 0; column < multipoles; column++) {
                int pivot = column;
                real pivot_norm = shear_abs2(
                    matrix[(size_t)column*multipoles + column]);
                int candidate;
                for (candidate = column + 1; candidate < multipoles;
                     candidate++) {
                    real candidate_norm = shear_abs2(
                        matrix[(size_t)candidate*multipoles + column]);
                    if (candidate_norm > pivot_norm) {
                        pivot = candidate;
                        pivot_norm = candidate_norm;
                    }
                }
                if (pivot_norm <= tolerance*tolerance) {
                    singular = TRUE;
                    break;
                }
                if (pivot != column) {
                    int k;
                    for (k = 0; k < multipoles; k++) {
                        shear_complex swap = matrix[
                            (size_t)column*multipoles + k];
                        matrix[(size_t)column*multipoles + k] = matrix[
                            (size_t)pivot*multipoles + k];
                        matrix[(size_t)pivot*multipoles + k] = swap;
                    }
                    for (k = 0; k < 4; k++) {
                        shear_complex swap = rhs[(size_t)column*4 + k];
                        rhs[(size_t)column*4 + k] = rhs[(size_t)pivot*4 + k];
                        rhs[(size_t)pivot*4 + k] = swap;
                    }
                }
                {
                    shear_complex divisor = matrix[
                        (size_t)column*multipoles + column];
                    int k;
                    for (k = 0; k < multipoles; k++)
                        matrix[(size_t)column*multipoles + k] = shear_div(
                            matrix[(size_t)column*multipoles + k], divisor);
                    for (k = 0; k < 4; k++)
                        rhs[(size_t)column*4 + k] = shear_div(
                            rhs[(size_t)column*4 + k], divisor);
                }
                for (row = 0; row < multipoles; row++) {
                    shear_complex factor;
                    int k;
                    if (row == column)
                        continue;
                    factor = matrix[(size_t)row*multipoles + column];
                    if (shear_abs2(factor) == 0.0)
                        continue;
                    for (k = 0; k < multipoles; k++)
                        matrix[(size_t)row*multipoles + k] = shear_sub(
                            matrix[(size_t)row*multipoles + k],
                            shear_mul(factor,
                                      matrix[(size_t)column*multipoles + k]));
                    for (k = 0; k < 4; k++)
                        rhs[(size_t)row*4 + k] = shear_sub(
                            rhs[(size_t)row*4 + k],
                            shear_mul(factor, rhs[(size_t)column*4 + k]));
                }
            }
            if (!singular) {
                int component;
                for (row = 0; row < multipoles; row++)
                    for (component = 0; component < 4; component++) {
                        size_t index = shear_gamma_index(
                            component, row, bin1, bin2, multipoles, bins);
                        gd->histShearGammaMultipoleRe[index] =
                            rhs[(size_t)row*4 + component].re;
                        gd->histShearGammaMultipoleIm[index] =
                            rhs[(size_t)row*4 + component].im;
                    }
            }
        }
    }

    free(matrix);
    free(rhs);
    return SUCCESS;
}

static void shear_reconstruct_angular(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      int bins, int nmax)
{
    int multipoles = 2*nmax + 1;
    real delta_phi = TWOPI/(real)cmd->sizeHistPhi;
    int component;
    int phi_bin;
    int bin1;
    int bin2;

    for (component = 0; component < 4; component++)
        for (phi_bin = 0; phi_bin < cmd->sizeHistPhi; phi_bin++) {
            real phi = -PI + (phi_bin + 0.5)*delta_phi;
            for (bin1 = 0; bin1 < bins; bin1++)
                for (bin2 = 0; bin2 < bins; bin2++) {
                    shear_complex sum = shear_make(0.0, 0.0);
                    int order;
                    for (order = -nmax; order <= nmax; order++) {
                        real window = order == 0
                            ? 1.0
                            : sin(0.5*order*delta_phi)
                              /(0.5*order*delta_phi);
                        shear_complex phase = shear_make(cos(order*phi),
                                                         sin(order*phi));
                        size_t index = shear_gamma_index(
                            component, order + nmax, bin1, bin2,
                            multipoles, bins);
                        shear_complex coefficient = shear_make(
                            gd->histShearGammaMultipoleRe[index],
                            gd->histShearGammaMultipoleIm[index]);
                        sum = shear_add(sum,
                                        shear_scale(shear_mul(coefficient, phase),
                                                    window/TWOPI));
                    }
                    {
                        size_t index = shear_angular_index(
                            component, phi_bin, bin1, bin2,
                            cmd->sizeHistPhi, bins);
                        gd->histShearGammaRe[index] = sum.re;
                        gd->histShearGammaIm[index] = sum.im;
                    }
                }
        }
}

static real shear_radial_center(struct cmdline_data *cmd,
                                struct global_data *gd, int bin)
{
    if (!cmd->useLogHist)
        return cmd->rminHist + (bin + 0.5)*gd->deltaR;
    if (cmd->rminHist == 0.0)
        return rpow(10.0,
                    ((bin + 0.5 - cmd->sizeHistN)/(real)cmd->logHistBinsPD)
                    + rlog10(cmd->rangeN));
    return cmd->rminHist*rpow(10.0, (bin + 0.5)/gd->i_deltaR);
}

static int shear_close_output(struct cmdline_data *cmd, FILE *stream,
                              const char *path)
{
    int failed = ferror(stream);

    if (fclose(stream) != 0)
        failed = TRUE;

    if (failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: cannot write %s", path);
        return FAILURE;
    }
    return SUCCESS;
}

static int shear_write_outputs(struct cmdline_data *cmd,
                               struct global_data *gd, int bins, int nmax)
{
    char path[MAXLENGTHOFFILES];
    FILE *stream;
    int bin;
    int component;
    int order;
    int bin1;
    int bin2;
    int phi_bin;
    int multipoles = 2*nmax + 1;

    if (format_checked(path, sizeof(path), "shear 2PCF output path",
                       "%s/histShearXi%s%s", cmd->rootDir,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        goto path_error;
    stream = fopen(path, "w");
    if (stream == NULL)
        goto open_error;
    fprintf(stream, "# convention: right-handed Cartesian flat sky; phi is "
                    "counterclockwise from +x\n");
    fprintf(stream, "# r xi_plus_re xi_plus_im xi_minus_re xi_minus_im weight\n");
    for (bin = 0; bin < bins; bin++)
        fprintf(stream, "%.17g %.17g %.17g %.17g %.17g %.17g\n",
                shear_radial_center(cmd, gd, bin),
                gd->histShearXiPlusRe[bin], gd->histShearXiPlusIm[bin],
                gd->histShearXiMinusRe[bin], gd->histShearXiMinusIm[bin],
                gd->histShearXiWeight[bin]);
    if (shear_close_output(cmd, stream, path) == FAILURE)
        return FAILURE;

    if (format_checked(path, sizeof(path), "shear multipole output path",
                       "%s/histShearGammaMultipoles%s%s", cmd->rootDir,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        goto path_error;
    stream = fopen(path, "w");
    if (stream == NULL)
        goto open_error;
    fprintf(stream, "# projection: Porth x-projection (Gamma^x)\n");
    fprintf(stream, "# catalogs: Z1=pivot Z2=radial_bin_1 Z3=radial_bin_2\n");
    fprintf(stream, "# component n radial_bin_1 radial_bin_2 gamma_re gamma_im upsilon_re upsilon_im\n");
    for (component = 0; component < 4; component++)
        for (order = -nmax; order <= nmax; order++)
            for (bin1 = 0; bin1 < bins; bin1++)
                for (bin2 = 0; bin2 < bins; bin2++) {
                    size_t index = shear_gamma_index(
                        component, order + nmax, bin1, bin2,
                        multipoles, bins);
                    fprintf(stream,
                            "%d %d %d %d %.17g %.17g %.17g %.17g\n",
                            component, order, bin1, bin2,
                            gd->histShearGammaMultipoleRe[index],
                            gd->histShearGammaMultipoleIm[index],
                            gd->histShearGammaNumeratorRe[index],
                            gd->histShearGammaNumeratorIm[index]);
                }
    if (shear_close_output(cmd, stream, path) == FAILURE)
        return FAILURE;

    if (format_checked(path, sizeof(path), "shear angular output path",
                       "%s/histShearGamma%s%s", cmd->rootDir,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        goto path_error;
    stream = fopen(path, "w");
    if (stream == NULL)
        goto open_error;
    fprintf(stream, "# projection: Porth x-projection (Gamma^x)\n");
    fprintf(stream, "# catalogs: Z1=pivot Z2=radial_bin_1 Z3=radial_bin_2\n");
    fprintf(stream, "# component phi_bin phi radial_bin_1 radial_bin_2 gamma_re gamma_im\n");
    for (component = 0; component < 4; component++)
        for (phi_bin = 0; phi_bin < cmd->sizeHistPhi; phi_bin++) {
            real phi = -PI + (phi_bin + 0.5)*TWOPI/cmd->sizeHistPhi;
            for (bin1 = 0; bin1 < bins; bin1++)
                for (bin2 = 0; bin2 < bins; bin2++) {
                    size_t index = shear_angular_index(
                        component, phi_bin, bin1, bin2,
                        cmd->sizeHistPhi, bins);
                    fprintf(stream, "%d %d %.17g %d %d %.17g %.17g\n",
                            component, phi_bin, phi, bin1, bin2,
                            gd->histShearGammaRe[index],
                            gd->histShearGammaIm[index]);
                }
        }
    return shear_close_output(cmd, stream, path);

path_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-shear-omp: output path is too long");
    return FAILURE;
open_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "octree-shear-omp: cannot open %s: %s", path, strerror(errno));
    return FAILURE;
}

global int searchcalc_octree_shear_omp(struct cmdline_data *cmd,
                                       struct global_data *gd,
                                       bodyptr *btable, INTEGER *nbody,
                                       INTEGER ipmin, INTEGER *ipmax,
                                       int cat1, int cat2, int cat3)
{
    real cpu_start = CPUTIME;
    int bins;
    int nmax;
    int multipoles;
    int denominator_multipoles;
    int ring_max;
    int threads = 1;
    size_t bins_squared;
    size_t gamma_count;
    size_t denominator_count;
    size_t angular_count;
    size_t ring_count;
    size_t threaded_ring_count;
    size_t threaded_bin_count;
    shear_complex *g_ring_first_all = NULL;
    shear_complex *w_ring_first_all = NULL;
    shear_complex *g_ring_second_all = NULL;
    shear_complex *w_ring_second_all = NULL;
    shear_complex *diag_g6_all = NULL;
    shear_complex *diag_g2_all = NULL;
    shear_complex *diag_abs2_all = NULL;
    real *diag_w2_all = NULL;
    shear_complex *xi_plus_all = NULL;
    shear_complex *xi_minus_all = NULL;
    real *xi_weight_all = NULL;
    INTEGER ip;

    if (cmd == NULL)
        return FAILURE;
    cmd->error_message[0] = '\0';
    if (shear_validate(cmd, gd, btable, nbody, ipmin, ipmax,
                       cat1, cat2, cat3)
        == FAILURE)
        return FAILURE;

    bins = cmd->sizeHistN;
    nmax = cmd->mChebyshev;
    multipoles = 2*nmax + 1;
    denominator_multipoles = 4*nmax + 1;
    ring_max = 2*nmax > nmax + 3 ? 2*nmax : nmax + 3;
#ifdef OPENMPCODE
    threads = omp_get_max_threads();
#endif
    if (threads < 1)
        threads = 1;

    if (shear_size_mul((size_t)bins, (size_t)bins, &bins_squared) == FAILURE
        || shear_size_mul(4*(size_t)multipoles, bins_squared,
                          &gamma_count) == FAILURE
        || shear_size_mul((size_t)denominator_multipoles, bins_squared,
                          &denominator_count) == FAILURE
        || shear_size_mul(4*(size_t)cmd->sizeHistPhi, bins_squared,
                          &angular_count) == FAILURE
        || shear_size_mul((size_t)(2*ring_max + 1), (size_t)bins,
                          &ring_count) == FAILURE
        || shear_size_mul((size_t)threads, ring_count,
                          &threaded_ring_count) == FAILURE
        || shear_size_mul((size_t)threads, (size_t)bins,
                          &threaded_bin_count) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "octree-shear-omp: histogram size overflow");
        return FAILURE;
    }
    if (shear_allocate_results(cmd, gd, bins, gamma_count,
                               denominator_count, angular_count) == FAILURE)
        return FAILURE;

#define SHEAR_ALLOC_WORK(pointer, count, label)                          \
    do {                                                                 \
        if (shear_calloc(cmd, (void **)&(pointer), (count),              \
                         sizeof(*(pointer)), (label)) == FAILURE)        \
            goto fail;                                                   \
    } while (0)
    SHEAR_ALLOC_WORK(g_ring_first_all, threaded_ring_count,
                     "thread first G rings");
    SHEAR_ALLOC_WORK(w_ring_first_all, threaded_ring_count,
                     "thread first W rings");
    SHEAR_ALLOC_WORK(g_ring_second_all, threaded_ring_count,
                     "thread second G rings");
    SHEAR_ALLOC_WORK(w_ring_second_all, threaded_ring_count,
                     "thread second W rings");
    SHEAR_ALLOC_WORK(diag_g6_all, threaded_bin_count, "thread G6 diagonal");
    SHEAR_ALLOC_WORK(diag_g2_all, threaded_bin_count, "thread G2 diagonal");
    SHEAR_ALLOC_WORK(diag_abs2_all, threaded_bin_count,
                     "thread absolute-shear diagonal");
    SHEAR_ALLOC_WORK(diag_w2_all, threaded_bin_count, "thread weight diagonal");
    SHEAR_ALLOC_WORK(xi_plus_all, threaded_bin_count, "thread xi-plus");
    SHEAR_ALLOC_WORK(xi_minus_all, threaded_bin_count, "thread xi-minus");
    SHEAR_ALLOC_WORK(xi_weight_all, threaded_bin_count, "thread xi weight");
#undef SHEAR_ALLOC_WORK

#pragma omp parallel private(ip)
    {
        int thread_id = 0;
        shear_pivot_workspace work;
#ifdef OPENMPCODE
        thread_id = omp_get_thread_num();
#endif
        work.cmd = cmd;
        work.gd = gd;
        work.bins = bins;
        work.ring_max = ring_max;
        work.g_ring_first = g_ring_first_all + (size_t)thread_id*ring_count;
        work.w_ring_first = w_ring_first_all + (size_t)thread_id*ring_count;
        work.g_ring_second = g_ring_second_all + (size_t)thread_id*ring_count;
        work.w_ring_second = w_ring_second_all + (size_t)thread_id*ring_count;
        work.diag_g6 = diag_g6_all + (size_t)thread_id*bins;
        work.diag_g2 = diag_g2_all + (size_t)thread_id*bins;
        work.diag_abs2 = diag_abs2_all + (size_t)thread_id*bins;
        work.diag_w2 = diag_w2_all + (size_t)thread_id*bins;
        work.xi_plus = xi_plus_all + (size_t)thread_id*bins;
        work.xi_minus = xi_minus_all + (size_t)thread_id*bins;
        work.xi_weight = xi_weight_all + (size_t)thread_id*bins;
        work.same_neighbor_catalog = cat2 == cat3;

#pragma omp for schedule(static,1) ordered
        for (ip = ipmin - 1; ip < ipmax[cat1]; ip++) {
            work.pivot = btable[cat1] + ip;
            memset(work.g_ring_first, 0,
                   ring_count*sizeof(*work.g_ring_first));
            memset(work.w_ring_first, 0,
                   ring_count*sizeof(*work.w_ring_first));
            memset(work.g_ring_second, 0,
                   ring_count*sizeof(*work.g_ring_second));
            memset(work.w_ring_second, 0,
                   ring_count*sizeof(*work.w_ring_second));
            memset(work.diag_g6, 0, (size_t)bins*sizeof(*work.diag_g6));
            memset(work.diag_g2, 0, (size_t)bins*sizeof(*work.diag_g2));
            memset(work.diag_abs2, 0,
                   (size_t)bins*sizeof(*work.diag_abs2));
            memset(work.diag_w2, 0, (size_t)bins*sizeof(*work.diag_w2));
            memset(work.xi_plus, 0, (size_t)bins*sizeof(*work.xi_plus));
            memset(work.xi_minus, 0, (size_t)bins*sizeof(*work.xi_minus));
            memset(work.xi_weight, 0, (size_t)bins*sizeof(*work.xi_weight));

            if (Update(work.pivot)
                && (!cballs_opt_read_mask(cmd)
                    || Mask(work.pivot) == MASK_NODE_VALID)) {
                work.g_ring_active = work.g_ring_first;
                work.w_ring_active = work.w_ring_first;
                work.collect_first_leg = TRUE;
                shear_walk_tree(&work, (nodeptr)roottable[cat2],
                                gd->rSizeTable[cat2]);
                if (cat3 == cat2) {
                    memcpy(work.g_ring_second, work.g_ring_first,
                           ring_count*sizeof(*work.g_ring_second));
                    memcpy(work.w_ring_second, work.w_ring_first,
                           ring_count*sizeof(*work.w_ring_second));
                } else {
                    work.g_ring_active = work.g_ring_second;
                    work.w_ring_active = work.w_ring_second;
                    work.collect_first_leg = FALSE;
                    shear_walk_tree(&work, (nodeptr)roottable[cat3],
                                    gd->rSizeTable[cat3]);
                }
            }

#pragma omp ordered
            shear_reduce_pivot(cmd, gd, &work, nmax);
        }
    }

    for (ip = 0; ip < bins; ip++) {
        real denominator = gd->histShearXiWeight[ip];
        gd->histShearXiPlusRe[ip] = cballs_normalize_or_zero(
            gd->histShearXiPlusRe[ip], denominator);
        gd->histShearXiPlusIm[ip] = cballs_normalize_or_zero(
            gd->histShearXiPlusIm[ip], denominator);
        gd->histShearXiMinusRe[ip] = cballs_normalize_or_zero(
            gd->histShearXiMinusRe[ip], denominator);
        gd->histShearXiMinusIm[ip] = cballs_normalize_or_zero(
            gd->histShearXiMinusIm[ip], denominator);
    }
    if (shear_solve_mode_coupling(cmd, gd, bins, nmax) == FAILURE)
        goto fail;
    shear_reconstruct_angular(cmd, gd, bins, nmax);
    gd->shearMultipoleMax = nmax;
    gd->shearAngularBins = cmd->sizeHistPhi;
    gd->cpusearch = CPUTIME - cpu_start;

    free(g_ring_first_all);
    free(w_ring_first_all);
    free(g_ring_second_all);
    free(w_ring_second_all);
    free(diag_g6_all);
    free(diag_g2_all);
    free(diag_abs2_all);
    free(diag_w2_all);
    free(xi_plus_all);
    free(xi_minus_all);
    free(xi_weight_all);
    if (!cballs_opt_no_out_hist(cmd)
        && shear_write_outputs(cmd, gd, bins, nmax) == FAILURE)
        return FAILURE;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "octree-shear-omp: completed in %g %s\n",
                        gd->cpusearch, PRNUNITOFTIMEUSED);
    return SUCCESS;

fail:
    free(g_ring_first_all);
    free(w_ring_first_all);
    free(g_ring_second_all);
    free(w_ring_second_all);
    free(diag_g6_all);
    free(diag_g2_all);
    free(diag_abs2_all);
    free(diag_w2_all);
    free(xi_plus_all);
    free(xi_minus_all);
    free(xi_weight_all);
    shear_clear_results(gd);
    return FAILURE;
}

#endif /* NDIM >= 2 */
