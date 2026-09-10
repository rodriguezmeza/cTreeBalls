/*
 * Weak-lensing shear 2PCF and 3PCF estimator.
 *
 * The exact path descends accepted cells to bodies.  The production path may
 * instead use a cell's weighted spin-2 first and second moments when its full
 * radial extent occupies one bin and its bearing error satisfies the same
 * order-dependent bound used by dual-node-style LogMultipole walks.
 *
 * OCTREE_SHEAR_SPHERICAL specializes the geometry for unit-vector catalogs.
 * Neighbor shears are parallel transported into the pivot's east/north frame
 * before the flat estimator algebra is applied.  Its cells store moments in
 * the tangent basis at their normalized centers, so accepted-node work never
 * combines spin-2 values expressed in different frames.
 */

#include "globaldefs.h"
#ifdef BALLS4SCANLEV
#include "octree_scan_frontier.h"
#endif
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdarg.h>

#ifndef SHEAR_ENGINE_NAME
#ifdef OCTREE_SHEAR_SPHERICAL
#define SHEAR_ENGINE_NAME "octree-shear-sphere-omp"
#define SHEAR_REQUIRED_DIMENSION 3
#else
#define SHEAR_ENGINE_NAME "octree-shear-omp"
#define SHEAR_REQUIRED_DIMENSION 2
#endif
#elif defined(OCTREE_SHEAR_SPHERICAL)
#define SHEAR_REQUIRED_DIMENSION 3
#else
#define SHEAR_REQUIRED_DIMENSION 2
#endif

#if NDIM < SHEAR_REQUIRED_DIMENSION

global int prepare_octree_shear_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btable, INTEGER *nbody)
{
    (void)gd;
    (void)btable;
    (void)nbody;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             SHEAR_ENGINE_NAME " requires NDIM >= %d",
             SHEAR_REQUIRED_DIMENSION);
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
             SHEAR_ENGINE_NAME " requires NDIM >= %d",
             SHEAR_REQUIRED_DIMENSION);
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
    bool allow_cells;
#ifdef OCTREE_SHEAR_SPHERICAL
    compute_vector pivot_unit;
    compute_vector pivot_east;
    compute_vector pivot_north;
#endif
} shear_pivot_workspace;

typedef struct {
    real *xi_plus_re;
    real *xi_plus_im;
    real *xi_minus_re;
    real *xi_minus_im;
    real *xi_weight;
    real *gamma_re;
    real *gamma_im;
    real *denominator_re;
    real *denominator_im;
} shear_result_accumulator;

#ifndef SHEAR_OMP_PIVOT_BLOCK_SIZE
#define SHEAR_OMP_PIVOT_BLOCK_SIZE 32
#endif
#if SHEAR_OMP_PIVOT_BLOCK_SIZE < 1
#error SHEAR_OMP_PIVOT_BLOCK_SIZE must be positive
#endif

/* Every two members claimed by one representative can be separated by as
 * much as 2*rsmooth.  Keep that complete group below the measured domain so
 * the representative moment never contains a resolved self/pair term. */
#define SHEAR_SMOOTH_MAX_RMIN_FRACTION 0.5

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

#ifdef BALLS4SCANLEV
static int shear_size_add(size_t a, size_t b, size_t *result)
{
    if (b > SIZE_MAX-a)
        return FAILURE;
    *result = a+b;
    return SUCCESS;
}
#endif

static int shear_calloc(struct cmdline_data *cmd, void **pointer,
                        size_t count, size_t item_size, const char *label)
{
    size_t bytes;

    *pointer = NULL;
    if (shear_size_mul(count, item_size, &bytes) == FAILURE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": allocation size overflow for %s", label);
        return FAILURE;
    }
    if (bytes == 0)
        bytes = 1;
    *pointer = calloc(1, bytes);
    if (*pointer == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": cannot allocate %zu bytes for %s",
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

static shear_complex shear_pivot_weighted_gamma(
        const shear_pivot_workspace *work)
{
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(work->cmd))
        return shear_make(Gamma1Rmin(work->pivot), Gamma2Rmin(work->pivot));
#endif
    return shear_scale(shear_make(Gamma1(work->pivot), Gamma2(work->pivot)),
                       Weight(work->pivot));
}

static real shear_pivot_weight(const shear_pivot_workspace *work)
{
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(work->cmd))
        return WeightRmin(work->pivot);
#endif
    return Weight(work->pivot);
}

#ifdef OCTREE_SHEAR_SPHERICAL
static real shear_dot3(const real *a, const real *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static bool shear_unit3(const real *position, compute_vector unit)
{
    real norm2 = shear_dot3(position, position);
    real inverse_norm;

    if (!(norm2 > 0.0) || !isfinite(norm2))
        return FALSE;
    inverse_norm = 1.0/rsqrt(norm2);
    unit[0] = position[0]*inverse_norm;
    unit[1] = position[1]*inverse_norm;
    unit[2] = position[2]*inverse_norm;
    return isfinite(unit[0]) && isfinite(unit[1]) && isfinite(unit[2]);
}

static bool shear_spherical_basis(const real *unit, compute_vector east,
                                  compute_vector north)
{
    real equatorial_norm = rsqrt(unit[0]*unit[0] + unit[1]*unit[1]);

    if (equatorial_norm > 64.0*DBL_EPSILON) {
        east[0] = -unit[1]/equatorial_norm;
        east[1] = unit[0]/equatorial_norm;
        east[2] = 0.0;
    } else {
        /* Longitude is undefined at a pole; select a stable local meridian. */
        east[0] = 1.0;
        east[1] = 0.0;
        east[2] = 0.0;
    }
    north[0] = unit[1]*east[2] - unit[2]*east[1];
    north[1] = unit[2]*east[0] - unit[0]*east[2];
    north[2] = unit[0]*east[1] - unit[1]*east[0];
    return isfinite(north[0]) && isfinite(north[1]) && isfinite(north[2]);
}

static bool shear_prepare_pivot_geometry(shear_pivot_workspace *work)
{
    return shear_unit3(Pos(work->pivot), work->pivot_unit)
        && shear_spherical_basis(work->pivot_unit, work->pivot_east,
                                 work->pivot_north);
}

static bool shear_spherical_phase(const shear_pivot_workspace *work,
                                  const real *position, shear_complex *phase)
{
    compute_vector tangent;
    real projection = shear_dot3(work->pivot_unit, position);
    real tangent_norm;

    tangent[0] = position[0] - projection*work->pivot_unit[0];
    tangent[1] = position[1] - projection*work->pivot_unit[1];
    tangent[2] = position[2] - projection*work->pivot_unit[2];
    tangent_norm = rsqrt(shear_dot3(tangent, tangent));
    if (!(tangent_norm > 64.0*DBL_EPSILON) || !isfinite(tangent_norm))
        return FALSE;
    phase->re = shear_dot3(tangent, work->pivot_east)/tangent_norm;
    phase->im = shear_dot3(tangent, work->pivot_north)/tangent_norm;
    return isfinite(phase->re) && isfinite(phase->im);
}

static bool shear_transport_rotation_to_pivot(
        const real *pivot_unit, const real *pivot_east,
        const real *pivot_north, const real *neighbor_position,
        shear_complex *rotation)
{
    compute_vector q_unit;
    compute_vector q_east;
    compute_vector q_north;
    compute_vector transported_east;
    compute_vector transported_north;
    real dot;
    real denominator;
    real c;
    real s;
    real orientation_norm;
    int axis;

    if (!shear_unit3(neighbor_position, q_unit)
        || !shear_spherical_basis(q_unit, q_east, q_north))
        return FALSE;
    (void)q_north;
    dot = shear_dot3(pivot_unit, q_unit);
    denominator = 1.0 + dot;
    if (!(denominator > 64.0*DBL_EPSILON) || !isfinite(denominator))
        return FALSE;
    for (axis = 0; axis < 3; axis++) {
        transported_east[axis] = pivot_east[axis]
            - shear_dot3(pivot_east, q_unit)/denominator
              *(pivot_unit[axis] + q_unit[axis]);
        transported_north[axis] = pivot_north[axis]
            - shear_dot3(pivot_north, q_unit)/denominator
              *(pivot_unit[axis] + q_unit[axis]);
    }
    c = shear_dot3(q_east, transported_east);
    s = shear_dot3(q_east, transported_north);
    orientation_norm = rsqrt(c*c + s*s);
    if (!(orientation_norm > 64.0*DBL_EPSILON)
        || !isfinite(orientation_norm))
        return FALSE;
    c /= orientation_norm;
    s /= orientation_norm;
    *rotation = shear_make(c*c - s*s, 2.0*c*s);
    return isfinite(rotation->re) && isfinite(rotation->im);
}

static bool shear_transport_to_pivot(const shear_pivot_workspace *work,
                                     nodeptr q, shear_complex input,
                                     shear_complex *transported)
{
    shear_complex rotation;

    if (!shear_transport_rotation_to_pivot(
            work->pivot_unit, work->pivot_east, work->pivot_north,
            Pos(q), &rotation))
        return FALSE;
    *transported = shear_mul(input, rotation);
    return isfinite(transported->re) && isfinite(transported->im);
}

#ifdef SMOOTHPIVOT
typedef struct {
    bodyptr pivot;
    bool frame_ready;
    compute_vector unit;
    compute_vector east;
    compute_vector north;
} shear_spherical_smooth_context;

static int shear_spherical_smooth_claim(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr pivot, bodyptr claimed, void *opaque_context)
{
    shear_spherical_smooth_context *context = opaque_context;
    shear_complex rotation;
    shear_complex weighted_gamma;
    shear_complex transported;

    (void)gd;
    if (context == NULL || pivot == NULL || claimed == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid spherical smoothing claim");
        return FAILURE;
    }
    if (!context->frame_ready || context->pivot != pivot) {
        context->pivot = pivot;
        context->frame_ready = shear_unit3(Pos(pivot), context->unit)
            && shear_spherical_basis(context->unit, context->east,
                                     context->north);
        if (!context->frame_ready) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME
                     ": cannot construct tangent frame for a smooth pivot");
            return FAILURE;
        }
    }
    if (!shear_transport_rotation_to_pivot(
            context->unit, context->east, context->north,
            Pos(claimed), &rotation)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME
                 ": spherical smooth-pivot transport is undefined");
        return FAILURE;
    }
    weighted_gamma = shear_scale(
        shear_make(Gamma1(claimed), Gamma2(claimed)), Weight(claimed));
    transported = shear_mul(weighted_gamma, rotation);
    if (!isfinite(transported.re) || !isfinite(transported.im)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME
                 ": non-finite spherical smooth-pivot shear");
        return FAILURE;
    }
    Gamma1Rmin(pivot) += transported.re;
    Gamma2Rmin(pivot) += transported.im;
    return SUCCESS;
}
#endif
#else
static bool shear_prepare_pivot_geometry(shear_pivot_workspace *work)
{
    (void)work;
    return TRUE;
}
#endif

static bool shear_center_geometry(const shear_pivot_workspace *work, nodeptr q,
                                  real *distance, shear_complex *phase,
                                  int *bin)
{
    compute_vector dr;
    real distance2;

    if (q == NULL || q == (nodeptr)work->pivot)
        return FALSE;
    DOTPSUBV(distance2, dr, Pos(work->pivot), Pos(q));
    if (!(distance2 > 0.0) || !isfinite(distance2))
        return FALSE;
    *distance = rsqrt(distance2);
    *bin = shear_radial_bin(work->cmd, work->gd, *distance);
    if (*bin < 0)
        return FALSE;
#ifdef OCTREE_SHEAR_SPHERICAL
    if (!shear_spherical_phase(work, Pos(q), phase))
        return FALSE;
#else
    /* dr is pivot-neighbor; the flat-sky phase is neighbor-pivot. */
    *phase = shear_make(-dr[0]/(*distance), -dr[1]/(*distance));
#endif
    return isfinite(phase->re) && isfinite(phase->im);
}

static real shear_physical_cell_radius(const shear_pivot_workspace *work,
                                       nodeptr q, real qsize)
{
    real radius;

#ifdef OCTREE_SHEAR_SPHERICAL
    if (work->cmd->theta > 0.0
        && isfinite((real)Radius(q)) && Radius(q) >= 0.0)
        radius = (real)Radius(q)*work->cmd->theta;
    else
        radius = 0.5*qsize*rsqrt((real)NDIM);
#else
    if (work->cmd->theta > 0.0 && isfinite((real)Radius(q)))
        radius = (real)Radius(q)*work->cmd->theta;
    else
        radius = 0.5*qsize*rsqrt((real)NDIM);
#endif
    (void)work;
    return isfinite(radius) && radius >= 0.0 ? radius : INFINITY;
}

#ifdef OCTREE_SHEAR_SPHERICAL
static bool shear_spherical_node_frame(
        struct cmdline_data *cmd, nodeptr q, real qsize,
        compute_vector unit, compute_vector east, compute_vector north,
        real *effective_radius, real *angular_radius)
{
    shear_pivot_workspace radius_work;
    real norm2;
    real norm;
    real radius = 0.0;

    if (q == NULL || !shear_unit3(Pos(q), unit)
        || !shear_spherical_basis(unit, east, north))
        return FALSE;
    norm2 = shear_dot3(Pos(q), Pos(q));
    norm = rsqrt(norm2);
    if (Type(q) == CELL) {
        memset(&radius_work, 0, sizeof(radius_work));
        radius_work.cmd = cmd;
        radius = shear_physical_cell_radius(&radius_work, q, qsize)
            + rabs(norm - 1.0);
    }
    if (!isfinite(radius) || radius < 0.0)
        return FALSE;
    if (effective_radius != NULL)
        *effective_radius = radius;
    if (angular_radius != NULL)
        *angular_radius = 2.0*rasin(MIN(1.0, 0.5*radius));
    return TRUE;
}

/* Re-express every spherical cell's spin moments in the tangent basis at
 * that cell's normalized center. Child moments are transported upward once;
 * ShearTransportError records a conservative spin-phase bound for the
 * hierarchical transport path and is enforced by accepted-node searches. */
static int shear_prepare_spherical_cell_moments(
        struct cmdline_data *cmd, nodeptr q, real qsize)
{
    compute_vector parent_unit;
    compute_vector parent_east;
    compute_vector parent_north;
    shear_complex gamma_sum = shear_make(0.0, 0.0);
    shear_complex gamma2_sum = shear_make(0.0, 0.0);
    real weight_sum = 0.0;
    real weight2_sum = 0.0;
    real gamma_abs2_sum = 0.0;
    real worst_error = 0.0;
    int child_index;

    if (q == NULL)
        return SUCCESS;
    if (Type(q) != CELL) {
        ShearTransportError(q) = 0.0;
        return SUCCESS;
    }
    for (child_index = 0; child_index < NSUB; child_index++) {
        nodeptr child = Subp(q)[child_index];

        if (child != NULL && Type(child) == CELL
            && shear_prepare_spherical_cell_moments(
                   cmd, child, 0.5*qsize) == FAILURE)
            return FAILURE;
    }
    if (!shear_spherical_node_frame(
            cmd, q, qsize, parent_unit, parent_east, parent_north,
            NULL, NULL)) {
        ShearWeightSum(q) = 0.0;
        Gamma1(q) = Gamma2(q) = 0.0;
        ShearGamma2Re(q) = ShearGamma2Im(q) = 0.0;
        ShearGammaAbs2(q) = ShearWeight2(q) = 0.0;
        ShearTransportError(q) = INFINITY;
        return SUCCESS;
    }

    for (child_index = 0; child_index < NSUB; child_index++) {
        nodeptr child = Subp(q)[child_index];
        compute_vector child_unit;
        compute_vector child_east;
        compute_vector child_north;
        shear_complex rotation;
        shear_complex gamma;
        shear_complex gamma2;
        real weight;
        real weight2;
        real gamma_abs2;
        real child_radius = 0.0;
        real child_angle = 0.0;
        real child_error = 0.0;
        real separation;
        real dot;

        if (child == NULL
            || (cballs_opt_read_mask(cmd)
                && Mask(child) == MASK_NODE_MASKED))
            continue;
        if (!shear_spherical_node_frame(
                cmd, child, 0.5*qsize, child_unit, child_east, child_north,
                &child_radius, &child_angle)) {
            ShearTransportError(q) = INFINITY;
            ShearWeightSum(q) = 0.0;
            return SUCCESS;
        }
        (void)child_east;
        (void)child_north;
        (void)child_radius;
        if (Type(child) == CELL) {
            if (!isfinite(ShearTransportError(child))) {
                ShearTransportError(q) = INFINITY;
                ShearWeightSum(q) = 0.0;
                return SUCCESS;
            }
            weight = ShearWeightSum(child);
            gamma = shear_make(Gamma1(child), Gamma2(child));
            gamma2 = shear_make(ShearGamma2Re(child),
                                ShearGamma2Im(child));
            weight2 = ShearWeight2(child);
            gamma_abs2 = ShearGammaAbs2(child);
            child_error = ShearTransportError(child);
        } else {
            weight = Weight(child);
            gamma = shear_make(Gamma1(child), Gamma2(child));
            gamma2 = shear_scale(shear_mul(gamma, gamma), weight*weight);
            weight2 = weight*weight;
            gamma_abs2 = weight2*shear_abs2(gamma);
        }
        if (!(weight >= 0.0) || !isfinite(weight)
            || !shear_transport_rotation_to_pivot(
                   parent_unit, parent_east, parent_north,
                   Pos(child), &rotation)) {
            ShearTransportError(q) = INFINITY;
            ShearWeightSum(q) = 0.0;
            return SUCCESS;
        }
        gamma_sum = shear_add(
            gamma_sum, shear_scale(shear_mul(gamma, rotation), weight));
        gamma2_sum = shear_add(
            gamma2_sum,
            shear_mul(gamma2, shear_mul(rotation, rotation)));
        weight_sum += weight;
        weight2_sum += weight2;
        gamma_abs2_sum += gamma_abs2;
        dot = MAX(-1.0, MIN(1.0, shear_dot3(parent_unit, child_unit)));
        separation = racos(dot);
        worst_error = MAX(worst_error,
                          child_error + 2.0*child_angle*separation);
    }
    ShearWeightSum(q) = weight_sum;
    if (weight_sum > 0.0) {
        Gamma1(q) = gamma_sum.re/weight_sum;
        Gamma2(q) = gamma_sum.im/weight_sum;
    } else {
        Gamma1(q) = Gamma2(q) = 0.0;
    }
    ShearGamma2Re(q) = gamma2_sum.re;
    ShearGamma2Im(q) = gamma2_sum.im;
    ShearWeight2(q) = weight2_sum;
    ShearGammaAbs2(q) = gamma_abs2_sum;
    ShearTransportError(q) = worst_error;
    return SUCCESS;
}
#endif

static bool shear_cell_geometry(shear_pivot_workspace *work, nodeptr q,
                                real qsize, real *distance,
                                shear_complex *phase, int *bin,
                                real *radius)
{
    compute_vector dr;
    real distance2;
    real lower;
    real upper;
    real angular_tolerance;
    real max_ratio;
    int lower_bin;
    int upper_bin;

    *radius = shear_physical_cell_radius(work, q, qsize);
#ifdef OCTREE_SHEAR_SPHERICAL
    {
        compute_vector center_unit;
        real center_norm2 = shear_dot3(Pos(q), Pos(q));
        real center_norm;

        if (!shear_unit3(Pos(q), center_unit)) {
            *distance = 0.0;
            *bin = -1;
            *phase = shear_make(0.0, 0.0);
            return FALSE;
        }
        center_norm = rsqrt(center_norm2);
        *radius += rabs(center_norm - 1.0);
        DOTPSUBV(distance2, dr, work->pivot_unit, center_unit);
    }
#else
    DOTPSUBV(distance2, dr, Pos(work->pivot), Pos(q));
#endif
    if (!(distance2 >= 0.0) || !isfinite(distance2)) {
        *distance = 0.0;
        *bin = -1;
        *phase = shear_make(0.0, 0.0);
        return FALSE;
    }
    *distance = rsqrt(distance2);
    *bin = shear_radial_bin(work->cmd, work->gd, *distance);
#ifdef OCTREE_SHEAR_SPHERICAL
    if (*distance > 0.0 && shear_spherical_phase(work, Pos(q), phase)) {
        /* phase was set in the pivot tangent basis */
    } else {
        *phase = shear_make(0.0, 0.0);
    }
#else
    *phase = *distance > 0.0
        ? shear_make(-dr[0]/(*distance), -dr[1]/(*distance))
        : shear_make(0.0, 0.0);
#endif
    if (*distance + *radius <= work->cmd->rminHist
        || *distance - *radius >= work->cmd->rangeN)
        return FALSE;
    if (*bin < 0 || !isfinite(phase->re) || !isfinite(phase->im))
        return FALSE;
    if (!work->allow_cells || !(work->cmd->theta > 0.0)
        || Nb(q) <= 0 || *distance <= *radius
        || (cballs_opt_read_mask(work->cmd)
            && Mask(q) != MASK_NODE_VALID))
        return FALSE;

    lower = *distance - *radius;
    upper = *distance + *radius;
    if (!(lower > work->cmd->rminHist && upper < work->cmd->rangeN))
        return FALSE;
    lower_bin = shear_radial_bin(work->cmd, work->gd, lower);
    upper_bin = shear_radial_bin(work->cmd, work->gd, upper);
    if (lower_bin != *bin || upper_bin != *bin)
        return FALSE;

    angular_tolerance = MIN(0.5*PI,
        work->cmd->theta*PI/(2.0*(real)work->ring_max + 1.0));
    max_ratio = angular_tolerance >= 0.5*PI
        ? 1.0 : rsin(MAX(0.0, angular_tolerance));
#ifdef OCTREE_SHEAR_SPHERICAL
    {
        const real cell_angle = 2.0*rasin(MIN(1.0, 0.5*(*radius)));
        const real separation = 2.0*rasin(MIN(1.0, 0.5*(*distance)));
        const real transport_error = ShearTransportError(q)
            + 2.0*cell_angle*separation;

        if (!isfinite(ShearTransportError(q))
            || transport_error > angular_tolerance)
            return FALSE;
    }
#endif
    return *radius/(*distance) <= max_ratio;
}

static void shear_accumulate_2pcf_sample(shear_pivot_workspace *work, int bin,
                                         shear_complex phase,
                                         shear_complex weighted_gamma,
                                         real weight)
{
    shear_complex z2 = shear_mul(phase, phase);
    shear_complex z4 = shear_mul(z2, z2);
    shear_complex pivot_weighted_gamma = shear_pivot_weighted_gamma(work);

    work->xi_plus[bin] = shear_add(
        work->xi_plus[bin],
        shear_mul(pivot_weighted_gamma, shear_conj(weighted_gamma)));
    work->xi_minus[bin] = shear_add(
        work->xi_minus[bin],
        shear_mul(shear_mul(pivot_weighted_gamma, weighted_gamma),
                  shear_conj(z4)));
    work->xi_weight[bin] += shear_pivot_weight(work)*weight;
}

static void shear_accumulate_3pcf_sample(shear_pivot_workspace *work, int bin,
                                         shear_complex phase,
                                         shear_complex weighted_gamma,
                                         real weight,
                                         shear_complex gamma2_sum,
                                         real gamma_abs2_sum,
                                         real weight2_sum)
{
    shear_complex phase_power = shear_make(1.0, 0.0);
    int order;
    size_t zero = shear_ring_index(0, bin, work->ring_max, work->bins);

    work->g_ring_active[zero] = shear_add(work->g_ring_active[zero],
                                         weighted_gamma);
    work->w_ring_active[zero].re += weight;
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

    if (work->collect_first_leg) {
        shear_complex z2 = shear_mul(phase, phase);
        shear_complex z6 = shear_mul(shear_mul(z2, z2), z2);
        work->diag_g6[bin] = shear_add(
            work->diag_g6[bin], shear_mul(gamma2_sum, shear_conj(z6)));
        work->diag_g2[bin] = shear_add(
            work->diag_g2[bin], shear_mul(gamma2_sum, shear_conj(z2)));
        work->diag_abs2[bin] = shear_add(
            work->diag_abs2[bin],
            shear_scale(shear_conj(z2), gamma_abs2_sum));
        work->diag_w2[bin] += weight2_sum;
    }
}

static void shear_accumulate_body_2pcf(shear_pivot_workspace *work, nodeptr q)
{
    real distance;
    real weight;
    int bin;
    shear_complex phase;
    shear_complex gamma;

    if (cballs_opt_read_mask(work->cmd) && Mask(q) != MASK_NODE_VALID)
        return;
    if (!shear_center_geometry(work, q, &distance, &phase, &bin))
        return;
    (void)distance;
    weight = Weight(q);
    gamma = shear_make(Gamma1(q), Gamma2(q));
#ifdef OCTREE_SHEAR_SPHERICAL
    if (!shear_transport_to_pivot(work, q, gamma, &gamma))
        return;
#endif
    shear_accumulate_2pcf_sample(work, bin, phase,
                                 shear_scale(gamma, weight), weight);
}

static void shear_accumulate_body_3pcf(shear_pivot_workspace *work, nodeptr q)
{
    real distance;
    real weight;
    int bin;
    shear_complex phase;
    shear_complex gamma;
    shear_complex weighted_gamma;

    if (cballs_opt_read_mask(work->cmd) && Mask(q) != MASK_NODE_VALID)
        return;
    if (!shear_center_geometry(work, q, &distance, &phase, &bin))
        return;
    (void)distance;
    weight = Weight(q);
    gamma = shear_make(Gamma1(q), Gamma2(q));
#ifdef OCTREE_SHEAR_SPHERICAL
    if (!shear_transport_to_pivot(work, q, gamma, &gamma))
        return;
#endif
    weighted_gamma = shear_scale(gamma, weight);
    shear_accumulate_3pcf_sample(
        work, bin, phase, weighted_gamma, weight,
        shear_mul(weighted_gamma, weighted_gamma),
        weight*weight*shear_abs2(gamma), weight*weight);
}

static void shear_accumulate_cell_2pcf(shear_pivot_workspace *work, nodeptr q,
                                       int bin, shear_complex phase)
{
    real weight = ShearWeightSum(q);
    shear_complex gamma = shear_make(Gamma1(q), Gamma2(q));
    shear_complex weighted_gamma;

#ifdef OCTREE_SHEAR_SPHERICAL
    if (!shear_transport_to_pivot(work, q, gamma, &gamma))
        return;
#endif
    weighted_gamma = shear_scale(gamma, weight);
    shear_accumulate_2pcf_sample(work, bin, phase, weighted_gamma, weight);
}

static void shear_accumulate_cell_3pcf(shear_pivot_workspace *work, nodeptr q,
                                       int bin, shear_complex phase)
{
    real weight = ShearWeightSum(q);
    shear_complex gamma = shear_make(Gamma1(q), Gamma2(q));
    shear_complex gamma2_sum = shear_make(
        ShearGamma2Re(q), ShearGamma2Im(q));
    shear_complex weighted_gamma;

#ifdef OCTREE_SHEAR_SPHERICAL
    shear_complex rotation;

    if (!shear_transport_rotation_to_pivot(
            work->pivot_unit, work->pivot_east, work->pivot_north,
            Pos(q), &rotation))
        return;
    gamma = shear_mul(gamma, rotation);
    gamma2_sum = shear_mul(
        gamma2_sum, shear_mul(rotation, rotation));
#endif
    weighted_gamma = shear_scale(gamma, weight);
    shear_accumulate_3pcf_sample(
        work, bin, phase, weighted_gamma, weight,
        gamma2_sum, ShearGammaAbs2(q), ShearWeight2(q));
}

/* Dedicated pair-only walk: no ring multipoles, diagonal moments, or 3PCF
 * storage are touched when options=only-2pcf is selected. */
static void shear_walk_tree_2pcf(shear_pivot_workspace *work, nodeptr q,
                                 real qsize)
{
    nodeptr child;

    if (q == NULL || q == (nodeptr)work->pivot)
        return;
    if (Type(q) == CELL) {
        real distance;
        real radius;
        int bin;
        shear_complex phase;
        bool accept;
        if (cballs_opt_read_mask(work->cmd)
            && Mask(q) == MASK_NODE_MASKED)
            return;
        accept = shear_cell_geometry(work, q, qsize, &distance, &phase,
                                          &bin, &radius);
        if (distance + radius <= work->cmd->rminHist
            || distance - radius >= work->cmd->rangeN)
            return;
        if (accept) {
            shear_accumulate_cell_2pcf(work, q, bin, phase);
            return;
        }
        for (child = More(q); child != Next(q); child = Next(child))
            shear_walk_tree_2pcf(work, child, qsize/2.0);
        return;
    }
    shear_accumulate_body_2pcf(work, q);
}

static void shear_walk_tree_3pcf(shear_pivot_workspace *work, nodeptr q,
                                 real qsize)
{
    nodeptr child;

    if (q == NULL || q == (nodeptr)work->pivot)
        return;
    if (Type(q) == CELL) {
        real distance;
        real radius;
        int bin;
        shear_complex phase;
        bool accept;
        if (cballs_opt_read_mask(work->cmd)
            && Mask(q) == MASK_NODE_MASKED)
            return;
        accept = shear_cell_geometry(work, q, qsize, &distance, &phase,
                                          &bin, &radius);
        if (distance + radius <= work->cmd->rminHist
            || distance - radius >= work->cmd->rangeN)
            return;
        if (accept) {
            shear_accumulate_cell_3pcf(work, q, bin, phase);
            return;
        }
        for (child = More(q); child != Next(q); child = Next(child))
            shear_walk_tree_3pcf(work, child, qsize/2.0);
        return;
    }
    shear_accumulate_body_3pcf(work, q);
}

static void shear_walk_tree_both(shear_pivot_workspace *work, nodeptr q,
                                 real qsize)
{
    nodeptr child;

    if (q == NULL || q == (nodeptr)work->pivot)
        return;
    if (Type(q) == CELL) {
        real distance;
        real radius;
        int bin;
        shear_complex phase;
        bool accept;
        if (cballs_opt_read_mask(work->cmd)
            && Mask(q) == MASK_NODE_MASKED)
            return;
        accept = shear_cell_geometry(work, q, qsize, &distance, &phase,
                                     &bin, &radius);
        if (distance + radius <= work->cmd->rminHist
            || distance - radius >= work->cmd->rangeN)
            return;
        if (accept) {
            shear_accumulate_cell_2pcf(work, q, bin, phase);
            shear_accumulate_cell_3pcf(work, q, bin, phase);
            return;
        }
        for (child = More(q); child != Next(q); child = Next(child))
            shear_walk_tree_both(work, child, qsize/2.0);
        return;
    }
    shear_accumulate_body_2pcf(work, q);
    shear_accumulate_body_3pcf(work, q);
}

static int shear_validate_catalog(struct cmdline_data *cmd, bodyptr table,
                                  INTEGER count, const char *role)
{
    INTEGER i;
    int axis;

    if (table == NULL || count <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": empty %s catalog", role);
        return FAILURE;
    }
    for (i = 0; i < count; i++) {
        bodyptr body_i = table + i;
        if (!isfinite(Pos(body_i)[0]) || !isfinite(Pos(body_i)[1])
            || !isfinite(Gamma1(body_i)) || !isfinite(Gamma2(body_i))
            || !isfinite(Weight(body_i)) || Weight(body_i) < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": invalid %s catalog row %" INTEGER_FMT,
                     role, i);
            return FAILURE;
        }
#ifdef OCTREE_SHEAR_SPHERICAL
        if (!isfinite(Pos(body_i)[2])
            || !(Pos(body_i)[0]*Pos(body_i)[0]
                 + Pos(body_i)[1]*Pos(body_i)[1]
                 + Pos(body_i)[2]*Pos(body_i)[2] > 0.0)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": invalid spherical position in %s "
                     "catalog row %" INTEGER_FMT, role, i);
            return FAILURE;
        }
#else
        for (axis = 2; axis < NDIM; axis++) {
            real reference = Pos(table)[axis];
            real tolerance = 64.0*DBL_EPSILON
                *(1.0 + rabs(reference) + rabs(Pos(body_i)[axis]));
            if (!isfinite(Pos(body_i)[axis])
                || rabs(Pos(body_i)[axis] - reference) > tolerance) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME ": %s catalog is not flat in axis %d",
                         role, axis);
                return FAILURE;
            }
        }
#endif
    }
    (void)axis;
    return SUCCESS;
}

static int shear_validate_common_plane(struct cmdline_data *cmd,
                                       bodyptr first, bodyptr second,
                                       const char *first_role,
                                       const char *second_role)
{
#ifdef OCTREE_SHEAR_SPHERICAL
    (void)cmd;
    (void)first;
    (void)second;
    (void)first_role;
    (void)second_role;
    return SUCCESS;
#else
    int axis;

    for (axis = 2; axis < NDIM; axis++) {
        real first_coordinate = Pos(first)[axis];
        real second_coordinate = Pos(second)[axis];
        real tolerance = 64.0*DBL_EPSILON
            *(1.0 + rabs(first_coordinate) + rabs(second_coordinate));

        if (rabs(first_coordinate - second_coordinate) > tolerance) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": %s and %s catalogs do not share "
                     "a tangent plane in axis %d",
                     first_role, second_role, axis);
            return FAILURE;
        }
    }
    return SUCCESS;
#endif
}

global int prepare_octree_shear_catalogs(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         bodyptr *btable, INTEGER *nbody)
{
    bodyptr p;
    int cat1;
    int cat2;
    int cat3;
    int ifile;
    int axis;
#ifndef OCTREE_SHEAR_SPHERICAL
    compute_vector minimum;
    compute_vector maximum;
    compute_vector center;
#endif

    if (cmd == NULL)
        return FAILURE;
    if (gd == NULL || btable == NULL || nbody == NULL || gd->ninfiles <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid catalog preparation state");
        return FAILURE;
    }

    cat1 = gd->iCatalogs[0];
    cat2 = gd->ninfiles >= 2 ? gd->iCatalogs[1] : cat1;
    cat3 = gd->ninfiles >= 3 ? gd->iCatalogs[2] : cat2;
    if (cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles
        || cat3 < 0 || cat3 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid catalog selection %d:%d:%d",
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

#if defined(SHEAR_SPHERE_BINARY_TWO_BALLS) && defined(SMOOTHPIVOT)
    /* Native octree builds normally convert/derive this value in treeload.c.
     * Independent binary-tree engines have no native-octree build stage, so
     * publish the same chord-space contract before the smoothing prepass. */
    if (!cballs_opt_smooth_pivot(cmd)) {
        gd->rsmooth[0] = 0.0;
    } else if (!strnull(cmd->rsmooth)) {
        double arcminutes;

        if (parse_double_checked(cmd->rsmooth, &arcminutes,
                                 cmd->error_message, _ERRORMSGSIZE_,
                                 "rsmooth") == FAILURE)
            return FAILURE;
        gd->rsmooth[0] = 2.0*rsin(
            0.5*(real)arcminutes*PI/(180.0*60.0));
    } else {
        gd->rsmooth[0] = MIN(0.01*cmd->rangeN, 0.5*cmd->rminHist);
    }
#endif

    for (ifile = 0; ifile < gd->ninfiles; ifile++) {
        if (btable[ifile] == NULL || nbody[ifile] <= 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME ": empty input catalog %d", ifile);
            return FAILURE;
        }
#ifdef OCTREE_SHEAR_SPHERICAL
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile]) {
            compute_vector unit;
            if (!shear_unit3(Pos(p), unit)) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME ": cannot normalize input catalog %d",
                         ifile);
                return FAILURE;
            }
            Pos(p)[0] = unit[0];
            Pos(p)[1] = unit[1];
            Pos(p)[2] = unit[2];
        }
#else
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
#endif
    }
#ifndef OCTREE_SHEAR_SPHERICAL
    DO_COORD(axis)
        center[axis] = 0.5*(minimum[axis] + maximum[axis]);
    for (ifile = 0; ifile < gd->ninfiles; ifile++)
        DO_BODY(p, btable[ifile], btable[ifile] + nbody[ifile])
            DO_COORD(axis)
                Pos(p)[axis] -= center[axis];
#endif

    (void)axis;
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
                 SHEAR_ENGINE_NAME ": invalid search state");
        return FAILURE;
    }
    if (cballs_opt_only_2pcf(cmd) && cballs_opt_only_3pcf(cmd)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": only-2pcf and only-3pcf are mutually exclusive");
        return FAILURE;
    }
    if (!isfinite(cmd->theta) || cmd->theta < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": theta must be finite and nonnegative");
        return FAILURE;
    }
    if (cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles
        || cat3 < 0 || cat3 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid catalog selection %d:%d:%d",
                 cat1, cat2, cat3);
        return FAILURE;
    }
    if (cmd->usePeriodic) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": periodic geometry is not supported");
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
                 SHEAR_ENGINE_NAME ": invalid histogram domain");
        return FAILURE;
    }
#ifdef OCTREE_SHEAR_SPHERICAL
    if (cmd->rangeN > 2.0 + 64.0*DBL_EPSILON) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": chord-distance rangeN cannot exceed 2");
        return FAILURE;
    }
#endif
    if (cmd->useLogHist && cmd->rminHist == 0.0
        && cmd->logHistBinsPD <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": logHistBinsPD must be positive");
        return FAILURE;
    }
    if (ipmin < 1 || ipmax[cat1] < ipmin || ipmax[cat1] > nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid pivot interval");
        return FAILURE;
    }
#ifndef SHEAR_SPHERE_BINARY_TWO_BALLS
    if (roottable[cat2] == NULL || roottable[cat3] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": a neighbor tree is not available");
        return FAILURE;
    }
#endif
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

static int shear_configure_smooth_radius(struct cmdline_data *cmd,
                                         struct global_data *gd)
{
#ifdef SMOOTHPIVOT
    real maximum;
    bool explicit_radius;

    if (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0)
        return SUCCESS;
    if (!isfinite(gd->rsmooth[0]) || gd->rsmooth[0] < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": invalid smooth-pivot radius %g",
                 gd->rsmooth[0]);
        return FAILURE;
    }

    explicit_radius = cmd->rsmooth != NULL && cmd->rsmooth[0] != '\0';
    if (!(cmd->rminHist > 0.0)) {
        if (explicit_radius) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     SHEAR_ENGINE_NAME
                     ": positive rsmooth requires rminHist > 0; "
                     "use no-smooth-pivot or raise rminHist");
            return FAILURE;
        }
        gd->rsmooth[0] = 0.0;
        return SUCCESS;
    }

    maximum = SHEAR_SMOOTH_MAX_RMIN_FRACTION*cmd->rminHist;
    if (gd->rsmooth[0] <= maximum)
        return SUCCESS;
    if (explicit_radius) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME
                 ": rsmooth=%g exceeds the safe limit 0.5*rminHist=%g; "
                 "the group would contain resolved pairs",
                 gd->rsmooth[0], maximum);
        return FAILURE;
    }

    verb_print_normal_info(
        cmd->verbose, cmd->verbose_log, gd->outlog,
        SHEAR_ENGINE_NAME
        ": limiting automatic rsmooth from %g to %g "
        "so 2*rsmooth <= rminHist\n",
        gd->rsmooth[0], maximum);
    gd->rsmooth[0] = maximum;
#else
    (void)cmd;
    (void)gd;
#endif
    return SUCCESS;
}

static int shear_allocate_results(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  size_t bins, size_t gamma_count,
                                  size_t denominator_count,
                                  size_t angular_count,
                                  bool run_2pcf, bool run_3pcf)
{
#define SHEAR_ALLOC_RESULT(member, count)                                \
    do {                                                                 \
        if (shear_calloc(cmd, (void **)&gd->member, (count),             \
                         sizeof(*gd->member), #member) == FAILURE)       \
            goto fail;                                                   \
    } while (0)

    shear_clear_results(gd);
    if (run_2pcf) {
        SHEAR_ALLOC_RESULT(histShearXiPlusRe, bins);
        SHEAR_ALLOC_RESULT(histShearXiPlusIm, bins);
        SHEAR_ALLOC_RESULT(histShearXiMinusRe, bins);
        SHEAR_ALLOC_RESULT(histShearXiMinusIm, bins);
        SHEAR_ALLOC_RESULT(histShearXiWeight, bins);
    }
    if (run_3pcf) {
        SHEAR_ALLOC_RESULT(histShearGammaNumeratorRe, gamma_count);
        SHEAR_ALLOC_RESULT(histShearGammaNumeratorIm, gamma_count);
        SHEAR_ALLOC_RESULT(histShearGammaMultipoleRe, gamma_count);
        SHEAR_ALLOC_RESULT(histShearGammaMultipoleIm, gamma_count);
        SHEAR_ALLOC_RESULT(histShearDenominatorRe, denominator_count);
        SHEAR_ALLOC_RESULT(histShearDenominatorIm, denominator_count);
        SHEAR_ALLOC_RESULT(histShearGammaRe, angular_count);
        SHEAR_ALLOC_RESULT(histShearGammaIm, angular_count);
    }
#undef SHEAR_ALLOC_RESULT
    return SUCCESS;

fail:
#undef SHEAR_ALLOC_RESULT
    shear_clear_results(gd);
    return FAILURE;
}

static void shear_bind_global_accumulator(
        struct global_data *gd, shear_result_accumulator *result)
{
    result->xi_plus_re = gd->histShearXiPlusRe;
    result->xi_plus_im = gd->histShearXiPlusIm;
    result->xi_minus_re = gd->histShearXiMinusRe;
    result->xi_minus_im = gd->histShearXiMinusIm;
    result->xi_weight = gd->histShearXiWeight;
    result->gamma_re = gd->histShearGammaNumeratorRe;
    result->gamma_im = gd->histShearGammaNumeratorIm;
    result->denominator_re = gd->histShearDenominatorRe;
    result->denominator_im = gd->histShearDenominatorIm;
}

#ifdef BALLS4SCANLEV
static void shear_bind_local_accumulator(
        real *storage, size_t bins, size_t gamma_count,
        size_t denominator_count, bool run_2pcf, bool run_3pcf,
        shear_result_accumulator *result)
{
    real *cursor = storage;

    memset(result, 0, sizeof(*result));
    if (run_2pcf) {
        result->xi_plus_re = cursor; cursor += bins;
        result->xi_plus_im = cursor; cursor += bins;
        result->xi_minus_re = cursor; cursor += bins;
        result->xi_minus_im = cursor; cursor += bins;
        result->xi_weight = cursor; cursor += bins;
    }
    if (run_3pcf) {
        result->gamma_re = cursor; cursor += gamma_count;
        result->gamma_im = cursor; cursor += gamma_count;
        result->denominator_re = cursor; cursor += denominator_count;
        result->denominator_im = cursor;
    }
}

static void shear_zero_accumulator(shear_result_accumulator *result,
                                   size_t bins, size_t gamma_count,
                                   size_t denominator_count)
{
    if (result->xi_plus_re != NULL) {
        memset(result->xi_plus_re, 0, bins*sizeof(*result->xi_plus_re));
        memset(result->xi_plus_im, 0, bins*sizeof(*result->xi_plus_im));
        memset(result->xi_minus_re, 0, bins*sizeof(*result->xi_minus_re));
        memset(result->xi_minus_im, 0, bins*sizeof(*result->xi_minus_im));
        memset(result->xi_weight, 0, bins*sizeof(*result->xi_weight));
    }
    if (result->gamma_re != NULL) {
        memset(result->gamma_re, 0, gamma_count*sizeof(*result->gamma_re));
        memset(result->gamma_im, 0, gamma_count*sizeof(*result->gamma_im));
        memset(result->denominator_re, 0,
               denominator_count*sizeof(*result->denominator_re));
        memset(result->denominator_im, 0,
               denominator_count*sizeof(*result->denominator_im));
    }
}

static void shear_merge_accumulator(shear_result_accumulator *target,
                                    const shear_result_accumulator *source,
                                    size_t bins, size_t gamma_count,
                                    size_t denominator_count)
{
    size_t i;

    if (source->xi_plus_re != NULL) {
        for (i = 0; i < bins; i++) {
            target->xi_plus_re[i] += source->xi_plus_re[i];
            target->xi_plus_im[i] += source->xi_plus_im[i];
            target->xi_minus_re[i] += source->xi_minus_re[i];
            target->xi_minus_im[i] += source->xi_minus_im[i];
            target->xi_weight[i] += source->xi_weight[i];
        }
    }
    if (source->gamma_re != NULL) {
        for (i = 0; i < gamma_count; i++) {
            target->gamma_re[i] += source->gamma_re[i];
            target->gamma_im[i] += source->gamma_im[i];
        }
        for (i = 0; i < denominator_count; i++) {
            target->denominator_re[i] += source->denominator_re[i];
            target->denominator_im[i] += source->denominator_im[i];
        }
    }
}
#endif

static void shear_reduce_pivot(shear_result_accumulator *result,
                               shear_pivot_workspace *work,
                               int nmax, bool run_2pcf, bool run_3pcf)
{
    int bin;
    int bin1;
    int bin2;
    int order;
    int multipoles = 2*nmax + 1;
    int denominator_offset = 2*nmax;
    shear_complex pivot_weighted_gamma = shear_pivot_weighted_gamma(work);
    real pivot_weight = shear_pivot_weight(work);

    if (run_2pcf) {
        for (bin = 0; bin < work->bins; bin++) {
            result->xi_plus_re[bin] += work->xi_plus[bin].re;
            result->xi_plus_im[bin] += work->xi_plus[bin].im;
            result->xi_minus_re[bin] += work->xi_minus[bin].re;
            result->xi_minus_im[bin] += work->xi_minus[bin].im;
            result->xi_weight[bin] += work->xi_weight[bin];
        }
    }
    if (!run_3pcf)
        return;

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
                value = shear_scale(value, pivot_weight);
                result->denominator_re[index] += value.re;
                result->denominator_im[index] += value.im;
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
                    result->gamma_re[index] += value.re;
                    result->gamma_im[index] += value.im;
                }
            }
        }
    }
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
                 SHEAR_ENGINE_NAME ": mode-coupling workspace overflow");
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
                 SHEAR_ENGINE_NAME ": cannot write %s", path);
        return FAILURE;
    }
    return SUCCESS;
}

static int shear_write_outputs(struct cmdline_data *cmd,
                               struct global_data *gd, int bins, int nmax,
                               bool run_2pcf, bool run_3pcf)
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

    if (run_2pcf) {
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
    }

    if (!run_3pcf)
        return SUCCESS;

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
             SHEAR_ENGINE_NAME ": output path is too long");
    return FAILURE;
open_error:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             SHEAR_ENGINE_NAME ": cannot open %s: %s", path, strerror(errno));
    return FAILURE;
}

static void shear_normalize_2pcf(struct global_data *gd, int bins)
{
    int bin;

    for (bin = 0; bin < bins; bin++) {
        real denominator = gd->histShearXiWeight[bin];

        gd->histShearXiPlusRe[bin] = cballs_normalize_or_zero(
            gd->histShearXiPlusRe[bin], denominator);
        gd->histShearXiPlusIm[bin] = cballs_normalize_or_zero(
            gd->histShearXiPlusIm[bin], denominator);
        gd->histShearXiMinusRe[bin] = cballs_normalize_or_zero(
            gd->histShearXiMinusRe[bin], denominator);
        gd->histShearXiMinusIm[bin] = cballs_normalize_or_zero(
            gd->histShearXiMinusIm[bin], denominator);
    }
}

#ifdef OCTREE_SHEAR_SPHERICAL
#include "../octree_shear_sphere_omp/shear_sphere_dual_tree.h"
#endif
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
#include "../shear_sphere_binary_2balls/shear_sphere_binary_scan.h"
#endif

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
    bool run_2pcf;
    bool run_3pcf;
    bool scan_2pcf;
    size_t bins_squared = 0;
    size_t gamma_count = 0;
    size_t denominator_count = 0;
    size_t angular_count = 0;
    size_t ring_count = 0;
    size_t threaded_ring_count = 0;
    size_t threaded_bin_count = 0;
#ifdef BALLS4SCANLEV
    size_t accumulator_stride = 0;
    size_t threaded_accumulator_count;
#endif
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
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    fcfc_balltreeptr kd_pivot_tree = NULL;
    fcfc_balltreeptr kd_first_tree = NULL;
    fcfc_balltreeptr kd_second_tree = NULL;
#endif
#ifdef BALLS4SCANLEV
    real *block_accumulator_all = NULL;
    shear_result_accumulator global_accumulator;
    bodyptr *pivot_order = NULL;
    INTEGER pivot_count = 0;
    INTEGER block_count;
#endif

    if (cmd == NULL)
        return FAILURE;
    cmd->error_message[0] = '\0';
#ifdef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    if (cballs_opt_legacy_one_ball(cmd)) {
        verb_print(cmd->verbose,
                   SHEAR_ENGINE_NAME ": dispatching to the "
                   "octree-shear-sphere-omp compatibility kernel\n");
        return searchcalc_octree_shear_sphere_omp(
            cmd, gd, btable, nbody, ipmin, ipmax, cat1, cat2, cat3);
    }
#endif
    if (shear_validate(cmd, gd, btable, nbody, ipmin, ipmax,
                       cat1, cat2, cat3)
        == FAILURE)
        return FAILURE;
    if (shear_configure_smooth_radius(cmd, gd) == FAILURE)
        return FAILURE;
    run_2pcf = !cballs_opt_only_3pcf(cmd);
    run_3pcf = !cballs_opt_only_2pcf(cmd);
    scan_2pcf = run_2pcf;
#ifdef OCTREE_SHEAR_SPHERICAL
#ifndef SHEAR_SPHERE_BINARY_TWO_BALLS
    {
        const int catalogs[3] = {cat1, cat2, cat3};
        bool prepared[MAXITEMS] = {FALSE};
        int selected;

        for (selected = 0; selected < 3; selected++) {
            int catalog = catalogs[selected];

            if (prepared[catalog])
                continue;
            if (roottable[catalog] == NULL) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         SHEAR_ENGINE_NAME
                         ": tree is not available for catalog %d", catalog);
                return FAILURE;
            }
            if (shear_prepare_spherical_cell_moments(
                    cmd, (nodeptr)roottable[catalog],
                    gd->rSizeTable[catalog]) == FAILURE)
                return FAILURE;
            prepared[catalog] = TRUE;
        }
    }
#endif
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
        shear_spherical_smooth_context smooth_context;

        memset(&smooth_context, 0, sizeof(smooth_context));
        if (prepare_smooth_pivots_with_accumulator(
                cmd, gd, btable, nbody, ipmin, ipmax, cat1, cat1,
                shear_spherical_smooth_claim, &smooth_context) == FAILURE)
            return FAILURE;
    }
#endif
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    if (SHEAR_SPHERE_BINARY_TREE_BUILD(
            cmd, gd, btable[cat1], nbody[cat1], cmd->nsmooth,
            TRUE, &kd_pivot_tree) == FAILURE
        || SHEAR_SPHERE_BINARY_TREE_BUILD(
            cmd, gd, btable[cat2], nbody[cat2], cmd->nsmooth,
            FALSE, &kd_first_tree) == FAILURE) {
        kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
        return FAILURE;
    }
    if (cat3 == cat2) {
        kd_second_tree = kd_first_tree;
    } else if (SHEAR_SPHERE_BINARY_TREE_BUILD(
                   cmd, gd, btable[cat3], nbody[cat3], cmd->nsmooth,
                   FALSE, &kd_second_tree) == FAILURE) {
        kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
        return FAILURE;
    }
    gd->ncellTable[cat1] = kd_pivot_tree->nnode;
    gd->ncellTable[cat2] = kd_first_tree->nnode;
    gd->ncellTable[cat3] = kd_second_tree->nnode;
#endif
    if (run_2pcf && !run_3pcf
        && (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0)) {
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
        bool dual_eligible = ipmin == 1 && ipmax[cat1] == nbody[cat1];
#else
        shear_sphere_pair_context context;

        context.cmd = cmd;
        context.gd = gd;
        context.bins = cmd->sizeHistN;
        bool dual_eligible = shear_sphere_dual_tree_eligible(
            &context, btable[cat1], nbody[cat1], ipmin, ipmax[cat1]);
#endif
        if (dual_eligible) {
            bins = cmd->sizeHistN;
            if (shear_allocate_results(cmd, gd, (size_t)bins, 0, 0, 0,
                                       TRUE, FALSE) == FAILURE) {
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
                kd_shear_release_trees(
                    kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
                return FAILURE;
            }
            gd->nbbcalc = 0;
            gd->nbccalc = 0;
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
            if (kd_shear_dual_tree_2pcf(
                    cmd, gd, kd_pivot_tree, kd_first_tree,
                    cat1 == cat2, bins) == FAILURE) {
#else
            if (shear_sphere_dual_tree_2pcf(
                    cmd, gd, (nodeptr)roottable[cat1],
                    (nodeptr)roottable[cat2], cat1 == cat2, bins)
                == FAILURE) {
#endif
                shear_clear_results(gd);
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
                kd_shear_release_trees(
                    kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
                return FAILURE;
            }
            shear_normalize_2pcf(gd, bins);
            gd->cpusearch = CPUTIME - cpu_start;
            if (!cballs_opt_no_out_hist(cmd)
                && shear_write_outputs(cmd, gd, bins, cmd->mChebyshev,
                                       TRUE, FALSE) == FAILURE) {
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
                kd_shear_release_trees(
                    kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
                return FAILURE;
            }
            verb_print_min_info(
                cmd->verbose, cmd->verbose_log, gd->outlog,
                SHEAR_ENGINE_NAME ": completed dual-tree 2PCF in %g %s\n",
                gd->cpusearch, PRNUNITOFTIMEUSED);
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
            kd_shear_release_trees(
                kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
            return SUCCESS;
        }
    }
#endif
#if defined(SMOOTHPIVOT) && !defined(OCTREE_SHEAR_SPHERICAL)
    if (cballs_opt_smooth_pivot(cmd)
        && prepare_smooth_pivots(cmd, gd, btable, nbody,
                                 ipmin, ipmax, cat1, cat1) == FAILURE)
        return FAILURE;
#endif
#ifdef BALLS4SCANLEV
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    pivot_count = kd_pivot_tree->npoint;
    if ((size_t)pivot_count > SIZE_MAX/sizeof(*pivot_order)
        || (pivot_order = malloc((size_t)pivot_count*sizeof(*pivot_order)))
             == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": pivot-order allocation failed");
        kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
        return FAILURE;
    }
    memcpy(pivot_order, kd_pivot_tree->bptr,
           (size_t)pivot_count*sizeof(*pivot_order));
#else
    if (cballs_octree_scan_body_order(
            cmd, gd, btable[cat1] + ipmin - 1,
            btable[cat1] + ipmax[cat1], cat1,
            &pivot_order, &pivot_count) == FAILURE)
        return FAILURE;
#ifdef SMOOTHPIVOT
    if (cballs_opt_smooth_pivot(cmd)) {
        INTEGER active_count = 0;
        for (INTEGER pivot_index = 0; pivot_index < pivot_count; pivot_index++)
            if (Update(pivot_order[pivot_index]))
                pivot_order[active_count++] = pivot_order[pivot_index];
        pivot_count = active_count;
    }
#endif
#endif
    verb_print(cmd->verbose,
               SHEAR_ENGINE_NAME ": BALLS4SCANLEV spatially ordered %"
               INTEGER_FMT " active body pivots\n", pivot_count);
#endif

    bins = cmd->sizeHistN;
    nmax = cmd->mChebyshev;
    multipoles = 2*nmax + 1;
    denominator_multipoles = 4*nmax + 1;
    ring_max = run_3pcf
        ? (2*nmax > nmax + 3 ? 2*nmax : nmax + 3)
        : 4;
#ifdef OPENMPCODE
    threads = omp_get_max_threads();
#endif
    if (threads < 1)
        threads = 1;

    if (shear_size_mul((size_t)threads, (size_t)bins,
                       &threaded_bin_count) == FAILURE
        || (run_3pcf
            && (shear_size_mul((size_t)bins, (size_t)bins,
                               &bins_squared) == FAILURE
                || shear_size_mul(4*(size_t)multipoles, bins_squared,
                                  &gamma_count) == FAILURE
                || shear_size_mul((size_t)denominator_multipoles, bins_squared,
                                  &denominator_count) == FAILURE
                || shear_size_mul(4*(size_t)cmd->sizeHistPhi, bins_squared,
                                  &angular_count) == FAILURE
                || shear_size_mul((size_t)(2*ring_max + 1), (size_t)bins,
                                  &ring_count) == FAILURE
                || shear_size_mul((size_t)threads, ring_count,
                                  &threaded_ring_count) == FAILURE))
#ifdef BALLS4SCANLEV
        || (run_2pcf
            && shear_size_mul(5, (size_t)bins, &accumulator_stride) == FAILURE)
        || (run_3pcf && gamma_count > SIZE_MAX/2)
        || (run_3pcf
            && shear_size_add(accumulator_stride, 2*gamma_count,
                              &accumulator_stride) == FAILURE)
        || (run_3pcf && denominator_count > SIZE_MAX/2)
        || (run_3pcf
            && shear_size_add(accumulator_stride, 2*denominator_count,
                              &accumulator_stride) == FAILURE)
        || shear_size_mul((size_t)threads, accumulator_stride,
                          &threaded_accumulator_count) == FAILURE
#endif
        ) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 SHEAR_ENGINE_NAME ": histogram size overflow");
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
        kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
#ifdef BALLS4SCANLEV
        free(pivot_order);
#endif
        return FAILURE;
    }
    if (shear_allocate_results(cmd, gd, bins, gamma_count,
                               denominator_count, angular_count,
                               run_2pcf, run_3pcf) == FAILURE) {
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
        kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
#ifdef BALLS4SCANLEV
        free(pivot_order);
#endif
        return FAILURE;
    }
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    if (run_2pcf && run_3pcf
        && (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0)
        && ipmin == 1 && ipmax[cat1] == nbody[cat1]) {
        if (kd_shear_dual_tree_2pcf(
                cmd, gd, kd_pivot_tree, kd_first_tree,
                cat1 == cat2, bins) == FAILURE)
            goto fail;
        scan_2pcf = FALSE;
    }
#endif
#ifdef OCTREE_SHEAR_SPHERICAL_TWO_BALLS
    if (run_2pcf && run_3pcf
        && (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0)) {
        shear_sphere_pair_context context;

        context.cmd = cmd;
        context.gd = gd;
        context.bins = bins;
        if (shear_sphere_dual_tree_eligible(
                &context, btable[cat1], nbody[cat1], ipmin, ipmax[cat1])) {
            if (shear_sphere_dual_tree_2pcf(
                    cmd, gd, (nodeptr)roottable[cat1],
                    (nodeptr)roottable[cat2], cat1 == cat2, bins)
                == FAILURE)
                goto fail;
            scan_2pcf = FALSE;
        }
    }
#endif

#define SHEAR_ALLOC_WORK(pointer, count, label)                          \
    do {                                                                 \
        if (shear_calloc(cmd, (void **)&(pointer), (count),              \
                         sizeof(*(pointer)), (label)) == FAILURE)        \
            goto fail;                                                   \
    } while (0)
    if (run_3pcf) {
        SHEAR_ALLOC_WORK(g_ring_first_all, threaded_ring_count,
                         "thread first G rings");
        SHEAR_ALLOC_WORK(w_ring_first_all, threaded_ring_count,
                         "thread first W rings");
        SHEAR_ALLOC_WORK(g_ring_second_all, threaded_ring_count,
                         "thread second G rings");
        SHEAR_ALLOC_WORK(w_ring_second_all, threaded_ring_count,
                         "thread second W rings");
        SHEAR_ALLOC_WORK(diag_g6_all, threaded_bin_count,
                         "thread G6 diagonal");
        SHEAR_ALLOC_WORK(diag_g2_all, threaded_bin_count,
                         "thread G2 diagonal");
        SHEAR_ALLOC_WORK(diag_abs2_all, threaded_bin_count,
                         "thread absolute-shear diagonal");
        SHEAR_ALLOC_WORK(diag_w2_all, threaded_bin_count,
                         "thread weight diagonal");
    }
    if (scan_2pcf) {
        SHEAR_ALLOC_WORK(xi_plus_all, threaded_bin_count, "thread xi-plus");
        SHEAR_ALLOC_WORK(xi_minus_all, threaded_bin_count, "thread xi-minus");
        SHEAR_ALLOC_WORK(xi_weight_all, threaded_bin_count, "thread xi weight");
    }
#ifdef BALLS4SCANLEV
    SHEAR_ALLOC_WORK(block_accumulator_all, threaded_accumulator_count,
                     "thread block result accumulators");
    shear_bind_global_accumulator(gd, &global_accumulator);
    block_count = pivot_count > 0
        ? 1 + (pivot_count - 1)/SHEAR_OMP_PIVOT_BLOCK_SIZE : 0;
#endif
#undef SHEAR_ALLOC_WORK

#ifdef BALLS4SCANLEV
#define SHEAR_SCAN_SHARED shared(pivot_order,pivot_count,block_count,         \
                                 block_accumulator_all,global_accumulator,    \
                                 accumulator_stride,gamma_count,             \
                                 denominator_count)
#else
#define SHEAR_SCAN_SHARED
#endif
#pragma omp parallel private(ip) SHEAR_SCAN_SHARED
    {
        int thread_id = 0;
        shear_pivot_workspace work;
        shear_result_accumulator accumulator;
#ifdef OPENMPCODE
        thread_id = omp_get_thread_num();
#endif
        work.cmd = cmd;
        work.gd = gd;
        work.bins = bins;
        work.ring_max = ring_max;
        work.g_ring_first = run_3pcf
            ? g_ring_first_all + (size_t)thread_id*ring_count : NULL;
        work.w_ring_first = run_3pcf
            ? w_ring_first_all + (size_t)thread_id*ring_count : NULL;
        work.g_ring_second = run_3pcf
            ? g_ring_second_all + (size_t)thread_id*ring_count : NULL;
        work.w_ring_second = run_3pcf
            ? w_ring_second_all + (size_t)thread_id*ring_count : NULL;
        work.diag_g6 = run_3pcf
            ? diag_g6_all + (size_t)thread_id*bins : NULL;
        work.diag_g2 = run_3pcf
            ? diag_g2_all + (size_t)thread_id*bins : NULL;
        work.diag_abs2 = run_3pcf
            ? diag_abs2_all + (size_t)thread_id*bins : NULL;
        work.diag_w2 = run_3pcf
            ? diag_w2_all + (size_t)thread_id*bins : NULL;
        work.xi_plus = scan_2pcf
            ? xi_plus_all + (size_t)thread_id*bins : NULL;
        work.xi_minus = scan_2pcf
            ? xi_minus_all + (size_t)thread_id*bins : NULL;
        work.xi_weight = scan_2pcf
            ? xi_weight_all + (size_t)thread_id*bins : NULL;
        work.same_neighbor_catalog = cat2 == cat3;
        work.allow_cells = !cballs_opt_no_one_ball(cmd) && cmd->theta > 0.0;

#ifdef BALLS4SCANLEV
        shear_bind_local_accumulator(
            block_accumulator_all + (size_t)thread_id*accumulator_stride,
            (size_t)bins, gamma_count, denominator_count,
            scan_2pcf, run_3pcf, &accumulator);
#pragma omp for schedule(dynamic,1) ordered
        for (INTEGER block = 0; block < block_count; block++) {
            const INTEGER first = block*SHEAR_OMP_PIVOT_BLOCK_SIZE;
            const INTEGER count = MIN((INTEGER)SHEAR_OMP_PIVOT_BLOCK_SIZE,
                                      pivot_count - first);
            shear_zero_accumulator(&accumulator, (size_t)bins, gamma_count,
                                   denominator_count);
            for (INTEGER offset = 0; offset < count; offset++) {
                work.pivot = pivot_order[first + offset];
#else
        shear_bind_global_accumulator(gd, &accumulator);
#pragma omp for schedule(static,1) ordered
        for (ip = ipmin - 1; ip < ipmax[cat1]; ip++) {
            work.pivot = btable[cat1] + ip;
#endif
            if (run_3pcf) {
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
            }
            if (scan_2pcf) {
                memset(work.xi_plus, 0, (size_t)bins*sizeof(*work.xi_plus));
                memset(work.xi_minus, 0, (size_t)bins*sizeof(*work.xi_minus));
                memset(work.xi_weight, 0,
                       (size_t)bins*sizeof(*work.xi_weight));
            }

            if (Update(work.pivot)
                && shear_prepare_pivot_geometry(&work)
                && (!cballs_opt_read_mask(cmd)
                    || Mask(work.pivot) == MASK_NODE_VALID)) {
                work.g_ring_active = work.g_ring_first;
                work.w_ring_active = work.w_ring_first;
                work.collect_first_leg = TRUE;
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
                kd_shear_walk(&work, kd_first_tree, 0,
                              scan_2pcf, run_3pcf);
#else
                if (scan_2pcf && run_3pcf)
                    shear_walk_tree_both(&work, (nodeptr)roottable[cat2],
                                         gd->rSizeTable[cat2]);
                else if (scan_2pcf)
                    shear_walk_tree_2pcf(&work, (nodeptr)roottable[cat2],
                                         gd->rSizeTable[cat2]);
                else
                    shear_walk_tree_3pcf(&work, (nodeptr)roottable[cat2],
                                         gd->rSizeTable[cat2]);
#endif
                if (run_3pcf && cat3 == cat2) {
                    memcpy(work.g_ring_second, work.g_ring_first,
                           ring_count*sizeof(*work.g_ring_second));
                    memcpy(work.w_ring_second, work.w_ring_first,
                           ring_count*sizeof(*work.w_ring_second));
                } else if (run_3pcf) {
                    work.g_ring_active = work.g_ring_second;
                    work.w_ring_active = work.w_ring_second;
                    work.collect_first_leg = FALSE;
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
                    kd_shear_walk(&work, kd_second_tree, 0, FALSE, TRUE);
#else
                    shear_walk_tree_3pcf(&work, (nodeptr)roottable[cat3],
                                         gd->rSizeTable[cat3]);
#endif
                }
            }

#ifdef BALLS4SCANLEV
                shear_reduce_pivot(&accumulator, &work, nmax,
                                   scan_2pcf, run_3pcf);
            }
#pragma omp ordered
            shear_merge_accumulator(&global_accumulator, &accumulator,
                                    (size_t)bins, gamma_count,
                                    denominator_count);
#else
#pragma omp ordered
            shear_reduce_pivot(&accumulator, &work, nmax,
                               scan_2pcf, run_3pcf);
#endif
        }
    }
#undef SHEAR_SCAN_SHARED

    if (run_2pcf)
        shear_normalize_2pcf(gd, bins);
    if (run_3pcf) {
        if (shear_solve_mode_coupling(cmd, gd, bins, nmax) == FAILURE)
            goto fail;
        shear_reconstruct_angular(cmd, gd, bins, nmax);
        gd->shearMultipoleMax = nmax;
        gd->shearAngularBins = cmd->sizeHistPhi;
    }
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
#ifdef BALLS4SCANLEV
    free(block_accumulator_all);
    free(pivot_order);
#endif
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
    if (!cballs_opt_no_out_hist(cmd)
        && shear_write_outputs(cmd, gd, bins, nmax,
                               run_2pcf, run_3pcf) == FAILURE)
        return FAILURE;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        SHEAR_ENGINE_NAME ": completed in %g %s\n",
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
#ifdef BALLS4SCANLEV
    free(block_accumulator_all);
    free(pivot_order);
#endif
#ifdef SHEAR_SPHERE_BINARY_TWO_BALLS
    kd_shear_release_trees(kd_pivot_tree, kd_first_tree, kd_second_tree);
#endif
    shear_clear_results(gd);
    return FAILURE;
}

#endif /* NDIM >= 2 */
