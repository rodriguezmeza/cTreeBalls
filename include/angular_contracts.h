#ifndef CBALLS_ANGULAR_CONTRACTS_H
#define CBALLS_ANGULAR_CONTRACTS_H

/* Scalar Fourier multipoles in 3D use observer-relative tangent directions,
 * not the interior angle between three-dimensional chords. */
static inline bool cballs_observer_frame(const struct cmdline_data *cmd)
{
#if NDIM == 3
    const char *method = cmd->searchMethod;
    return !cmd->usePeriodic && method != NULL
        && (!strcmp(method, "octree-sincos-omp")
            || !strncmp(method, "octree-ggg-", 11)
            || !strncmp(method, "octree-balls4-", 14)
            || !strncmp(method, "octree-2balls-", 14)
            || !strcmp(method, "kdtree-omp")
            || !strncmp(method, "balltree-", 9));
#else
    (void)cmd;
    return FALSE;

#endif
}

#if NDIM == 3
static inline bool cballs_angular_basis(const cballs_storage_real *pivot,
        real normal[3], real first[3], real second[3])
{
    real norm = hypot(hypot((real)pivot[0], (real)pivot[1]), (real)pivot[2]);
    if (!(norm > 0.0) || !isfinite(norm)) return FALSE;
    int axis = 0;
    for (int k = 0; k < 3; k++) normal[k] = (real)pivot[k] / norm;
    for (int k = 1; k < 3; k++)
        if (fabs(normal[k]) < fabs(normal[axis])) axis = k;
    for (int k = 0; k < 3; k++)
        first[k] = (k == axis ? 1.0 : 0.0) - normal[axis]*normal[k];
    norm = sqrt(1.0 - normal[axis]*normal[axis]);
    for (int k = 0; k < 3; k++) first[k] /= norm;
    CROSSVP(second, first, normal);
    return TRUE;
}
#endif

/* dr = pivot - neighbor. The handedness matches TreeCorr's positive spherical
 * LogMultipole modes; basis orientation cancels between the two legs. */
static inline bool cballs_angular_phase(const cballs_storage_real *pivot,
                                         const compute_vector dr,
                                         real *cosine, real *sine)
{
    real x, y;
#if NDIM == 3
    real normal[3], first[3], second[3];
    if (!cballs_angular_basis(pivot, normal, first, second)) return FALSE;
    x = y = 0.0;
    for (int k = 0; k < 3; k++) {
        x -= dr[k]*first[k];
        y -= dr[k]*second[k];
    }
    /* Coincident/radial/antipodal legs have no bearing. Exclude them from
     * angular multipoles, independently of the ordinary 2PCF pair count. */
    real norm = hypot(x, y);
    if (!(norm > 32.0*DBL_EPSILON*hypot(hypot(dr[0], dr[1]), dr[2])))
        return FALSE;
#else
    (void)pivot;
    x = -dr[0]; y = -dr[1];
    real norm = hypot(x, y);
    if (!(norm > 0.0)) return FALSE;
#endif
    *cosine = x/norm;
    *sine = y/norm;
    return isfinite(*cosine) && isfinite(*sine);
}

/* Conservative variation of a projected bearing, including pivot motion. */
static inline real cballs_angular_extent(const cballs_storage_real *pivot,
                                         const compute_vector dr,
                                         real pivot_radius, real neighbor_radius)
{
    real distance = 0.0;
    for (int k = 0; k < NDIM; k++) distance += dr[k]*dr[k];
    distance = sqrt(distance);
#if NDIM == 3
    real norm = hypot(hypot((real)pivot[0], (real)pivot[1]), (real)pivot[2]);
    if (!(norm > pivot_radius)) return 1.0;
    real cross[3], normal[3];
    for (int k = 0; k < 3; k++) normal[k] = (real)pivot[k]/norm;
    CROSSVP(cross, normal, dr);
    real transverse = hypot(hypot(cross[0], cross[1]), cross[2]);
    real error = pivot_radius + neighbor_radius
        + distance*pivot_radius/(norm-pivot_radius);
    return transverse > error ? error/transverse : 1.0;
#else
    (void)pivot;
    real error = pivot_radius + neighbor_radius;
    return distance > error ? error/distance : 1.0;
#endif
}

static inline bool cballs_angular_cell_ok(const struct cmdline_data *cmd,
        const cballs_storage_real *pivot, const compute_vector dr,
        real pivot_radius, real neighbor_radius)
{
#ifdef TPCF
    if (cballs_opt_only_2pcf(cmd)) return TRUE;
    const real extent = cballs_angular_extent(pivot, dr, pivot_radius, neighbor_radius);
    const real orders = cballs_opt_edge_corrections(cmd) ? 2.0*cmd->mChebyshev : cmd->mChebyshev;
    const real tolerance = MIN(0.5*PI, MAX(0.0, cmd->theta)*PI/(2.0*orders+1.0));
    return extent < 1.0 && extent <= sin(tolerance);
#else
    (void)cmd; (void)pivot; (void)dr; (void)pivot_radius; (void)neighbor_radius;
    return TRUE;
#endif
}

#endif
