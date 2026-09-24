/* ==============================================================================
 MODULE: cballsutils.c			    [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "tree_contracts.h"

#include <limits.h>
#include <stdint.h>
#include <string.h>

#include <errno.h>
#include <sys/wait.h>

global int cballs_compiled_ndim(void)
{
    return NDIM;
}

global int cballs_max_memory_catalogs(void)
{
    return MAXITEMS;
}

global int cballs_load_memory_catalog(struct cmdline_data *cmd,
                                      struct global_data *gd,
                                      int ifile,
                                      const double *positions,
                                      size_t nbody,
                                      int ndim,
                                      const double *kappa,
                                      const double *weights,
                                      const unsigned char *mask,
                                      const double *gamma1,
                                      const double *gamma2)
{
    bodyptr catalog;
    size_t i;
    int k;

    if (cmd == NULL || gd == NULL) return FAILURE;
    if (ifile < 0 || ifile >= MAXITEMS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog index %d is outside [0, %d)",
                 ifile, MAXITEMS);
        return FAILURE;
    }
    if (ifile == 0) {
        for (k = 0; k < MAXITEMS; k++) {
            if (bodytable[k] != NULL) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "cannot replace a live C catalog; clean the run first");
                return FAILURE;
            }
        }
        gd->ninfiles = 0;
    }
    if (ifile != gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalogs must be loaded contiguously from index 0");
        return FAILURE;
    }
    if (positions == NULL || ndim != NDIM) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory positions must have compiled dimension %d", NDIM);
        return FAILURE;
    }
    if (nbody < 3) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d needs at least 3 bodies", ifile);
        return FAILURE;
    }
#ifdef LONGINT
    if (nbody > (size_t)LONG_MAX) {
#else
    if (nbody > (size_t)INT_MAX) {
#endif
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d is too large for INTEGER", ifile);
        return FAILURE;
    }
    if (nbody > SIZE_MAX / sizeof(body)
        || nbody > (size_t)(LONG_MAX / (long)sizeof(body))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "in-memory catalog %d byte size is not representable", ifile);
        return FAILURE;
    }
    if ((gamma1 == NULL) != (gamma2 == NULL)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "gamma1 and gamma2 must either both be supplied or both omitted");
        return FAILURE;
    }
#ifndef THREEPCFSHEAR
    if (gamma1 != NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "gamma arrays require a THREEPCFSHEAR build");
        return FAILURE;
    }
#endif

    for (i = 0; i < nbody; i++) {
        for (k = 0; k < NDIM; k++) {
            if (!isfinite(positions[i * (size_t)NDIM + (size_t)k])) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "in-memory catalog %d contains a non-finite position",
                         ifile);
                return FAILURE;
            }
        }
        if (kappa != NULL && !isfinite(kappa[i])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite kappa", ifile);
            return FAILURE;
        }
        if (weights != NULL && !isfinite(weights[i])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite weights", ifile);
            return FAILURE;
        }
        if (mask != NULL && mask[i] > 1) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d mask values must be 0 or 1", ifile);
            return FAILURE;
        }
#ifdef THREEPCFSHEAR
        if (gamma1 != NULL
            && (!isfinite(gamma1[i]) || !isfinite(gamma2[i]))) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "in-memory catalog %d contains non-finite shear", ifile);
            return FAILURE;
        }
#endif
    }

    if (cballs_calloc_checked((void **)&catalog,nbody,sizeof(body),
                             "in-memory catalog",cmd->error_message,_ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    if (catalog == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "not enough memory for in-memory catalog %d", ifile);
        return FAILURE;
    }

    for (i = 0; i < nbody; i++) {
        bodyptr p = catalog + i;
        Type(p) = BODY;
        Update(p) = TRUE;
        Update2(p) = TRUE;
        Mask(p) = mask == NULL || mask[i] ? MASK_NODE_VALID
                                          : MASK_NODE_MASKED;
        Mass(p) = 1.0;
        Kappa(p) = kappa == NULL ? 1.0 : (REAL)kappa[i];
        Weight(p) = weights == NULL ? 1.0 : (REAL)weights[i];
        for (k = 0; k < NDIM; k++)
            Pos(p)[k] = (REAL)positions[i * (size_t)NDIM + (size_t)k];
#ifdef THREEPCFSHEAR
        Gamma1(p) = gamma1 == NULL ? 1.0 : (REAL)gamma1[i];
        Gamma2(p) = gamma2 == NULL ? 1.0 : (REAL)gamma2[i];
#endif
        Id(p) = (INTEGER)i + 1;
    }

    bodytable[ifile] = catalog;
    gd->nbodyTable[ifile] = (INTEGER)nbody;
    gd->ninfiles = ifile + 1;
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += (INTEGER)(nbody * sizeof(body));
    cmd->nbody = (INTEGER)nbody;
    for (k = 0; k < NDIM; k++)
        gd->Box[k] = cmd->lengthBox;
    gd->input_comment = "Python in-memory catalog";
    return SUCCESS;
}

global int cballs_set_memory_forest_ids(struct cmdline_data *cmd,
                                        struct global_data *gd, int ifile,
                                        const int64_t *forest_ids, size_t nbody)
{
#if (defined(LYAFORESTOMP) || defined(LYAFORESTMPI)) && NDIM == 3
    size_t i;
    int k;
    REAL minimum[NDIM], maximum[NDIM];
    if (cmd == NULL || gd == NULL) return FAILURE;
    if (ifile < 0 || ifile >= gd->ninfiles || forest_ids == NULL
        || bodytable[ifile] == NULL || nbody != (size_t)gd->nbodyTable[ifile]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "invalid in-memory forest catalog");
        return FAILURE;
    }
    for (k = 0; k < NDIM; k++)
        minimum[k] = maximum[k] = Pos(bodytable[ifile])[k];
    for (i = 0; i < nbody; i++) {
        bodyptr p = bodytable[ifile] + i;
        REAL distance = hypot(hypot(Pos(p)[0], Pos(p)[1]), Pos(p)[2]);
#ifdef LONGINT
        const int64_t id_min = LONG_MIN, id_max = LONG_MAX;
#else
        const int64_t id_min = INT_MIN, id_max = INT_MAX;
#endif
        if (forest_ids[i] < id_min || forest_ids[i] > id_max
            || !isfinite(distance) || distance <= 0.0
            || !isfinite(Weight(p)) || Weight(p) < 0.0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "forest pixels require representable INTEGER IDs, "
                     "positive observer distance and non-negative weights");
            return FAILURE;
        }
        LyaForestId(p) = (INTEGER)forest_ids[i];
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
        Octree3pcf3dLosId(p) = (INTEGER)forest_ids[i];
#endif
        LyaDistance(p) = distance;
        for (k = 0; k < NDIM; k++) {
            LyaLOS(p)[k] = Pos(p)[k] / distance;
            minimum[k] = MIN(minimum[k], Pos(p)[k]);
            maximum[k] = MAX(maximum[k], Pos(p)[k]);
        }
    }
    for (k = 0; k < NDIM; k++) gd->Box[k] = maximum[k] - minimum[k];
#if defined(OCTREE3PCF3DOMP) || defined(OCTREE3PCF3DMPI)
    gd->octree3pcf3d_los_ids[ifile] = TRUE;
#endif
    gd->input_comment = "Python in-memory Lyman-alpha forest catalog";
    return SUCCESS;
#else
    (void)gd; (void)ifile; (void)forest_ids; (void)nbody;
    if (cmd != NULL)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "forest catalogs require LYAFORESTOMP or LYAFORESTMPI and NDIM=3");
    return FAILURE;
#endif
}
