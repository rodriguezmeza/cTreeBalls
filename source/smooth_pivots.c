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

#ifdef SMOOTHPIVOT
local size_t smooth_hash_cell(const INTEGER cell[NDIM], size_t mask)
{
    unsigned long long hash = 1469598103934665603ULL;
    int k;

    DO_COORD(k) {
        hash ^= (unsigned long long)cell[k];
        hash *= 1099511628211ULL;
    }

    return (size_t)hash & mask;
}

local int smooth_body_cell(struct cmdline_data* cmd,
                           struct global_data* gd,
                           bodyptr p,
                           const real cell_width[NDIM],
                           const INTEGER periodic_cells[NDIM],
                           INTEGER cell[NDIM])
{
    int k;

    DO_COORD(k) {
        real coordinate;

        if (!isfinite(Pos(p)[k])) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "prepare_smooth_pivots: non-finite body coordinate[%d]=%g",
                     k, Pos(p)[k]);
            return FAILURE;
        }

        if (cmd->usePeriodic) {
            real box = gd->Box[k];

            if (!isfinite(box) || box <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "prepare_smooth_pivots: invalid periodic box[%d]=%g",
                         k, box);
                return FAILURE;
            }

            coordinate = Pos(p)[k] + 0.5*box;
            coordinate -= floor(coordinate/box)*box;
            cell[k] = (INTEGER)floor(coordinate/cell_width[k]);
            if (cell[k] >= periodic_cells[k])
                cell[k] = periodic_cells[k] - 1;
        } else {
            cell[k] = (INTEGER)floor(Pos(p)[k]/cell_width[k]);
        }
    }

    return SUCCESS;
}

local int smooth_neighbor_cells(INTEGER center, INTEGER periodic_cells,
                                bool periodic, INTEGER neighbors[3])
{
    int count = 0;
    int offset;

    for (offset = -1; offset <= 1; offset++) {
        INTEGER candidate = center + offset;
        int i;

        if (periodic) {
            candidate %= periodic_cells;
            if (candidate < 0)
                candidate += periodic_cells;
        }

        for (i = 0; i < count; i++)
            if (neighbors[i] == candidate)
                break;
        if (i == count)
            neighbors[count++] = candidate;
    }

    return count;
}

local int smooth_claim_cell(struct cmdline_data* cmd,
                            struct global_data* gd,
                            bodyptr p, bodyptr scan_table,
                            const INTEGER target[NDIM],
                            const INTEGER *body_cells,
                            const INTEGER *bucket_heads,
                            const INTEGER *bucket_next,
                            size_t bucket_mask,
                            cballs_smooth_claim_accumulator accumulator,
                            void *accumulator_context)
{
    INTEGER iq;
    size_t bucket = smooth_hash_cell(target, bucket_mask);

    for (iq = bucket_heads[bucket]; iq >= 0; iq = bucket_next[iq]) {
        bodyptr q;
        real dr1;
        compute_vector dr;
        int k;

        DO_COORD(k)
            if (body_cells[(size_t)iq*NDIM + k] != target[k])
                break;
        if (k != NDIM)
            continue;

        q = scan_table + iq;
        if (p == q)
            continue;
        if (cballs_opt_read_mask(cmd)
            && Mask(q) == MASK_NODE_MASKED)
            continue;
        if (!accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)
            || dr1 > gd->rsmooth[0])
            continue;

        if (Update(q) == TRUE) {
            if (accumulator != NULL
                && accumulator(cmd, gd, p, q, accumulator_context) == FAILURE)
                return FAILURE;
            Update(q) = FALSE;
            NbRmin(p) += 1;
#ifndef NOWKAvg
            KappaRmin(p) += Weight(q)*Kappa(q);
#else
            KappaRmin(p) += Kappa(q);
#endif
            WeightRmin(p) += Weight(q);
#ifdef THREEPCFSHEAR
            if (accumulator == NULL) {
                Gamma1Rmin(p) += Weight(q)*Gamma1(q);
                Gamma2Rmin(p) += Weight(q)*Gamma2(q);
            }
#endif
        } else {
            NbRminOverlap(p) += 1;
        }
    }

    return SUCCESS;
}

/*
 * Build smoothing groups in stable pivot order before an OpenMP search starts.
 * Search workers may read the resulting body fields but must not mutate them.
 */
global int prepare_smooth_pivots_with_accumulator(
                                 struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 bodyptr *btable, INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2,
                                 cballs_smooth_claim_accumulator accumulator,
                                 void *accumulator_context)
{
    INTEGER first, last, npivot, nscan;
    INTEGER *body_cells = NULL;
    INTEGER *bucket_heads = NULL;
    INTEGER *bucket_next = NULL;
    INTEGER periodic_cells[NDIM];
    real cell_width[NDIM];
    size_t bucket_count = 16;
    size_t bucket_mask;
    bodyptr p, q;
    int k;
    int status = FAILURE;

    if (btable == NULL || nbody == NULL || ipmax == NULL
        || cat1 < 0 || cat1 >= gd->ninfiles
        || cat2 < 0 || cat2 >= gd->ninfiles) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid catalog arguments");
        return FAILURE;
    }

    first = ipmin - 1;
    last = ipmax[cat1];
    npivot = nbody[cat1];
    nscan = nbody[cat2];
    if (first < 0 || last < first || last > npivot || nscan < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid pivot or scan bounds");
        return FAILURE;
    }
    if (!isfinite(gd->rsmooth[0]) || gd->rsmooth[0] < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: invalid smoothing radius %g",
                 gd->rsmooth[0]);
        return FAILURE;
    }
    if (!cballs_opt_smooth_pivot(cmd) || gd->rsmooth[0] == 0.0) {
        DO_BODY(q, btable[cat2], btable[cat2] + nscan)
            Update(q) = TRUE;
        if (cat1 != cat2)
            DO_BODY(p, btable[cat1], btable[cat1] + npivot)
                Update(p) = TRUE;
        DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            #ifndef NOWKAvg
            KappaRmin(p) = Weight(p)*Kappa(p);
#else
            KappaRmin(p) = Kappa(p);
#endif
            WeightRmin(p) = Weight(p);
#ifdef THREEPCFSHEAR
            Gamma1Rmin(p) = Weight(p)*Gamma1(p);
            Gamma2Rmin(p) = Weight(p)*Gamma2(p);
#endif
        }
        return SUCCESS;
    }
    if ((size_t)nscan > ((size_t)-1)/(NDIM*sizeof(INTEGER))
        || (size_t)nscan > ((size_t)-1)/2) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: catalog is too large");
        return FAILURE;
    }

    while (bucket_count < (size_t)nscan*2) {
        if (bucket_count > ((size_t)-1)/2) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "prepare_smooth_pivots: hash table size overflow");
            return FAILURE;
        }
        bucket_count <<= 1;
    }
    if (bucket_count > ((size_t)-1)/sizeof(*bucket_heads)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: hash table allocation overflow");
        return FAILURE;
    }
    bucket_mask = bucket_count - 1;

    size_t map_cells, map_bytes, map_heads;
    if (!cballs_size_mul((size_t)nscan,NDIM+1,&map_cells)
        || !cballs_size_mul(map_cells,sizeof(INTEGER),&map_bytes)
        || !cballs_size_mul(bucket_count,sizeof(INTEGER),&map_heads)
        || !cballs_size_add(map_bytes,map_heads,&map_bytes)) {
        snprintf(cmd->error_message,_ERRORMSGSIZE_,"smoothing ownership map dimension overflow");
        goto cleanup;
    }
    if (cballs_memory_preflight(map_bytes,"complete smoothing ownership map",
                               cmd->error_message,_ERRORMSGSIZE_) == FAILURE) goto cleanup;
    if (cballs_malloc_checked((void **)&body_cells,(size_t)nscan*NDIM,sizeof(*body_cells),
            "smoothing cells",cmd->error_message,_ERRORMSGSIZE_) == FAILURE
        || cballs_malloc_checked((void **)&bucket_next,(size_t)nscan,sizeof(*bucket_next),
            "smoothing links",cmd->error_message,_ERRORMSGSIZE_) == FAILURE
        || cballs_malloc_checked((void **)&bucket_heads,bucket_count,sizeof(*bucket_heads),
            "smoothing heads",cmd->error_message,_ERRORMSGSIZE_) == FAILURE) goto cleanup;
    if (body_cells == NULL || bucket_next == NULL || bucket_heads == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "prepare_smooth_pivots: unable to allocate temporary ownership map");
        goto cleanup;
    }

    DO_COORD(k) {
        if (cmd->usePeriodic) {
            real box = gd->Box[k];
            if (!isfinite(box) || box <= 0.0) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "prepare_smooth_pivots: invalid periodic box[%d]=%g",
                         k, box);
                goto cleanup;
            }
            periodic_cells[k] = (INTEGER)floor(box/gd->rsmooth[0]);
            if (periodic_cells[k] < 1)
                periodic_cells[k] = 1;
            cell_width[k] = box/(real)periodic_cells[k];
        } else {
            periodic_cells[k] = 0;
            cell_width[k] = gd->rsmooth[0];
        }
    }

    for (size_t ibucket = 0; ibucket < bucket_count; ibucket++)
        bucket_heads[ibucket] = -1;

    DO_BODY(q, btable[cat2], btable[cat2] + nscan) {
        INTEGER iq = q - btable[cat2];
        INTEGER *cell = body_cells + (size_t)iq*NDIM;
        size_t bucket;

        if (smooth_body_cell(cmd, gd, q, cell_width, periodic_cells, cell)
            == FAILURE)
            goto cleanup;
        bucket = smooth_hash_cell(cell, bucket_mask);
        bucket_next[iq] = bucket_heads[bucket];
        bucket_heads[bucket] = iq;
        Update(q) = TRUE;
    }

    if (cat1 != cat2)
        DO_BODY(p, btable[cat1], btable[cat1] + npivot)
            Update(p) = TRUE;

    DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
        NbRmin(p) = 1;
        NbRminOverlap(p) = 0;
        #ifndef NOWKAvg
            KappaRmin(p) = Weight(p)*Kappa(p);
#else
            KappaRmin(p) = Kappa(p);
#endif
        WeightRmin(p) = Weight(p);
#ifdef THREEPCFSHEAR
        Gamma1Rmin(p) = Weight(p)*Gamma1(p);
        Gamma2Rmin(p) = Weight(p)*Gamma2(p);
#endif
    }

    DO_BODY(p, btable[cat1] + first, btable[cat1] + last) {
        INTEGER center[NDIM];
        INTEGER neighbors[NDIM][3];
        int neighbor_count[NDIM];

        if ((cballs_opt_read_mask(cmd) && Mask(p) == MASK_NODE_MASKED)
            || Update(p) == FALSE)
            continue;
        if (smooth_body_cell(cmd, gd, p, cell_width, periodic_cells, center)
            == FAILURE)
            goto cleanup;
        DO_COORD(k)
            neighbor_count[k] = smooth_neighbor_cells(center[k],
                                                       periodic_cells[k],
                                                       cmd->usePeriodic,
                                                       neighbors[k]);

#if NDIM == 3
        for (int i0 = 0; i0 < neighbor_count[0]; i0++)
            for (int i1 = 0; i1 < neighbor_count[1]; i1++)
                for (int i2 = 0; i2 < neighbor_count[2]; i2++) {
                    INTEGER target[NDIM] = {
                        neighbors[0][i0], neighbors[1][i1], neighbors[2][i2]
                    };
                    if (smooth_claim_cell(
                            cmd, gd, p, btable[cat2], target,
                            body_cells, bucket_heads, bucket_next, bucket_mask,
                            accumulator, accumulator_context) == FAILURE)
                        goto cleanup;
                }
#elif NDIM == 2
        for (int i0 = 0; i0 < neighbor_count[0]; i0++)
            for (int i1 = 0; i1 < neighbor_count[1]; i1++) {
                INTEGER target[NDIM] = {
                    neighbors[0][i0], neighbors[1][i1]
                };
                if (smooth_claim_cell(
                        cmd, gd, p, btable[cat2], target,
                        body_cells, bucket_heads, bucket_next, bucket_mask,
                        accumulator, accumulator_context) == FAILURE)
                    goto cleanup;
            }
#else
#error prepare_smooth_pivots supports only NDIM=2 or NDIM=3
#endif
    }

    status = SUCCESS;

cleanup:
    free(bucket_heads);
    free(bucket_next);
    free(body_cells);
    return status;
}

global int prepare_smooth_pivots(struct cmdline_data* cmd,
                                 struct global_data* gd,
                                 bodyptr *btable, INTEGER *nbody,
                                 INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2)
{
    return prepare_smooth_pivots_with_accumulator(
        cmd, gd, btable, nbody, ipmin, ipmax, cat1, cat2, NULL, NULL);
}
#endif
