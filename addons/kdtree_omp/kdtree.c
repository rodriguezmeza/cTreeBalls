/* ==============================================================================
 MODULE: kdtree.c        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Based on: zeno lib
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: kd = init_kdtree(cmd, gd, btab, nbody, nbody);
      build_kddtree(cmd, gd, nbucket);
      finish_kdtree(kd);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"
#include "kdtree.h"

local void set_radius(struct cmdline_data*, ballnode *, bodyptr *, int, int);
local void set_cofm(struct cmdline_data *, ballnode *, bodyptr *, int, int);
local void set_bounds(bound *, bodyptr *, INTEGER, INTEGER);
local INTEGER  median_index(bodyptr *, int, INTEGER, INTEGER);
local void combine_nodes(struct cmdline_data*,
                         ballnode *, bodyptr *, ballnode *, ballnode *);
local void upward_pass(struct cmdline_data*, ballxptr, int);

// Alloc ball tree memory, set bounds for root
//  check that kd->bptr will consume additional
//  memory to btab->Pos...
//  in such a case change kd-bptr[i] -> Pos(p)
ballxptr init_kdtree(struct cmdline_data* cmd,
                   struct  global_data* gd,
                   bodyptr btab, INTEGER nbody)
{
    ballxptr kd;
    INTEGER i, j;

    if (nbody <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "init_kdtree: catalog is empty");
        return NULL;
    }
    kd = calloc(1, sizeof(*kd));
    if (kd == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "init_kdtree: unable to allocate tree context");
        return NULL;
    }
    kd->npoint = nbody;
    kd->body_base = btab;
    if (cballs_malloc_checked((void **)&kd->bptr, (size_t)nbody,
            sizeof(*kd->bptr), "KDTREE body pointers", cmd->error_message,
            _ERRORMSGSIZE_) == FAILURE) {
        finish_kdtree(kd);
        return NULL;
    }
    gd->bytes_tot += nbody*sizeof(bodyptr);
    verb_print(cmd->verbose,
        "Allocated %g MByte for particle storage in kd tree structure.\n",
               nbody*sizeof(bodyptr)*INMB);
    for (i = j = 0; i < nbody; i++)
        kd->bptr[j++] = nthBody(btab, i);           // ( btab + i )
    if (j != kd->npoint) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "init_kdtree: body pointer publication failed");
        finish_kdtree(kd);
        return NULL;
    }
    set_bounds(&kd->bnd, kd->bptr, 0, kd->npoint - 1);

    return kd;
}

// Free kd tree memory
void finish_kdtree(ballxptr kd)
{
    if (kd == NULL) return;
#ifdef SINGLEP
    free(kd->packed_points);
#endif
    free(kd->bptr);
    free(kd->body_order);
    free(kd->ntab);
    free(kd);
}

// https://www.geeksforgeeks.org/bitwise-operators-in-c-cpp/
// Build ball tree, compute number of nbodies, nodes, and nsplit
//  alloc memory for ball tree nodes and set the main tree walking loop
int build_kdtree(struct cmdline_data* cmd,
                 struct  global_data* gd,
                 ballxptr kd, int nbucket)
{
    INTEGER n, m, j;
    int k, i, d, ct;
    ballnode *ntab;

    //B Find number of nodes and number of times to split
    if (kd == NULL || kd->npoint <= 0 || nbucket <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "build_kdtree: invalid tree or leaf capacity");
        return FAILURE;
    }
    n = kd->npoint;
    k = 1;
    while (n > nbucket) {
        n = n>>1;                                   // a>>b = a/pow(2,b)
        k = k<<1;                                   // a<<b = a*pow(2,b)
    }
    kd->nnode = k<<1;
    kd->nsplit = k;
    //E

    if (cballs_calloc_checked((void **)&kd->ntab, (size_t)kd->nnode,
            sizeof(*kd->ntab), "KDTREE nodes", cmd->error_message,
            _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    ntab = kd->ntab;
    gd->bytes_tot += kd->nnode*sizeof(ballnode);
    verb_print(cmd->verbose,
               "Number of nbodies, nodes, and nsplit: %" INTEGER_FMT " %d %d\n",
               kd->npoint, kd->nnode, kd->nsplit);
    verb_print(cmd->verbose,
               "Allocated %g MByte for particle storage in kd node tab.\n\n",
               kd->nnode*sizeof(ballnode)*INMB);

    //B Initialize root node
    ntab[KDROOT].first = 0;			                // index of first body in root
    ntab[KDROOT].last = kd->npoint-1;               // index of last body in root
    ntab[KDROOT].bnd = kd->bnd;
    i = KDROOT;
    ct = KDROOT;
    SetNext(ct);
    //E

    //B loop splitting nodes starting with KDROOT (i=1)
    //      and set its next (ct) to stop.
    for ( ; ; ) {
        if (i < kd->nsplit) {
            //B find longest dimension
            d = 0;
            DO_COORD(k) {
                if (ntab[i].bnd.maxb[k]-ntab[i].bnd.minb[k] >
                    ntab[i].bnd.maxb[d]-ntab[i].bnd.minb[d])
                    d = k;
            }
            //E

            //B Splitting using median
            m = median_index(kd->bptr, d, ntab[i].first, ntab[i].last);
            //E
            ntab[i].dim = d;
            ntab[i].split = Pos(kd->bptr[m])[d];
            ntab[Lower(i)].bnd = ntab[i].bnd;
            ntab[Lower(i)].bnd.maxb[d] = ntab[i].split;
            ntab[Lower(i)].first = ntab[i].first;
            ntab[Lower(i)].last = m-1;

            ntab[Upper(i)].bnd = ntab[i].bnd;
            ntab[Upper(i)].bnd.minb[d] = ntab[i].split;
            ntab[Upper(i)].first = m;
            ntab[Upper(i)].last = ntab[i].last;

            i = Lower(i);
        } else {
            ntab[i].dim = -1;
            SetNext(i);
            if (i == ct) break;
        }
    }
    //E

    upward_pass(cmd, kd, KDROOT);

    if (cballs_malloc_checked((void **)&kd->body_order,
            (size_t)kd->npoint, sizeof(*kd->body_order),
            "KDTREE body order", cmd->error_message,
            _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    for (j = 0; j < kd->npoint; j++)
        kd->body_order[kd->bptr[j] - kd->body_base] = j;
    gd->bytes_tot += kd->npoint * sizeof(*kd->body_order);

#ifdef SINGLEP
    if (cballs_malloc_checked((void **)&kd->packed_points,
            (size_t)kd->npoint, sizeof(*kd->packed_points),
            "KDTREEOMP packed leaf points", cmd->error_message,
            _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
    for (j = 0; j < kd->npoint; j++) {
        SETV(kd->packed_points[j].pos, Pos(kd->bptr[j]));
        kd->packed_points[j].kappa = Kappa(kd->bptr[j]);
        kd->packed_points[j].weighted_kappa = Weight(kd->bptr[j])*Kappa(kd->bptr[j]);
    }
    gd->bytes_tot += kd->npoint * sizeof(*kd->packed_points);
#endif

    return SUCCESS;
}

//  Compute cell radius
local void set_radius(struct cmdline_data* cmd,
                      ballnode *kd, bodyptr *bptr, int lo, int hi)
{
    int i, k;
    real d, dmax;

    dmax = 0.0;
    DO_COORD(k)
        dmax += rsqr((real)Pos(bptr[lo])[k] - (real)kd->cmpos[k]);

    for (i = lo + 1; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr((real)Pos(bptr[i])[k] - (real)kd->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    kd->bnd.geometric_radius = cballs_store_upper_bound(rsqrt(dmax));
    if (cmd->theta == 0.0) {
        kd->bnd.radius = cballs_store_upper_bound(MAX_REAL_NUMBER);
    } else {
        const real radius = rsqrt(dmax)/cmd->theta;
        kd->bnd.radius = cballs_store_search_bound(radius);
    }
}

//  Compute cell inertia tensor and deformation factor
local void set_inertia(struct cmdline_data *cmd, ballnode *kd,
                       bodyptr *bptr, int lo, int hi)
{
    int i, k;
    int l;

    CLRM(kd->Ixy);

    for (i = lo; i <= hi; ++i) {
        if (cballs_opt_read_mask(cmd)
            && Mask(bptr[i]) == MASK_NODE_MASKED) continue;
        DO_COORD(k) {
            DO_COORD(l) {
                kd->Ixy[k][l] +=
                Mass(bptr[i])*Pos(bptr[i])[k]*Pos(bptr[i])[l];
            }
        }
    }

    //B Computation of deformation factor
    real etap;
    real etax;

    //B etaxy
    real Ixx = kd->Ixy[0][0];
    real Iyy = kd->Ixy[1][1];
    real Ixy = kd->Ixy[0][1];
    const real dxy = Ixx + Iyy;
    etap = dxy != 0.0 ? (Ixx - Iyy)/dxy : 0.0;
    etax = dxy != 0.0 ? 2.0*Ixy/dxy : 0.0;
    kd->etaxy = rsqrt( rsqr(etap) + rsqr(etax) );
#if NDIM == 3
    //B etaxz
    real Izz = kd->Ixy[2][2];
    real Ixz = kd->Ixy[0][2];
    const real dxz = Ixx + Izz;
    etap = dxz != 0.0 ? (Ixx - Izz)/dxz : 0.0;
    etax = dxz != 0.0 ? 2.0*Ixz/dxz : 0.0;
    kd->etaxz = rsqrt( rsqr(etap) + rsqr(etax) );
    //B etayz
    real Iyz = kd->Ixy[1][2];
    const real dyz = Iyy + Izz;
    etap = dyz != 0.0 ? (Iyy - Izz)/dyz : 0.0;
    etax = dyz != 0.0 ? 2.0*Iyz/dyz : 0.0;
    kd->etayz = rsqrt( rsqr(etap) + rsqr(etax) );
#endif
    //E
}

//  Computes cell center of mass and averages scalar fields
local void set_cofm(struct cmdline_data *cmd, ballnode *kd,
                    bodyptr *bptr, int lo, int hi)
{
    compute_vector cmpos_sum;
    int i;
    real KappaAvg = 0.0;

    kd->weight = 0.0;
    kd->weighted_kappa_sum = 0.0;
    kd->weighted_kappa_sq_sum = 0.0;
    kd->kappa_sum = 0.0;
    kd->kappa_sq_sum = 0.0;
    kd->weight_sum = 0.0;
    kd->weight_sq_sum = 0.0;
    kd->valid_count = 0;
    CLRV(cmpos_sum);

    for (i = lo; i <= hi; ++i) {
        if (cballs_opt_read_mask(cmd)
            && Mask(bptr[i]) == MASK_NODE_MASKED) continue;
#ifdef KappaAvgON
        KappaAvg += KappaAvg(bptr[i]);
#else
        KappaAvg += Kappa(bptr[i]);
#endif
        kd->valid_count++;
        kd->kappa_sum += Kappa(bptr[i]);
        kd->kappa_sq_sum += Kappa(bptr[i])*Kappa(bptr[i]);
        kd->weight_sum += Weight(bptr[i]);
        kd->weight_sq_sum += Weight(bptr[i])*Weight(bptr[i]);
        kd->weight += Mass(bptr[i]);
        const real field = Weight(bptr[i])*Kappa(bptr[i]);
        kd->weighted_kappa_sum += field;
        kd->weighted_kappa_sq_sum += field*field;
        int k;
        DO_COORD(k)
            cmpos_sum[k] += Mass(bptr[i]) * (real)Pos(bptr[i])[k];
    }
    if (kd->weight > 0.0) {
        int k;
        DO_COORD(k)
            kd->cmpos[k] = (cballs_storage_real)
                (cmpos_sum[k] / kd->weight);
    } else {
        SETV(kd->cmpos, kd->bnd.center);
    }
    
    kd->kappa = kd->valid_count > 0
        ? KappaAvg/(real)kd->valid_count : 0.0;
}

//  Compute bounds from body pointers in specified range
local void set_bounds(bound *bndptr, bodyptr *bptr, INTEGER lo, INTEGER hi)
{
    int k, i;
    bound bnd;
    DO_COORD(k)                                         // initialize bounds
        bnd.maxb[k] =  bnd.minb[k] = Pos(bptr[lo])[k];

    for (i = lo + 1; i <= hi; ++i) {                        // find actual bounds
        DO_COORD(k) {
        if (bnd.minb[k] > Pos(bptr[i])[k])
            bnd.minb[k] = Pos(bptr[i])[k];
        else if (bnd.maxb[k] < Pos(bptr[i])[k])
            bnd.maxb[k] = Pos(bptr[i])[k];
      }
    }

    DO_COORD(k) {
        bnd.width[k] = bnd.maxb[k] - bnd.minb[k];
        bnd.center[k] = 0.5*(bnd.maxb[k] + bnd.minb[k]);
    }

    *bndptr = bnd;				                        // store actual bounds
}

// Partly sort body pointers in a specified range,
//      and return index of median, using JST's median algorithm
#define SwapBody(b1,b2)  { bodyptr _tmp; _tmp = b1; b1 = b2; b2 = _tmp; }

local INTEGER median_index(bodyptr *p, int d, INTEGER lo, INTEGER hi)
{
    INTEGER i, j, m;
    real f;

    m = j = (lo + hi) / 2;
    while (lo < hi) {
        m = (lo + hi) / 2;
        f = Pos(p[m])[d];
        SwapBody(p[m], p[hi]);
        i = hi - 1;
        m = lo;
        while (Pos(p[m])[d] < f)
            ++m;
        while (m < i) {
            while (Pos(p[i])[d] >= f)
                if (--i == m)
                    break;
            SwapBody(p[m], p[i]);
            --i;
            while (Pos(p[m])[d] < f)
                ++m;
        }
        SwapBody(p[m], p[hi]);
        if (j <= m)
            hi = m - 1;
        if (j >= m)
            lo = m + 1;
    }

  return m;
}

#undef SwapBody


//  Adjust bounds of each node to fit bodies exactly
local void upward_pass(struct cmdline_data* cmd, ballxptr kd, int cell)
{
    ballnode *ntab = kd->ntab;
//    int d;
//    int k;
//    real radius;
    bodyptr *bptr = kd->bptr;

    if (ntab[cell].dim != -1) {                        // not a terminal node?
        upward_pass(cmd, kd, Lower(cell));
        upward_pass(cmd, kd, Upper(cell));
        combine_nodes(cmd,
                      &ntab[cell], bptr, &ntab[Lower(cell)], &ntab[Upper(cell)]);
    } else {                                        // scan bodies in node
        set_bounds(&ntab[cell].bnd, kd->bptr,
                   ntab[cell].first, ntab[cell].last);
        set_cofm(cmd, &ntab[cell], kd->bptr,
                 ntab[cell].first, ntab[cell].last);
        set_radius(cmd, &ntab[cell], kd->bptr,
                   ntab[cell].first, ntab[cell].last);
        set_inertia(cmd, &ntab[cell], kd->bptr,
                   ntab[cell].first, ntab[cell].last);
    }
}

//  Combine two nodes: bounding, center of mass pos, and radius
local void combine_nodes(struct cmdline_data* cmd,
                         ballnode *pout, bodyptr *bptr,
                         ballnode *p1, ballnode *p2)
{
    int k;

    // AABB (Minimum Axis-aligned Bounding Box)
    DO_COORD(k) {
        pout->bnd.minb[k] = MIN(p2->bnd.minb[k], p1->bnd.minb[k]);
        pout->bnd.maxb[k] = MAX(p2->bnd.maxb[k], p1->bnd.maxb[k]);
    }

    // width and geometric center
    DO_COORD(k) {
        pout->bnd.width[k] = (pout->bnd.maxb[k]-pout->bnd.minb[k]);
        pout->bnd.center[k] = (pout->bnd.maxb[k]+pout->bnd.minb[k])*0.5;
    }

    // cmpos
    compute_vector cmpos_sum;
    pout->weight = 0.0;
    CLRV(cmpos_sum);

    pout->weight += p1->weight;
    DO_COORD(k)
        cmpos_sum[k] += (real)p1->cmpos[k] * p1->weight;

    pout->weight += p2->weight;
    DO_COORD(k)
        cmpos_sum[k] += (real)p2->cmpos[k] * p2->weight;

    if (pout->weight > 0.0) {
        DO_COORD(k)
            pout->cmpos[k] = (cballs_storage_real)
                (cmpos_sum[k] / pout->weight);
    } else {
        SETV(pout->cmpos, pout->bnd.center);
    }

    // KappaAvg
    pout->valid_count = p1->valid_count + p2->valid_count;
    pout->kappa = pout->valid_count > 0
        ? (p1->kappa_sum+p2->kappa_sum)/(real)pout->valid_count : 0.0;
    pout->weighted_kappa_sum = p1->weighted_kappa_sum + p2->weighted_kappa_sum;
    pout->weighted_kappa_sq_sum = p1->weighted_kappa_sq_sum + p2->weighted_kappa_sq_sum;
    pout->kappa_sum = p1->kappa_sum + p2->kappa_sum;
    pout->kappa_sq_sum = p1->kappa_sq_sum + p2->kappa_sq_sum;
    pout->weight_sum = p1->weight_sum + p2->weight_sum;
    pout->weight_sq_sum = p1->weight_sq_sum + p2->weight_sq_sum;

    // radius
    int i;
    int lo, hi;
    real d, dmax;

    lo = p1->first;
    hi = p1->last;

    dmax = 0.0;
    DO_COORD(k)
        dmax += rsqr((real)Pos(bptr[lo])[k] - (real)pout->cmpos[k]);

    for (i = lo + 1; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr((real)Pos(bptr[i])[k] - (real)pout->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    lo = p2->first;
    hi = p2->last;

    for (i = lo; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr((real)Pos(bptr[i])[k] - (real)pout->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    pout->bnd.geometric_radius = cballs_store_upper_bound(rsqrt(dmax));
    if (cmd->theta == 0.0) {
        pout->bnd.radius = cballs_store_upper_bound(MAX_REAL_NUMBER);
    } else {
        const real radius = rsqrt(dmax)/cmd->theta;
        pout->bnd.radius = cballs_store_search_bound(radius);
    }
    //B Computation of inertia and deformation factor
    int l;

    CLRM(pout->Ixy);

    lo = p1->first;
    hi = p1->last;
    for (i = lo; i <= hi; ++i) {
        DO_COORD(k) {
            DO_COORD(l) {
                if (!cballs_opt_read_mask(cmd)
                    || Mask(bptr[i]) != MASK_NODE_MASKED)
                pout->Ixy[k][l] +=
                            Mass(bptr[i])*Pos(bptr[i])[k]*Pos(bptr[i])[l];
            }
        }
    }

    lo = p2->first;
    hi = p2->last;
    for (i = lo; i <= hi; ++i) {
        DO_COORD(k) {
            DO_COORD(l) {
                if (!cballs_opt_read_mask(cmd)
                    || Mask(bptr[i]) != MASK_NODE_MASKED)
                pout->Ixy[k][l] +=
                            Mass(bptr[i])*Pos(bptr[i])[k]*Pos(bptr[i])[l];
            }
        }
    }

    //B Computation of deformation factor
    real etap;
    real etax;

    //B etaxy
    real Ixx = pout->Ixy[0][0];
    real Iyy = pout->Ixy[1][1];
    real Ixy = pout->Ixy[0][1];
    const real dxy = Ixx + Iyy;
    etap = dxy != 0.0 ? (Ixx - Iyy)/dxy : 0.0;
    etax = dxy != 0.0 ? 2.0*Ixy/dxy : 0.0;
    pout->etaxy = rsqrt( rsqr(etap) + rsqr(etax) );
#if NDIM == 3
    //B etaxz
    real Izz = pout->Ixy[2][2];
    real Ixz = pout->Ixy[0][2];
    const real dxz = Ixx + Izz;
    etap = dxz != 0.0 ? (Ixx - Izz)/dxz : 0.0;
    etax = dxz != 0.0 ? 2.0*Ixz/dxz : 0.0;
    pout->etaxz = rsqrt( rsqr(etap) + rsqr(etax) );
    //B etayz
    real Iyz = pout->Ixy[1][2];
    const real dyz = Iyy + Izz;
    etap = dyz != 0.0 ? (Iyy - Izz)/dyz : 0.0;
    etax = dyz != 0.0 ? 2.0*Iyz/dyz : 0.0;
    pout->etayz = rsqrt( rsqr(etap) + rsqr(etax) );
#endif
    //E
    //E
}
