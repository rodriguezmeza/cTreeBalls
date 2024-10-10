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

#include <assert.h>
#include "globaldefs.h"
#include "kdtree.h"

local void set_radius(struct cmdline_data*, ballnode *, bodyptr *, int, int);
local void set_cofm(ballnode *, bodyptr *, int, int);
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

    kd = (ballxptr) allocate(sizeof(ballcontext));
    kd->npoint = nbody;
    kd->bptr = (bodyptr *) allocate(nbody * sizeof(bodyptr));
    gd->bytes_tot += nbody*sizeof(bodyptr);
    verb_print(cmd->verbose,
        "Allocated %g MByte for particle storage in kd tree structure.\n",
               nbody*sizeof(bodyptr)*INMB);
    for (i = j = 0; i < nbody; i++)
        kd->bptr[j++] = nthBody(btab, i);           // ( btab + i )
    assert(j == kd->npoint);                        // Check counting all
                                                    //  Terminates run if not
    set_bounds(&kd->bnd, kd->bptr, 0, kd->npoint - 1);

    return kd;
}

// Free kd tree memory
void finish_kdtree(ballxptr kd)
{
    free(kd->bptr);
    free(kd->ntab);
    free(kd);
}

// https://www.geeksforgeeks.org/bitwise-operators-in-c-cpp/
// Build ball tree, compute number of nbodies, nodes, and nsplit
//  alloc memory for ball tree nodes and set the main tree walking loop
void build_kdtree(struct cmdline_data* cmd,
                  struct  global_data* gd,
                  ballxptr kd, int nbucket)
{
    INTEGER n, m;
    int k, i, d, ct;
    ballnode *ntab;

    //B Find number of nodes and number of times to split
    n = kd->npoint;
    k = 1;
    while (n > nbucket) {
        n = n>>1;                                   // a>>b = a/pow(2,b)
        k = k<<1;                                   // a<<b = a*pow(2,b)
    }
    kd->nnode = k<<1;
    kd->nsplit = k;
    //E

    ntab = kd->ntab = (ballnode *) allocate(kd->nnode * sizeof(ballnode));
    gd->bytes_tot += kd->nnode*sizeof(ballnode);
    verb_print(cmd->verbose,
               "Number of nbodies, nodes, and nsplit: %d %d %d\n",
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
}

//  Compute cell radius
local void set_radius(struct cmdline_data* cmd,
                      ballnode *kd, bodyptr *bptr, int lo, int hi)
{
    int i, k;
    real d, dmax;

    dmax = 0.0;
    DO_COORD(k)
        dmax += rsqr(Pos(bptr[lo])[k] - kd->cmpos[k]);

    for (i = lo + 1; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr(Pos(bptr[i])[k] - kd->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    kd->bnd.radius = rsqrt(dmax)/cmd->theta;        // set radius
                                                    //  and rescale it
}

//  Compute cell inertia tensor and deformation factor
local void set_inertia(ballnode *kd, bodyptr *bptr, int lo, int hi)
{
    int i, k;
    int l;

    CLRM(kd->Ixy);

    for (i = lo; i <= hi; ++i) {
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
    etap = (Ixx - Iyy)/(Ixx + Iyy);
    etax = 2.0*Ixy/(Ixx + Iyy);
    kd->etaxy = rsqrt( rsqr(etap) + rsqr(etax) );
#if NDIM == 3
    //B etaxz
    real Izz = kd->Ixy[2][2];
    real Ixz = kd->Ixy[0][2];
    etap = (Ixx - Izz)/(Ixx + Izz);
    etax = 2.0*Ixz/(Ixx + Izz);
    kd->etaxz = rsqrt( rsqr(etap) + rsqr(etax) );
    //B etayz
    real Iyz = kd->Ixy[1][2];
    etap = (Iyy - Izz)/(Iyy + Izz);
    etax = 2.0*Iyz/(Iyy + Izz);
    kd->etayz = rsqrt( rsqr(etap) + rsqr(etax) );
#endif
    //E
}

//  Computes cell center of mass and averages scalar fields
local void set_cofm(ballnode *kd, bodyptr *bptr, int lo, int hi)
{
    vector tmpv;
    int i;
    real KappaAvg = 0.0;

    kd->weight = 0.0;
    CLRV(kd->cmpos);

    for (i = lo; i <= hi; ++i) {
        KappaAvg += KappaAvg(bptr[i]);
        kd->weight += Mass(bptr[i]);
        MULVS(tmpv, Pos(bptr[i]), Mass(bptr[i]));
        ADDV(kd->cmpos, kd->cmpos, tmpv);
    }
    if (kd->weight > 0.0) {
        DIVVS(kd->cmpos, kd->cmpos, kd->weight);
    } else {
        SETV(kd->cmpos, kd->bnd.center);
    }
    
    kd->kappa = KappaAvg/((real)(hi-lo+1));
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
        set_cofm(&ntab[cell], kd->bptr,
                 ntab[cell].first, ntab[cell].last);
        set_radius(cmd, &ntab[cell], kd->bptr,
                   ntab[cell].first, ntab[cell].last);
        set_inertia(&ntab[cell], kd->bptr,
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
    vector tmpv;
    pout->weight = 0.0;
    CLRV(pout->cmpos);

    pout->weight += p1->weight;
    MULVS(tmpv, p1->cmpos, p1->weight);
    ADDV(pout->cmpos, pout->cmpos, tmpv);

    pout->weight += p2->weight;
    MULVS(tmpv, p2->cmpos, p2->weight);
    ADDV(pout->cmpos, pout->cmpos, tmpv);

    if (pout->weight > 0.0) {
        DIVVS(pout->cmpos, pout->cmpos, pout->weight);
    } else {
        SETV(pout->cmpos, pout->bnd.center);
    }

    // KappaAvg
    pout->kappa = 0.5*(p1->kappa + p2->kappa);

    // radius
    int i;
    int lo, hi;
    real d, dmax;

    lo = p1->first;
    hi = p1->last;

    dmax = 0.0;
    DO_COORD(k)
        dmax += rsqr(Pos(bptr[lo])[k] - pout->cmpos[k]);

    for (i = lo + 1; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr(Pos(bptr[i])[k] - pout->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    lo = p2->first;
    hi = p2->last;

    for (i = lo; i <= hi; ++i) {
        d = 0.0;
        DO_COORD(k)
            d += rsqr(Pos(bptr[i])[k] - pout->cmpos[k]);
        if (d > dmax)
            dmax = d;
    }

    pout->bnd.radius = rsqrt(dmax)/cmd->theta;              // set radius
                                                            //  and rescale it
    //B Computation of inertia and deformation factor
    int l;

    CLRM(pout->Ixy);

    lo = p1->first;
    hi = p1->last;
    for (i = lo; i <= hi; ++i) {
        DO_COORD(k) {
            DO_COORD(l) {
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
    etap = (Ixx - Iyy)/(Ixx + Iyy);
    etax = 2.0*Ixy/(Ixx + Iyy);
    pout->etaxy = rsqrt( rsqr(etap) + rsqr(etax) );
#if NDIM == 3
    //B etaxz
    real Izz = pout->Ixy[2][2];
    real Ixz = pout->Ixy[0][2];
    etap = (Ixx - Izz)/(Ixx + Izz);
    etax = 2.0*Ixz/(Ixx + Izz);
    pout->etaxz = rsqrt( rsqr(etap) + rsqr(etax) );
    //B etayz
    real Iyz = pout->Ixy[1][2];
    etap = (Iyy - Izz)/(Iyy + Izz);
    etax = 2.0*Iyz/(Iyy + Izz);
    pout->etayz = rsqrt( rsqr(etap) + rsqr(etax) );
#endif
    //E
    //E
}
