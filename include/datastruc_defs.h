/*==============================================================================
 HEADER: data_struc_defs.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definition of N-Body data structure
 Language: C
 Use: '#include "...."
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#ifndef _data_struc_defs_h
#define _data_struc_defs_h

#include "mask_utils.h"

typedef struct _node {
    short type;
    bool update;
    bool update2;
    short mask;                                     // MASK_NODE_* state below
//B look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    bool updatepivot;
#endif
//E

    real mass;                                      // to weight counts...
    real kappa;                                     // scalar field of interest

#ifdef THREEPCFSHEAR
    real gamma1;
    real gamma2;
#endif

    real weight;                                    // to weight fields...

    vector pos;

    struct _node *next;

#ifdef bhistON
    int *bhistNsub;
    real *bhistXi2pcfsub;
    real **histXi;
    bool histON;
#endif

    bool selected;                                  // To see the bodies
                                                    //  belonging to a cell
#ifdef BODY3ON
    INTEGER nbb;                                    // BODY3:
                                                    // If comes from smoothing
                                                    //  gives number of
                                                    //  smoothingbodies
                                                    // Body will be tagged NBODY3
#endif

//B BALLS
    int lev;
    int idxscanlev;
#ifdef IdNodeON
    INTEGER Id;
#endif
#ifdef DEBUG
    bool hit;
#endif

    INTEGER nb;
    cballs_storage_real radius;

#ifdef SMOOTHPIVOT
    real kapparmin;                                 // Sum_p-in-rmin kappa_p
    real weightrmin;                                // to weight fields...
    INTEGER nbrmin;
    INTEGER nbrmin_overlap;
#endif
//E

//B socket:
#ifdef ADDONS
#include "datastruc_defs_include_00.h"
#endif
//E

} node, *nodeptr;

#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Update2(x) (((nodeptr) (x))->update2)
#define Mask(x) (((nodeptr) (x))->mask)

//B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
#define UpdatePivot(x) (((nodeptr) (x))->updatepivot)
#endif
//E

#define Mass(x)   (((nodeptr) (x))->mass)
#define Kappa(x)    (((nodeptr) (x))->kappa)

#ifdef THREEPCFSHEAR
//B shear
#define Gamma1(x)    (((nodeptr) (x))->gamma1)
#define Gamma2(x)    (((nodeptr) (x))->gamma2)
//E
#endif

#define Weight(x)   (((nodeptr) (x))->weight)


#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)
#define bNsub(x)    (((nodeptr) (x))->bhistNsub)
#define bXi2pcfsub(x)    (((nodeptr) (x))->bhistXi2pcfsub)
#define hXi(x)    (((nodeptr) (x))->histXi)
#define hON(x)    (((nodeptr) (x))->histON)

#ifdef BODY3ON
#define Nbb(x)    (((nodeptr) (x))->nbb)            // BODY3
#endif

//B BALLS
#ifdef IdNodeON
#define IdNode(x)    (((nodeptr) (x))->Id)
#endif
#define Level(x)    (((nodeptr) (x))->lev)
#define IDXSCAN(x)    (((nodeptr) (x))->idxscanlev)
#ifdef DEBUG
#define HIT(x)    (((nodeptr) (x))->hit)
#endif
//E

//B To see the bodies belonging to a cell:
#define Selected(x)    (((nodeptr) (x))->selected)
//E

//B Balls-correction.
#define Nb(x) (((nodeptr) (x))->nb)
#define Radius(x) (((nodeptr) (x))->radius)

#ifdef SMOOTHPIVOT
#define KappaRmin(x)   (((nodeptr) (x))->kapparmin) // Sum_p-in-rmin kappa_p
#define WeightRmin(x)  (((nodeptr) (x))->weightrmin)// Sum_p-in-rmin weight_p
#define NbRmin(x) (((nodeptr) (x))->nbrmin)
#define NbRminOverlap(x) (((nodeptr) (x))->nbrmin_overlap)
#endif
//E

#define BODY 00
#define BODY3 03                                    // Smooth
#define CELL 02
// BALLS
#define NODEBODY 04
#define NODECELL 05
//E


typedef struct {
    node bodynode;
    INTEGER Id;

//B socket:
#ifdef ADDONS
#include "datastruc_defs_include_01.h"
#endif
//E

} body, *bodyptr;

#define Id(x)    (((bodyptr) (x))->Id)
#define nthBody(bp,n)  ((bp) + (n))

#define NSUB (1 << NDIM)

// The meaning of the structure and its components can be changed
typedef struct {
    node cellnode;
    cballs_storage_real size;
    nodeptr more;
    nodeptr subp[NSUB];
#ifdef IdCellON
    INTEGER Id;
#endif
    bool inside;
} cell, *cellptr;
 
//B To debug cells:
#define Size(x) (((cellptr) (x))->size)
//E
#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->subp)
#ifdef IdNodeON
#define IdCell(x)   (((cellptr) (x))->Id)
#endif
#define Inside(x)   (((cellptr) (x))->inside)


//B socket:
#ifdef ADDONS
#include "datastruc_defs_include_02.h"
#endif
//E

//B Tree search

// Alternative definition for VWrap. (0,0,..) center of the box

#define VWrap1(v, t)                                            \
   if (v[t] >= 0.5 * gd->Box[t])      v[t] -= gd->Box[t];       \
   else if (v[t] < -0.5 * gd->Box[t]) v[t] += gd->Box[t]

// Alternative definition for VWrap. (0,0,..) lower edge of the box

#define VWrap2(v, t)                                            \
   if (v[t] >= gd->Box[t])      v[t] -= gd->Box[t];             \
   else if (v[t] < 0.0) v[t] += gd->Box[t]

//Chose one PBC::
// By now it is only working with boxes centered at (0,0,...)
// then choose VWrap1

#define VWrap       VWrap1
//#define VWrap       VWrap2

#if NDIM == 2
#define VWrapAll(v)                                             \
   {VWrap (v, 0);                                               \
   VWrap (v, 1);}
#endif

#if NDIM == 3
#define VWrapAll(v)                                             \
   {VWrap (v, 0);                                               \
   VWrap (v, 1);                                                \
   VWrap (v, 2);}
#endif

//E ! Tree search




//B Macro useful to compute chebyshev polynomials
#define CHEBYSHEVTUOMPSINCOS                                            \
{REAL xicosmphi, xisinmphi; int m, mmax;                                \
    mmax = cmd->mChebyshev + 1;                                         \
    if (mmax >= 1) {                                                    \
        hist->ChebsT[1] = 1.0;                                          \
        hist->histXithreadcos[1][n] += xi * hist->ChebsT[1];            \
        hist->ChebsU[1] = 0.0;                                          \
        hist->histXithreadsin[1][n] += xi * hist->ChebsU[1] * sinphi;   \
    }                                                                   \
    if (mmax >= 2) {                                                    \
        hist->ChebsT[2] = cosphi;                                       \
        hist->histXithreadcos[2][n] += xi * hist->ChebsT[2];            \
        hist->ChebsU[2] = 1.0;                                          \
        hist->histXithreadsin[2][n] += xi * hist->ChebsU[2] * sinphi;   \
    }                                                                   \
    if (mmax >= 3) {                                                    \
        hist->ChebsT[3] = 2.0*cosphi*cosphi - 1.0;                      \
        hist->histXithreadcos[3][n] += xi * hist->ChebsT[3];            \
        hist->ChebsU[3] = 2.0*cosphi;                                   \
        hist->histXithreadsin[3][n] += xi * hist->ChebsU[3] * sinphi;   \
    }                                                                   \
    for (m = 4; m <= mmax; m++) {                                       \
        hist->ChebsT[m] = 2.0*cosphi*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        hist->histXithreadcos[m][n] += xi * hist->ChebsT[m];            \
        hist->ChebsU[m] = 2.0*cosphi*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        hist->histXithreadsin[m][n] += xi * hist->ChebsU[m] * sinphi;   \
    }}

#define CHEBYSHEVTUOMPSINCOSANY CHEBYSHEVTUOMPSINCOS
//E


#define DO_BODY(p,start,finish)     for (p = start; p < finish; p++)
#define DO_DESCENDENTS(q,p)         for (q = More(p); q != Next(p); q = Next(q))
#define DO_COORD(k)                 for ((k)=0; (k)<NDIM; (k)++)

#define DO_SUNS(k)                 for (k=0; k < ((int)rpow(2,NDIM)); k++)

//B Input file formats:
#define INCOLUMNS       0
#define INCOLUMNSALL    10
#define INCOLUMNSBIN    2
#define INCOLUMNSBINALL 11
#define INNULL          1
//B version 1.0.1
#define INTAKAHASHI      3                           // Takahashi
//E

//B Output file formats:
#define OUTCOLUMNS       0
#define OUTCOLUMNSALL    3
#define OUTCOLUMNSBIN    2
#define OUTCOLUMNSBINALL 4
#define OUTNULL          1
//E


//B Search options:
#define SEARCHNULL              0
#define OCTREESINCOSOMPMETHOD     24
//E

//Rotation angle in radians. To use in a sphere (3D case)
#define ROTANGLE                0.01

//B socket:
#ifdef ADDONS
#include "datastruc_defs_include_03.h"
#endif
//E

#endif // ! _data_struc_defs_h

