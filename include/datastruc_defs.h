/*==============================================================================
 HEADER: data_struc_defs.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definition of N-Body data structure
 Language: C
 Use: '#include "...."
 ==============================================================================*/
//        1          2          3          4          5          6          7
 
#ifndef _data_struc_defs_h
#define _data_struc_defs_h

#define IPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&(param);												\
  id[nt++]=INT;}

#define RPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)										\
  {strcpy(tag[nt],paramtext);										\
  addr[nt]=&param;													\
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)									\
  {strcpy(tag[nt],paramtext);										\
  param=(string) malloc(n);											\
  addr[nt]=param;													\
  id[nt++]=STRING;}


typedef struct _node {
    short type;
    bool update;
    REAL weight;                                    // The weight can be other
                                                    //  scalar field of interest
    REAL kappa;
    vector pos;
    struct _node *next;
#ifdef bhistON
    int *bhistNsub;
    REAL *bhistXi2pcfsub;
    REAL **histXi;
    bool histON;
#endif
    
//B kd-tree
    byte flags;                                     // status in b-b calculation
    byte newlevel;                                  //  level set by latest
                                                    //  force calc.
//E
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
    REAL radius;

//    INTEGER nbrmin;
//    INTEGER nbrmin_overlap;
//E

#ifdef KappaAvgON
    REAL kappaavg;
#endif

} node, *nodeptr;

#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Weight(x)   (((nodeptr) (x))->weight)
#define Kappa(x)    (((nodeptr) (x))->kappa)

#ifdef KappaAvgON
#define KappaAvg(x)    (((nodeptr) (x))->kappaavg)
#endif

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

//#define NbRmin(x) (((nodeptr) (x))->nbrmin)
//#define NbRminOverlap(x) (((nodeptr) (x))->nbrmin_overlap)
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
//B kd-tree
    REAL smooth;                                    // smoothing length
//E
} body, *bodyptr;

#define Id(x)    (((bodyptr) (x))->Id)

//B kd-tree
#define Flags(x)     (((nodeptr) (x))->flags)
#define NewLevel(x)  (((nodeptr) (x))->newlevel)
#define NthBody(bp,n)  ((bp) + (n))
#define INQUE    0x04            // body listed in current priority que
#define DONE     0x08            // smoothing operation is complete
#define Smooth(x)    (((bodyptr) (x))->smooth)
#define SetFlag(x,f)  (Flags(x) |= (f))
#define ClrFlag(x,f)  (Flags(x) &= ~(f))
#define InQue(x)    ((Flags(x) & INQUE) != 0)
#define Done(x)     ((Flags(x) & DONE) != 0)
//E

#define NSUB (1 << NDIM)

// The meaning of the structure and its components can be changed
typedef struct {
    node cellnode;
    REAL size;
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

#if !defined(global)                                // global def question must
#  define global extern                             //  be here
#endif

//B Tree search

// Alternative definition for VWrap. (0,0,..) center of the box

#define VWrap1(v, t)                                         \
   if (v[t] >= 0.5 * gd.Box[t])      v[t] -= gd.Box[t];     \
   else if (v[t] < -0.5 * gd.Box[t]) v[t] += gd.Box[t]

// Alternative definition for VWrap. (0,0,..) lower edge of the box

#define VWrap2(v, t)                                        \
   if (v[t] >= gd.Box[t])      v[t] -= gd.Box[t];           \
   else if (v[t] < 0.0) v[t] += gd.Box[t]

//Chose one PBC::
// By now it is only working with boxes centered at (0,0,...)
// then choose VWrap1

#define VWrap       VWrap1
//#define VWrap       VWrap2

#if NDIM == 2
#define VWrapAll(v)                                         \
   {VWrap (v, 0);                                           \
   VWrap (v, 1);}
#endif

#if NDIM == 3
#define VWrapAll(v)                                         \
   {VWrap (v, 0);                                           \
   VWrap (v, 1);                                            \
   VWrap (v, 2);}
#endif

//E ! Tree search




//B Macros useful to compute chebyshev polynomials

//B BALLS4
#define CHEBYSHEVOMPCC                                      \
{real xicosmphi; int m;                                               \
      hist->Chebs[1] = 1.0;                                    \
   xicosmphi = xi * hist->Chebs[1];                    \
   hist->histXithread[1][n] += xicosmphi;                   \
   hist->Chebs[2] = cosphi;                                 \
   xicosmphi = xi * hist->Chebs[2];                    \
   hist->histXithread[2][n] += xicosmphi;                   \
   hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xi * hist->Chebs[3];                    \
   hist->histXithread[3][n] += xicosmphi;                   \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
       xicosmphi = xi * hist->Chebs[m];                \
       hist->histXithread[m][n] += xicosmphi;               \
   }}

#define CHEBYSHEVOMPCCBACKUPDELETE                                      \
  {hist->Chebs[1] = 1.0;                                    \
   xicosmphi = xj * xi * hist->Chebs[1];                    \
   hist->histXithread[1][n] += xicosmphi;                   \
   hist->Chebs[2] = cosphi;                                 \
   xicosmphi = xj * xi * hist->Chebs[2];                    \
   hist->histXithread[2][n] += xicosmphi;                   \
   hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xj * xi * hist->Chebs[3];                    \
   hist->histXithread[3][n] += xicosmphi;                   \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
       xicosmphi = xj * xi * hist->Chebs[m];                \
       hist->histXithread[m][n] += xicosmphi;               \
   }}

//E


#define CHEBYSHEVMPIOMP                                     \
  {hist_omp.Chebs[1] = 1.0;                                 \
   xicosmphi = xi * hist_omp.Chebs[1];                      \
   hist_omp.histXithread[1][n] += xicosmphi;                \
   hist_omp.Chebs[2] = cosphi;                              \
   xicosmphi = xi * hist_omp.Chebs[2];                      \
   hist_omp.histXithread[2][n] += xicosmphi;                \
   hist_omp.Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);       \
   xicosmphi = xi * hist_omp.Chebs[3];                      \
   hist_omp.histXithread[3][n] += xicosmphi;                \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist_omp.Chebs[m] = 2.0*(cosphi)*hist_omp.Chebs[m-1] - hist_omp.Chebs[m-2];  \
       xicosmphi = xi * hist_omp.Chebs[m];                  \
       hist_omp.histXithread[m][n] += xicosmphi;            \
   }}


// Same as above but to use in search=tree-omp-sincos only
#ifdef MANUALCHEBYSHEV
#define CHEBYSHEVTUOMPSINCOS                                      \
{REAL xicosmphi,xisinmphi; int m;                                 \
    REAL cosphi2, cosphi3, cosphi4;                               \
    cosphi2 = cosphi*cosphi; cosphi3 = cosphi2*cosphi;            \
    cosphi4 = cosphi3*cosphi;                                     \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                     \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                     \
    hist->ChebsT[3] = 2.0*cosphi2 - (1.0);                        \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                     \
    hist->ChebsT[4] = -3.0*cosphi + 4.0*cosphi3;                  \
    xicosmphi = xi * hist->ChebsT[4];                             \
    hist->histXithreadcos[4][n] += xicosmphi;                     \
    hist->ChebsT[5] = 1.0 - 8.0*cosphi2 + 8.0*cosphi4;            \
    xicosmphi = xi * hist->ChebsT[5];                             \
    hist->histXithreadcos[5][n] += xicosmphi;                     \
    hist->ChebsT[6] = 5.0*cosphi - 20.0*cosphi3 + 16.0*cosphi4*cosphi; \
    xicosmphi = xi * hist->ChebsT[6];                             \
    hist->histXithreadcos[6][n] += xicosmphi;                     \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                     \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                     \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                     \
    hist->ChebsU[4] = -1.0 + 4.0*cosphi2;                         \
    xisinmphi = xi * hist->ChebsU[4] * sinphi;                    \
    hist->histXithreadsin[4][n] += xisinmphi;                     \
    hist->ChebsU[5] = -4.0*cosphi + 8.0*cosphi3;                  \
    xisinmphi = xi * hist->ChebsU[5] * sinphi;                    \
    hist->histXithreadsin[5][n] += xisinmphi;                     \
    hist->ChebsU[6] = 1.0 -12.0*cosphi2 + 16.0*cosphi4;           \
    xisinmphi = xi * hist->ChebsU[6] * sinphi;                    \
    hist->histXithreadsin[6][n] += xisinmphi;                     \
    for (m=7; m<=cmd.mchebyshev+1; m++){                          \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                 \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#else
#define CHEBYSHEVTUOMPSINCOS                                      \
{real xicosmphi,xisinmphi; int m;                                 \
    hist->ChebsT[1] = 1.0;                                        \
    xicosmphi = xi * hist->ChebsT[1];                             \
    hist->histXithreadcos[1][n] += xicosmphi;                     \
    hist->ChebsT[2] = cosphi;                                     \
    xicosmphi = xi * hist->ChebsT[2];                             \
    hist->histXithreadcos[2][n] += xicosmphi;                     \
    hist->ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);              \
    xicosmphi = xi * hist->ChebsT[3];                             \
    hist->histXithreadcos[3][n] += xicosmphi;                     \
    hist->ChebsU[1] = 0.0;                                        \
    xisinmphi = xi * hist->ChebsU[1] * sinphi;                    \
    hist->histXithreadsin[1][n] += xisinmphi;                     \
    hist->ChebsU[2] = 1.0;                                        \
    xisinmphi = xi * hist->ChebsU[2] * sinphi;                    \
    hist->histXithreadsin[2][n] += xisinmphi;                     \
    hist->ChebsU[3] = 2.0*cosphi;                                 \
    xisinmphi = xi * hist->ChebsU[3] * sinphi;                    \
    hist->histXithreadsin[3][n] += xisinmphi;                     \
    for (m=4; m<=cmd.mchebyshev+1; m++){                          \
        hist->ChebsT[m] = 2.0*(cosphi)*hist->ChebsT[m-1] - hist->ChebsT[m-2]; \
        xicosmphi = xi * hist->ChebsT[m];                         \
        hist->histXithreadcos[m][n] += xicosmphi;                 \
        hist->ChebsU[m] = 2.0*(cosphi)*hist->ChebsU[m-1] - hist->ChebsU[m-2]; \
        xisinmphi = xi * hist->ChebsU[m] * sinphi;                \
        hist->histXithreadsin[m][n] += xisinmphi;                 \
    }}
#endif

#define CHEBYSHEVOMPBALLS                                   \
  {hist->Chebs[1] = 1.0;                                    \
   xicosmphi = xi * hist->Chebs[1];                         \
   hist->histXithread[1][n] += xj*xicosmphi;                \
   hist->Chebs[2] = cosphi;                                 \
   xicosmphi = xi * hist->Chebs[2];                         \
   hist->histXithread[2][n] += xj*xicosmphi;                \
   hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xi * hist->Chebs[3];                         \
   hist->histXithread[3][n] += xj*xicosmphi;                \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
       xicosmphi = xi * hist->Chebs[m];                     \
       hist->histXithread[m][n] += xj*xicosmphi;            \
   }}

#define CHEBYSHEVDIRECTOMP                                  \
  {hist.Chebs[1] = 1.0;                                     \
   xicosmphi = xi * hist.Chebs[1];                          \
   hist.histXithread[1][n] += xicosmphi;                    \
   hist.Chebs[2] = cosphi;                                  \
   xicosmphi = xi * hist.Chebs[2];                          \
   hist.histXithread[2][n] += xicosmphi;                    \
   hist.Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);           \
   xicosmphi = xi * hist.Chebs[3];                          \
   hist.histXithread[3][n] += xicosmphi;                    \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist.Chebs[m] = 2.0*(cosphi)*hist.Chebs[m-1] - hist.Chebs[m-2];  \
       xicosmphi = xi * hist.Chebs[m];                      \
       hist.histXithread[m][n] += xicosmphi;                \
   }}


#define NOCHEBYSHEVOMP                                      \
    {real Cheb;                                             \
    for (m=1; m<=cmd.mchebyshev+1; m++){                    \
      Cheb = rcos((real)(m-1) * racos(cosphi));             \
        xicosmphi = xi * Cheb;                              \
        hist->histXithread[m][n] += xicosmphi;              \
    }}

#define NOCHEBYSHEVDIRECTOMP                                      \
    {real Cheb;                                             \
    for (m=1; m<=cmd.mchebyshev+1; m++){                    \
      Cheb = rcos((real)(m-1) * racos(cosphi));             \
        xicosmphi = xi * Cheb;                              \
        hist.histXithread[m][n] += xicosmphi;              \
    }}

// Same as above. Unify
#define NOCHEBYSHEV                                         \
{real Cheb;for (m=1; m<=cmd.mchebyshev+1; m++){             \
        Cheb = rcos((real)(m-1) * racos(cosphi));           \
        xicosmphi = xi * Cheb;                              \
        gd.histXi[m][n] += xicosmphi;                       \
    }}

#define CHEBYSHEV                                           \
{Chebs[1] = 1.0;                                  \
   xicosmphi = xi * Chebs[1];                               \
   gd.histXi[1][n] += xicosmphi;                            \
   Chebs[2] = cosphi;                                       \
   xicosmphi = xi * Chebs[2];                               \
   gd.histXi[2][n] += xicosmphi;                            \
   Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);                \
   xicosmphi = xi * Chebs[3];                               \
   gd.histXi[3][n] += xicosmphi;                            \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
          Chebs[m] = 2.0*(cosphi)*Chebs[m-1] - Chebs[m-2];  \
          xicosmphi = xi * Chebs[m];                        \
          gd.histXi[m][n] += xicosmphi;                     \
   }}

#define CHEBYSHEVNEW                                        \
{hist.Chebs[1] = 1.0;                                       \
   xicosmphi = xi * hist.Chebs[1];                          \
   gd.histXi[1][n] += xicosmphi;                            \
   hist.Chebs[2] = cosphi;                                  \
   xicosmphi = xi * hist.Chebs[2];                          \
   gd.histXi[2][n] += xicosmphi;                            \
   hist.Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);           \
   xicosmphi = xi * hist.Chebs[3];                          \
   gd.histXi[3][n] += xicosmphi;                            \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist.Chebs[m] = 2.0*(cosphi)*hist.Chebs[m-1] - hist.Chebs[m-2];  \
       xicosmphi = xi * hist.Chebs[m];                      \
          gd.histXi[m][n] += xicosmphi;                     \
   }}

#define CHEBYSHEVNEWNEW                                     \
{hist->Chebs[1] = 1.0;                                       \
   xicosmphi = xi * hist->Chebs[1];                          \
   gd.histXi[1][n] += xicosmphi;                            \
   hist->Chebs[2] = cosphi;                                  \
   xicosmphi = xi * hist->Chebs[2];                          \
   gd.histXi[2][n] += xicosmphi;                            \
   hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);           \
   xicosmphi = xi * hist->Chebs[3];                          \
   gd.histXi[3][n] += xicosmphi;                            \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
       xicosmphi = xi * hist->Chebs[m];                      \
          gd.histXi[m][n] += xicosmphi;                     \
   }}

#define CHEBYSHEVOMP1                                       \
  {Chebs[1] = 1.0;                                          \
   xicosmphi = xi * Chebs[1];                               \
   histXithread[1][n] += xicosmphi;                         \
   Chebs[2] = cosphi;                                       \
   xicosmphi = xi * Chebs[2];                               \
   histXithread[2][n] += xicosmphi;                         \
   Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);                \
   xicosmphi = xi * Chebs[3];                               \
   histXithread[3][n] += xicosmphi;                         \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
          Chebs[m] = 2.0*(cosphi)*Chebs[m-1] - Chebs[m-2];  \
          xicosmphi = xi * Chebs[m];                        \
          histXithread[m][n] += xicosmphi;                  \
   }}

//E End chebyshev definitions


#define DO_BODY(p,start,finish)     for (p = start; p < finish; p++)
#define DO_DESCENDENTS(q,p)         for (q = More(p); q != Next(p); q = Next(q))
#define DO_COORD(k)                 for ((k)=0; (k)<NDIM; (k)++)

#define DO_SUNS(k)                 for (k=0; k < ((int)rpow(2,NDIM)); k++)

//B Input file formats:
#define INCOLUMNS       0
#define INCOLUMNSBIN    2
#define INNULL          1
// takahasi
#define INTAKAHASI      3
#define INCOLUMNS2DTO3D 4
//E

//B Output file formats:
#define OUTCOLUMNS       0
#define OUTCOLUMNSBIN    2
#define OUTNULL          1
//E


//B Search options:
#define SEARCHNULL              0
#define TREEOMPMETHODSINCOS     24
//E

//Rotation angle in radians. To use in a sphere (3D case)
#define ROTANGLE                0.01

#ifdef ADDONS
#include "datastruc_defs_include.h"
#endif

#endif // ! _data_struc_defs_h

