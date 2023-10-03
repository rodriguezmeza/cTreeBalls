/*==============================================================================
 HEADER: data_struc_defs.h		[tpcf]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definition of N-Body data structure
 Language: C
 Use: '#include "...."
 ==============================================================================*/

 
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


// The weight can be -> other scalar field of interest
typedef struct _node {
    short type;
    bool update;
    real weight;
    real kappa;
    vector pos;
    struct _node *next;
    int *bhistNsub;
    real *bhistXi2pcfsub;
    real **histXi;
    bool histON;
    
    byte flags;                // status in b-b calculation
    byte newlevel;            // level set by latest force calc.

//B To see the bodies belonging to a cell:
    bool selected;
//E
//BODY3
    INTEGER nbb;     // If comes from smoothing gives number of smoothingbodies
                    // Body will be tagged NBODY3

// BALLS
    int lev;
    int idxscanlev;
    INTEGER Id;
#ifdef DEBUG
    bool hit;
#endif
//

} node, *nodeptr;

#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Weight(x)   (((nodeptr) (x))->weight)
#define Kappa(x)    (((nodeptr) (x))->kappa)
#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)
#define bNsub(x)    (((nodeptr) (x))->bhistNsub)
#define bXi2pcfsub(x)    (((nodeptr) (x))->bhistXi2pcfsub)
#define hXi(x)    (((nodeptr) (x))->histXi)
#define hON(x)    (((nodeptr) (x))->histON)
// BODY3
#define Nbb(x)    (((nodeptr) (x))->nbb)
// BALLS
#define IdNode(x)    (((nodeptr) (x))->Id)
#define Level(x)    (((nodeptr) (x))->lev)
#define IDXSCAN(x)    (((nodeptr) (x))->idxscanlev)
#ifdef DEBUG
#define HIT(x)    (((nodeptr) (x))->hit)
#endif
//

//B To see the bodies belonging to a cell:
#define Selected(x)    (((nodeptr) (x))->selected)
//E

#define BODY 00
#define BODY3 03    // Smooth
#define CELL 02
// BALLS
#define NODEBODY 04
#define NODECELL 05
//

// The meaning of the structure and its components can be changed
typedef struct {
    node bodynode;
    INTEGER Id;

    real smooth;                // smoothing length
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
    INTEGER nb;
    real radius;
//B To debug cells:
    real size;
//
    nodeptr more;
//    union {
    nodeptr subp[NSUB];
    INTEGER Id;
//B NOLSST:
    bool inside;
//    nodeptr up;
//E
} cell, *cellptr;
 
#define Radius(x) (((cellptr) (x))->radius)
//B To debug cells:
#define Size(x) (((cellptr) (x))->size)
//E
#define Nb(x) (((cellptr) (x))->nb)

#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->subp)
#define IdCell(x)   (((cellptr) (x))->Id)
//B NOLSST:
#define Inside(x)   (((cellptr) (x))->inside)
//#define Up(x)   (((cellptr) (x))->up)
//E

#if !defined(global)            // global def question must be here
#  define global extern
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

#define CHEBYSHEVOMP                                        \
  {hist->Chebs[1] = 1.0;                                    \
   xicosmphi = xi * hist->Chebs[1];                         \
   hist->histXithread[1][n] += xicosmphi;                   \
   hist->Chebs[2] = cosphi;                                 \
   xicosmphi = xi * hist->Chebs[2];                         \
   hist->histXithread[2][n] += xicosmphi;                   \
   hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
   xicosmphi = xi * hist->Chebs[3];                         \
   hist->histXithread[3][n] += xicosmphi;                   \
   for (m=4; m<=cmd.mchebyshev+1; m++){                     \
       hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
       xicosmphi = xi * hist->Chebs[m];                     \
       hist->histXithread[m][n] += xicosmphi;               \
   }}

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

// This approximation gives wrong answer. Check!
//#define ONETHIRD 0.33333333333333333
//#define FOUR45       0.0888888888888888888888
//#define XCOS        (1.0-cosphi)
//Cheb = rcos((real)(m-1) * rsqrt(2.0*XCOS+ONETHIRD*rsqr(XCOS)+FOUR45*rpow(XCOS,3.0)) );

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
#define TREEOMPMETHOD           9
#define DIRECT3PCFOMP           25
#define TREE3PCFBFOMPMETHOD     45
#define SEARCHNULL              0
#define TREEOMPMETHODSINCOS     24
#define BALLSOMPMETHOD         46
//E

//Rotation angle in radians. To use in a sphere (3D case)
#define ROTANGLE                0.01



#endif // ! _data_struc_defs_h

