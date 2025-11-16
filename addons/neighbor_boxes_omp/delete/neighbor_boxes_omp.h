/* ==============================================================================
 MODULE: neighbor_boxes.h        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza.
 Based on: zeno lib
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: include "kdtree_box.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _neighbor_boxes_omp_h
#define _neighbor_boxes_omp_h

/*
//B from define.h
#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef int lint;
#endif //_LONGIDS

extern char _fnameData[128];
extern char _fnameOut[128];
extern int _input_format;
extern int _use_tree;
extern int _max_tree_order;
extern int _max_tree_nparts;
extern int _use_pm;
extern lint _n_objects;
extern float _l_box;
extern float _l_box_half;
extern int _n_grid;
//                MACROS            
// Other possible macros
//_DEBUG, _VERBOSE, _TRUE_ACOS, _LOGBIN

//#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
//#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
//#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)

#ifndef N_LOGINT
#define N_LOGINT 20 //# bins per decade for logarithmic binning
#endif

//Monopole
//#ifndef _I_R_MAX
//real _I_R_MAX;
//#define _I_R_MAX 0.005 //1/r_max
//#endif
//#ifndef _LOG_R_MAX
//real _LOG_R_MAX;
//#define _LOG_R_MAX 2.30103 //log10(r_max)
//#endif
//#ifndef _NB_R
//int _NB_R;
//#define _NB_R 256 //# bins in r for monopole correlation
//#endif
//E

//B from define.h
typedef struct {
  lint np;          //#objects in the catalog
  double *pos;
} CatalogBalls;         //Catalog (double precision)

typedef struct branch {
  float x_lo[3];
  float x_hi[3];
  char leaf;
  lint np;
  void *sons; //sons
} branch; //Tree node/branch

typedef struct {
  int np;
  double *pos;
} NeighborBox; //Neighbor box
//E

//B from common.h
//void timer(int i);

lint _linecount(FILE *f);

int _optimal_nside(double lb,double rmax,lint np);

void _free_catalog(CatalogBalls cat);

//void _error_mem_out(void);

//void _error_open_file(char *fname);

//void _error_read_line(char *fname,lint nlin);

#ifdef _DEBUG
void write_cat(CatalogBalls cat,char *fn);

void write_grid(double *grid,char *fn);

void write_tree(branch *tree,char *fn);
#endif //_DEBUG
//E

//B from neighbours.h
void _free_boxes(int nside,NeighborBox *boxes);
NeighborBox *_catalog_to_boxes(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               int n_box_side,CatalogBalls cat);
//E

*/

#endif  /* ! _neighbor_boxes_omp_h */
