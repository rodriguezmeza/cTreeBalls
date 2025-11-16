/* ==============================================================================
 MODULE: neighbor_boxes_omp.c        [cTreeBalls]
 Written by: Mario A. Rodriguez-Meza.
 Based on: cute_box code
 Starting date:    april 2023
 Purpose: 2-point correlation functions on boxes computation
 Language: C
 Use: kd = init_kdtree(cmd, gd, btab, nbody, nbody);
      build_kddtree(cmd, gd, nbucket);
      finish_kdtree(kd);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "neighbor_boxes_omp.h"

/*
//B from define.c

////////// Input parameters //////////
///
//File names
//char _fnameData[128]="default";   //Data catalog filename
char _fnameOut[128]="default";    //Output filename

//File format
int _input_format=-1;

//Parameters
float _l_box=-1; //Box size
float _l_box_half=-1;

//Tree stuff
int _use_tree=-1;
int _max_tree_order=-1;
int _max_tree_nparts=-1;

//PM stuff
int _use_pm=-1;        //Should I use PM?
int _n_grid=-1;        //# cells per side (CUTE_box)
///
//////////////////////////////////////

////////// Internal variables //////////
///
lint _n_objects=-1;          //# objects to read from the files
///
////////////////////////////////////////

//E
*/

/*
//B from common.c

#ifndef FRACTION_AR
#define FRACTION_AR 8.0
#endif //FRACTION_AR

lint _linecount(FILE *f)
{
  // Returns #lines in f
  lint i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

int _optimal_nside(double lb,double rmax,lint np)
{
  // Estimates a good candidate for the size of
  // a set of nearest-neighbor searching boxes
//    verb_print_debug(1, "\nAqui voy (2): %g %g %ld %g\n", lb, rmax, np, FRACTION_AR);
  int nside1=(int)(FRACTION_AR*lb/rmax);
//    verb_print_debug(1, "\nAqui voy (3): %d\n", nside1);
  int nside2=(int)(pow(0.5*np,0.3333333));
//    verb_print_debug(1, "\nAqui voy (3): %d\n", nside2);

  return MIN(nside1,nside2);
}

void _free_catalog(CatalogBalls cat)
{
  // Frees position arrays in catalog
  if(cat.np>0)
    free(cat.pos);
}

#ifdef _DEBUG_
void write_cat(CatalogBalls cat,char *fn)
{
  // Writes catalog into file fn.
  // Only used for debugging
  FILE *fr;
  lint ii;
  fr=fopen(fn,"w");
  if(fr==NULL) _error_open_file(fn);
  for(ii=0;ii<cat.np;ii++) {
    fprintf(fr,"%lf %lf %lf\n",cat.pos[3*ii],
        cat.pos[3*ii+1],cat.pos[3*ii+2]);
  }
  fclose(fr);
}

void write_grid(double *grid,char *fn)
{
  //////
  // Writes grid into file fn.
  // Only used for debugging
  FILE *fr;
  lint ii;
  double agrid=_l_box/_n_grid;
  fr=fopen(fn,"w");
  if(fr==NULL) _error_open_file(fn);
  for(ii=0;ii<_n_grid;ii++) {
    lint jj;
    double x=(ii+0.5)*agrid;
    for(jj=0;jj<_n_grid;jj++) {
      lint kk;
      double y=(jj+0.5)*agrid;
      for(kk=0;kk<_n_grid;kk++) {
    double z=(kk+0.5)*agrid;
    lint index=kk+_n_grid*(jj+_n_grid*ii);
    fprintf(fr,"%lf %lf %lf %lf\n",x,y,z,grid[index]);
      }
    }
  }
  fclose(fr);
}

static void print_branch(branch *br,FILE *fil)
{
  //////
  // Recursive branch writer
  int ii;
  if(br->leaf) {
    for(ii=0;ii<br->np;ii++) {
      fprintf(fil,"%lf %lf %lf \n",
          ((double *)(br->sons))[3*ii],
          ((double *)(br->sons))[3*ii+1],
          ((double *)(br->sons))[3*ii+2]);
    }
  }
  else {
    for(ii=0;ii<8;ii++)
      print_branch(((branch **)(br->sons))[ii],fil);
  }
}

void write_tree(branch *tree,char *fn)
{
  //////
  // Writes all particles in tree into file fn.
  // Only used for debugging
  FILE *fr;
  fr=fopen(fn,"w");
  print_branch(tree,fr);
  fclose(fr);
}
#endif //_DEBUG

//E


void _free_boxes(int nside,NeighborBox *boxes)
{
  //////
  // Frees all memory associated with a box
  // set of size nside
  int ii;

  for(ii=0;ii<nside*nside*nside;ii++) {
    if(boxes[ii].np>0)
      free(boxes[ii].pos);
  }
  
  free(boxes);
}

// creates boxes for nearest-neighbor searching
NeighborBox *_catalog_to_boxes(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               int n_box_side,CatalogBalls cat)
{
    string routineName = "_catalog_to_boxes";
  lint ii;
  int nside;
  NeighborBox *boxes;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: Building neighbor boxes \n",
                           routineName, gd->bytes_tot*INMB);
  nside=n_box_side;
  printf("  There will be %d boxes per side with a size of %lf \n",
     nside,_l_box/nside);

   boxes=(NeighborBox *)malloc(nside*nside*nside*sizeof(NeighborBox));
    if(boxes==NULL) error("CUTE: Out of memory!!\n");;
  for(ii=0;ii<nside*nside*nside;ii++)
    boxes[ii].np=0;

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz;

    ix=(int)(cat.pos[3*ii]/_l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/_l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/_l_box*nside);

    (boxes[ix+nside*(iy+nside*iz)].np)++;
  }

  for(ii=0;ii<nside*nside*nside;ii++) {
    int npar=boxes[ii].np;
    if(npar>0) {
      boxes[ii].pos=(double *)malloc(3*npar*sizeof(double));
        if(boxes[ii].pos==NULL) error("CUTE: Out of memory!!\n");
      boxes[ii].np=0;
    }
  }

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz,index,offset;

    ix=(int)(cat.pos[3*ii]/_l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/_l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/_l_box*nside);
    index=ix+nside*(iy+nside*iz);
    offset=3*boxes[index].np;
    (boxes[index].pos)[offset]=cat.pos[3*ii];
    (boxes[index].pos)[offset+1]=cat.pos[3*ii+1];
    (boxes[index].pos)[offset+2]=cat.pos[3*ii+2];
    (boxes[index].np)++;
  }

  printf("\n");
  
  return boxes;
}

*/


