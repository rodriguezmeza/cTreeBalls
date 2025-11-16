/* ==============================================================================
 MODULE: search_neighbor_boxes_omp.c        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Based on: cute_box by David Alonso
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: kd = searchcalc_neighbor_boxes_omp(cmd, gd, btab, nbody,
                                    ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#ifdef _LONGIDS_
typedef long lint;
#else
typedef int lint;
#endif

#ifndef _N_LOGINT_
#define _N_LOGINT_ 20                               // bins per decade for
#endif                                              //  logarithmic binning

// parameters
float lbox;                                         // Box size

local int NB_R_;
local real LOG_R_MAX_;
//B this line causes results cballs vs cute-box be different in a +- 40%
local real I_R_MAX_;
// these lines instead gives very good agreement (0% rel error)
//      if use these lines make a similar correction in addons/Makefile_addons_settings...
//#ifndef I_R_MAX_
//#define I_R_MAX_ 0.005 //1/r_max
//#endif
//E

static double I_DR;
static double R2_MAX;

typedef struct {
  int np;
  double *pos;
} NeighborBox;

int _optimal_nside(double lb,double rmax,lint np);
void _free_boxes(int nside,NeighborBox *boxes);
NeighborBox *_catalog_to_boxes(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               int n_box_side);
local void make_CF(unsigned long long DD[], int nD,
                   double corr[], double ercorr[]);
local void _write_CF(char *fname,double *corr,double *ercorr,
                     unsigned long long *DD);
local void pass_run_params(struct  cmdline_data* cmd,
                           struct  global_data* gd);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);
local void _run_monopole_corr_neighbors(struct  cmdline_data* cmd,
                                        struct  global_data* gd, int cat1, int cat2);

/*
 Search omp parallel/serial routine using normal tree method:

 To be called using: search=neighbor-boxes-omp

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histNN,
    *                                    histNNSubXi2pcf, histNNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int searchcalc_neighbor_boxes_omp(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr *btable, INTEGER *nbody,
                                    INTEGER ipmin, INTEGER *ipmax,
                                    int cat1, int cat2)
{
    string routineName = "searchcalc_octree_box_omp";
    double cpustart;
    int nn;
    
    cpustart = CPUTIME;
    print_info(cmd, gd);

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_hist_sincos(cmd, gd);

/*
    INTEGER ipfalse;
    ipfalse=0;
//B kappa Avg Rmin
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
//E
*/
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);
    pass_run_params(cmd, gd);

    _run_monopole_corr_neighbors(cmd, gd, cat1, cat2);

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nGoing out: CPU time = %lf %s\n",
                        CPUTIME-cpustart, PRNUNITOFTIMEUSED);

  return 0;
}

// pass cBalls parameters to CB
void pass_run_params(struct  cmdline_data* cmd,
                     struct  global_data* gd)
{
    string routineName = "pass_run_params";
  
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: cBalls run parameters -> CB\n", routineName);

    lbox=cmd->lengthBox;

//B this line causes results cballs vs cute-box be different in a +- 40%
    I_R_MAX_ = 1.0/cmd->rangeN;
//E
    NB_R_ = cmd->sizeHistN;
    LOG_R_MAX_ = log10(cmd->rangeN);
    I_DR=I_R_MAX_*NB_R_;
    R2_MAX=1./(I_R_MAX_*I_R_MAX_);
}

// correlator for monopole in the periodic-box case
//  using neighbor boxes
void _corr_mono_box_neighbors(struct  cmdline_data* cmd,
                              struct  global_data* gd,
                              int nside,NeighborBox *boxes,
                 lint np,
                 unsigned long long hh[])
{
    string routineName = "_corr_mono_box_neighbors";
    double agrid=lbox/nside;
    double r_max=1/I_R_MAX_;
    int index_max=(int)(r_max/agrid)+1;
    int i;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: boxes will be correlated up to %d box sizes \n",
                           routineName, index_max);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");


  for(i=0;i<NB_R_;i++)
    hh[i]=0;                                        // clear shared histogram

#pragma omp parallel default(none)                  \
  shared(index_max,nside,boxes,hh,lbox,np,    \
    agrid,NB_R_,R2_MAX,I_DR,LOG_R_MAX_,bodytable,cmd,gd)
  {
      lint ii;
      double a2grid=agrid*agrid;
      unsigned long long hthread[NB_R_];            // histogram for each thread

      bodyptr p;

      for(ii=0;ii<NB_R_;ii++)
          hthread[ii]=0;                            // clear private histogram

#pragma omp for nowait schedule(dynamic)
      for(ii=0;ii<np;ii++) {
          int ix0,iy0,iz0;
          double x0,y0,z0;
          int idz;
          p = bodytable[0]+ii;

          x0=Pos(p)[0];
          y0=Pos(p)[1];
          z0=Pos(p)[2];

      ix0=(int)(x0/lbox*nside);
      iy0=(int)(y0/lbox*nside);
      iz0=(int)(z0/lbox*nside);

      for(idz=-index_max;idz<=index_max;idz++) {
    int idy,idz_dist2;
    int iwrapz=0;
    int iz1=iz0+idz;
    if(iz1<0) {
      iz1+=nside;
      iwrapz=1;
    } else if(iz1>=nside) {
      iz1-=nside;
      iwrapz=1;
    }
    idz_dist2=MAX(0,abs(idz)-1);
    idz_dist2=idz_dist2*idz_dist2;
    for(idy=-index_max;idy<=index_max;idy++) {
      int idx,idy_dist2;
      int iwrapy=0;
      int iy1=iy0+idy;
      if(iy1<0) {
        iy1+=nside;
        iwrapy=1;
      } else if(iy1>=nside) {
        iy1-=nside;
        iwrapy=1;
      }
      idy_dist2=MAX(0,abs(idy)-1);
      idy_dist2=idy_dist2*idy_dist2;
      for(idx=-index_max;idx<=index_max;idx++) {
        int ibox,idx_dist;
        int iwrapx=0;
        int ix1=ix0+idx;
        double d2max;
        int jj;
        if(ix1<0) {
          ix1+=nside;
          iwrapx=1;
        } else if(ix1>=nside) {
          ix1-=nside;
          iwrapx=1;
        }
        ibox=ix1+nside*(iy1+nside*iz1);
        idx_dist=MAX(0,abs(idx)-1);
        d2max=a2grid*(idx_dist*idx_dist+idy_dist2+idz_dist2);
        if(d2max>R2_MAX) continue;
        for(jj=0;jj<boxes[ibox].np;jj++) {
          double xr[3];
          double r2;
          int ir;
          xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj]);
          xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);
          xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);
          if(iwrapx) xr[0]=lbox-xr[0];
          if(iwrapy) xr[1]=lbox-xr[1];
          if(iwrapz) xr[2]=lbox-xr[2];
          r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
          if(r2>R2_MAX) continue;
#ifdef _LOGBIN_
          if(r2>0) {
        ir=(int)(_N_LOGINT_*(0.5*log10(r2)-LOG_R_MAX_)+NB_R_);
        if((ir<NB_R_)&&(ir>=0))
          (hthread[ir])++;
          }
#else // ! _LOGBIN_
            ir=(int)(sqrt(r2)*I_DR);
          if(ir<NB_R_) //Check bound
        (hthread[ir])++;
#endif // ! _LOGBIN_
        }
      }
    }
      }
    } // ! end pragma omp for

#pragma omp critical
    {
      for(ii=0;ii<NB_R_;ii++) //Check bound
    hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}

// main routine for monopole using neighbor boxes
local void _run_monopole_corr_neighbors(struct  cmdline_data* cmd,
                                  struct  global_data* gd, int cat1, int cat2)
{
    string routineName = "_run_monopole_corr_neighbors";
  lint n_dat;
  int nside;
  NeighborBox *boxes;
  unsigned long long DD[NB_R_];
  double corr[NB_R_],ercorr[NB_R_];

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: Correlation function parameters:\n", routineName);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: Range: %.3lf < r < %.3lf (Mpc/h)\n",
                           routineName, 0., 1./(I_R_MAX_));
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: #bins: %d\n",routineName, NB_R_);
#ifdef _LOGBIN_
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: using logarithmic binning with %d bins per decade \n",
                           routineName, _N_LOGINT_);
#else
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: resolution: Dr = %.3lf (Mpc/h)\n",
                           routineName, 1./(I_R_MAX_*NB_R_));
#endif // !_LOGBIN_

    //B cBalls structure...
    bodyptr p;
    int k;
    int ii;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: Box: %g %g %g\n",
                           routineName, gd->Box[0], gd->Box[1], gd->Box[2]);
    for(ii=0;ii<gd->nbodyTable[cat1];ii++) {
        p = bodytable[cat1]+ii;
        if (scanopt(cmd->options, "cute-box-fmt"))
            DO_COORD(k)
                Pos(p)[k] += 0.5*gd->Box[k];
    }
    //E

    nside=_optimal_nside(lbox,1./I_R_MAX_,gd->nbodyTable[cat1]);
    boxes=_catalog_to_boxes(cmd, gd, nside);

#ifdef _DEBUG_
    write_cat(cat_dat,"output_cute_box/debug_DatCat.dat");
#endif

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: correlating...\n", routineName);
    _corr_mono_box_neighbors(cmd, gd, nside, boxes,
                             gd->nbodyTable[cat1], DD);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: writing output...\n", routineName);
    make_CF(DD,gd->nbodyTable[cat1],corr,ercorr);
    _write_CF(gd->fpfnamehistCFFileName,corr,ercorr,DD);
    gd->flagPrint = FALSE;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: cleaning up...\n", routineName);
    _free_boxes(nside,boxes);
}


#ifndef FRACTION_AR
#define FRACTION_AR 8.0
#endif

// estimates a good candidate for the size of
//  a set of nearest-neighbor searching boxes
int _optimal_nside(double lb,double rmax,lint np)
{
  int nside1=(int)(FRACTION_AR*lb/rmax);
  int nside2=(int)(pow(0.5*np,0.3333333));

  return MIN(nside1,nside2);
}

#ifdef _DEBUG_
// writes catalog
//  only used for debugging
void write_cat(CatalogBalls cat,char *fn)
{
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
#endif // ! _DEBUG_


// frees all memory associated with a box
//  set of size nside
void _free_boxes(int nside,NeighborBox *boxes)
{
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
                               int n_box_side)
{
    string routineName = "_catalog_to_boxes";
    lint ii;
    int nside;
    NeighborBox *boxes;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: Building neighbor boxes \n",
                           routineName, gd->bytes_tot*INMB);
    nside=n_box_side;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "%s: there will be %d boxes per side with a size of %lf\n",
                           routineName, nside,lbox/nside);

    boxes=(NeighborBox *)malloc(nside*nside*nside*sizeof(NeighborBox));
    if(boxes==NULL) error("%s: out of memory!!\n", routineName);;
    for(ii=0;ii<nside*nside*nside;ii++)
    boxes[ii].np=0;

    bodyptr p;
    for(ii=0;ii<gd->nbodyTable[0];ii++) {
        int ix,iy,iz;
        p = bodytable[0]+ii;

        ix=(int)(Pos(p)[0]/lbox*nside);
        iy=(int)(Pos(p)[1]/lbox*nside);
        iz=(int)(Pos(p)[2]/lbox*nside);

        (boxes[ix+nside*(iy+nside*iz)].np)++;
    }

    for(ii=0;ii<nside*nside*nside;ii++) {
        int npar=boxes[ii].np;
        if(npar>0) {
            boxes[ii].pos=(double *)malloc(3*npar*sizeof(double));
            if(boxes[ii].pos==NULL) error("%s: out of memory!!\n", routineName);
            boxes[ii].np=0;
        }
    }

    for(ii=0;ii<gd->nbodyTable[0];ii++) {
        int ix,iy,iz,index,offset;
        p = bodytable[0]+ii;

        ix=(int)(Pos(p)[0]/lbox*nside);
        iy=(int)(Pos(p)[1]/lbox*nside);
        iz=(int)(Pos(p)[2]/lbox*nside);

        index=ix+nside*(iy+nside*iz);
        offset=3*boxes[index].np;

        (boxes[index].pos)[offset]=Pos(p)[0];
        (boxes[index].pos)[offset+1]=Pos(p)[1];
        (boxes[index].pos)[offset+2]=Pos(p)[2];

        (boxes[index].np)++;
  }

  return boxes;
}

// creates correlation function and poisson errors
//  from pair counts DD, DR and RR
// The correlation function is estimated as:
//      xi=(V/v(r))*(DD(r)/N^2)
//      where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
//      total # particles.
//  Note that, since in this case we have simple periodic boundary conditions,
//      no random catalogs are needed.
void make_CF(unsigned long long DD[],int nD,
         double corr[],double ercorr[])
{
  double *edd;
  double rho_av=nD/(lbox*lbox*lbox);
  int ii;

  edd=(double *)malloc(sizeof(double)*NB_R_);
  if(edd==NULL)
    error("CUTE: Out of memory!!\n");

#ifndef _LOGBIN_
  DD[0]-=nD;                                        // substract diagonal
#endif
  for(ii=0;ii<NB_R_;ii++)
    edd[ii]=1./sqrt((double)DD[ii]);

  for(ii=0;ii<NB_R_;ii++) {
    if(DD[ii]==0) {
      corr[ii]=0;
      ercorr[ii]=0;
    }
    else {
      double r0,r1,vr,rho_r;
#ifdef _LOGBIN_
      r0=pow(10.,(((double)ii-NB_R_)/_N_LOGINT_)+LOG_R_MAX_);
      r1=pow(10.,(((double)ii+1-NB_R_)/_N_LOGINT_)+LOG_R_MAX_);
#else
      r0=ii/(I_R_MAX_*NB_R_);
      r1=(ii+1)/(I_R_MAX_*NB_R_);
#endif // !_LOGBIN_
      vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      rho_r=DD[ii]/(nD*vr);
      corr[ii]=rho_r/rho_av-1;
      ercorr[ii]=(1+corr[ii])*edd[ii];
    }
  }

  free(edd);
}

// writes correlation function
local void _write_CF(char *fname,double *corr,double *ercorr,
          unsigned long long *DD)
{
  FILE *fo;
  int ii;

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"CUTE: Error opening output file %s",fname);
    fprintf(stderr,", using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error("CUTE: Couldn't open file %s \n",oname);
  }

  for(ii=0;ii<NB_R_;ii++) {
    double rr;
#ifdef _LOGBIN_
    rr=pow(10,((ii+0.5)-NB_R_)/_N_LOGINT_+LOG_R_MAX_);
#else
    rr=(ii+0.5)/(NB_R_*I_R_MAX_);
#endif // ! _LOGBIN_
    fprintf(fo,"%lE %lE %lE %llu \n",
        rr,corr[ii],ercorr[ii],DD[ii]);
  }
  fclose(fo);
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose, "Search: Running ... (neighbor-boxes-omp) \n");

    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    verb_print(cmd->verbose, "computing only 2pcf... \n");

    return SUCCESS;
}

