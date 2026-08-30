/* ==============================================================================
 MODULE: search_neighbor_boxes_omp.c        [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Modified version of some routines from cute_box by David Alonso
 See also Rapaport, The art of molecular dynamics simulations, 2nd edition (2004)
 Starting date:    april 2023
 Purpose: 2/3-point correlation functions computation
 Language: C
 Use: kd = searchcalc_neighbor_boxes_omp(cmd, gd, btab, nbody,
                                    ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include <pthread.h>

//B parameters
float lbox;                                         // Box size
local int NB_R_;                                    // number of radial bins
local real LOG_R_MAX_;                              // log of radial range
//B this line causes results cballs
//  vs cute-box be different in a +- 40%
local real I_R_MAX_;                                // inverse of radial range
// these lines instead gives very good
//      agreement (0% rel error)
//      if use these lines make a similar correction
//      in addons/Makefile_addons_settings...
//#ifndef I_R_MAX_
//#define I_R_MAX_ 0.005 //1/r_max
//#endif
//E
static double I_DR;
static double R2_MAX;
//E

//B definition of structures
typedef struct {
    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;
} gdl_hist_omp, *gdlptr_hist_omp;

typedef struct {
  int np;
  double *pos;
} NeighborBox;
//E

//B definition of functions
int _optimal_nside(double lb,double rmax,INTEGER np);
void _free_boxes(int nside,NeighborBox *boxes);

local int _catalog_to_boxes(struct cmdline_data *cmd,
                            struct global_data *gd,
                            int n_box_side,
                            NeighborBox **boxes_out);

local int make_CF(struct cmdline_data *cmd,
                  struct global_data *gd,
                  unsigned long long DD[], int nD,
                  double corr[], double ercorr[]);

local int _write_CF(struct  cmdline_data* cmd,
                     struct  global_data* gd,
                     char *fname,double *corr,double *ercorr,
                     unsigned long long *DD);
local int PrintHistNN(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd);
local void pass_run_params(struct  cmdline_data* cmd,
                           struct  global_data* gd);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);
local int _run_monopole_corr_neighbors(struct cmdline_data* cmd,
                                        struct global_data* gd,
                                        int cat1, int cat2);
local int _corr_mono_box_neighbors(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   int nside,NeighborBox *boxes,
                                   INTEGER np, unsigned long long hh[]);
#ifdef _DEBUG_
local int write_cat(struct cmdline_data* cmd,
                     struct global_data* gd, char *fn);
#endif

local int search_init_hist_omp(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               gdlptr_hist_omp hist);

//E

/*
 Search omp parallel/serial routine using neighbor boxes method:

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
    * Histograms outputs (in global gd): histNN,
    *                                    histXi2pcf, histCF,
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

    search_init_gd_hist_sincos(cmd, gd);

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: Total allocated %g MByte storage so far.\n",
                           routineName, gd->bytes_tot*INMB);
    pass_run_params(cmd, gd);

    if (_run_monopole_corr_neighbors(cmd, gd, cat1, cat2) == FAILURE)
    return FAILURE;

    gd->cpusearch = CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nGoing out: CPU time = %lf %s\n",
                        CPUTIME-cpustart, PRNUNITOFTIMEUSED);

  return SUCCESS;
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

// main routine for monopole using neighbor boxes
local int _run_monopole_corr_neighbors(struct  cmdline_data* cmd,
                                  struct  global_data* gd, int cat1, int cat2)
{
    string routineName = "_run_monopole_corr_neighbors";
    int status = FAILURE;
    INTEGER n_dat;
    int nside;
    NeighborBox *boxes = NULL;
    unsigned long long DD[NB_R_];
    double corr[NB_R_],ercorr[NB_R_];

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: Correlation function parameters:\n", routineName);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: Range: %.3lf < r < %.3lf (Mpc/h)\n",
                           routineName, 0., 1./(I_R_MAX_));
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: #bins: %d\n",routineName, NB_R_);
#ifdef LOGBINCBON
#ifdef _LOGBIN_
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: using logarithmic binning with %d bins per decade \n",
                               routineName, cmd->logHistBinsPD);
#else
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: resolution: Dr = %.3lf (Mpc/h)\n",
                               routineName, 1./(I_R_MAX_*NB_R_));
#endif // !_LOGBIN_
#else // ! LOGBINCBON
    if (cmd->useLogHist) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: using logarithmic binning with %d bins per decade \n",
                               routineName, cmd->logHistBinsPD);
    } else {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%s: resolution: Dr = %.3lf (Mpc/h)\n",
                               routineName, 1./(I_R_MAX_*NB_R_));
    }
#endif // ! LOGBINCBON

    //B cBalls structure...
    bodyptr p;
    int k;
    int ii;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: Box: %g %g %g\n",
                           routineName, gd->Box[0], gd->Box[1], gd->Box[2]);
    if (cballs_opt_cute_box_fmt(cmd))
    for(ii=0;ii<gd->nbodyTable[cat1];ii++) {
        p = bodytable[cat1]+ii;
        DO_COORD(k)
            Pos(p)[k] += 0.5*gd->Box[k];
    }
    //E

    nside=_optimal_nside(lbox,1./I_R_MAX_,gd->nbodyTable[cat1]);
    if (_catalog_to_boxes(cmd, gd, nside, &boxes) == FAILURE)
        return FAILURE;

#ifdef SMOOTHPIVOT
    if (prepare_smooth_pivots(cmd, gd, bodytable, gd->nbodyTable,
                              1, gd->nbodyTable, cat1, cat2) == FAILURE)
        goto cleanup;
#endif

#ifdef _DEBUG_
    char OutputFileName[32];
    if (format_checked(OutputFileName, sizeof(OutputFileName),
        "OutputFileName", "%s/%s%s",
                       cmd->rootDir,"debug_DatCat",EXTFILES) != 0)
        return FAILURE;
    write_cat(cmd, gd, OutputFileName);
#endif

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: correlating...\n", routineName);
    if (_corr_mono_box_neighbors(cmd, gd, nside, boxes,
                                 gd->nbodyTable[cat1], DD) == FAILURE)
        goto cleanup;
    
    // ===============================================
    //B Saving histograms section: case 2pCORRELATION:
    // ===============================================
    if (make_CF(cmd, gd, DD, gd->nbodyTable[cat1], corr, ercorr) == FAILURE) {
        _free_boxes(nside, boxes);
        return FAILURE;
    }
    
    if (cballs_opt_compute_histn(cmd)) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\n\t%s: printing neighbor-boxes-omp method...\n\n",
                               routineName);
        PRINT_OR_FAIL(PrintHistNN(cmd, gd));
        PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));
    } else {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "\n%s: writing output...\n", routineName);
        if (_write_CF(cmd, gd, gd->fpfnamehistCFFileName,corr,ercorr,DD) == FAILURE)
            goto cleanup;
    }
    gd->flagPrint = FALSE;
    // ===============================================
    //E Saving histograms section: case GGGCORRELATION
    // ===============================================

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: cleaning up...\n", routineName);
    
    status = SUCCESS;

    cleanup:
        if (boxes != NULL)
            _free_boxes(nside,boxes);

    return status;
}

// correlator for monopole in the periodic-box case
//  using neighbor boxes
local int _corr_mono_box_neighbors(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   int nside,NeighborBox *boxes,
                                   INTEGER np, unsigned long long hh[])
{
    string routineName = "_corr_mono_box_neighbors";
    double agrid=lbox/nside;
    double r_max=1/I_R_MAX_;
    int index_max=(int)(r_max/agrid)+1;
    int i;

    int catt1=0;
#ifdef OPENMPCODE
    ThreadCount(cmd, gd, gd->nbodyTable[catt1], catt1);
    pthread_t main_thread_id = pthread_self();
#endif
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: boxes will be correlated up to %d box sizes \n",
                           routineName, index_max);
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");


    for (i = 0; i < NB_R_; i++) {
        hh[i] = 0;                                  // clear shared histogram
    }
    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

#ifdef SMOOTHPIVOT
    bodyptr p;
    INTEGER ipfalse;
    ipfalse=0;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n%s: rsmooth = %e\n",
                           routineName, gd->rsmooth[0]);
#else
#ifdef BALLS4SCANLEV
    int ifile=0;
    bodyptr p;
    DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
    if (MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile)
        == FAILURE)
        return FAILURE;
    int k;
    int cat1=0;
    if (cballs_opt_cute_box_fmt(cmd))
    for(int ii=0;ii<gd->nbodyTable[cat1];ii++) {
        p = bodytable[cat1]+ii;
        DO_COORD(k)
            Pos(p)[k] += 0.5*gd->Box[k];
    }
#endif
#endif // ! SMOOTHPIVOT

#ifndef BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(none)                  \
shared(index_max,nside,boxes,hh,lbox,np,            \
agrid,NB_R_,R2_MAX,I_DR,LOG_R_MAX_,bodytable,       \
cmd,gd,main_thread_id,                              \
ipfalse)
#else
#pragma omp parallel default(none)                  \
  shared(index_max,nside,boxes,hh,lbox,np,main_thread_id, \
    agrid,NB_R_,R2_MAX,I_DR,LOG_R_MAX_,bodytable,cmd,gd)
#endif
#else // ! BALLS4SCANLEV
#ifdef SMOOTHPIVOT
#pragma omp parallel default(none)                  \
shared(index_max,nside,boxes,hh,lbox,np,            \
agrid,NB_R_,R2_MAX,I_DR,LOG_R_MAX_,bodytable,       \
cmd,gd,nodetablescanlevB4,main_thread_id,           \
ipfalse)
#else
#pragma omp parallel default(none)                  \
  shared(index_max,nside,boxes,hh,lbox,np,          \
    agrid,NB_R_,R2_MAX,I_DR,LOG_R_MAX_,bodytable,   \
    cmd,gd,nodetablescanlevB4,main_thread_id)
#endif
#endif // ! BALLS4SCANLEV
  {
      INTEGER ii;
      INTEGER ip;
      double a2grid=agrid*agrid;
      unsigned long long hthread[NB_R_];            // histogram for each thread

      pthread_t current_thread_id = pthread_self();

#ifndef BALLS4SCANLEV
      bodyptr p;
#else
      nodeptr p;
#endif

#ifdef SMOOTHPIVOT
      INTEGER ipfalsethreads;
      ipfalsethreads = 0;
#endif

      gdl_hist_omp hist;
      search_init_hist_omp(cmd, gd, &hist);

      for(ii=0;ii<NB_R_;ii++)
          hthread[ii]=0;                            // clear private histogram

      int cat1=0;

#pragma omp for nowait schedule(static,1)
#ifndef ORIGINALCB
#ifndef BALLS4SCANLEV
      for (p = bodytable[cat1]; p < bodytable[cat1] + np; p++) {
          ii = p - bodytable[cat1];
#else
      for (INTEGER ii=0; ii< gd->nnodescanlevTableB4[cat1]; ii++) {
          p = nodetablescanlevB4[cat1][ii];
#endif
#else
      for(ii=0;ii<np;ii++) {
          p = bodytable[0]+ii;
#endif
          int ix0,iy0,iz0;
          double x0,y0,z0;
          int idz;

#ifdef SMOOTHPIVOT
          if (Update(p) == FALSE) {
              ipfalsethreads++;
              continue;
          }
#endif

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
                          hist.nbbcalcthread += 1;
#ifdef LOGBINCBON
#ifdef _LOGBIN_
                              if(r2>0) {
                                  ir=(int)(cmd->logHistBinsPD
                                           *(0.5*log10(r2)-LOG_R_MAX_)+NB_R_);
                                  if((ir<NB_R_)&&(ir>=0))
                                      (hthread[ir])++;
                              }
#else // ! _LOGBIN_
                              ir=(int)(sqrt(r2)*I_DR);
                              if(ir<NB_R_)          // check radial range
                                  (hthread[ir])++;
#endif // ! _LOGBIN_
#else // ! LOGBINCBON
                          if (cmd->useLogHist) {
                              if(r2>0) {
                                  ir=(int)(cmd->logHistBinsPD
                                           *(0.5*log10(r2)-LOG_R_MAX_)+NB_R_);
                                    if((ir<NB_R_)&&(ir>=0))
                                        (hthread[ir])++;
                              }
                          } else {
                              ir=(int)(sqrt(r2)*I_DR);
                              if(ir<NB_R_)          // check radial range
                                  (hthread[ir])++;
                          }
#endif // ! LOGBINCBON
                      } // ! jj=0;jj<boxes[ibox].np
                  } // ! idx=-index_max;idx<=index_max
              } // ! idy=-index_max;idy<=index_max
          } // ! idz=-index_max;idz<=index_max

#ifndef BALLS4SCANLEV
          ip = p - bodytable[cat1] + 1;
#else
          ip = ii+1;
#endif
          if (ip%gd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, ".");
          if (ip%cmd->stepState == 0)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%d ", ip);

      } // ! end pragma omp for

      if (main_thread_id == current_thread_id)
          verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "\n");
      /* Publish in thread-ID order so floating-point sums are repeatable. */
#ifdef OPENMPCODE
      int thread_id = omp_get_thread_num();
      int thread_count = omp_get_num_threads();
#else
      int thread_id = 0;
      int thread_count = 1;
#endif
      for (int thread_turn = 0; thread_turn < thread_count; thread_turn++) {
#pragma omp barrier
      if (thread_id == thread_turn) {
        for(ii=0;ii<NB_R_;ii++) //Check bound
            hh[ii]+=hthread[ii]; //Add private histograms to shared one
#ifdef SMOOTHPIVOT
        ipfalse += ipfalsethreads;
#endif
        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
      }
      }
#pragma omp barrier

    } // end pragma omp parallel
    
#ifdef SMOOTHPIVOT
    real den, num;
#ifdef BALLS4SCANLEV
    num = (real)gd->nnodescanlevTableB4[cat1];
    den = (real)(gd->nnodescanlevTableB4[cat1]-ipfalse);
#else
    num = np;
    den = (real)(np-ipfalse);
#endif
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: p falses found = %ld and num = %e, den = %e\n",
                           routineName, ipfalse, num, den);

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: p falses found = %ld\n", routineName, ipfalse);

    INTEGER ifalsecount;
    ifalsecount = 0;
    INTEGER itruecount;
    itruecount = 0;

    for(int ii=0;ii<np;ii++) {
        p = bodytable[0]+ii;
        if (Update(p) == FALSE) {
            ifalsecount++;
        } else {
            itruecount++;
        }
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: p falses found (after re-count) = %ld\n",
                           routineName, ifalsecount);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: p true found (after re-count) = %ld\n",
                           routineName, itruecount);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: total = %ld\n",
                           routineName, itruecount+ifalsecount);
#endif

    return SUCCESS;
}

#ifndef FRACTION_AR
#define FRACTION_AR 8.0
#endif

// estimates a good candidate for the size of
//  a set of nearest-neighbor searching boxes
int _optimal_nside(double lb,double rmax,INTEGER np)
{
  int nside1=(int)(FRACTION_AR*lb/rmax);
  int nside2=(int)(pow(0.5*np,0.3333333));

  return MIN(nside1,nside2);
}

#ifdef _DEBUG_
// writes catalog
//  only used for debugging
local int write_cat(struct cmdline_data* cmd,
                     struct global_data* gd, char *fn)
{
    string routineName = "write_cat";
    FILE *fr;
    INTEGER ii;
    fr=fopen(fn,"w");
    if(fr==NULL) error("error_open_file",fn);
    bodyptr p;
    for(ii=0;ii<gd->nbodyTable[0];ii++) {
        p = bodytable[0]+ii;
        WRITE_OUTPUT_OR_FAIL(fr, fn,
                             "%lf %lf %lf\n",Pos(p)[0],Pos(p)[1],Pos(p)[2]);
    }
    CLOSE_OUTPUT_OR_FAIL(fr, fn);

    return SUCCESS;
}
#endif // ! _DEBUG_


// frees all memory associated with a box
//  set of size nside
void _free_boxes(int nside,NeighborBox *boxes)
{
  int ii;
    
    for (ii = 0; ii < nside*nside*nside; ii++) {
        if (boxes[ii].pos != NULL) {
            free(boxes[ii].pos);
            boxes[ii].pos = NULL;
        }
    }
    free(boxes);
}

//B
local int _neighbor_box_index_for_body(struct cmdline_data *cmd,
                                    const char *routineName,
                                    bodyptr p,
                                    INTEGER ibody,
                                    int nside,
                                    int *index_out)
{
    int ix, iy, iz;

    if (index_out == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                "%s: index_out is NULL", routineName);
        return FAILURE;
    }

    if (nside <= 0 || lbox <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                "%s: invalid box grid nside=%d lbox=%g",
                routineName, nside, lbox);
        return FAILURE;
    }

    if (Pos(p)[0] < 0.0 || Pos(p)[0] >= lbox ||
        Pos(p)[1] < 0.0 || Pos(p)[1] >= lbox ||
        Pos(p)[2] < 0.0 || Pos(p)[2] >= lbox) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                "%s: body %ld position outside [0,lbox): %g %g %g lbox=%g",
                routineName, (long)ibody,
                Pos(p)[0], Pos(p)[1], Pos(p)[2], lbox);
        return FAILURE;
    }

    ix = (int)(Pos(p)[0] / lbox * nside);
    iy = (int)(Pos(p)[1] / lbox * nside);
    iz = (int)(Pos(p)[2] / lbox * nside);

    if (ix < 0 || ix >= nside || iy < 0 || iy >= nside || iz < 0 || iz >= nside) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                "%s: body %ld invalid box index %d %d %d for nside=%d",
                routineName, (long)ibody, ix, iy, iz, nside);
        return FAILURE;
    }

    *index_out = ix + nside * (iy + nside * iz);
    return SUCCESS;
}
//E

// creates boxes for nearest-neighbor searching
local int _catalog_to_boxes(struct cmdline_data *cmd,
                                  struct global_data *gd,
                                  int n_box_side,
                                  NeighborBox **boxes_out)
{
    string routineName = "_catalog_to_boxes";
    INTEGER ii;
    int nside;
    NeighborBox *boxes = NULL;

    if (boxes_out == NULL)
        cBALLS_FAIL(cmd, "%s: boxes_out is NULL\n", routineName);

    *boxes_out = NULL;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: Building neighbor boxes \n",
                           routineName, gd->bytes_tot*INMB);
    nside=n_box_side;
    if (nside <= 0 || lbox <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid box grid nside=%d lbox=%g",
                 routineName, nside, lbox);
        return FAILURE;
    }
    
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "%s: there will be %d boxes per side with a size of %lf\n",
                           routineName, nside,lbox/nside);

    boxes=(NeighborBox *)malloc(nside*nside*nside*sizeof(NeighborBox));
    if (boxes == NULL)
        cBALLS_FAIL(cmd, "%s: out of memory!!\n", routineName);

    for (ii = 0; ii < nside*nside*nside; ii++) {
        boxes[ii].np = 0;
        boxes[ii].pos = NULL;
    }

    for(ii=0;ii<nside*nside*nside;ii++)
    boxes[ii].np=0;

    bodyptr p;
    for(ii=0;ii<gd->nbodyTable[0];ii++) {
        int index;
        p = bodytable[0]+ii;

        if (_neighbor_box_index_for_body(cmd, routineName, p, ii, nside, &index) == FAILURE)
            goto fail;

        boxes[index].np++;
    }

    //B
    for (ii = 0; ii < nside*nside*nside; ii++) {
        int npar = boxes[ii].np;
        if (npar > 0) {
            boxes[ii].pos = (double *)malloc(3*npar*sizeof(double));
            if (boxes[ii].pos == NULL) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "%s: out of memory!!\n", routineName);
                goto fail;
            }
            boxes[ii].np = 0;
        }
    }
    //E
    
    for(ii=0;ii<gd->nbodyTable[0];ii++) {
        int index, offset;
        p = bodytable[0] + ii;

        if (_neighbor_box_index_for_body(cmd, routineName,
                                         p, ii, nside, &index) == FAILURE)
            goto fail;

        offset = 3 * boxes[index].np;
        boxes[index].pos[offset] = Pos(p)[0];
        boxes[index].pos[offset + 1] = Pos(p)[1];
        boxes[index].pos[offset + 2] = Pos(p)[2];
        boxes[index].np++;
  }


    *boxes_out = boxes;
    return SUCCESS;

    fail:
        if (boxes != NULL)
            _free_boxes(nside, boxes);
        *boxes_out = NULL;
        return FAILURE;
}

local int search_init_hist_omp(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdlptr_hist_omp hist)
{
    string routineName = "search_init_hist_omp";
    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;
    return SUCCESS;
}

// creates correlation function and poisson errors
//  from pair counts DD, DR and RR
// The correlation function is estimated as:
//      xi=(V/v(r))*(DD(r)/N^2)
//      where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
//      total # particles.
//  Note that, since in this case we have simple periodic boundary conditions,
//      no random catalogs are needed.
local int make_CF(struct cmdline_data *cmd,
                        struct global_data *gd,
                        unsigned long long DD[], int nD,
                        double corr[], double ercorr[])
{
    string routineName = "make_CF";
    double *edd;
    double rho_av=nD/(lbox*lbox*lbox);
    int ii;

    edd=(double *)malloc(sizeof(double)*NB_R_);
    if(edd==NULL)
        cBALLS_FAIL(cmd, "%s: Out of memory!!\n", routineName);

#ifdef LOGBINCBON
#ifndef _LOGBIN_
    DD[0]-=nD;                                          // substract diagonal
#endif
#else // ! LOGBINCBON
    if (!cmd->useLogHist) {
        DD[0]-=nD;                                      // substract diagonal
    }
#endif // ! LOGBINCBON

    for(ii=0;ii<NB_R_;ii++)
        edd[ii]=1./sqrt((double)DD[ii]);

    for(ii=0;ii<NB_R_;ii++) {
        if(DD[ii]==0) {
            corr[ii]=0;
            ercorr[ii]=0;
        }
        else {
            double r0,r1,vr,rho_r;
#ifdef LOGBINCBON
#ifdef _LOGBIN_
            r0=pow(10.,(((double)ii-NB_R_)/cmd->logHistBinsPD)+LOG_R_MAX_);
            r1=pow(10.,(((double)ii+1-NB_R_)/cmd->logHistBinsPD)+LOG_R_MAX_);
#else
            r0=ii/(I_R_MAX_*NB_R_);
            r1=(ii+1)/(I_R_MAX_*NB_R_);
#endif // !_LOGBIN_
#else // ! LOGBINCBON
            if (cmd->useLogHist) {
                r0=pow(10.,(((double)ii-NB_R_)/cmd->logHistBinsPD)+LOG_R_MAX_);
                r1=pow(10.,(((double)ii+1-NB_R_)/cmd->logHistBinsPD)+LOG_R_MAX_);
            } else {
                r0=ii/(I_R_MAX_*NB_R_);
                r1=(ii+1)/(I_R_MAX_*NB_R_);
            }
#endif // ! LOGBINCBON
            vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
            rho_r=DD[ii]/(nD*vr);
            corr[ii]=rho_r/rho_av-1;
            ercorr[ii]=(1+corr[ii])*edd[ii];
            gd->histNN[ii+1] = DD[ii];
            gd->histCF[ii+1] = corr[ii];
        }
    }

  free(edd);

  return SUCCESS;
}

// writes correlation function
local int _write_CF(struct  cmdline_data* cmd,
                     struct  global_data* gd,
                     char *fname,double *corr,double *ercorr,
                     unsigned long long *DD)
{
    string routineName = "_write_CF";
    FILE *fo;
    int ii;

    if (stropen_checked(fname, "w!", &fo,
                        cmd->error_message, _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;

    for(ii=0;ii<NB_R_;ii++) {
        double rr;
#ifdef LOGBINCBON
#ifdef _LOGBIN_
        rr=pow(10,( (ii+0.5)-NB_R_)/(cmd->logHistBinsPD) + LOG_R_MAX_);
#else
        rr=(ii+0.5)/(NB_R_*I_R_MAX_);
#endif // ! _LOGBIN_
#else // ! LOGBINCBON
        if (cmd->useLogHist) {
            rr=pow(10,( (ii+0.5)-NB_R_)/(cmd->logHistBinsPD) + LOG_R_MAX_);
        } else {
            rr=(ii+0.5)/(NB_R_*I_R_MAX_);
        }
#endif // ! LOGBINCBON
        WRITE_OUTPUT_OR_FAIL(fo, fname,
                             "%lE %lE %lE %llu \n",
                             rr,corr[ii],ercorr[ii],DD[ii]);
    }

    CLOSE_OUTPUT_OR_FAIL(fo, fname);

    return SUCCESS;
}

local int PrintHistNN(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "PrintHistNN";
    real rBin, rbinlog;
    int n;
    stream outstr;

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName, "w!");

    verb_print_q(2, cmd->verbose,
                "Printing : to a file %s ...\n",gd->fpfnamehistNNFileName);

    for (n=1; n<=cmd->sizeHistN; n++) {
#ifdef LOGBINCBON
#ifdef _LOGBIN_
            rbinlog = ((real)(n-0.5-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            rBin=rpow(10.0,rbinlog);
#else // ! _LOGBIN_
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
#endif // ! _LOGBIN_
#else // ! LOGBINCBON
        if (cmd->useLogHist) {
            rbinlog = ((real)(n-0.5-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
#endif // ! LOGBINCBON
        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName,
                             "%16.8e %16.8e\n",rBin,gd->histNN[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistNNFileName);

    if (cballs_opt_and_cf(cmd))
        PRINT_OR_FAIL(PrintHistCF(cmd, gd));

    return SUCCESS;
}

local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "PrintHistCF";
    real rBin, rbinlog;
    int n;
    stream outstr;
    //B correct cute-box-rmin
    real deltaR;
    if ((cballs_opt_cute_box_rmin(cmd)))
        deltaR = cmd->rangeN/cmd->sizeHistN;
    //E

    OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName, "w!");

    verb_print_q(2, cmd->verbose,
                "Printing : to a file %s ...\n",gd->fpfnamehistCFFileName);

    for (n=1; n<=cmd->sizeHistN; n++) {
#ifdef LOGBINCBON
#ifdef _LOGBIN_
        //B correct cute-box-rmin
        rbinlog = ((real)(n-0.5-cmd->sizeHistN))/cmd->logHistBinsPD
                                + rlog10(cmd->rangeN);
        rBin=rpow(10.0,rbinlog);
        //E
#else // ! _LOGBIN_
        //B correct cute-box-rmin
        if ((cballs_opt_cute_box_rmin(cmd))) {
            deltaR = cmd->rangeN/cmd->sizeHistN;
            rBin = ((real)n-0.5)*deltaR;
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        //E
#endif // ! _LOGBIN_
#else // ! LOGBINCBON
        if (cmd->useLogHist) {
            //B correct cute-box-rmin
            rbinlog = ((real)(n-0.5-cmd->sizeHistN))/cmd->logHistBinsPD
                                    + rlog10(cmd->rangeN);
            rBin=rpow(10.0,rbinlog);
            //E
        } else {
            //B correct cute-box-rmin
            if ((cballs_opt_cute_box_rmin(cmd))) {
                deltaR = cmd->rangeN/cmd->sizeHistN;
                rBin = ((real)n-0.5)*deltaR;
            } else {
                rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
            }
            //E
        }
#endif // ! LOGBINCBON

        WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName,
                             "%16.8e %16.8e\n",rBin,gd->histCF[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnamehistCFFileName);

    return SUCCESS;
}

local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "PrintHistXi2pcf";
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    if (format_checked(namebuf, sizeof(namebuf),
        "namebuf", "%s%s%s", gd->fpfnamehistXi2pcfFileName,
                       cmd->suffixOutFiles, EXTFILES) != 0)
        return FAILURE;
    verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                             "%16.8e %16.8e\n",rBin,gd->histXi2pcf[n]);
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);

    return SUCCESS;
}

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    string routineName = "print_info";

    verb_print(cmd->verbose, "Search: Running ... (neighbor-boxes-omp) \n");

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "%s: warning!! fix lengthBox accordingly to catalog box you give %s\n",
            routineName,
            "otherwise you may get wrong results or even segmentation fault...");
    if (!cballs_opt_only_pos(cmd))
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "%s: warning!! if catalog doesn´t have convergence field consider using 'only-pos' in options %s\n",
            routineName,
            "otherwise you may get wrong results or even segmentation fault...");

    if (cballs_opt_smooth_pivot(cmd))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    verb_print(cmd->verbose, "computing only 2pcf... \n");
#ifdef LOGBINCBON
    verb_print(cmd->verbose, "activated faster internal loop... \n");
#endif
#ifdef _LOGBIN_
    verb_print(cmd->verbose, "activated radial logscale... \n");
#endif
#ifdef SMOOTHPIVOT
    verb_print(cmd->verbose, "activated SMOOTHPIVOT... \n");
#endif
#ifdef ORIGINALCB
    verb_print(cmd->verbose, "activated ORIGINALCB... \n");
#endif
#ifdef _DEBUG_
    verb_print(cmd->verbose, "saving some useful debuging info... \n");
#endif

    return SUCCESS;
}
