/* ==============================================================================
 MODULE: search_octree_kk_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 2-point correlation function computation
 Language: C
 Use: searchcalc_octree_kk_omp(cmd, gd, btable, nbody,
                                           ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

#if THREEDIMCODE

typedef struct {
    realptr histNthread;
    // 2pcf
        realptr histNNSubXi2pcfthread;
    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
        realptr histNNSubXi2pcfthreadp;
        realptr histNNSubXi2pcfthreadtotal;
#endif
    //E

    real *histXi2pcfthread;
    real *histXi2pcfthreadsub;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

} gdhist_sincos_omp_kk, *gdhistptr_sincos_omp_kk;


local void normal_walktree_sincos(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr, nodeptr, real,
                                  gdhistptr_sincos_omp_kk, int *, int *);
local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr, cellptr, cellptr,
                          gdhistptr_sincos_omp_kk, int *, int *);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_kk hist,
                               int *nbList, int *intList);

local int search_init_gd_sincos_omp_kk(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
local int search_init_sincos_omp_kk(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist);
local int search_free_sincos_omp_kk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist);
local int computeBodyProperties_sincos_kk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_kk hist);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
#ifdef DEBUG
local char pivotsfilePath[MAXLENGTHOFFILES];
local FILE *outpivots;
#endif

#ifdef ADDPIVOTNEIGHBOURS
#define FACTIVENB  10
#define FACTIVEINT  10
 
local int actlenNb;
local int *activeNb;

local int actlenInt;
local int *activeInt;

local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_kk hist, int);

#endif
#endif
//E

/*
 Search routine using octree method:

 To be called using: search=octree-kk-omp

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `btable`: Input: point table array
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histN,
    *                                    histNNSubXi2pcf, histNNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int searchcalc_octree_kk_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax, 
                                             int cat1, int cat2)
{
    double cpustart;

    cpustart = CPUTIME;
    print_info(cmd, gd);

//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
#ifdef DEBUG
    sprintf(pivotsfilePath,"%s/pivot_info%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outpivots=fopen(pivotsfilePath, "w")))
        error("\nsearchcalc_tc_kk_omp: error opening file '%s' \n",
              pivotsfilePath);
#endif
#endif
//E

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_sincos_omp_kk(cmd, gd);

    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
#endif

#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
//    if (cballs_opt_smooth_pivot(cmd)) {
        //B Alloc memory for neighbour lists
        int Nrsmooth;
        Nrsmooth = 0.25*nbody[cat2]*rsqr(gd->rsmooth[0]);
        actlenNb = FACTIVENB * Nrsmooth;
        verb_print(cmd->verbose, "\n- Nrsmooth and actlenNb: %d %d\n",
                   Nrsmooth, actlenNb);
        verb_log_print(cmd->verbose,gd->outlog,
                       "searchcalc: actlenNb = %ld\n",actlenNb);
        activeNb = (int *) allocate(actlenNb * sizeof(int));
        //E
        //B Alloc memory for intercation lists
        int NrangeN;
        NrangeN = 0.25*nbody[cat2]*rsqr(cmd->rangeN);
        actlenInt = FACTIVEINT * NrangeN;
        verb_print(cmd->verbose, "- NrangeN and actlenInt: %d %d\n",
                   NrangeN, actlenInt);
        verb_log_print(cmd->verbose,gd->outlog,
                       "searchcalc: actlenInt = %ld\n",actlenInt);
        activeInt = (int *) allocate(actlenInt * sizeof(int));
        //E
//    }
#endif
#endif
    //E

    verb_print(cmd->verbose,
        "\nsearchcalc_tc_kk_omp: Total allocated %g MByte storage so far.\n",
        gd->bytes_tot/(1024.0*1024.0));

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose,
                   "\nRunning...\n - Completed pivot node:\n");

#ifdef DEBUG
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,outpivots,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#endif
#else // ! DEBUG
#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           actlenNb,activeNb,actlenInt,activeInt,                           \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#else
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           ipmin,ipmax,cat1,cat2,ipfalse,icountNbRmin,icountNbRminOverlap)
#endif
#else // ! SMOOTHPIVOT
#pragma omp parallel default(none)                                          \
    shared(cmd,gd,btable,nbody,roottable,                         \
           ipmin,ipmax,cat1,cat2)
#endif
#endif // ! DEBUG
  {
    bodyptr p;
    bodyptr q;
    int n, m, ip;
    int i;

    //B init:
    gdhist_sincos_omp_kk hist;
    search_init_sincos_omp_kk(cmd, gd, &hist);
    //E

    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    INTEGER ipfalsethreads;
    ipfalsethreads = 0;
    INTEGER icountNbRminthread;
    icountNbRminthread=0;
    INTEGER icountNbRminOverlapthread;
    icountNbRminOverlapthread=0;
#endif

//#ifdef ADDPIVOTNEIGHBOURS
    int nbList;
    int intList;
//#endif
    //E

#pragma omp for nowait schedule(dynamic)
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
          //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
          NbRmin(p) = 1;
          NbRminOverlap(p) = 0;
          KappaRmin(p) = Kappa(p);
#endif
//#ifdef ADDPIVOTNEIGHBOURS
          nbList=0;
          intList=0;
//#endif
#ifdef SMOOTHPIVOT
//          if (cballs_opt_smooth_pivot(cmd)) {
              if (Update(p) == FALSE) {
                  ipfalsethreads++;
                  continue;
              }
//          }
#endif
          //E

//B segment to be included below...
          //B Set histograms to zero for the pivot
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histXi2pcfthreadsub[n] = 0.0;  // Affects only 2pcf
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
              hist.histNNSubXi2pcfthreadp[n] = 0.;// Affects only 2pcf
#endif
//E
          }
          //E
//E

          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist, &nbList, &intList);
          computeBodyProperties_sincos_kk(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
//          if (cballs_opt_smooth_pivot(cmd)) {
              if (cmd->verbose_log>=3)
                  verb_log_print(cmd->verbose_log, gd->outlog,
                                 " - Summing nbList: %ld\n", nbList);
              for (i = 0; i < nbList; i++) {        // loop over neighbours
                  if (cmd->verbose_log>=3)
                      verb_log_print(cmd->verbose_log, gd->outlog,
                                     " - Summing intList: %ld\n", intList);
                  q = btable[cat1] + activeNb[i];

                  sumnode_nblist_omp(cmd, gd, btable, ipmin, ipmax, cat1, cat2,
                                     q, &hist, intList);
                  computeBodyProperties_sincos_kk(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
              } // ! end i loop
//          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS
#endif

          ip = p - btable[cat1] + 1;
          //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
          icountNbRminthread += NbRmin(p);
          icountNbRminOverlapthread += NbRminOverlap(p);
#endif
#ifdef DEBUG
#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
          fprintf(outpivots,"%ld \t%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p), intList,
                  KappaRmin(p)/NbRmin(p));
#else
          fprintf(outpivots,"%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p),
                  KappaRmin(p)/NbRmin(p));
#endif
#endif
#endif
          //E
          if (cmd->verbose >= VERBOSENORMALINFO) {
              if (ip%cmd->stepState == 0) {
                  verb_print(cmd->verbose, "%d\n", ip);
              }
          } else
              if (ip%cmd->stepState == 0) {
                  verb_log_print(cmd->verbose_log, gd->outlog,
                                 " - Completed pivot: %ld\n", ip);
              }
      } // end do body p

#pragma omp critical
    {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gd->histNN[n] += hist.histNthread[n];
            gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
            gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
#endif
//E
            gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
        }

        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
        //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
#endif
        //E
    } // ! critical

    search_free_sincos_omp_kk(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose, "\n\n");             // end of completed pivot

    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    real xi, den, num;
    int mm;
//    if (cballs_opt_smooth_pivot(cmd)) {
        num = (real)nbody[cat1];
        den = (real)(nbody[cat1]-ipfalse);
#ifdef ADDPIVOTNEIGHBOURS
        xi = 1.0;
#else
        xi = num/den;
#endif
        if (cmd->verbose>=VERBOSENORMALINFO)
            verb_print(cmd->verbose,
                       "octree-kkk-omp: p falses found = %ld and %e %e %e\n",
                       ipfalse, num, den, xi);
//    }
#endif
    //E

    int nn;

    if (!cballs_opt_asymmetric(cmd)) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
                printf("%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
#else
            if (cmd->verbose>3)
                printf("%d %e\n", nn,
                   gd->histNNSubXi2pcf[nn]);
#endif
            gd->histXi2pcf[nn] /= 2.0;
            gd->histNNSubXi2pcf[nn] /= 2.0;
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
            gd->histNNSubXi2pcftotal[nn] /= 2.0;
//            if (cballs_opt_smooth_pivot(cmd)) {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
//            } else {
#else
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
//            }
#endif
//E
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
#ifdef SMOOTHPIVOT
            if (cmd->verbose>3)
            printf(0,"%d %e %e\n", nn,
                   gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
//            if (cballs_opt_smooth_pivot(cmd)) {
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
//            } else {
#else
            if (cmd->verbose>3)
            printf(0,"%d %e\n", nn,
                   gd->histNNSubXi2pcf[nn]);
                gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
//            }
#endif
        }
    }

    if (cballs_opt_compute_histn(cmd)) {
#ifdef SMOOTHPIVOT
//        if (cballs_opt_smooth_pivot(cmd)) {
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
//        } else {
#else
            search_compute_HistN(cmd, gd, nbody[cat1]);
//        }
#endif
    }

#ifdef SMOOTHPIVOT
//    if (cballs_opt_smooth_pivot(cmd)) {
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "octree-kkk-omp: p falses found = %ld\n",ipfalse);
            //B kappa Avg Rmin
            verb_print(cmd->verbose,
                       "octree-kkk-omp: count NbRmin found = %ld\n",
                       icountNbRmin);
            verb_print(cmd->verbose,
                       "octree-kkk-omp: count overlap found = %ld\n",
                       icountNbRminOverlap);
        }
        
        bodyptr pp;
        INTEGER ifalsecount;
        ifalsecount = 0;
        INTEGER itruecount;
        itruecount = 0;
        for (pp = btable[cat1] + ipmin -1; pp < btable[cat1] + ipmax[cat1]; pp++) {
            if (Update(pp) == FALSE) {
                ifalsecount++;
            } else {
                itruecount++;
            }
        }
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose, "octree-kkk-omp: p falses found = %ld\n",
                       ifalsecount);
            verb_print(cmd->verbose, "octree-kkk-omp: p true found = %ld\n",itruecount);
            verb_print(cmd->verbose, "octree-kkk-omp: total = %ld\n",
                       itruecount+ifalsecount);
        }
        //E
//    }
#endif

#ifdef DEBUG
    fclose(outpivots);                              // Close file to debug pivots
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nGoing out: CPU time = %lf %s\n",
               CPUTIME-cpustart, PRNUNITOFTIMEUSED);

    return SUCCESS;
}

local void normal_walktree_sincos(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr p, nodeptr q, real qsize,
                                  gdhistptr_sincos_omp_kk hist,
                                  int *nbList, int *intList)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    compute_vector dr;
#endif

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                    if (!cballs_opt_no_one_ball(cmd)) {
                        accept_body(cmd, gd,
                                    p, (nodeptr)q, &dr1, dr);
                        if ( (Radius(p)+Radius(q))/(dr1) < gd->deltaR)
                            sumnode_sincos_cell(cmd, gd, btable, cat2, p,
                                                ((cellptr) q),
                                                ((cellptr) q+1), hist,
                                                nbList, intList);
                        else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(cmd, gd, btable, cat2,
                                                       p,l,qsize/2, hist,
                                                       nbList, intList);
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(cmd, gd, btable, cat2,
                                                   p,l,qsize/2, hist,
                                                   nbList, intList);
                    }
                }
            } else { // ! Type(q) == CELL
                sumnode_sincos(cmd, gd, btable, cat2,
                               p, ((cellptr)q), ((cellptr)q+1), hist,
                               nbList, intList);
            } // ! Type(q) == CELL
        } // ! p != q
    }
}

#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
//B kappa Avg Rmin
local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_kk hist, int intList)
{
    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    compute_vector dr;
#endif
    int i;
    int n;
    real xi;

    for (i = 0; i < intList; i++) {
        q = btable[cat2] + activeInt[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd->rminHist) {
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                  + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubXi2pcfthreadp[n] =
                hist->histNNSubXi2pcfthreadp[n] + 1.;
                //E
                xi = Weight(q)*Kappa(q);
                hist->histXi2pcfthreadsub[n] += xi;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! end loop i
}
//E
#endif
#endif

local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr p, cellptr start, cellptr finish,
                          gdhistptr_sincos_omp_kk hist,
                          int *nbList, int *intList)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    compute_vector dr;
#endif
    int n;
    real xi;
    int iq;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
//        if (cballs_opt_smooth_pivot(cmd)) {
            if (dr1<=gd->rsmooth[0]) {
                if (Update(q)==TRUE) {
#ifdef ADDPIVOTNEIGHBOURS
                    iq = (bodyptr)q-btable[cat2];
                    activeNb[*nbList]=iq;
                    *nbList +=1;
                    if (*nbList > actlenNb)
                        error("nbList: too many neighbors, %d %d\n",
                              *nbList, actlenNb);
#endif
                    Update(q) = FALSE;
                    NbRmin(p) += 1;
                    KappaRmin(p) += Weight(q)*Kappa(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
//        }
#endif
        //E
        if(dr1>cmd->rminHist) {
#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
//            if (cballs_opt_smooth_pivot(cmd)) {
                iq = (bodyptr)q-btable[cat2];
                activeInt[*intList]=iq;
                *intList +=1;
                if (*intList > actlenInt)
                    error("intList: too many neighbors\n");
//            }
#endif
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNNSubXi2pcfthread[n] =
                hist->histNNSubXi2pcfthread[n] + 1.;
                //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
                hist->histNNSubXi2pcfthreadp[n] =
                hist->histNNSubXi2pcfthreadp[n] + 1.;
#endif
                //E
                xi = Weight(q)*Kappa(q);
                hist->histXi2pcfthreadsub[n] += xi;
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist
    } // ! accept_body
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_kk hist,
                               int *nbList, int *intList)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    compute_vector dr;
#endif
    int n;
    real xi;
    REAL cosphi,sinphi;

    q = start;
    if (cballs_opt_read_mask(cmd))
        if (Mask(q) != MASK_NODE_VALID) return;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            INTEGER iq;
            iq = (bodyptr)q-btable[cat2];
            activeInt[*intList]=iq;
            *intList +=1;
            if (*intList > actlenInt)
                error("intList: too many neighbors\n");
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histNNSubXi2pcfthread[n] =
                hist->histNNSubXi2pcfthread[n] + Nb(q);
                #ifndef NOWKAvg
                xi = Weight(q)*Kappa(q);
#else
                xi = Nb(q)*Kappa(q);
#endif
                hist->histXi2pcfthreadsub[n] += xi;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! accept_body
}

//B Routines as in cballsutils

local int search_init_gd_sincos_omp_kk(struct  cmdline_data* cmd,
                                        struct  global_data* gd)
{
    int n;
    int m;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNNSub[n] = 0.0;

    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
//E
        gd->histXi2pcf[n] = 0.0;
    }

    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_init_sincos_omp_kk(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_kk hist)
{
    int n;
    
    hist->histNthread = dvector(1,cmd->sizeHistN);
// 2pcf
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
//E
//
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_kk(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_kk hist)
{
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
//B kappa Avg Rmin
#ifdef SMOOTHPIVOT
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
//E
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);

    return SUCCESS;
}

local int computeBodyProperties_sincos_kk(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_kk hist)
{
    int n;
    int m;
    real xi_2p;

// check Weight factor... must be an average of Weights
#ifdef ADDPIVOTNEIGHBOURS
    xi_2p = Weight(p)*Kappa(p);
#else
    xi_2p = Weight(p)*Kappa(p);
    //B kappa Avg Rmin
#ifdef SMOOTHPIVOT
//    if (cballs_opt_smooth_pivot(cmd)) {
        xi_2p = (KappaRmin(p)/NbRmin(p));
//    }
#endif
    //E
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
    }

    return SUCCESS;
}


//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc: Using octree-kk-omp... \n");

    if (cmd->usePeriodic==TRUE)
        error("CheckParameters: can´t have periodic boundaries and OCTREEKKKOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        error("CheckParameters: can´t have normal scale hist and OCTREEKKKOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
            cmd->useLogHist);
#ifndef TPCF
        error("CheckParameters: can´t have TPCF OFF and OCTREEKKKOMP definition \nSet TPCF ON\n");
#endif
#if NDIM == 2
    if (cmd->computeTPCF) {
        verb_print(cmd->verbose,
                   "searchcalc: Using octree-kkk-omp... \n");
        verb_print(cmd->verbose,
                   "and with computing TPCF... \n");
        verb_print(cmd->verbose,
                   "CheckParameters: \n");
        error("OCTREEKKKOMP definition works only in a 3D unit sphere")
    }
#endif
#ifdef POLARAXIS
    verb_print(cmd->verbose, "with POLARAXIS... \n");
#endif
#ifdef SMOOTHPIVOT
#ifdef ADDPIVOTNEIGHBOURS
//    if (cballs_opt_smooth_pivot(cmd))
        verb_print(cmd->verbose, "with ADDPIVOTNEIGHBOURS... \n");
#endif
#endif
    if (cballs_opt_no_one_ball(cmd))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
#ifdef SMOOTHPIVOT
//    if (cballs_opt_smooth_pivot(cmd))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
#endif
    if (cballs_opt_default_rsmooth(cmd))
        verb_print(cmd->verbose, "with option default-rsmooth... \n");
    if (cballs_opt_fix_rsmooth(cmd))
        verb_print(cmd->verbose, "with option fix-rsmooth... \n");

    return SUCCESS;
}


#ifdef ADDPIVOTNEIGHBOURS
#undef FACTIVENB
#undef FACTIVEINT
#endif

#endif // ! THREEDIMCODE
