/* ==============================================================================
 MODULE: search_octree_nn_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_octree_nn_omp(cmd, gd, btable, nbody,
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
        realptr histNNSubN2pcfthread;
    //B kappa Avg Rmin
        realptr histNNSubN2pcfthreadp;
        realptr histNNSubN2pcfthreadtotal;
    //E

    real *histN2pcfthread;
    real *histN2pcfthreadsub;

    INTEGER nbbcalcthread;
    INTEGER nbccalcthread;

} gdhist_sincos_omp_nn, *gdhistptr_sincos_omp_nn;


local void normal_walktree_sincos(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  bodyptr *btable, int cat2,
                                  bodyptr, nodeptr, real,
                                  gdhistptr_sincos_omp_nn, int *, int *);
local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr, cellptr, cellptr,
                          gdhistptr_sincos_omp_nn, int *, int *);
local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_nn hist,
                               int *nbList, int *intList);

local int search_init_gd_sincos_omp_nn(struct  cmdline_data* cmd,
                                        struct  global_data* gd);
local int search_init_sincos_omp_nn(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                     gdhistptr_sincos_omp_nn hist);
local int search_free_sincos_omp_nn(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                     gdhistptr_sincos_omp_nn hist);
local int computeBodyProperties_sincos_nn(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                       gdhistptr_sincos_omp_nn hist);
local int print_info(struct cmdline_data* cmd,
                     struct  global_data* gd);

//B kappa Avg Rmin
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
                              gdhistptr_sincos_omp_nn hist, int);

#endif

//E

global int nncorrelation_octree_nn_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax,
                                             int cat1, int cat2)
{
    int n;
    real *histDD;
    real *histDR;
    real *histRD;
    real *histRR;
    real sumDD;
    real sumDR;
    real sumRD;
    real sumRR;
    real histRRfactor;
    real histDRfactor;
    real histRDfactor;

    verb_print(cmd->verbose,
               "nncorrelation: Using octree-nn-omp... \n");

    histDD = dvector(1,cmd->sizeHistN);
    histDR = dvector(1,cmd->sizeHistN);
    histRD = dvector(1,cmd->sizeHistN);
    histRR = dvector(1,cmd->sizeHistN);

    //B first DD
    verb_print(cmd->verbose,
               "nncorrelation: computing DD... \n");
    searchcalc_octree_nn_omp(cmd, gd,
                             bodytable, gd->nbodyTable,
                             1, gd->nbodyTable,
                             cat1, cat1);
    sumDD = 0.0;
    for (n=1; n<=cmd->sizeHistN; n++) {
        histDD[n] = gd->histN2pcf[n];
        sumDD += gd->histNN[n];
    }
    //E
    //B second DR
    if (!scanopt(cmd->options, "NNStandard")) {
        verb_print(cmd->verbose,
                   "nncorrelation: computing DR... \n");
        searchcalc_octree_nn_omp(cmd, gd,
                                 bodytable, gd->nbodyTable,
                                 1, gd->nbodyTable,
                                 cat1, cat2);
        sumDR = 0.0;
        for (n=1; n<=cmd->sizeHistN; n++) {
            histDR[n] = gd->histN2pcf[n];
            sumDR += gd->histNN[n];
        }
    }
    //E
    //B second RD
    if (!scanopt(cmd->options, "NNStandard")) {
        verb_print(cmd->verbose,
                   "nncorrelation: computing RD... \n");
        searchcalc_octree_nn_omp(cmd, gd,
                                 bodytable, gd->nbodyTable,
                                 1, gd->nbodyTable,
                                 cat2, cat1);
        sumRD = 0.0;
        for (n=1; n<=cmd->sizeHistN; n++) {
            histRD[n] = gd->histN2pcf[n];
            sumRD += gd->histNN[n];
        }
    }
    //E
    //B third RR
    verb_print(cmd->verbose,
               "nncorrelation: computing RR... \n");
    searchcalc_octree_nn_omp(cmd, gd,
                             bodytable, gd->nbodyTable,
                             1, gd->nbodyTable,
                             cat2, cat2);
    sumRR = 0.0;
    for (n=1; n<=cmd->sizeHistN; n++) {
        histRR[n] = gd->histN2pcf[n];
        sumRR += gd->histNN[n];
    }
    //E

    for (n=1; n<=cmd->sizeHistN; n++) {
        histRR[n] = MAX(histRR[n],1.0);
        histDR[n] = MAX(histDR[n],1.0);
        histRD[n] = MAX(histRD[n],1.0);
    }

    histRRfactor = sumDD/sumRR;
    histDRfactor = sumDD/sumDR;
    histRDfactor = sumDD/sumRD;

    //B computation of NN estimator
    //  consider in options the several stimators
    //B Landy-Szalay1: Xi = (DD - 2DR + RR)/RR
    if (scanopt(cmd->options, "NNLandySzalay1")) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            gd->histN2pcf[n] = (
                                histDD[n]
                                - 2.0*histDR[n]*histDRfactor
                                + histRR[n]*histRRfactor
                                )/(histRR[n]*histRRfactor);
        }
    }
    //B Landy-Szalay2: Xi = (DD - DR - RD + RR)/RR
    if (scanopt(cmd->options, "NNLandySzalay2")) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            gd->histN2pcf[n] = (
                                histDD[n] - histDR[n]*histDRfactor
                                - histRD[n]*histRDfactor
                                + histRR[n]*histRRfactor
                                )/(histRR[n]*histRRfactor);
        }
    }
    //B standard: Xi = DD/RR - 1
    if (scanopt(cmd->options, "NNStandard")) {
        for (n=1; n<=cmd->sizeHistN; n++) {
            gd->histN2pcf[n] = (
                                histDD[n]
                                - histRR[n]*histRRfactor
                                )/(histRR[n]*histRRfactor);
        }
    }
    //E
    //E

    //B some statistics
    double aveDD, adevDD, sdevDD, varDD, skewDD, curtDD;
    moment(histDD, cmd->sizeHistN, &aveDD, &adevDD, &sdevDD,
           &varDD, &skewDD, &curtDD);
    double aveDR, adevDR, sdevDR, varDR, skewDR, curtDR;
    moment(histDR, cmd->sizeHistN, &aveDR, &adevDR, &sdevDR,
           &varDR, &skewDR, &curtDR);
    double aveRD, adevRD, sdevRD, varRD, skewRD, curtRD;
    moment(histRD, cmd->sizeHistN, &aveRD, &adevRD, &sdevRD,
           &varRD, &skewRD, &curtRD);
    double aveRR, adevRR, sdevRR, varRR, skewRR, curtRR;
    moment(histRR, cmd->sizeHistN, &aveRR, &adevRR, &sdevRR,
           &varRR, &skewRR, &curtRR);

    real varXifactor;

    if (scanopt(cmd->options, "NNStandard"))
        varXifactor = 1.0
        + histRRfactor*aveRR/aveDD;
    if (scanopt(cmd->options, "NNLandySzalay1"))
        varXifactor = 1.0
        + 2.0*histDRfactor*aveDR/aveDD
        + histRRfactor*aveRR/aveDD;
    if (scanopt(cmd->options, "NNLandySzalay1"))
        varXifactor = 1.0
        + histDRfactor*aveDR/aveDD
        + histRDfactor*aveRD/aveDD
        + histRRfactor*aveRR/aveDD;

    real varDDnum;
    real *varRRweight;
    varRRweight = dvector(1,cmd->sizeHistN);
    varDDnum = aveDD * rsqr(varXifactor);
    for (n=1; n<=cmd->sizeHistN; n++)
        varRRweight[n] = histRR[n] * histRRfactor;

    real **covarRR;
    covarRR=dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);



    free_dmatrix(covarRR,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(varRRweight,1,cmd->sizeHistN);
    //E

    free_dvector(histRR,1,cmd->sizeHistN);
    free_dvector(histRD,1,cmd->sizeHistN);
    free_dvector(histDR,1,cmd->sizeHistN);
    free_dvector(histDD,1,cmd->sizeHistN);

    return SUCCESS;
}


/*
 Search routine using octree method:

 To be called using: search=octree-nn-omp

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
    *                                    histNNSubN2pcf, histNNSubN2pcftotal,
    *                                    histN2pcf, histN,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int searchcalc_octree_nn_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                                             INTEGER ipmin, INTEGER *ipmax, 
                                             int cat1, int cat2)
{
    double cpustart;

    cpustart = CPUTIME;
    print_info(cmd, gd);

//B kappa Avg Rmin
#ifdef DEBUG
    sprintf(pivotsfilePath,"%s/pivot_info%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outpivots=fopen(pivotsfilePath, "w")))
        error("\nsearchcalc_tc_nn_omp: error opening file '%s' \n",
              pivotsfilePath);
#endif
//E

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    search_init_gd_sincos_omp_nn(cmd, gd);

    //B kappa Avg Rmin
    INTEGER ipfalse;
    ipfalse=0;
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;

#ifdef ADDPIVOTNEIGHBOURS
    if (scanopt(cmd->options, "smooth-pivot")) {
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
    }
#endif
    //E

    verb_print(cmd->verbose,
        "\nsearchcalc_tc_nn_omp: Total allocated %g MByte storage so far.\n",
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
#endif // ! DEBUG
  {
    bodyptr p;
    bodyptr q;
    int n, m, ip;
    int i;

    //B init:
    gdhist_sincos_omp_nn hist;
    search_init_sincos_omp_nn(cmd, gd, &hist);
    //E

    //B kappa Avg Rmin
    INTEGER ipfalsethreads;
    ipfalsethreads = 0;
    INTEGER icountNbRminthread;
    icountNbRminthread=0;
    INTEGER icountNbRminOverlapthread;
    icountNbRminOverlapthread=0;

//#ifdef ADDPIVOTNEIGHBOURS
    int nbList;
    int intList;
//#endif
    //E

#pragma omp for nowait schedule(dynamic)
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
          //B kappa Avg Rmin
          NbRmin(p) = 1;
          NbRminOverlap(p) = 0;
          WeightRmin(p) = Weight(p);
//#ifdef ADDPIVOTNEIGHBOURS
          nbList=0;
          intList=0;
//#endif
          if (scanopt(cmd->options, "smooth-pivot")) {
              if (Update(p) == FALSE) {
                  ipfalsethreads++;
                  continue;
              }
          }
          //E

//B segment to be included below...
          //B Set histograms to zero for the pivot
          for (n = 1; n <= cmd->sizeHistN; n++) {
              hist.histN2pcfthreadsub[n] = 0.0;  // Affects only 2pcf
//B kappa Avg Rmin
              hist.histNNSubN2pcfthreadp[n] = 0.;// Affects only 2pcf
//E
          }
          //E
//E

          normal_walktree_sincos(cmd, gd, btable, cat2,
                                 p, ((nodeptr) roottable[cat2]),
                                 gd->rSizeTable[cat2], &hist, &nbList, &intList);
          computeBodyProperties_sincos_nn(cmd, gd, p,
                                           ipmax[cat1]-ipmin+1, &hist);
#ifdef ADDPIVOTNEIGHBOURS
          if (scanopt(cmd->options, "smooth-pivot")) {
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
                  computeBodyProperties_sincos_nn(cmd, gd, q,
                                                   ipmax[cat1]-ipmin+1, &hist);
              } // ! end i loop
          } // ! scanoption smooth-pivot
#endif // ! ADDPIVOTNEIGHBOURS

          ip = p - btable[cat1] + 1;
          //B kappa Avg Rmin
          icountNbRminthread += NbRmin(p);
          icountNbRminOverlapthread += NbRminOverlap(p);
#ifdef DEBUG
#ifdef ADDPIVOTNEIGHBOURS
          fprintf(outpivots,"%ld \t%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p), intList,
                  WeightRmin(p)/NbRmin(p));
#else
          fprintf(outpivots,"%ld \t%ld \t%ld \t\t%g\n",
                  ip, NbRmin(p), NbRminOverlap(p),
                  WeightRmin(p)/NbRmin(p));
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
            gd->histNNSubN2pcf[n] += hist.histNNSubN2pcfthread[n];
//B kappa Avg Rmin
            gd->histNNSubN2pcftotal[n] += hist.histNNSubN2pcfthreadtotal[n];
//E
            gd->histN2pcf[n] += hist.histN2pcfthread[n];
        }

        gd->nbbcalc += hist.nbbcalcthread;
        gd->nbccalc += hist.nbccalcthread;
        //B kappa Avg Rmin
        ipfalse += ipfalsethreads;
        icountNbRmin += icountNbRminthread;
        icountNbRminOverlap += icountNbRminOverlapthread;
        //E
    } // ! critical

    search_free_sincos_omp_nn(cmd, gd, &hist);     // free memory
  } // end pragma omp parallel

    if (cmd->verbose >= VERBOSENORMALINFO)
        verb_print(cmd->verbose, "\n\n");             // end of completed pivot

    //B kappa Avg Rmin
    real xi, den, num;
    int mm;
    if (scanopt(cmd->options, "smooth-pivot")) {
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
    }
    //E

    int nn;

    if (!scanopt(cmd->options, "asymmetric")) {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
                printf("%d %e %e\n", nn,
                   gd->histNNSubN2pcf[nn], gd->histNNSubN2pcftotal[nn]);
            gd->histN2pcf[nn] /= 2.0;
            gd->histNNSubN2pcf[nn] /= 2.0;
//B kappa Avg Rmin
            gd->histNNSubN2pcftotal[nn] /= 2.0;
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histN2pcf[nn] /= MAX(gd->histNNSubN2pcftotal[nn],1.0);
            } else {
                gd->histN2pcf[nn] /= MAX(gd->histNNSubN2pcf[nn],1.0);
            }
//E
        }
    } else {
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            if (cmd->verbose>3)
            printf(0,"%d %e %e\n", nn,
                   gd->histNNSubN2pcf[nn], gd->histNNSubN2pcftotal[nn]);
            if (scanopt(cmd->options, "smooth-pivot")) {
                gd->histN2pcf[nn] /= MAX(gd->histNNSubN2pcftotal[nn],1.0);
            } else {
                gd->histN2pcf[nn] /= MAX(gd->histNNSubN2pcf[nn],1.0);
            }
        }
    }

    if (scanopt(cmd->options, "compute-HistN")) {
        if (scanopt(cmd->options, "smooth-pivot")) {
            search_compute_HistN(cmd, gd, nbody[cat1]-ipfalse);
        } else {
            search_compute_HistN(cmd, gd, nbody[cat1]);
        }
    }

    if (scanopt(cmd->options, "smooth-pivot")) {
        if (cmd->verbose>=VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "octree-nn-omp: p falses found = %ld\n",ipfalse);
            //B kappa Avg Rmin
            verb_print(cmd->verbose,
                       "octree-nn-omp: count NbRmin found = %ld\n",
                       icountNbRmin);
            verb_print(cmd->verbose,
                       "octree-nn-omp: count overlap found = %ld\n",
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
            verb_print(cmd->verbose, "octree-nn-omp: p falses found = %ld\n",
                       ifalsecount);
            verb_print(cmd->verbose, "octree-nn-omp: p true found = %ld\n",itruecount);
            verb_print(cmd->verbose, "octree-nn-omp: total = %ld\n",
                       itruecount+ifalsecount);
        }
        //E
    }

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
                                  gdhistptr_sincos_omp_nn hist,
                                  int *nbList, int *intList)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd->options, "no-one-ball")) {
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

#ifdef ADDPIVOTNEIGHBOURS
//B kappa Avg Rmin
local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr *btable,
                              INTEGER ipmin, INTEGER *ipmax,
                              int cat1, int cat2,
                              bodyptr p,
                              gdhistptr_sincos_omp_nn hist, int intList)
{
    bodyptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
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
                hist->histNNSubN2pcfthread[n] =
                hist->histNNSubN2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubN2pcfthreadp[n] =
                hist->histNNSubN2pcfthreadp[n] + 1.;
                //E
                xi = Weight(q);
                hist->histN2pcfthreadsub[n] += xi;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! end loop i
}
//E
#endif

local void sumnode_sincos(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr *btable, int cat2,
                          bodyptr p, cellptr start, cellptr finish,
                          gdhistptr_sincos_omp_nn hist,
                          int *nbList, int *intList)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xi;
    int iq;

    q = start;
    if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
        //B kappa Avg Rmin
        if (scanopt(cmd->options, "smooth-pivot")) {
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
                    WeightRmin(p) += Weight(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
        }
        //E
        if(dr1>cmd->rminHist) {
#ifdef ADDPIVOTNEIGHBOURS
            if (scanopt(cmd->options, "smooth-pivot")) {
                iq = (bodyptr)q-btable[cat2];
                activeInt[*intList]=iq;
                *intList +=1;
                if (*intList > actlenInt)
                    error("intList: too many neighbors\n");
            }
#endif
            if (cmd->rminHist==0)
                n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                    - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNNSubN2pcfthread[n] =
                hist->histNNSubN2pcfthread[n] + 1.;
                //B kappa Avg Rmin
                hist->histNNSubN2pcfthreadp[n] =
                hist->histNNSubN2pcfthreadp[n] + 1.;
                //E
                xi = Weight(q);
                hist->histN2pcfthreadsub[n] += xi;
                hist->nbbcalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1>cmd->rminHist
    } // ! accept_body
}

local void sumnode_sincos_cell(struct  cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr *btable, int cat2,
                               bodyptr p, cellptr start, cellptr finish,
                               gdhistptr_sincos_omp_nn hist,
                               int *nbList, int *intList)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;
    real xi;
    REAL cosphi,sinphi;

    q = start;
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
                hist->histNNSubN2pcfthread[n] =
                hist->histNNSubN2pcfthread[n] + 1.0;
                xi = Weight(q);
                hist->histN2pcfthreadsub[n] += xi;
                hist->nbccalcthread += 1;
            } // ! 1 < n < sizeHistN
        } // ! dr1 > rminHist
    } // ! accept_body
}

//B Routines as in cballsutils

local int search_init_gd_sincos_omp_nn(struct  cmdline_data* cmd,
                                        struct  global_data* gd)
{
    int n;
    int m;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNNSub[n] = 0.0;

    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubN2pcf[n] = 0.0;
//B kappa Avg Rmin
        gd->histNNSubN2pcftotal[n] = 0.0;
//E
        gd->histN2pcf[n] = 0.0;
    }

    gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

local int search_init_sincos_omp_nn(struct  cmdline_data* cmd,
                                     struct  global_data* gd,
                                     gdhistptr_sincos_omp_nn hist)
{
    int n;
    
    hist->histNthread = dvector(1,cmd->sizeHistN);
// 2pcf
    hist->histNNSubN2pcfthread = dvector(1,cmd->sizeHistN);
//B kappa Avg Rmin
    hist->histNNSubN2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubN2pcfthreadtotal = dvector(1,cmd->sizeHistN);
//E
//
    hist->histN2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histN2pcfthreadsub = dvector(1,cmd->sizeHistN);

    hist->nbbcalcthread = 0;
    hist->nbccalcthread = 0;

    return SUCCESS;
}

local int search_free_sincos_omp_nn(struct  cmdline_data* cmd,
                                         struct  global_data* gd,
                                      gdhistptr_sincos_omp_nn hist)
{
    free_dvector(hist->histN2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histN2pcfthread,1,cmd->sizeHistN);
//B kappa Avg Rmin
    free_dvector(hist->histNNSubN2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubN2pcfthreadp,1,cmd->sizeHistN);
//E
    free_dvector(hist->histNNSubN2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);

    return SUCCESS;
}

local int computeBodyProperties_sincos_nn(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp_nn hist)
{
    int n;
    int m;
    real xi_2p;

// check Weight factor... must be an average of Weights
#ifdef ADDPIVOTNEIGHBOURS
    xi_2p = Weight(p);
#else
    xi_2p = Weight(p);
    //B kappa Avg Rmin
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi_2p = (WeightRmin(p)/NbRmin(p));
    }
    //E
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histN2pcfthread[n] += xi_2p*hist->histN2pcfthreadsub[n];
    }

    return SUCCESS;
}


//E Routines as in cballsutils

local int print_info(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    verb_print(cmd->verbose,
               "searchcalc: Using octree-nn-omp... \n");

    if (cmd->usePeriodic==TRUE)
        error("CheckParameters: can´t have periodic boundaries and OCTREEKKKOMP definition (usePeriodic=%d)\nSet usePeriodic=false\n",
            cmd->usePeriodic);
    if (cmd->useLogHist==FALSE)
        error("CheckParameters: can´t have normal scale hist and OCTREEKKKOMP definition (useLogHist=%d)\nSet useLogHist=true\n",
            cmd->useLogHist);
    if (cmd->computeTPCF==FALSE)
        error("CheckParameters: can´t have computeTPCF=false and OCTREEKKKOMP definition (computeTPCF=%d)\nSet computeTPCF=true\n",
            cmd->computeTPCF);
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
#ifdef ADDPIVOTNEIGHBOURS
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose, "with ADDPIVOTNEIGHBOURS... \n");
#endif
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (scanopt(cmd->options, "smooth-pivot"))
        verb_print(cmd->verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
    if (scanopt(cmd->options, "default-rsmooth"))
        verb_print(cmd->verbose, "with option default-rsmooth... \n");
    if (scanopt(cmd->options, "fix-rsmooth"))
        verb_print(cmd->verbose, "with option fix-rsmooth... \n");

    return SUCCESS;
}


#ifdef ADDPIVOTNEIGHBOURS
#undef FACTIVENB
#undef FACTIVEINT
#endif

#endif // ! THREEDIMCODE
