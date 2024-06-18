/* ==============================================================================
!	MODULE: search_tree_3pcf_direct_omp.c		[cTreeBalls]
!	Written by: M.A. Rodriguez-Meza
!    Starting date:    april 2023
!    Purpose: 3-point correlation function computation
!	Language: C
!	Use: searchcalc_xxx_omp();
!	Major revisions:
!==============================================================================*/
//        1          2          3          4          5          6          7


// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"


//B COMIENZA METODO NORMAL (3pcf - brute force) DE BUSQUEDA

local void normal_walktree_nblist_omp(struct cmdline_data* cmd,
                                      struct  global_data* gd,
                                      bodyptr, nodeptr, real);
local void find_nblist_omp(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           bodyptr, cellptr, cellptr);
local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr p, INTEGER *, gdhistptr_omp_3pcfbf hist);

local int search_init_omp_3pcfbf(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_3pcfbf hist);
local int search_free_omp_3pcfbf(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_3pcfbf hist);
local int computeBodyProperties_omp_3pcfbf(struct cmdline_data* cmd,
                                           struct  global_data* gd,
                                           bodyptr, int, gdhist_omp_3pcfbf);
local int search_compute_HistN_3pcfbf(struct cmdline_data* cmd,
                                      struct  global_data* gd,
                                      int nbody);

#if !defined(FACTIVE)
#  define FACTIVE  0.75
#endif
 
local int actlen;
local int *activenb;
local int nblist;


/*
 Search serial routine using tree brute force direct method:

 To be called using: search=tree-3pcf-direct-omp

 Arguments:
    * `btable`: Input: point table array
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histN,
    *                                    histNSubXi2pcf, histNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return:
    void
 */
global int searchcalc_normal_3pcf_direct_omp(struct cmdline_data* cmd,
                                             struct  global_data* gd,
                                             bodyptr *btable, INTEGER *nbody,
                            INTEGER ipmin, INTEGER *ipmax, int cat1, int cat2)
{
    double cpustart;
    int nn;

    cpustart = CPUTIME;
    verb_print(cmd->verbose, 
               "searchcalc_normal: Running 3pcf-direct ... (tree-omp) \n");

    ThreadCount(cmd, gd, nbody[cat1], cat1);

//B Init:
    search_init_gd_hist(cmd, gd);
    if (cmd->computeTPCF) {
        int mm, ll;
        for (nn = 1; nn <= cmd->sizeHistN; nn++) {
            gd->histNNN[nn] = 0.0;
            for (mm = 1; mm <= cmd->sizeHistN; mm++)
                for (ll = 1; ll <= cmd->sizeHistTheta; ll++) {
                    gd->histXi3pcf[nn][mm][ll] = 0.0;
                    gd->histNNNSub[nn][mm][ll] = 0.0;
                }
        }
    }
//E
    
    actlen = FACTIVE * 216 * 512 * gd->tdepth;
    verb_log_print(cmd->verbose,gd->outlog,"searchcalc: actlen = %ld",actlen);
    activenb = (int *) allocate(actlen * sizeof(int));

#pragma omp parallel default(none)   shared(cmd,gd,btable,nbody,root,roottable,nblist,actlen,activenb,ipmin,ipmax,cat1,cat2)
  {
    bodyptr p;
    int n, ip;
    int m, l;
    gdhist_omp_3pcfbf hist;

    INTEGER nbbcalcthread = 0;
    search_init_omp_3pcfbf(cmd, gd, &hist);

#pragma omp for nowait schedule(dynamic)
      for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
//B Set histograms to zero for the pivot
        hist.ipcount = 0;
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;           // Affect only 3pcf evaluation
            hist.histXi2pcfthreadsub[n] = 0.0;
        }
          if (cmd->computeTPCF) {
              for (n = 1; n <= cmd->sizeHistN; n++)
                  for (m = 1; m <= cmd->sizeHistN; m++)
                      for (l = 1; l <= cmd->sizeHistTheta; l++)
                          hist.histNNNSubthread[n][m][l] = 0.0;
          }
//E
        nblist=0;

        normal_walktree_nblist_omp(cmd, gd, p, ((nodeptr) roottable[cat2]),
                                   gd->rSizeTable[cat2]);
        int_piksrt(nblist, activenb);
        verb_log_print(cmd->verbose_log, gd->outlog,
                       " - Summing nblist: %ld\n", nblist);
        sumnode_nblist_omp(cmd, gd, p,&nbbcalcthread, &hist);
        computeBodyProperties_omp_3pcfbf(cmd, gd, p, ipmax[cat1]-ipmin+1, hist);
          if (cmd->computeTPCF) {
              for (n = 1; n <= cmd->sizeHistN; n++)
                  for (m = 1; m <= cmd->sizeHistN; m++)
                      for (l = 1; l <= cmd->sizeHistTheta; l++)
                          hist.histXi3pcfthread[n][m][l] = 
                                    Kappa(p)*MAX(hist.histNNNSubthread[n][m][l],1.0);
          }
        ip = p - btable[cat1] + 1;
        if (ip%cmd->stepState == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog, " - Completed pivot: %ld\n", ip);
        }
    } // end do body p

#pragma omp critical
    {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gd->histN[n] += hist.histNthread[n];
            gd->histNSub[n] += hist.histNSubthread[n];
            gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
        }
        if (cmd->computeTPCF) {
            //        int m, l;
            for (n = 1; n <= cmd->sizeHistN; n++) {
                gd->histNNN[n] += hist.histNNNthread[n];
                for (m = 1; m <= cmd->sizeHistN; m++)
                    for (l = 1; l <= cmd->sizeHistTheta; l++) {
                        gd->histNNNSub[n][m][l] += hist.histNNNSubthread[n][m][l];
                        gd->histXi3pcf[n][m][l] += hist.histXi3pcfthread[n][m][l];
                    }
            }
        }
        gd->nbbcalc += nbbcalcthread;
    }

    search_free_omp_3pcfbf(cmd, gd, &hist);

  } // end pragma omp parallel

    for (nn = 1; nn <= cmd->sizeHistN; nn++) {
        gd->histXi2pcf[nn] /= 2.0;
        gd->histNSub[nn] = gd->histN[nn];
        gd->histNSub[nn] /= 2.0;
        gd->histXi2pcf[nn] /= MAX(gd->histNSub[nn],1.0);
    }

    if (scanopt(cmd->options, "compute-HistN"))
        search_compute_HistN_3pcfbf(cmd, gd, ipmax[cat1]-ipmin+1);

//B Computation of tpcf:
    // Make version not using GSL
#ifdef TPCF
    if (cmd->computeTPCF) {
        int i;
        double data[cmd->sizeHistTheta];
        
        gsl_fft_real_wavetable * real;
        gsl_fft_real_workspace * work;
        
        work = gsl_fft_real_workspace_alloc (cmd->sizeHistTheta);
        real = gsl_fft_real_wavetable_alloc (cmd->sizeHistTheta);
        
        for (nn = 1; nn <= cmd->sizeHistN; nn++)
            for (i = 1; i <= cmd->sizeHistN; i++) {
                for (ll = 0; ll < cmd->sizeHistTheta; ll++)
                    data[ll] = gd->histXi3pcf[nn][i][ll+1];
                gsl_fft_real_transform (data, 1, cmd->sizeHistTheta, real, work);
                for (mm = 1; mm <= cmd->mchebyshev+1; mm++) {
                    gd->histZetaM[mm][nn][i] = data[mm-1];
                }
            }
        
        gsl_fft_real_wavetable_free (real);
    }
#endif
//E

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return _SUCCESS_;
}

local void normal_walktree_nblist_omp(struct cmdline_data* cmd,
                                      struct  global_data* gd,
                                      bodyptr p, nodeptr q, real qsize)
{
    nodeptr l;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (reject_cell(cmd, gd, (nodeptr)p, q, qsize)) {
                } else {
                    for (l = More(q); l != Next(q); l = Next(l)) {
                        normal_walktree_nblist_omp(cmd, gd, p,l,qsize/2);
                    }
                }
            } else
                find_nblist_omp(cmd, gd, p, ((cellptr) q),( (cellptr) q+1));
        }
    }
}

local void find_nblist_omp(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           bodyptr p, cellptr start, cellptr finish)
{
    cellptr q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int iq;

    for (q = start; q < finish; q++) {
        if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
            gd->nbbcalc += 1;
            iq = (bodyptr)q-bodytab;
            activenb[nblist]=iq;
            nblist +=1;
            if (nblist > actlen)
                error("find_nblist: too many neighbors\n");
        }
    }
}

local void sumnode_nblist_omp(struct cmdline_data* cmd,
                              struct  global_data* gd,
                              bodyptr p, INTEGER *nbbcalcthread, gdhistptr_omp_3pcfbf hist)
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
    real theta1;

    for (i = 0; i < nblist-1; i++) {
        q = bodytab + activenb[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd->rminHist) {
            theta1 = angle_dxdy(dr[0], dr[1]);
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                hist->histXi2pcfthreadsub[n] += Kappa(q);
                if (cmd->computeTPCF) {
                    hist->histNNNthread[n] = hist->histNNNthread[n] + 1.;
                }
            }
        }

        if (cmd->computeTPCF) {
            bodyptr h;
            int m, j, l;
#ifdef SINGLEP
            float dr1_h;
            float dr_h[NDIM];
#else
            real dr1_h;
            vector dr_h;
#endif
            real theta, theta2;
            for (j = i+1; j < nblist; j++) {
                h = bodytab + activenb[j];
                accept_body(cmd, gd, p, (nodeptr)h, &dr1_h, dr_h);
                if(dr1_h>cmd->rminHist) {
                    theta2 = angle_dxdy(dr_h[0], dr_h[1]);
                    theta = rabs(theta2-theta1);
                    m = (int) ((dr1_h-cmd->rminHist) * gd->i_deltaR) + 1;
                    l = (int) (theta / gd->deltaTheta) + 1;
                    if ( (n<=cmd->sizeHistN && n>=1)
                        && (m<=cmd->sizeHistN && m>=1) && (l<=cmd->sizeHistTheta && l>=1)) {
                        hist->histNNNSubthread[n][m][l] = hist->histNNNSubthread[n][m][l] + 1.;
                        hist->histXi3pcfthread[n][m][l] += Kappa(h) * Kappa(q);
                        *nbbcalcthread += 1;
                    }
                }
            } // ! end loop j
        } // ! TPCF
    } // ! end loop i

//
    i = nblist;
        q = bodytab + activenb[i];
        accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd->rminHist) {
            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
            if (n<=cmd->sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                hist->histXi2pcfthreadsub[n] += Kappa(q);
            }
        }

}

//B Routines as in treeutils
local int search_init_omp_3pcfbf(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_3pcfbf hist)
{
    int n;
    int m;
    int l;

    hist->Chebs = dvector(1,cmd->mchebyshev+1);
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNSubthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNthread[n] = 0.0;
        hist->histNSubthread[n] = 0.0;
        hist->histXi2pcfthread[n] = 0.0;
        hist->histXi2pcfthreadsub[n] = 0.0;
    }

    hist->histNNNthread = dvector(1,cmd->sizeHistN);
    hist->histNNNSubthread = dmatrix3D(1,cmd->sizeHistN,1,cmd->sizeHistN,1,cmd->sizeHistTheta);
    hist->histXi3pcfthread = dmatrix3D(1,cmd->sizeHistN,1,cmd->sizeHistN,1,cmd->sizeHistTheta);
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNNNthread[n] = 0.0;
        for (m = 1; m <= cmd->sizeHistN; m++)
            for (l = 1; l <= cmd->sizeHistTheta; l++)
                hist->histXi3pcfthread[n][m][l] = 0.0;
    }

    return _SUCCESS_;
}

local int search_free_omp_3pcfbf(struct cmdline_data* cmd,
                                 struct  global_data* gd,
                                 gdhistptr_omp_3pcfbf hist)
{
    free_dmatrix3D(hist->histNNNSubthread,1,cmd->sizeHistN,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(hist->histXi3pcfthread,1,cmd->sizeHistN,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(hist->histNNNthread,1,cmd->sizeHistN);

    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNSubthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
    free_dvector(hist->Chebs,1,cmd->mchebyshev+1);

    return _SUCCESS_;
}

local int computeBodyProperties_omp_3pcfbf(struct cmdline_data* cmd,
                                           struct  global_data* gd,
                                           bodyptr p, int nbody, gdhist_omp_3pcfbf hist)
{
    int n;
    real xi, xi_2p;

// BODY3
    if (Type(p) == BODY) {
        xi = Kappa(p)/nbody;
        xi_2p = Kappa(p);
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
//
    for (n=1; n<=cmd->sizeHistN; n++) {
        hist.histXi2pcfthread[n] += xi_2p*hist.histXi2pcfthreadsub[n];
    }

    return _SUCCESS_;
}

local int search_compute_3pcfbf_Xi(struct cmdline_data* cmd,
                                   struct  global_data* gd,
                                   int nbody)
{
    int k;
    int n;
    real normFac;
    real Vol;

    Vol = 1.0;
    DO_COORD(k)
        Vol = Vol*gd->Box[k];

// Check logarithmic scale here:
    if (NDIM == 3)
        normFac = Vol / (2.0 * PI * rpow(gd->deltaR, 3.0) * nbody * nbody);
    else if (NDIM == 2)
        normFac = Vol / (PI * rpow(gd->deltaR, 2.0) * nbody * nbody);
    else error("\n\nWrong NDIM!\n\n");

    normFac /= 2.0;
    for (n = 1; n <= cmd->sizeHistN; n++)
        if (NDIM == 3)
            gd->histCF[n] = gd->histN[n] * normFac / rsqr(n-0.5);
        else if (NDIM == 2)
            gd->histCF[n] = gd->histN[n] * normFac / ((int)n-0.5);
        else error("\n\nWrong NDIM!\n\n");

    return _SUCCESS_;
}

local int search_compute_HistN_3pcfbf(struct cmdline_data* cmd,
                                      struct  global_data* gd,
                                      int nbody)
{
    int n;
    real normFac;

    normFac = 0.5;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histN[n] *= normFac;

    if (scanopt(cmd->options, "and-CF"))
        search_compute_3pcfbf_Xi(cmd, gd, nbody);

    return _SUCCESS_;
}
//E

//E TERMINA METODO NORMAL (3pcf - brute force) DE BUSQUEDA
