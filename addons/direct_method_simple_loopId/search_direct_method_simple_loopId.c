/* ==============================================================================
 MODULE: search_direct_simple_loopId.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_direct_simple_sincos_loopId(cmd, gd,
                                      btable, nbody, ipmin, ipmax, cat1, cat2);
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

//B chebyshev polynomial definition
#define CHEBYSHEVTUSINCOS                                   \
{real xicosmphi,xisinmphi; int m;                           \
    hist.ChebsT[1] = 1.0;                                   \
    xicosmphi = xi * hist.ChebsT[1];                        \
    gd->histXicos[1][n] += xicosmphi;                       \
    hist.ChebsT[2] = cosphi;                                \
    xicosmphi = xi * hist.ChebsT[2];                        \
    gd->histXicos[2][n] += xicosmphi;                       \
    hist.ChebsT[3] = 2.0*(cosphi)*(cosphi) - (1.0);         \
    xicosmphi = xi * hist.ChebsT[3];                        \
    gd->histXicos[3][n] += xicosmphi;                       \
    hist.ChebsU[1] = 0.0;                                   \
    xisinmphi = xi * hist.ChebsU[1] * sinphi;               \
    gd->histXisin[1][n] += xisinmphi;                       \
    hist.ChebsU[2] = 1.0;                                   \
    xisinmphi = xi * hist.ChebsU[2] * sinphi;               \
    gd->histXisin[2][n] += xisinmphi;                       \
    hist.ChebsU[3] = 2.0*cosphi;                            \
    xisinmphi = xi * hist.ChebsU[3] * sinphi;               \
    gd->histXisin[3][n] += xisinmphi;                       \
    for (m=4; m<=cmd->mChebyshev+1; m++){                   \
        hist.ChebsT[m] = 2.0*(cosphi)*hist.ChebsT[m-1] - hist.ChebsT[m-2]; \
        xicosmphi = xi * hist.ChebsT[m];                    \
        gd->histXicos[m][n] += xicosmphi;                   \
        hist.ChebsU[m] = 2.0*(cosphi)*hist.ChebsU[m-1] - hist.ChebsU[m-2]; \
        xisinmphi = xi * hist.ChebsU[m] * sinphi;           \
        gd->histXisin[m][n] += xisinmphi;                   \
    }}
//E

//B Define a structure:
typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real **xiOUTVPcossin;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real **histZetaMtmpcossin;
    real *ChebsT;
    real *ChebsU;
    real *histXi2pcfsub;
    //B weights-norm
    realptr histW;
    realptr histWW;
    //E

    compute_vector q0;
    real drpq2, drpq;
    compute_vector dr0;
} gdhist_sincos_direct, *gdhistptr_sincos_direct;
//E

local int search_init_sincos_gd(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct);
local int search_init_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct);
local int search_free_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct);
local int computeBodyProperties_sincos_direct(struct cmdline_data* cmd,
                                       struct  global_data* gd,
                                       bodyptr, int, gdhistptr_sincos_direct);

/*
 Search serial routine using direct method:

 To be called using: search=direct-simple-sincos

 Arguments:
    * `btable`: Input: point table array
    * `nbody`: Input: number of points in table array
    * `ipmin`: Input: minimum point in table array to analyse
    * `ipmax`: Input: maximum point in table array to analyse
    * `cat1`: Input: catalog tag to act as pivot catalog
    * `cat2`: Input: catalog tag to act as a scanning catalog
 Global tructures used: gd, cmd
 Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histN,
    *                                    histNNSubXi2pcf, histNNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
 Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return:
    void
 */
global int searchcalc_direct_simple_sincos_loopId(struct cmdline_data* cmd,
                                           struct  global_data* gd,
                                           bodyptr *btable, INTEGER *nbody,
                                           INTEGER ipmin, INTEGER *ipmax,
                                           int cat1, int cat2)
{
    string routineName = "searchcalc_direct_simple_sincos_loopId";
    bodyptr p, q;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    compute_vector dr;
#endif
    int k, n;
    double cpustart;
    real cosphi, sinphi;
    gdhist_sincos_direct hist;

    real xi;
    real s, sy;
    compute_vector pr0;

    cpustart = CPUTIME;
    verb_print(cmd->verbose,
               "Search: Running... (direct-simple-sincos-loopId) \n");

    debug_tracking_s("001",routineName);

    search_init_sincos_gd(cmd, gd, &hist);
    search_init_sincos(cmd, gd, &hist);

    INTEGER ipfalse;
    ipfalse = 0;
    INTEGER iptrue;
    iptrue = 0;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "\n\tTotalID: %ld\n", gd->totalID);

    INTEGER ip;

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nRunning...\n - Completed pivot node:\n");

    debug_tracking("002");

    INTEGER jp=0;
    INTEGER jId=0;

    for (p = btable[cat1]; p < btable[cat1]; p++) {
        jp = p - btable[cat1] + 1;
        if (jp%gd->totalID) jId += 1;

        if (Update(p) == FALSE) {
            ipfalse++;
            continue;
        } else {
            iptrue++;
        }

        for (n = 1; n <= cmd->sizeHistN; n++) {
#ifdef THREEPCFCONVERGENCE
                gd->histNNSub[n] = 0.0;             // Affects only 3pcf
#endif
            hist.histW[n] = 0.0;
            hist.histXi2pcfsub[n] = 0.0;            // Affects only 2pcf
        }
#ifdef THREEPCFCONVERGENCE
            CLRM_ext_ext(gd->histXicos, cmd->mChebyshev+1, cmd->sizeHistN);
            CLRM_ext_ext(gd->histXisin, cmd->mChebyshev+1, cmd->sizeHistN);
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
            hist.drpq = rsqrt(hist.drpq2);
#endif // ! NDIM
#endif

        for (INTEGER i = 1; i <= gd->sizeID; i++) {

            q = btable[cat1] + jId + i - 1;
            if (Id(p) != Id(q)) continue;
            if (Update(q) == FALSE) continue;
            if (p == q) continue;

            if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
                if (cmd->useLogHist) {
                    if(dr1>cmd->rminHist) {
                        if (cmd->rminHist==0)
                            n = (int)(cmd->logHistBinsPD*(rlog10(dr1)
                                - rlog10(cmd->rangeN))
                                + cmd->sizeHistN) + 1;
                        else
                            n = (int)(rlog10(dr1/cmd->rminHist)
                                        * gd->i_deltaR) + 1;
                        if (n<=cmd->sizeHistN && n>=1) {
                            gd->histNN[n] = gd->histNN[n] + 1.;
                            hist.histW[n] = hist.histW[n] + Weight(q);
                            gd->histNNSubXi2pcf[n] = gd->histNNSubXi2pcf[n]
                                                    + 1.;
                            gd->histNNSub[n] = gd->histNNSub[n] + 1.;
                            xi = Weight(q)*Kappa(q);
#ifdef THREEPCFCONVERGENCE
                            REAL cosphi,sinphi;
#if NDIM == 3
                            DOTVP(s, dr, hist.dr0);
                            cosphi = s/(dr1*hist.drpq);
                            CROSSVP(pr0,hist.dr0,Pos(p));
                            DOTVP(sy, dr, pr0);
                            sinphi = rsqrt(1.0 - rsqr(cosphi));
                            if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                            cosphi = -dr[0]/dr1;
                            sinphi = -dr[1]/dr1;
#endif // ! NDIM
                            if (rabs(cosphi)>1.0)
                                verb_log_print(cmd->verbose, gd->outlog,
                                        "Warning!... %s: %g\n",
                                        "cossphi must be in (-1,1)", cosphi);
                                CHEBYSHEVTUSINCOS;
#endif
                            hist.histXi2pcfsub[n] += xi;
                            gd->nbbcalc += 1;
                        } // ! 1 <= n <= sizeHistN
                    } // ! dr1 > rminHist
                } else { // !useLogHist
                    if(dr1>cmd->rminHist) {
                        n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                        if (n<=cmd->sizeHistN && n>=1) {
                            gd->histNN[n] = gd->histNN[n] + 1.;
                            hist.histW[n] = hist.histW[n] + Weight(q);
                            gd->histNNSubXi2pcf[n] = gd->histNNSubXi2pcf[n]
                                                    + 1.;
                            gd->histNNSub[n] = gd->histNNSub[n] + 1.;
                            xi = Weight(q)*Kappa(q);
#ifdef THREEPCFCONVERGENCE
                            real cosphi,sinphi;
#if NDIM == 3
                            DOTVP(s, dr, hist.dr0);
                            cosphi = s/(dr1*hist.drpq);
                            CROSSVP(pr0,hist.dr0,Pos(p));
                            DOTVP(sy, dr, pr0);
                            sinphi = rsqrt(1.0 - rsqr(cosphi));;
                            if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                            cosphi = -dr[0]/dr1;
                            sinphi = -dr[1]/dr1;
#endif // ! NDIM
                            if (rabs(cosphi)>1.0)
                                verb_log_print(cmd->verbose, gd->outlog,
                                    "Warning!... %s: %g\n",
                                    "cossphi must be in (-1,1)", cosphi);
                            CHEBYSHEVTUSINCOS;
#endif
                            hist.histXi2pcfsub[n] += xi;
                            gd->nbbcalc += 1;   // B-B encounter
                        } // 1 < n < sizeHistN
                    } // dr1 > rmin
                } // !useLogHist
            } // accept body
//            } // p != q
        } // end do body q

        computeBodyProperties_sincos_direct(cmd, gd, p, nbody[cat1], &hist);

        ip = p - btable[cat1] + 1;
        if (ip%gd->stepState == 0)
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, ".");
        if (ip%cmd->stepState == 0)
            verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "%d ", ip);
    } // end do body p

    debug_tracking("003");

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n");

    verb_print(cmd->verbose,
                   "\t%s: true and false kappas  (loop-Id only) : %ld %ld\n",
                   routineName, iptrue, ipfalse);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] /= 2.0;
        gd->histNNSubXi2pcf[n] /= 2.0;
        if (cballs_opt_weights_norm(cmd)) {
            gd->histXi2pcf[n] /= MAX(hist.histWW[n],1.0);
        } else {
            gd->histXi2pcf[n] /= MAX(gd->histNNSubXi2pcf[n],1.0);
        }
        gd->histXi2pcf[n] /= gd->totalID;
    }

    if (cballs_opt_compute_histn(cmd))
        search_compute_HistN(cmd, gd, nbody[cat1]);

    search_free_sincos(cmd, gd, &hist);

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    debug_tracking_s("004 ...final", routineName);

    return SUCCESS;
}

//B Several routines like the ones in treeutils module

// Set to zero global histogram arrays that keep the binning information
//  to be save in the output files. Also some counters are set to zero
local int search_init_sincos_gd(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct hist)
{
    int n;
    int m;

#ifdef THREEPCFCONVERGENCE
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNNSubXi2pcf[n] = 0.0;
        gd->histNN[n] = 0.0;
    }
    gd->nbbcalc = 0;                                // B-B encounter

    return SUCCESS;
}

// Allocate local histogram arrays that keep temporal binning information
local int search_init_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct hist)
{
    int n;
    int m;

    hist->histW = dvector(1,cmd->sizeHistN);
    hist->histWW = dvector(1,cmd->sizeHistN);

#ifdef THREEPCFCONVERGENCE
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
#endif

//#ifdef THREEPCFCONVERGENCE
        hist->histXi2pcfsub = dvector(1,cmd->sizeHistN);
//#endif

#ifdef THREEPCFCONVERGENCE
        hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->xiOUTVPcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMtmpcossin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif

    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histW[n] = 0.0;
        hist->histWW[n] = 0.0;
        hist->histXi2pcfsub[n] = 0.0;
    }

    return SUCCESS;
}

// Free allocation memory of local histogram arrays
local int search_free_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos_direct hist)
{
#ifdef THREEPCFCONVERGENCE
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        free_dmatrix(hist->histZetaMtmpcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        free_dmatrix(hist->xiOUTVPcossin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
#endif

    free_dvector(hist->histXi2pcfsub,1,cmd->sizeHistN);

#ifdef THREEPCFCONVERGENCE
        free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
#endif

    free_dvector(hist->histWW,1,cmd->sizeHistN);
    free_dvector(hist->histW,1,cmd->sizeHistN);

    return SUCCESS;
}

// Compute properties and fill temporal and global histogram arrays
//  that keep the binning information to be save in the output files
local int computeBodyProperties_sincos_direct(struct cmdline_data* cmd,
                                       struct  global_data* gd, bodyptr p,
                                       int nbody, gdhistptr_sincos_direct hist)
{
    int n;
    int m;
    real xi, xi_2p;
    real wi;

    xi = Kappa(p)/nbody;
    xi_2p = Weight(p)*Kappa(p);
    wi = Weight(p);

#ifdef THREEPCFCONVERGENCE
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++) {
                gd->histXicos[m][n] /= MAX(gd->histNNSub[n],1.0);
                gd->histXisin[m][n] /= MAX(gd->histNNSub[n],1.0);
            }
        
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos, gd->histXicos[m],
                      gd->histXicos[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin, gd->histXisin[m],
                      gd->histXisin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos, gd->histXisin[m],
                      gd->histXicos[m],cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            OUTVP_ext(hist->xiOUTVPcossin,
                      gd->histXicos[m], gd->histXisin[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMtmpcossin,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,
                      hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,
                      hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            MULMS_ext(hist->histZetaMtmpcossin,
                      hist->xiOUTVPcossin,xi,cmd->sizeHistN);
            ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                     hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                     hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                     hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            ADDM_ext(gd->histZetaMcossin[m],gd->histZetaMcossin[m],
                     hist->histZetaMtmpcossin,cmd->sizeHistN);
        }
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] += xi_2p*hist->histXi2pcfsub[n];
        hist->histWW[n] += wi*hist->histW[n];
    }

    return SUCCESS;
}

//E Several routines like the ones in treeutils module


#undef CHEBYSHEVTUSINCOS
