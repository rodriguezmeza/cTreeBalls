/* ==============================================================================
!	MODULE: direct_simple.c		[cTreeBalls]
!	Written by: M.A. Rodriguez-Meza
!    Starting date:    april 2023
!    Purpose: 3-point correlation function computation
!	Language: C
!	Use: direct_simple();
!	Major revisions:
!==============================================================================*/
//        1          2          3          4          5          6          7


#include "globaldefs.h"

//B Some macros and definitions
#define NOCHEBYSHEVTU                                       \
{real ChebT, ChebU;for (m=1; m<=cmd->mchebyshev+1; m++){    \
    ChebT = rcos((real)(m-1) * racos(cosphi));              \
    ChebU = sinphi * rsin((real)(m-1) * racos(cosphi))/rsin(racos(cosphi)); \
    xicosmphi = xi * ChebT;                                 \
    xisinmphi = xi * ChebU;                                 \
    gd->histXicos[m][n] += xicosmphi;                       \
    gd->histXisin[m][n] += xisinmphi;                       \
    }}

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
    for (m=4; m<=cmd->mchebyshev+1; m++){                   \
        hist.ChebsT[m] = 2.0*(cosphi)*hist.ChebsT[m-1] - hist.ChebsT[m-2]; \
        xicosmphi = xi * hist.ChebsT[m];                    \
        gd->histXicos[m][n] += xicosmphi;                   \
        hist.ChebsU[m] = 2.0*(cosphi)*hist.ChebsU[m-1] - hist.ChebsU[m-2]; \
        xisinmphi = xi * hist.ChebsU[m] * sinphi;           \
        gd->histXisin[m][n] += xisinmphi;                   \
    }}

//B Define a structure:
typedef struct {
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real *ChebsT;
    real *ChebsU;
    real *histXi2pcfsub;

    vector q0;
    real drpq2, drpq;
    vector dr0;
} gdhist_sincos, *gdhistptr_sincos;
//E

//B Tree utilities protodefs
local int search_init(struct cmdline_data* cmd,
                      struct  global_data* gd, gdhistptr hist);
local int search_free(struct cmdline_data* cmd,
                      struct  global_data* gd, gdhistptr hist);
local int computeBodyProperties(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                bodyptr p, int nbody, gdhistptr hist);

local int search_init_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos);
local int search_free_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos);
local int computeBodyProperties_sincos(struct cmdline_data* cmd,
                                       struct  global_data* gd,
                                       bodyptr, int, gdhistptr_sincos);
//E

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
    * Global tructures used: gd, cmd
    * Histograms outputs (in global gd): histZetaMcos, histZetaMsin,
    *                                    histZetaMsincos, histN,
    *                                    histNSubXi2pcf, histNSubXi2pcftotal,
    *                                    histXi2pcf, histXi,
    * Counting encounters (in global gd): nbbcalc, nbccalc, ncccalc
 Return:
    void
 */
global int searchcalc_direct_simple_sincos(struct cmdline_data* cmd, 
                                           struct  global_data* gd,
                                           bodyptr *btable, INTEGER *nbody,
                            INTEGER ipmin, INTEGER *ipmax, int cat1, int cat2)
{
    bodyptr p, q;
    real dr1;
    int k, n;
    vector dr;
    double cpustart;
    real cosphi, sinphi;
    gdhist_sincos hist;

    cpustart = CPUTIME;
    verb_print(cmd->verbose, "Search: Running ... (direct-simple-sincos) \n");

    search_init_sincos(cmd, gd, &hist);

    real xi;
    real s, sy;
    vector pr0;
    for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
        for (n = 1; n <= cmd->sizeHistN; n++) {
            if (cmd->computeTPCF) {
                gd->histNSub[n] = 0.0;              // Affects only 3pcf
            }
            hist.histXi2pcfsub[n] = 0.0;            // Affects only 2pcf
        }
        if (cmd->computeTPCF) {
            CLRM_ext_ext(gd->histXicos, cmd->mchebyshev+1, cmd->sizeHistN);
            CLRM_ext_ext(gd->histXisin, cmd->mchebyshev+1, cmd->sizeHistN);
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
            hist.drpq = rsqrt(hist.drpq2);
#endif // ! NDIM
        }

        for (q = btable[cat2] + ipmin -1; q < btable[cat2] + ipmax[cat2]; q++) {
            if (p != q) {
                if (accept_body(cmd, gd, p, (nodeptr)q, &dr1, dr)) {
                    if (cmd->useLogHist) {
                        if(dr1>cmd->rminHist) {
                            if (cmd->rminHist==0)
                                n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                                          + cmd->sizeHistN) + 1;
                            else
                                n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                            if (n<=cmd->sizeHistN && n>=1) {
                                gd->histN[n] = gd->histN[n] + 1.;
                                gd->histNSubXi2pcf[n] = gd->histNSubXi2pcf[n] + 1.;
                                gd->histNSub[n] = gd->histNSub[n] + 1.;
                                xi = Kappa(q);
                                if (cmd->computeTPCF) {
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
                                                       "Warning!... cossphi must be in (-1,1): %g\n",
                                                       cosphi);
                                    CHEBYSHEVTUSINCOS;
                                }
                                hist.histXi2pcfsub[n] += xi;
                                gd->nbbcalc += 1;
                            }
                        }
                    } else { // !useLogHist
                        if(dr1>cmd->rminHist) {
                            n = (int) ( (dr1-cmd->rminHist) * gd->i_deltaR) + 1;
                            if (n<=cmd->sizeHistN && n>=1) {
                                gd->histN[n] = gd->histN[n] + 1.;
                                gd->histNSubXi2pcf[n] = gd->histNSubXi2pcf[n] + 1.;
                                gd->histNSub[n] = gd->histNSub[n] + 1.;
                                xi = Kappa(q);
                                if (cmd->computeTPCF) {
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
                                                       "Warning!... cossphi must be in (-1,1): %g\n",
                                                       cosphi);
                                    CHEBYSHEVTUSINCOS;
                                }
                                hist.histXi2pcfsub[n] += xi;
                                gd->nbbcalc += 1;   // B-B encounter
                            } // 1 < n < sizeHistN
                        } // dr1 > rmin
                    } // !useLogHist
                } // accept body
            } // p != Q
        } // end do body q

        computeBodyProperties_sincos(cmd, gd, p, nbody[cat1], &hist);
//E
        INTEGER ip;
        ip = p - btable[cat1] + 1;
        if (ip%cmd->stepState == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog,
                           " - Completed pivot: %ld\n", ip);
        }
    } // end do body p

    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] /= 2.0;
        gd->histNSubXi2pcf[n] /= 2.0;
        gd->histXi2pcf[n] /= MAX(gd->histNSubXi2pcf[n],1.0);
    }

    if (scanopt(cmd->options, "compute-HistN"))
        search_compute_HistN(cmd, gd, nbody[cat1]);

    search_free_sincos(cmd, gd, &hist);

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return _SUCCESS_;
}

//B Several routines like the ones in treeutils module
local int search_init(struct cmdline_data* cmd,
                      struct  global_data* gd, gdhistptr hist)
{
    int n;

    if (cmd->computeTPCF) {
        hist->Chebs = dvector(1,cmd->mchebyshev+1);
        hist->xiOUTVP = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    }
    hist->histXi2pcfsub = dvector(1,cmd->sizeHistN);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histXi2pcfsub[n] = 0.0;
    }

    return _SUCCESS_;
}

local int search_free(struct cmdline_data* cmd,
                      struct  global_data* gd, gdhistptr hist)
{
    free_dvector(hist->histXi2pcfsub,1,cmd->sizeHistN);
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVP,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dvector(hist->Chebs,1,cmd->mchebyshev+1);
    }

    return _SUCCESS_;
}

local int computeBodyProperties(struct cmdline_data* cmd,
                                struct  global_data* gd, bodyptr p, int nbody, gdhistptr hist)
{
    int n;
    int m;
    real xi, xi_2p;

    xi = Kappa(p)/nbody;
    xi_2p = Kappa(p);
//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mchebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++)
                gd->histXi[m][n] /= MAX(gd->histNSub[n],1.0);
        for (m=1; m<=cmd->mchebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVP, gd->histXi[m], gd->histXi[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmp,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd->sizeHistN);
            ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],hist->histZetaMtmp,cmd->sizeHistN);
        }
//#endif
    }

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] += xi_2p*hist->histXi2pcfsub[n];
    }

    return _SUCCESS_;
}

local int search_init_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos hist)
{
    int n;
    int m;

    hist->ChebsT = dvector(1,cmd->mchebyshev+1);
    hist->ChebsU = dvector(1,cmd->mchebyshev+1);
    hist->histXi2pcfsub = dvector(1,cmd->sizeHistN);
    hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

    for (m = 1; m <= cmd->mchebyshev+1; m++) {
        CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
        CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histXi2pcfsub[n] = 0.0;
        gd->histN[n] = 0.0;
        for (m = 1; m <= cmd->mchebyshev+1; m++)
            gd->histXi[m][n] = 0.0;
    }
    gd->nbbcalc = 0;                                // B-B encounter

    return _SUCCESS_;
}

local int search_free_sincos(struct cmdline_data* cmd,
                             struct  global_data* gd, gdhistptr_sincos hist)
{
    free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);

    free_dvector(hist->histXi2pcfsub,1,cmd->sizeHistN);

    free_dvector(hist->ChebsU,1,cmd->mchebyshev+1);
    free_dvector(hist->ChebsT,1,cmd->mchebyshev+1);

    return _SUCCESS_;
}

local int computeBodyProperties_sincos(struct cmdline_data* cmd,
                                       struct  global_data* gd, bodyptr p, int nbody,
                                            gdhistptr_sincos hist)
{
    int n;
    int m;
    real xi, xi_2p;

    if (Type(p) == BODY) {
        xi = Kappa(p)/nbody;
        xi_2p = Kappa(p);
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }

    for (m=1; m<=cmd->mchebyshev+1; m++)
        for (n=1; n<=cmd->sizeHistN; n++) {
            gd->histXicos[m][n] /= MAX(gd->histNSub[n],1.0);
            gd->histXisin[m][n] /= MAX(gd->histNSub[n],1.0);
        }

    for (m=1; m<=cmd->mchebyshev+1; m++){
        OUTVP_ext(hist->xiOUTVPcos, gd->histXicos[m], gd->histXicos[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsin, gd->histXisin[m], gd->histXisin[m],cmd->sizeHistN);
        OUTVP_ext(hist->xiOUTVPsincos, gd->histXisin[m],
                  gd->histXicos[m],cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
        CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
        MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
        ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                 hist->histZetaMtmpcos,cmd->sizeHistN);
        ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                 hist->histZetaMtmpsin,cmd->sizeHistN);
        ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                 hist->histZetaMtmpsincos,cmd->sizeHistN);
    }

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->histXi2pcf[n] += xi_2p*hist->histXi2pcfsub[n];
    }

    return _SUCCESS_;
}

//E Several routines like the ones in treeutils module


#undef NOCHEBYSHEVTU
#undef CHEBYSHEVTUSINCOS
