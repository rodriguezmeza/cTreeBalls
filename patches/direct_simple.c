/* ==============================================================================
!	MODULE: direct_simple.c		[tpcf]                                          !
!	Written by: M.A. Rodriguez-Meza.											!
!    Starting date:    april 2023                                               !
!    Purpose: 3-point correlation function computation                          !
!	Language: C																	!
!	Use: direct_simple();													    !
!	Major revisions:															!
!==============================================================================*/


//B Set of search using brute force:
//
// search=direct-simple :: evalHistograms_direct_simple()
//
//E

#ifndef _direct_simple_c
#define _direct_simple_c

#include "globaldefs.h"

// search=direct-simple
global int searchcalc_direct_simple(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    bodyptr p, q;
    int n;
    vector dr;
    double cpustart;
    real xi;
    real cosphi;
    real xicosmphi;
    int m;
    real dr1;
    gdhist hist;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "Search: Running ... (direct-simple) \n");

    search_init(&hist);
    search_init_gd_hist();

    DO_BODY(p, btab+ipmin-1, btab+ipmax) {
//B Set histograms to zero for the pivot
        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histNSub[n] = 0.0;               // Affect only 3pcf evaluation
            hist.histXi2pcfsub[n] = 0.0;
        }
#ifdef TPCF
        CLRM_ext_ext(gd.histXi, cmd.mchebyshev+1, cmd.sizeHistN);
#endif
//E

#ifdef TPCF
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
            hist.drpq = rsqrt(hist.drpq2);
#endif
#endif
        DO_BODY(q, btab+ipmin-1, btab+ipmax) {
            if (p != q) {
                if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
//                    if(dr1>0) {
                    if(dr1>cmd.rminHist) {
                        if (cmd.rminHist==0)
                            n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                        else
                            n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;

                        if (n<=cmd.sizeHistN && n>=1) {
                            gd.histN[n] = gd.histN[n] + 1.;
                            gd.histNSub[n] = gd.histNSub[n] + 1.;
                            xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                            real s;
                            DOTVP(s, dr, hist.dr0);
                            cosphi = -s/(dr1*hist.drpq);
#else
                            cosphi = -dr[1]/dr1;
#endif
                            if (rabs(cosphi)>1.0)
                                verb_log_print(cmd.verbose, gd.outlog,
                                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                                    cosphi);
//                            CHEBYSHEVNEW
                            real Cheb;
                            for (m=1; m<=cmd.mchebyshev+1; m++){
                                    Cheb = rcos((real)(m-1) * racos(cosphi));
                                    xicosmphi = xi * Cheb;
                                    gd.histXi[m][n] += xicosmphi;
                                }
#endif // ! TPCF
                    hist.histXi2pcfsub[n] += xi;
                    gd.nbbcalc += 1;                // B-B encounter
                        }
                    }
#else // ! LOGHIST
                    if(dr1>cmd.rminHist) {
                    n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                        if (n<=cmd.sizeHistN && n>=1) {
                    gd.histN[n] = gd.histN[n] + 1.;
                    gd.histNSub[n] = gd.histNSub[n] + 1.;
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                        real s;
                        DOTVP(s, dr, hist.dr0);
                        cosphi = -s/(dr1*hist.drpq);
#else
                        cosphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
//                    CHEBYSHEVNEW
                    real Cheb;
                    for (m=1; m<=cmd.mchebyshev+1; m++){
                            Cheb = rcos((real)(m-1) * racos(cosphi));
                            xicosmphi = xi * Cheb;
                            gd.histXi[m][n] += xicosmphi;
                        }

#endif // ! TPCF
                    hist.histXi2pcfsub[n] += xi;
                    gd.nbbcalc += 1;                // B-B encounter
                        }
                    }
#endif // ! LOGHIST
                }
            }
        } // end do body q

        computeBodyProperties(p, ipmax-ipmin+1, &hist);
    } // end do body p

    for (n = 1; n <= cmd.sizeHistN; n++) {
        gd.histXi2pcf[n] /= 2.0;
        gd.histNSub[n] = gd.histN[n];
        gd.histNSub[n] /= 2.0;
        gd.histXi2pcf[n] /= MAX(gd.histNSub[n],1.0);
    }

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN(ipmax-ipmin+1);

    search_free(&hist);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return _SUCCESS_;
}

#endif    // ! _direct_simple_c
