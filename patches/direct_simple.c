/* ==============================================================================
!	MODULE: direct_simple.c		[cTreeBalls]                                    !
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

// This routine computes only 3pcf. It works in 2D; Periodic=0; normal-scale; rminHist=0
// Set properly Makefile_settings
//
// Test: cballs nbody=10000 search=direct-simple-sincos lbox=2000 rangeN=200
//
// search=direct-simple-sincos
global int evalHistograms_direct_simple_sincos(void)
{
    bodyptr p, q;
    real rr, rrRange;
    int k, n;
    vector dr;
    double cpustart;
    real cosphi, sinphi;
    int m;
    real **xiOUTVPcos;
    real **xiOUTVPsin;
    real **xiOUTVPsincos;
    real **histZetaMtmpcos;
    real **histZetaMtmpsin;
    real **histZetaMtmpsincos;
    real *ChebsT;
    real *ChebsU;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "Search: Running ... (direct-simple-sincos) \n");

//B Allocating memory to arrays and setting the initial state to compute nPCF
    ChebsT = dvector(1,cmd.mchebyshev+1);
    ChebsU = dvector(1,cmd.mchebyshev+1);
    xiOUTVPcos = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    xiOUTVPsin = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    xiOUTVPsincos = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    histZetaMtmpcos = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    histZetaMtmpsin = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    histZetaMtmpsincos = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);

    for (m = 1; m <= cmd.mchebyshev+1; m++) {
        CLRM_ext(gd.histZetaMcos[m], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsin[m], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsincos[m], cmd.sizeHistN);
    }
    for (n = 1; n <= cmd.sizeHistN; n++) {
        gd.histN[n] = 0.0;
        for (m = 1; m <= cmd.mchebyshev+1; m++)
            gd.histXi[m][n] = 0.0;
    }
    gd.nbbcalc = 0;                             // B-B encounter
//E
    rrRange = rsqr(cmd.rangeN);

    real xi;
    real ChebT, ChebU;
    real xicosmphi, xisinmphi;

    DO_BODY(p, bodytab, bodytab+cmd.nbody) {

        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histNSub[n] = 0.0;
        }
        CLRM_ext_ext(gd.histXicos, cmd.mchebyshev+1, cmd.sizeHistN);
        CLRM_ext_ext(gd.histXisin, cmd.mchebyshev+1, cmd.sizeHistN);

        DO_BODY(q, bodytab, bodytab+cmd.nbody) {
            if (p != q) {
                DO_COORD(k) {
                    dr[k] = Pos(p)[k] - Pos(q)[k];
                }
                rr=0.0;
                DO_COORD(k) {
                    rr += rsqr(dr[k]);
                }
                if (rr < rrRange) {                 // Find a B, compute nPCF properties
                    n = (int) (rsqrt(rr) * gd.i_deltaR) + 1;
                    gd.histN[n] = gd.histN[n] + 1.;
                    gd.histNSub[n] = gd.histNSub[n] + 1.;
                    xi = Kappa(q);
                    cosphi = (Pos(q)[0] - Pos(p)[0])/rsqrt(rr);     // x, y
                    sinphi = (Pos(q)[1] - Pos(p)[1])/rsqrt(rr);     // x, y

                    if (rabs(cosphi)>1.0) error("\ncossphi must be in (-1,1): %g\n",cosphi);
//
                    CHEBYSHEVTU;
//                    NOCHEBYSHEVTU;
                    gd.nbbcalc += 1;                // B-B encounter
                }
            }
        } // end do body q

//B Computing properties for pivot body p:
        xi = Kappa(p)/cmd.nbody;
        for (m=1; m<=cmd.mchebyshev+1; m++)
            for (n=1; n<=cmd.sizeHistN; n++) {
                gd.histXicos[m][n] /= MAX(gd.histNSub[n],1.0);
                gd.histXisin[m][n] /= MAX(gd.histNSub[n],1.0);
            }

        for (m=1; m<=cmd.mchebyshev+1; m++){
            OUTVP_ext(xiOUTVPcos, gd.histXicos[m], gd.histXicos[m],cmd.sizeHistN);
            OUTVP_ext(xiOUTVPsin, gd.histXisin[m], gd.histXisin[m],cmd.sizeHistN);
            OUTVP_ext(xiOUTVPsincos, gd.histXisin[m], gd.histXicos[m],cmd.sizeHistN);
            CLRM_ext(histZetaMtmpcos,cmd.sizeHistN);
            CLRM_ext(histZetaMtmpsin,cmd.sizeHistN);
            CLRM_ext(histZetaMtmpsincos,cmd.sizeHistN);
            MULMS_ext(histZetaMtmpcos,xiOUTVPcos,xi,cmd.sizeHistN);
            MULMS_ext(histZetaMtmpsin,xiOUTVPsin,xi,cmd.sizeHistN);
            MULMS_ext(histZetaMtmpsincos,xiOUTVPsincos,xi,cmd.sizeHistN);
            ADDM_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],histZetaMtmpcos,cmd.sizeHistN);
            ADDM_ext(gd.histZetaMsin[m],gd.histZetaMsin[m],histZetaMtmpsin,cmd.sizeHistN);
            ADDM_ext(gd.histZetaMsincos[m],gd.histZetaMsincos[m],histZetaMtmpsincos,cmd.sizeHistN);
        }
//E
    } // end do body p

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN(cmd.nbody);

//B Freeing allocated memory
    free_dmatrix(histZetaMtmpsincos,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(histZetaMtmpsin,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(histZetaMtmpcos,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(xiOUTVPsincos,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(xiOUTVPsin,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(xiOUTVPcos,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dvector(ChebsU,1,cmd.mchebyshev+1);
    free_dvector(ChebsT,1,cmd.mchebyshev+1);
//E

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return _SUCCESS_;
}

#endif    // ! _direct_simple_c
