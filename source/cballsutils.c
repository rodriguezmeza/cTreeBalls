/* ==============================================================================
 MODULE: cballsutils.c			    [cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"


global int search_init_sincos_omp(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
    int n;
    int m;

#ifdef TPCF
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
#endif
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
//B 2pcf
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
#endif
//E
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
#ifdef TPCF
        hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        
        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mChebyshev+1,
                                             1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,
                      1,cmd->sizeHistN,1,cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        hist->histZetaMthreadcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);

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
       hist->histNthread[n] = 0.0;
       hist->histNNSubthread[n] = 0.0;
       hist->histNNSubXi2pcfthread[n] = 0.0;
#ifdef SMOOTHPIVOT
       hist->histNNSubXi2pcfthreadp[n] = 0.0;
       hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
#endif
       hist->histXi2pcfthread[n] = 0.0;
       hist->histXi2pcfthreadsub[n] = 0.0;
   }

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        }
#endif

    return SUCCESS;
}

global int search_free_sincos_omp(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
#ifdef TPCF
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
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        free_dmatrix3D(hist->histZetaMthreadcossin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsincos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
#ifdef SMOOTHPIVOT
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
#endif
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
#ifdef TPCF
        free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
#endif

    return SUCCESS;
}

global int computeBodyProperties_sincos(struct  cmdline_data* cmd,
                                            struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp hist)
{
    int n;
    int m;
    real xi, xi_2p;

    //B Normalization of histograms
    if (Type(p) == BODY) {
#ifdef NOSTANDARNORMHIST
        xi = Kappa(p);
        xi_2p = Kappa(p);
#ifdef SMOOTHPIVOT
            xi_2p = KappaRmin(p);
            xi = NbRmin(p)*xi_2p;
#endif
#else // ! NOSTANDARNORMHIST
        xi = Kappa(p)/nbody;
#ifdef BALLS4SCANLEV
        xi_2p = (Weight(p)/Nb(p))*Kappa(p);
#else
        xi_2p = Weight(p)*Kappa(p);
#endif
#ifdef SMOOTHPIVOT
#ifdef BALLS4SCANLEV
            xi_2p = KappaRmin(p);
#endif
            xi = NbRmin(p)*xi_2p/nbody;
#endif
#endif // ! NOSTANDARNORMHIST
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
    //E Normalization of histograms

#ifdef TPCF
        for (m=1; m<=cmd->mChebyshev+1; m++)
            //B Normalization of histograms
#ifdef NOSTANDARNORMHIST
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= 1.0;
                hist->histXithreadsin[m][n] /= 1.0;
            }
#else
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            }
#endif
            //E
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos,
                      hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin,
                      hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos,
                      hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            OUTVP_ext(hist->xiOUTVPcossin,
                      hist->histXithreadcos[m], hist->histXithreadsin[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMtmpcossin,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,
                      hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            MULMS_ext(hist->histZetaMtmpcossin,
                      hist->xiOUTVPcossin,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],
                     hist->histZetaMthreadcos[m],
                     hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],
                     hist->histZetaMthreadsin[m],
                     hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],
                     hist->histZetaMthreadsincos[m],
                     hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            ADDM_ext(hist->histZetaMthreadcossin[m],
                     hist->histZetaMthreadcossin[m],
                     hist->histZetaMtmpcossin,cmd->sizeHistN);
        }
#endif

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
    }

    return SUCCESS;
}

global int search_init_gd_hist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;
    int n1, n2, l;

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaM[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
    }
    
#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaGmRe[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaGmIm[m], cmd->sizeHistN);
        }
#endif

    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

global int search_init_gd_hist_sincos(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

#ifdef TPCF
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        }
#endif
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
#ifdef SMOOTHPIVOT
        gd->histNNSubXi2pcftotal[n] = 0.0;
#endif
        gd->histXi2pcf[n] = 0.0;
#ifdef TPCF
            for (m = 1; m <= cmd->mChebyshev+1; m++) {
                // HERE MUST BE gd->histXicos and gd->histXisin
            }
#endif
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

//B Computation of histogram of all B-B encounters
// The correlation function is estimated as:
//    xi=(V/v(r))*(DD(r)/N^2)
// where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
// total # particles.
//
// Note: only rminHistN = 0 works and agree with CUTE_BOX
//  you may try options=cute-box-rmin to correct a bit the results...
//      but for biger values of rminHistN differences grow...
//
local int search_compute_Xi(struct  cmdline_data* cmd,
                            struct  global_data* gd, int nbody)
{
    int k;
    int n;
    real normFac;
    real Vol;
    //B correct cute-box-rmin
    real deltaR;
    if ((scanopt(cmd->options, "cute-box-rmin")))
        deltaR = cmd->rangeN/cmd->sizeHistN;
    //E

    Vol = 1.0;
    DO_COORD(k)
        Vol = Vol*gd->Box[k];

if (!cmd->useLogHist) {
    if ((scanopt(cmd->options, "cute-box"))) {
        gd->histNN[1]-=nbody;
    }
}
    real *edd;
    real *corr;
    real *ercorr;
    edd = dvector(1,cmd->sizeHistN);
    corr = dvector(1,cmd->sizeHistN);
    ercorr = dvector(1,cmd->sizeHistN);
    real rho_av=(real)nbody/Vol;

    for (n = 1; n <= cmd->sizeHistN; n++)
        edd[n] = 1./rsqrt(gd->histNN[n]);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        if(gd->histNN[n]==0) {
            corr[n]=0;
            ercorr[n]=0;
        } else {
            double r0,r1,vr,rho_r;
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    r0 = rpow(10.0, ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                    r1 = rpow(10.0, ((real)(n+1-cmd->sizeHistN))/cmd->logHistBinsPD
                              + rlog10(cmd->rangeN) );
                } else {
                    r0 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n))*gd->deltaR );
                    r1 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n+1))*gd->deltaR );
                }
            } else {
                //B correct cute-box-rmin
                if ((scanopt(cmd->options, "cute-box-rmin"))) {
                    r0=(real)n*deltaR;
                    r1=(real)(n+1)*deltaR;
                } else {
                    r0=(real)n*gd->deltaR;
                    r1=(real)(n+1)*gd->deltaR;
                }
            }

#if (NDIM==3)
            if (scanopt(cmd->options, "cute-box")) {
                //B this version does not give same results as CB
                //      although the programming is the same...
                vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
                rho_r=gd->histNN[n]/((real)nbody*vr);
                corr[n]=rho_r/rho_av-1;             // Correlation function
                ercorr[n]=(1+corr[n])*edd[n];       // Poisson errors
                gd->histCF[n] = corr[n];            // Original line
                //E
            } else {
                if (cmd->useLogHist) {
                    vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
// rho_r/rho_av = ( histNN[n]/(nbody*vr) ) / (nbody/Vol)
                    normFac = Vol/(vr*((real)(nbody*nbody)));
                    gd->histCF[n] = gd->histNN[n] * normFac - 1.0;
                } else {
                    //B correct cute-box-rmin
                    if ((scanopt(cmd->options, "cute-box-rmin"))) {
                        normFac = Vol/(2.0*PI*rpow(deltaR,3.0)*nbody*nbody);
                    } else {
                        normFac = Vol/(2.0*PI*rpow(gd->deltaR,3.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5);
// This line gives results in agreement with CB:
                    }
                    gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5) -1.0;
                    //E
                }
            }
#else
            if (scanopt(cmd->options, "cute-box")) {
                // This should be CB version...
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            } else {
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5);
// This line gives results in agreement with CB:
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            }
#endif // ! NDIM
        }
    }

    free_dvector(ercorr,1,cmd->sizeHistN);
    free_dvector(corr,1,cmd->sizeHistN);
    free_dvector(edd,1,cmd->sizeHistN);

    return SUCCESS;
}


global int search_compute_HistN(struct  cmdline_data* cmd, 
                                struct  global_data* gd, int nbody)
{
    int n;
    real normFac;

//B Check this factor is correct...
// to agree with cute_box normalization commented out these lines
//B these does not work!!
//    normFac = 1.0;
//E
    normFac = 0.5;
    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNN[n] *= normFac;
//E
    if (scanopt(cmd->options, "and-CF"))
        search_compute_Xi(cmd, gd, nbody);

    return SUCCESS;
}


// No se usa qsize: update
global bool reject_cell(struct  cmdline_data* cmd, struct
                        global_data* gd, nodeptr p, nodeptr q, real qsize)
{
    real drpq, drpq2;
    vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}


//B 2023.11.22
global bool reject_balls(struct  cmdline_data* cmd,
                         struct  global_data* gd, nodeptr p, nodeptr q,
                         real *drpq, vector dr)
{
    real drpq2;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    *drpq = rsqrt(drpq2);

    if ( *drpq >= gd->Rcut + Radius(p) + Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}
//E

global bool reject_cell_balls(struct  cmdline_data* cmd,
                              struct  global_data* gd, nodeptr p, nodeptr q,
                              real *drpq, vector dr)
{
    real drpq2;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    *drpq = rsqrt(drpq2);

    if ( *drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

global bool reject_bodycell(struct  cmdline_data* cmd,
                            struct  global_data* gd, nodeptr p, nodeptr q)
{
    real drpq, drpq2;
    vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

global bool reject_cellcell(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q)
{
    real drpq, drpq2;
    vector dr;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    drpq = rsqrt(drpq2);

    if ( drpq >= gd->Rcut+Radius(p)+Radius(q) )
        return (TRUE);
    else
        return (FALSE);
}

#ifdef SINGLEP
global bool accept_body(struct  cmdline_data* cmd, struct  global_data* gd,
                        bodyptr p, nodeptr q, float *drpq, float *dr)
#else
global bool accept_body(struct  cmdline_data* cmd, struct  global_data* gd,
                        bodyptr p, nodeptr q, real *drpq, vector dr)
#endif
{
    real drpq2;

    DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
    }
    *drpq = rsqrt(drpq2);

    if ( *drpq<gd->Rcut )
        return (TRUE);
    else
        return (FALSE);
}


//B Calculate number of threads
// Use in the terminal:
// export OMP_NUM_THREADS=8
// to set the maximum number of threads that can be used in a run
global int ThreadCount(struct  cmdline_data* cmd, struct  global_data* gd,
                       INTEGER nbody, int cat1) {
    int nthreads=0;
#pragma omp parallel
    {
#pragma omp atomic
        nthreads++;
    }
    verb_print(cmd->verbose, "\nUsing %d threads \n",nthreads);

//B socket:
#ifdef ADDONS
#include "cballsutils_include_01.h"
#endif
//E

    return SUCCESS;
}
//E


//B coordinate transformations routines

//
// All angles and input distances are in radians.
// Will be useful to convert to degrees, hours, arcmin and arcsec.
// And write this info in log file.
//
// deg = rads * 180/PI // rads to degrees.
// deg = rads * 180/(PI*15.0) // rads to hours.
// deg = (rads * 180/PI) * 60 // rads to arcmin.
// deg = (rads * 180/PI) * 60 * 60 // rads to arcsec.
//
// Also consider wrapping angles to the range (-pi, pi]:
// if ( angle <= -PI ) angle += 2.0*PI
// if ( angle >   PI ) angle -= 2.0*PI
//
// Transform to x,y,z coordinates:
// ra = rhs.getX();
// dec = rhs.getY();
//    x = cosdec * cosra;
//    y = cosdec * sinra;
//    z = sindec;
//

#define ARFKEN              0
#define NOARFKEN            1
#define GALACTIC            2
#define ECLIPTIC            3
#define CELESTIAL           4                   // also known as equatorial

//local int coordinate_string_to_int(string, int *);
//local int transfmt_int;
local int spherical_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                  real theta, real phi, vector xyz);
local int galactic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                 real latitud, real longitud, vector xyz);
local int ecliptic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                 real DEC, real RA, vector xyz);
local int celestial_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                  real DEC, real RA, vector xyz);

global int coordinate_string_to_int(struct cmdline_data* cmd,
                                    struct  global_data* gd)
{
    gd->coordTag = -1;
    gd->coordTag = ARFKEN;
    if (scanopt(cmd->options, "arfken")) gd->coordTag = ARFKEN;
    if (scanopt(cmd->options, "no-arfken")) gd->coordTag = NOARFKEN;
    if (scanopt(cmd->options, "galactic")) gd->coordTag = GALACTIC;
    if (scanopt(cmd->options, "ecliptic")) gd->coordTag = ECLIPTIC;
    if (scanopt(cmd->options, "celestial")) gd->coordTag = CELESTIAL;

    return SUCCESS;
}

global int coordinate_transformation(struct cmdline_data* cmd, struct  global_data* gd,
                                    real theta, real phi, vector xyz)
{
    switch(gd->coordTag) {
        case ARFKEN:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case NOARFKEN:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case GALACTIC:
            galactic_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case ECLIPTIC:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
        case CELESTIAL:
            celestial_to_cartesians(cmd, gd, theta, phi, xyz); break;
        default:
            spherical_to_cartesians(cmd, gd, theta, phi, xyz); break;
    }

    return SUCCESS;
}

local int spherical_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real theta, real phi, vector xyz)
{
    real ra, dec;

    if (scanopt(cmd->options, "no-arfken")) {
        //B this would be a galactic->cartesian tranformation...
        //      if RA is theta and DEC is phi
        //      more preciselly: RA is the longitud (l)
        //                      DEC is the latitud (b)
        ra = theta;
        dec = PIO2 - phi;
        xyz[0] = rcos(dec)*rcos(ra);
        xyz[1] = rcos(dec)*rsin(ra);
        xyz[2] = rsin(dec);
        //E
    } else {
        ra = phi;
        dec = theta;
// Standard transformation. See Arfken
        xyz[0] = rsin(dec)*rcos(ra);
        xyz[1] = rsin(dec)*rsin(ra);
        xyz[2] = rcos(dec);
    }

    return SUCCESS;
}

// latitud (b): the angle of an object northward of the galactic equator
// longitud (l): the angular distance of an object eastward along the galactic equator
local int galactic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real latitud, real longitud, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\ngalactic_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    phi = longitud;
    theta = latitud;
    //B Standard transformation. See Arfken
    xyz[0] = rsin(theta)*rcos(phi);
    xyz[1] = rsin(theta)*rsin(phi);
    xyz[2] = rcos(theta);
    //E

    return SUCCESS;
}

local int ecliptic_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real DEC, real RA, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\necliptic_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    if (scanopt(cmd->options, "ra-reversed")) {
        phi = TWOPI - RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    } else {
        phi = RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    }
    
    return SUCCESS;
}

local int celestial_to_cartesians(struct  cmdline_data* cmd,
                                   struct  global_data* gd,
                                   real DEC, real RA, vector xyz)
{
    real phi, theta;

    verb_print(cmd->verbose, "\ncelestial_to_cartesians: coordTag: %d\n",
               gd->coordTag);

    if (scanopt(cmd->options, "ra-reversed")) {
        phi = TWOPI - RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    } else {
        phi = RA;
        theta = PIO2 - DEC;
        // Standard transformation. See Arfken
        xyz[0] = rsin(theta)*rcos(phi);
        xyz[1] = rsin(theta)*rsin(phi);
        xyz[2] = rcos(theta);
    }
    
    return SUCCESS;
}

global int spherical_periodic_condition(real *thetaL, real *thetaR,
                                        real *phiL, real *phiR)
{
    while(*thetaL<0) {
        *thetaL += 2.0*PI;
    }
    while(*thetaL>2.0*PI) {
        *thetaL -= 2.0*PI;
    }
    while(*thetaR<0) {
        *thetaR += 2.0*PI;
    }
    while(*thetaR>2.0*PI) {
        *thetaR -= 2.0*PI;
    }

    while(*phiL<0) {
        *phiL += PI;
    }
    while(*phiL>PI) {
        *phiL -= PI;
    }
    while(*phiR<0) {
        *phiR += PI;
    }
    while(*phiR>PI) {
        *phiR -= PI;
    }

    return SUCCESS;
}

#undef ARFKEN
#undef NOARFKEN
#undef GALACTIC
#undef ECLIPTIC
#undef CELESTIAL                                // also known as equatorial
#undef NULLTRANSFORM

//E coordinate transformations routines


//B section of several routines to do pre/post processing

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETASTDDEV \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

#define MHISTZETAHEADER \
"# [1] rBins; [2] diagonal; [3] theta2=Nbins/4.0; [4] theta2=2.0*Nbins/4.0; \
[5] theta2=3.0*Nbins/4.0; [6] theta2=4.0*Nbins/4.0 - 1.0\n"

#define MHISTZETAHEADERSTDDEV \
"# [1] rBins; [2] diagonal &SD; [3] theta2=Nbins/4.0 &SD; [4] theta2=2.0*Nbins/4.0 &SD; \
[5] theta2=3.0*Nbins/4.0 &SD; [6] theta2=4.0*Nbins/4.0 - 1.0 &SD\n"

// routine to compute covariance matriz of correlations of Takahashi simulations
//  test were done on:
//  run/Cosma/Takahasi_nres12_balls-omp/zs9_balls-omp/rxxx/Output files
global int statHistogram(struct cmdline_data* cmd, struct  global_data* gd)
{
    string routine_name = "statHistogram";
    char namebuf1[256];
    char namebuf2[256];
    char namebuf3[256];
    char namebuf4[256];
    char namebuf5[256];
    char namebuf6[256];
    struct stat buf;
    stream instr1;
    stream instr2;
    stream instr3;
    stream outstr1;
    stream outstr2;
    char rootDirPath[MAXLENGTHOFFILES];
    int nrealization = 0;
    int status1 = 1;
    int status2 = 1;
    int status3 = 1;
    int ifile = 0;
    int i;
    int j;
    int m;
    int n1;
    int n2;
    int npts1;
    int npts2;
    int npts3;
    real matElement;
    real matElement2;
    real matElementIm;
    real **mat1;
    real **mat2;
    real **mattmp;
    real ***mat3;
    real ***matAvg;
    real ***matStdDev;
    real ***matCovMat;
    int sizeHistN;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    real ZetaStdDev;
    real Zeta2StdDev;
    real Zeta3StdDev;
    real Zeta4StdDev;
    real Zeta5StdDev;
    int Nbins;
    real *rBin;
    int NR;

    sizeHistN = cmd->sizeHistN;
    mat1 = dmatrix(1,sizeHistN,1,sizeHistN);
    mat2 = dmatrix(1,sizeHistN,1,sizeHistN);

//B read one file to got needed info to go...
    verb_print(cmd->verbose,
               "\n%s: first read two file to get some needed info:",
               routine_name);
    //B we use rootDir string to construct input/output file names
    sprintf(rootDirPath,"%s/%s%03d/%s/",
            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
    //E
    m=1;
    //B file 1
    sprintf(namebuf1,"%s%s_%d%s",
            rootDirPath, gd->infilenames[ifile], m, EXTFILES);
    verb_print(cmd->verbose,
                "\n%s: opening file %s...",
               routine_name,namebuf1);
    if (stat(namebuf1, &buf) == 0)               // no input file exists?
        instr1 = stropen(namebuf1, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name,namebuf1);
        status1 = 0;
    }
    //E file 1
    //B file 2
    sprintf(namebuf2,"%s%s_%d%s",
            rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
    verb_print(cmd->verbose,
               "\n%s: opening file %s...",
               routine_name,namebuf2);
    if (stat(namebuf2, &buf) == 0)               // no input file exists?
        instr2 = stropen(namebuf2, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name,namebuf2);
        status2 = 0;
    }
    //E file 2

    if (status1 == 0 || status2 == 0) {
        if (status1 != 0) fclose(instr1);
        if (status2 != 0) fclose(instr2);
        error("\n%s: one of the input file does not exist: %s %s",
              routine_name, namebuf1, namebuf2);
    } else {
        //B processing input files
        inout_InputDataMatrix(namebuf1, mat1, &npts1,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        inout_InputDataMatrix(namebuf2, mat2, &npts2,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts1 != npts2) {
            error("\n%s: in realization r%03d: %s %d %d",
                  routine_name, nrealization, "npts are different:", npts1,npts2);
        }
        cmd->sizeHistN = npts1;
        fclose(instr1);
        fclose(instr2);
        //E
    } // ! status
    verb_print(cmd->verbose,
               "\n%s: mChebyshev and sizeHistN: %d %d\n",
               routine_name, cmd->mChebyshev, cmd->sizeHistN);
    free_dmatrix(mat2,1,sizeHistN,1,sizeHistN);
    free_dmatrix(mat1,1,sizeHistN,1,sizeHistN);
//E read one file

//B read rBin file
    rBin = dvector(1,sizeHistN);
    //B file 3
    verb_print(cmd->verbose,
               "\n%s: reading rbins...",
               routine_name, cmd->mChebyshev, cmd->sizeHistN);
    sprintf(namebuf5,"%s%s%s",
            rootDirPath, "rbins", EXTFILES);
    verb_print(cmd->verbose,
               "\n%s: opening file %s...",routine_name, namebuf5);
    if (stat(namebuf5, &buf) == 0)               // no input file exists?
        instr3 = stropen(namebuf5, "r");
    else {
        verb_print(cmd->verbose,
                   "\n%s: Input file does not exist: %s\n",
                   routine_name, namebuf5);
        status3 = 0;
    }
    //E file 3
    if (status3 == 0) {
        error("\n%s: the input file does not exist: %s",
              routine_name, namebuf5);
    } else {
        //B processing input file
        inout_InputDataVector(namebuf5, rBin, &npts3,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts3 != cmd->sizeHistN) {
            error("\n%s: in realization r%03d: %s %d %d",
                  routine_name, nrealization, "npts is not equal to sizeHistN:",
                  npts3,cmd->sizeHistN);
        }
        fclose(instr3);
        //E
    } // ! status3
    verb_print(cmd->verbose, "\ndone.\n");
//B read rBin file

    Nbins = cmd->sizeHistN;

    mat1 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat2 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat3 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    matAvg =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    matStdDev =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);


    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(matAvg[m], cmd->sizeHistN);
        CLRM_ext(matStdDev[m], cmd->sizeHistN);
    }

//B read realization files and compute mean values
    verb_print(cmd->verbose,
               "\n%s: mean value computation:", routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/%s%03d/%s/",
                cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            sprintf(namebuf1,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) == 0)               // no input file exists?
                instr1 = stropen(namebuf1, "r");
            else {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",
                       routine_name, namebuf2);
            if (stat(namebuf2, &buf) == 0)               // no input file exists?
                instr2 = stropen(namebuf2, "r");
            else {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
                if (status1 != 0) fclose(instr1);
                if (status2 != 0) fclose(instr2);
                break;
            } else {
                //B output file 1
                sprintf(namebuf3, "%s%s_%d%s",
                        rootDirPath, cmd->outfile, m, EXTFILES);
                verb_print(cmd->verbose,
                           "\n%s: opening file %s... to save statistics",
                           routine_name, namebuf3);
                outstr1 = stropen(namebuf3, "w!");
                //E output file 1
                //B output file 2
                sprintf(namebuf4, "%s%s%s_%d%s",
                        rootDirPath, "m", cmd->outfile, m, EXTFILES);
                verb_print(cmd->verbose,
                           "\n%s: opening file %s... to save statistics",
                           routine_name, namebuf4);
                outstr2 = stropen(namebuf4, "w!");
                //E output file 2

                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                error("\n%s: in realization r%03d: %s %d %d",
                      routine_name, nrealization,
                      "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        matElement = mat1[n1][n2]+mat2[n1][n2];
                        if (scanopt(cmd->options, "same-infiles"))
                            matElement *= 0.5;
                        mat3[m][n1][n2] = matElement;
                        matAvg[m][n1][n2] += matElement;
                        fprintf(outstr1,"%16.8e ",matElement);
                    }
                    fprintf(outstr1,"\n");
                }

                fprintf(outstr2,MHISTZETAHEADER);
                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    Zeta = mat3[m][n1][n1];
                    Zeta2 = mat3[m][n1][(int)(Nbins/4.0)];
                    Zeta3 = mat3[m][n1][(int)(2.0*Nbins/4.0)];
                    Zeta4 = mat3[m][n1][(int)(3.0*Nbins/4.0)];
                    Zeta5 = mat3[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                    fprintf(outstr2,MHISTZETA,
                            rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
                }

                fclose(instr1);
                fclose(instr2);
                fclose(outstr1);
                fclose(outstr2);
                //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            sprintf(namebuf3, "%s%s_%s_%d%s",
                    rootDirPath, cmd->outfile, "Avg", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            outstr1 = stropen(namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    matAvg[m][n1][n2] /= (nrealization-1);
                    fprintf(outstr1,"%16.8e ",matAvg[m][n1][n2]);
                }
                fprintf(outstr1,"\n");
            }
            fclose(outstr1);

            //B output file 2
            sprintf(namebuf6, "%s%s%s_%s_%d%s",
                    rootDirPath, "m", cmd->outfile, "Avg", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf6);
            outstr2 = stropen(namebuf6, "w!");
            //E output file 2
            fprintf(outstr2,MHISTZETAHEADER);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                Zeta = matAvg[m][n1][n1];
                Zeta2 = matAvg[m][n1][(int)(Nbins/4.0)];
                Zeta3 = matAvg[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = matAvg[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = matAvg[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr2,MHISTZETA,rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            fclose(outstr2);

        } // ! end m loop
    } // ! nrealization > 3
//E read realization files and compute mean values


//
//B read realization files and compute std values
//
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
               "\n\n%s: standard deviation computation:", routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/%s%03d/%s/",
                cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            sprintf(namebuf1,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                           routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                error("\n%s: in realization r%03d: %s %d %d",
                      routine_name, nrealization,
                      "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (scanopt(cmd->options, "same-infiles"))
                            matElement = rabs(
                                           0.5*(mat1[n1][n2]+mat2[n1][n2])
                                           -matAvg[m][n1][n2]
                                           );
                        else
                            matElement = rabs(
                                           mat1[n1][n2]+mat2[n1][n2]
                                           -matAvg[m][n1][n2]
                                           );

                        matStdDev[m][n1][n2] += matElement;
                    }
                }
                //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            sprintf(namebuf3, "%s%s_%s_%d%s",
                    rootDirPath, cmd->outfile, "StdDev", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            outstr1 = stropen(namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    matStdDev[m][n1][n2] /= (nrealization-1);
                    fprintf(outstr1,"%16.8e ",matStdDev[m][n1][n2]);
                }
                fprintf(outstr1,"\n");
            }
            fclose(outstr1);

            //B output file 2
            sprintf(namebuf6, "%s%s%s_%s_%d%s",
                    rootDirPath, "m", cmd->outfile, "StdDev", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf6);
            outstr2 = stropen(namebuf6, "w!");
            //E output file 2
            fprintf(outstr2,MHISTZETAHEADERSTDDEV);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                Zeta = matAvg[m][n1][n1];
                ZetaStdDev = matStdDev[m][n1][n1];
                Zeta2 = matAvg[m][n1][(int)(Nbins/4.0)];
                Zeta2StdDev = matStdDev[m][n1][(int)(Nbins/4.0)];
                Zeta3 = matAvg[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta3StdDev = matStdDev[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = matAvg[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta4StdDev = matStdDev[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = matAvg[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                Zeta5StdDev = matStdDev[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr2,MHISTZETASTDDEV,
                        rBin[n1],Zeta,ZetaStdDev,
                        Zeta2,Zeta2StdDev,
                        Zeta3,Zeta3StdDev,
                        Zeta4,Zeta4StdDev,
                        Zeta5,Zeta5StdDev
                        );
            }
            fclose(outstr2);

        } // ! end m loop
    } // ! nrealization > 3

//
//E read realization files and compute std values
//

//
//B read realization files and compute covariance matrices values
//
    Nbins = cmd->sizeHistN*cmd->sizeHistN;

    int Nv;
    int N=cmd->sizeHistN;
    real **Zv;
    real **ZvAvg;
    Nv = cmd->sizeHistN*cmd->sizeHistN;
    Zv = dmatrix(1,cmd->mChebyshev+1,1,Nv);
    ZvAvg = dmatrix(1,cmd->mChebyshev+1,1,Nv);

    mattmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    matCovMat =
            dmatrix3D(1,cmd->mChebyshev+1,1,Nv,1,Nv);

    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    CLRM_ext(mattmp, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(matCovMat[m], Nv);
    }
    CLRM_ext_ext(Zv, cmd->mChebyshev+1, Nv);
    CLRM_ext_ext(ZvAvg, cmd->mChebyshev+1, Nv);

    int n3, n4;
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
               "\n\n%s: covariance matrices computation:", routine_name);
    verb_print(cmd->verbose,
               "\n%s: reading files and computing covariance matrices...",
               routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/%s%03d/%s/",
                cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            sprintf(namebuf1,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                            "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2

            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    error("\n%s: in realization r%03d: %s %d %d",
                          routine_name, nrealization,
                          "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (scanopt(cmd->options, "same-infiles"))
                            Zv[m][N*(n1-1)+n2] += 0.5*(mat1[n1][n2]+mat2[n1][n2]);
                        else
                            Zv[m][N*(n1-1)+n2] += mat1[n1][n2]+mat2[n1][n2];
                    }
                }
                    //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    NR = nrealization-1;

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed: %d\n",
               routine_name, NR);
    
    verb_print(cmd->verbose,
               "\n%s: computing mean vectors...\n",
               routine_name);
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (i=1; i<=Nv; i++) {
            ZvAvg[m][i] = Zv[m][i]/((real)(NR-1));          // sample mean
        }
    }


    //B computing covariance matrices
    status1 = 1;
    status2 = 1;
    nrealization = 0;

    verb_print(cmd->verbose,
    "\n%s: reading files again and final computation of covariance matrices...",
               routine_name);
    do {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/%s%03d/%s/",
                cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B file 1
            sprintf(namebuf1,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                           "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s...",routine_name, namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
            {
                verb_print(cmd->verbose,
                            "\n%s: Input file does not exist: %s\n",
                            routine_name, namebuf2);
                status2 = 0;
            }
            //E file 2

            if (status1 == 0 || status2 == 0) {
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                      cmd->verbose, cmd->verbose_log, gd->outlog
                                      );
                if (npts1 != npts2) {
                    error("\n%s: in realization r%03d: %s %d %d",
                          routine_name, nrealization,
                          "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        if (scanopt(cmd->options, "same-infiles")) {
                            matElement=
                            0.5*(mat1[n1][n2]+mat2[n1][n2])-ZvAvg[m][N*(n1-1)+n2];
                            for (n3=1; n3<=cmd->sizeHistN; n3++) {
                                for (n4=1; n4<=cmd->sizeHistN; n4++) {
                                    matElement2=
                                    0.5*(mat1[n3][n4]+mat2[n3][n4])-ZvAvg[m][N*(n3-1)+n4];
                                    matCovMat[m][N*(n1-1)+n2][N*(n3-1)+n4] +=
                                    matElement*matElement2;
                                }
                            }
                        } else {
                            matElement=
                            mat1[n1][n2]+mat2[n1][n2]-ZvAvg[m][N*(n1-1)+n2];
                            for (n3=1; n3<=cmd->sizeHistN; n3++) {
                                for (n4=1; n4<=cmd->sizeHistN; n4++) {
                                    matElement2=
                                    mat1[n3][n4]+mat2[n3][n4]-ZvAvg[m][N*(n3-1)+n4];
                                    matCovMat[m][N*(n1-1)+n2][N*(n3-1)+n4] +=
                                    matElement*matElement2;
                                }
                            }

                        }
                    }
                }
                    //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    NR = nrealization-1;

    verb_print(cmd->verbose,
               "\n%s: number of realization analyzed (again): %d\n",
               routine_name, NR);

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        for (i=1; i<=Nv; i++) {
            for (int j=1; j<=Nv; j++) {
                matCovMat[m][i][j] /= ((real)(NR-1));          // sample mean
            }
        }
    }
    //E computing covariance matrices

    verb_print(cmd->verbose,
               "\n%s: saving covariance matrices...",
               routine_name);

    if (NR > 2) {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            sprintf(namebuf3, "%s%s_%s_%d%s",
                    rootDirPath, cmd->outfile, "CovMat", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\n%s: opening file %s... to save statistics",
                       routine_name, namebuf3);
            outstr1 = stropen(namebuf3, "w!");
            //E output file 1
            for (n1=1; n1<=Nv; n1++) {
                for (n2=1; n2<=Nv; n2++) {
                    fprintf(outstr1,"%16.8e ",matCovMat[m][n1][n2]);
                }
                fprintf(outstr1,"\n");
            }
            fclose(outstr1);
        } // ! end m loop
    } // ! nrealization > 3
    verb_print(cmd->verbose," done.\n");

//
//E read realization files and compute covariance matrices values
//

    free_dmatrix(Zv,1,cmd->mChebyshev+1,1,cmd->sizeHistN*cmd->sizeHistN);
    free_dmatrix(ZvAvg,1,cmd->mChebyshev+1,1,cmd->sizeHistN*cmd->sizeHistN);
    free_dmatrix3D(matCovMat,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mattmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(matStdDev,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(matAvg,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat3,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat2,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(rBin,1,cmd->sizeHistN);

    return SUCCESS;
}

// routine to compute edge corrections using two saved histZetaM histograms
global int computeEdgeCorrections(struct cmdline_data* cmd,
                                  struct  global_data* gd)
{
    char namebuf1[256];
    char namebuf2[256];
    char namebuf3[256];
    char namebuf4[256];
    char namebuf5[256];
    struct stat buf;
    stream instr1;
    stream instr2;
    stream instr3;
    stream outstr1;
    stream outstr2;
    int status1 = 1;
    int status2 = 1;
    int status3 = 1;
    int ifile = 0;
    int m;
    int n1;
    int n2;
    int npts1;
    int npts2;
    int npts3;
    real **mat1;
    real **mat2;
    real ***mat3;
    real ***mat4;
    real ***mat5;
    real *rBins;
    int mChebyshev;
    real rBin, rbinlog;

    char rootDirPath1[MAXLENGTHOFFILES];
    char preFileName1[MAXLENGTHOFFILES];
    char rootDirPath2[MAXLENGTHOFFILES];
    char preFileName2[MAXLENGTHOFFILES];

    extractInputRootDir(gd->infilenames[ifile], rootDirPath1, preFileName1, ifile,
                        cmd->verbose, cmd->verbose_log, gd->outlog
                        );
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: rootDir input file(%d) %s\n",
                ifile, rootDirPath1);
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: preFileName input file (%d) %s\n",
                 ifile, preFileName1);
    extractInputRootDir(gd->infilenames[ifile], rootDirPath2, preFileName2, ifile+1,
                        cmd->verbose, cmd->verbose_log, gd->outlog
                        );
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: rootDir input file(%d) %s\n",
                ifile+1, rootDirPath2);
    verb_print_q(2,cmd->verbose,
                "computeEdgeCorrections: preFileName input file (%d) %s\n",
                 ifile+1, preFileName2);
    
    //B read one file to get needed info...
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: first read one file to get some needed info:");
    m=1;
    //B file 1
    sprintf(namebuf1,"%s/%s_%s_%d%s",
            rootDirPath1, preFileName1, "cos", m, EXTFILES);
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
    if (stat(namebuf1, &buf) == 0)               // no input file exists?
        instr1 = stropen(namebuf1, "r");
    else {
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf1);
        status1 = 0;
    }
    //E file 1
    int nrow, ncol;
    if (status1 == 0) {
        fclose(instr1);
        error("\ncomputeEdgeCorrections: %s: %s",
              "input file does not exist",
              namebuf1);
    } else {
        //B processing input files
        inout_InputDataMatrix_info(namebuf1, &nrow, &ncol,
                                   cmd->verbose, cmd->verbose_log, gd->outlog
                                   );
        if (nrow != ncol) {
            error("\ncomputeEdgeCorrections: in input files: %s %d %d",
                  "nrow and ncol are different:", nrow, ncol);
        }
        fclose(instr1);
        //E
    } // ! status
    cmd->sizeHistN = nrow;
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: sizeHistN: %d\n",
               cmd->sizeHistN);
    //E read one file to get needed info...

    // should know the size of the matrix
    mat1 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    mat2 = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);

//B read rBin file
    rBins = dvector(1,cmd->sizeHistN);
    //B file 3
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: reading rbins...");
    sprintf(namebuf3,"%s/%s%s",
            rootDirPath1, "rbins", EXTFILES);
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: opening file %s...",namebuf3);
    if (stat(namebuf3, &buf) == 0)               // no input file exists?
        instr3 = stropen(namebuf3, "r");
    else {
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                    namebuf3);
        status3 = 0;
    }
    //E file 3
    if (status3 == 0) {
        error("\ncomputeEdgeCorrections: the input file does not exist: %s",
              namebuf3);
    } else {
        //B processing input file
        inout_InputDataVector(namebuf3, rBins, &npts3,
                              cmd->verbose, cmd->verbose_log, gd->outlog
                              );
        if (npts3 != cmd->sizeHistN) {
            error("\ncomputeEdgeCorrections: in rBin file: %s %d %d",
                  "npts is not equal to sizeHistN:",
                  npts3,cmd->sizeHistN);
        }
        fclose(instr3);
        //E
    } // ! status3
    verb_print(cmd->verbose, "\ndone.\n");
//E read rBin file

    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: touching multipole files:");
    mChebyshev = 1;
    status1 = 1;
    do {
        sprintf(namebuf1,"%s/%s_%s_%d%s",
                rootDirPath1, preFileName1, "cos", mChebyshev, EXTFILES);
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0) {            // no input file exists?
            instr1 = stropen(namebuf1, "r");
            fclose(instr1);
        } else {
            verb_print(cmd->verbose,
                       "\nstatHistograms: Input file does not exist: %s\n",
                       namebuf1);
            status1 = 0;
            mChebyshev--;
        }
        mChebyshev++;
    } while (status1);

    cmd->mChebyshev = --mChebyshev - 1;
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: found %d multipole files\n",
               cmd->mChebyshev+1);

    //B reading histograms files set 1
    verb_print(cmd->verbose,
               "\ncomputeEdgeCorrections: reading multipole files set 1:");
    mat3 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    CLRM_ext(mat1, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
    }
    //B input monopoles
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        //B input files cos type
        status1 = 1;
        sprintf(namebuf1,"%s/%s_%s_%d%s",
                rootDirPath1, preFileName1, "cos", m, EXTFILES);
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0)               // no input file exists?
            instr1 = stropen(namebuf1, "r");
        else {
            verb_print(cmd->verbose,
                       "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                        namebuf1);
            status1 = 0;
        }
        if (status1 == 0) {
            fclose(instr1);
            error("\ncomputeEdgeCorrections: input file does not exist: %s",
                  namebuf1);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts1 != cmd->sizeHistN) {
                error("\ncomputeEdgeCorrections: in file %s: %s %d %d",
                      namebuf1, "npts, sizeHistN are different:",
                      npts1, cmd->sizeHistN);
            }
            fclose(instr1);
            //E
        } // ! status
        //E input files cos type
        
        //B input files sin type
        status2 = 1;
        sprintf(namebuf2,"%s/%s_%s_%d%s",
                rootDirPath1, preFileName1, "sin", m, EXTFILES);
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf2);
        if (stat(namebuf2, &buf) == 0)               // no input file exists?
            instr2 = stropen(namebuf2, "r");
        else {
            verb_print(cmd->verbose,
                       "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                        namebuf2);
            status2 = 0;
        }
        if (status2 == 0) {
            fclose(instr2);
            error("\ncomputeEdgeCorrections: input file does not exist: %s",
                  namebuf2);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts2 != cmd->sizeHistN) {
                error("\ncomputeEdgeCorrections: in file %s: %s %d %d",
                      namebuf1, "npts, sizeHistN are different:",
                      npts2, cmd->sizeHistN);
            }
            fclose(instr2);
            //E
        } // ! status
        //E input files sin type

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                mat3[m][n1][n2] = mat1[n1][n2] + mat2[n1][n2];
            }
        }
    } // ! loop monopoles
    //E input monopoles
    //E reading histograms files set 1

    //B reading histograms files set 2 multipoles of N
    verb_print(cmd->verbose,
               "\n\ncomputeEdgeCorrections: reading multipole files set 2:");
    mat4 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    CLRM_ext(mat1, cmd->sizeHistN);
    CLRM_ext(mat2, cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat4[m], cmd->sizeHistN);
    }
    //B input monopoles
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        //B input files cos type
        status1 = 1;
        sprintf(namebuf1,"%s/%s_%s_%d%s",
                rootDirPath2, preFileName2, "cos_N", m, EXTFILES);
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf1);
        if (stat(namebuf1, &buf) == 0)               // no input file exists?
            instr1 = stropen(namebuf1, "r");
        else {
            verb_print(cmd->verbose,
                       "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                        namebuf1);
            status1 = 0;
        }
        if (status1 == 0) {
            fclose(instr1);
            error("\ncomputeEdgeCorrections: input file does not exist: %s",
                  namebuf1);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf1, mat1, &npts1,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts1 != cmd->sizeHistN) {
                error("\ncomputeEdgeCorrections: in file %s: %s %d %d",
                      namebuf1, "npts, sizeHistN are different:",
                      npts1, cmd->sizeHistN);
            }
            fclose(instr1);
            //E
        } // ! status
        //E input files cos type
        
        //B input files sin type
        status2 = 1;
        sprintf(namebuf2,"%s/%s_%s_%d%s",
                rootDirPath2, preFileName2, "sin_N", m, EXTFILES);
        verb_print(cmd->verbose,
                   "\ncomputeEdgeCorrections: opening file %s...",namebuf2);
        if (stat(namebuf2, &buf) == 0)               // no input file exists?
            instr2 = stropen(namebuf2, "r");
        else {
            verb_print(cmd->verbose,
                       "\ncomputeEdgeCorrections: Input file does not exist: %s\n",
                        namebuf2);
            status2 = 0;
        }
        if (status2 == 0) {
            fclose(instr2);
            error("\ncomputeEdgeCorrections: input file does not exist: %s",
                  namebuf2);
        } else {
            //B processing input file
            inout_InputDataMatrix(namebuf2, mat2, &npts2,
                                  cmd->verbose, cmd->verbose_log, gd->outlog
                                  );
            if (npts2 != cmd->sizeHistN) {
                error("\ncomputeEdgeCorrections: in file %s: %s %d %d",
                      namebuf1, "npts, sizeHistN are different:",
                      npts2, cmd->sizeHistN);
            }
            fclose(instr2);
            //E
        } // ! status
        //E input files sin type

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                mat4[m][n1][n2] = mat1[n1][n2] + mat2[n1][n2];
            }
        }
    } // ! loop monopoles
    //E input monopoles
    //E reading histograms files set 2

    mat5 = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat5[m], cmd->sizeHistN);
    }

    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            matrixClm(cmd, gd, mat3, mat4, n1, n2, mat5);
            if (cmd->verbose_log>=3)  {
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "\n\nhistZetaM elements again (%d, %d):\n\n",
                               n1, n2);
                for (m=1; m<=cmd->mChebyshev+1; m++) {
                        verb_log_print(cmd->verbose_log, gd->outlog,
                                       "%g\n", mat5[m][n1][n2]);
                }
            }
        }
    }

    verb_print_q(2, cmd->verbose, "\n\n");
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf4, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_EE", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf4);
        outstr1 = stropen(namebuf4, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr1,"%16.8e ",mat5[m][n1][n2]);
            }
            fprintf(outstr1,"\n");
        }
        fclose(outstr1);
    }

    //B  and saves matrix ZetaM for each m multipole at a set of theta2 angles
    if (scanopt(cmd->options, "out-m-HistZeta")) {
        real Zeta;
        real Zeta2;
        real Zeta3;
        real Zeta4;
        real Zeta5;
        int Nbins;
        
        verb_print_q(2, cmd->verbose, "Printing : logscale ... %d %g %g\n",
                     cmd->useLogHist, cmd->rminHist, cmd->rangeN);
        Nbins = cmd->sizeHistN;
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            sprintf(namebuf5, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                    "_EE", m, EXTFILES);
            verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",
                         namebuf5);
            outstr2 = stropen(namebuf5, "w!");
            fprintf(outstr2,MHISTZETAHEADER);
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                if (cmd->useLogHist) {
                    gd->deltaR = rlog10(cmd->rangeN/cmd->rminHist)/cmd->sizeHistN;
                    if (cmd->rminHist==0) {
                        rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                        + rlog10(cmd->rangeN);
                    } else {
                        rbinlog = rlog10(cmd->rminHist)
                                    + ((real)(n1)-0.5)*gd->deltaR;
                    }
                    rBin=rpow(10.0,rbinlog);
                } else {
                    gd->deltaR = (cmd->rangeN-cmd->rminHist)/cmd->sizeHistN;
                    rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
                }
                Zeta = mat5[m][n1][n1];
                Zeta2 = mat5[m][n1][(int)(Nbins/4.0)];
                Zeta3 = mat5[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = mat5[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = mat5[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr2,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
            }
            fclose(outstr2);
        }
    }
    //E

    free_dmatrix3D(mat5,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat4,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(mat3,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat2,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(mat1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(rBins,1,cmd->sizeHistN);

    return SUCCESS;
}

#ifdef USEGSL
global int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    double ***mat3, double ***mat4,
                    int n1, int n2, double ***mat5)
{
    // 0 <= l, m < 2*mChebyshev + 1
    // 1 <= n1, n2 <= sizeHistN
    int l, m;
    int lmx;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    real C1;
    int s;

    gsl_vector * bl = gsl_vector_alloc (neqs);
    gsl_matrix * Clm = gsl_matrix_alloc (neqs, neqs);
    gsl_matrix * ClmChk = gsl_matrix_alloc (neqs, neqs);
    gsl_vector *x = gsl_vector_alloc (neqs);
    gsl_permutation * p = gsl_permutation_alloc (neqs);

    gsl_vector *t = gsl_vector_alloc (neqs);
    gsl_matrix * u = gsl_matrix_alloc (neqs, neqs);
    real v;

    for (l=0; l<neqs; l++) {
        gsl_vector_set(bl, l, 0.0);
        for (m=0; m<neqs; m++) {
            gsl_matrix_set(Clm, l, m, 0.0);
            gsl_matrix_set(ClmChk, l, m, 0.0);
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nMatrix and b elements (%d, %d):\n\n", n1, n2);
    for (l=0; l<neqs; l++) {
        if (l<=mx)
            lmx = mx-(l+1)+2;
        else
            lmx = (l-1)+2-mx;
//B PIVOTLOOP
        gsl_vector_set(bl, l, mat3[lmx][n1][n2]/mat4[1][n1][n2]);
//E
        if (l<=mx) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), gsl_vector_get(bl, l), l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, gsl_vector_get(bl, l), l, lmx);
        }
        for (m=0; m<neqs; m++) {
            if (l-m>=-mx && l-m<0) {
//B PIVOTLOOP
                C1 = mat4[m-l+1][n1][n2]/mat4[1][n1][n2];
//E
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                //B
                if (scanopt(cmd->options, "full-sky")) {
                    if (l!=m) {
                        gsl_matrix_set(Clm, l, m, 0.0);
                        gsl_matrix_set( ClmChk, l, m, 0.0);
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            } // ! l-m >= -mx && l-m < 0
            if (l-m>=0 && l-m<=mx) {
                C1 = mat4[l-m+1][n1][n2]/mat4[1][n1][n2];
                gsl_matrix_set(Clm, l, m, C1);
                gsl_matrix_set( ClmChk, l, m, gsl_matrix_get(Clm, l, m) );
                //B
                if (scanopt(cmd->options, "full-sky")) {
                    if (l!=m) {
                        gsl_matrix_set(Clm, l, m, 0.0);
                        gsl_matrix_set( ClmChk, l, m, 0.0);
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", gsl_matrix_get(Clm, l, m));
                continue;
            } // ! l-m >= 0 && l-m <= mx
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", gsl_matrix_get(Clm, l, m));
        } // ! loop m
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    } // ! loop l

    gsl_linalg_LU_decomp (Clm, p, &s);
    gsl_linalg_LU_solve (Clm, p, bl, x);
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"x = \n");
        gsl_vector_fprintf (gd->outlog, x, "%g");
    }

    // check A x = b
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=0; l<neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ( gsl_matrix_get(ClmChk,l,m)*gsl_vector_get(x,m) );
            }
            gsl_vector_set(t, l, v);
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",
                           gsl_vector_get(bl,l),gsl_vector_get(t,l));
        }
    }

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nhistZetaM elements:\n\n");
    for (l=0; l<neqs; l++) {
        if (l>=mx) {
            lmx = l+1-mx;
            mat5[lmx][n1][n2] = gsl_vector_get(x, l);
            if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%g\n", mat5[lmx][n1][n2]);
        }
    }

    gsl_matrix_free (u);
    gsl_vector_free (t);
    gsl_permutation_free (p);
    gsl_vector_free (x);
    gsl_matrix_free (ClmChk);
    gsl_matrix_free (Clm);
    gsl_vector_free (bl);

    return SUCCESS;
}
#else // ! USEGSL
// routine to compute matriz elements Clm and its inverse
//  Be careful with the use of NONORMHISTON and NMultipoles switches:
//  NMultipolesON = 1 and NONORMHISTON = 1
global int matrixClm(struct cmdline_data* cmd, struct  global_data* gd,
                    double ***mat3, double ***mat4,
                    int n1, int n2, double ***mat5)
{
// 1 <= l, m <= 2*mChebyshev + 1
// 1 <= n1, n2 <= sizeHistN
    int l, m;
    int j;
    int *indx;
    real p;
    real **Clm;
    real **ClmChk;
    real **u;
    real *bl;
    real *blChk;
    real *t;
    real C1;

    int lm, lp;
    int neqs=2*cmd->mChebyshev+1;
    int mx=cmd->mChebyshev;
    int lmx;

    Clm = dmatrix(1,neqs,1,neqs);
    ClmChk = dmatrix(1,neqs,1,neqs);
    u = dmatrix(1,neqs,1,neqs);
    bl = dvector(1,neqs);
    blChk = dvector(1,neqs);
    t = dvector(1,neqs);
    indx = ivector(1,neqs);

    CLRM_ext(Clm,neqs);
    CLRM_ext(ClmChk,neqs);
    CLRV_ext(bl,neqs);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\n\nMatrix and b elements (%d, %d):\n\n", n1, n2);
    for (l=1; l<=neqs; l++) {
        if (l<=mx+1)
            lmx = (mx-(l+1)+2) +1;
        else
            lmx = (l-2)+2-mx;

        bl[l] = mat3[lmx][n1][n2]/mat4[1][n1][n2];
        blChk[l] = bl[l];

        if (l<=mx+1) {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               -(lmx-1), bl[l], l, lmx);
        } else {
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "b%d : %g :: %d %d\n",
                               lmx-1, bl[l], l, lmx);
        }
        for (m=1; m<=neqs; m++) {
            if (l-m>=-mx && l-m<0) {

                //B check indexs m-l+1 !! m start at 1 and GSL version starts at 0
                C1 = mat4[m-l+1][n1][n2]/mat4[1][n1][n2];

                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];

                //B
                if (scanopt(cmd->options, "full-sky")) {
                    if (l!=m) {
                        Clm[l][m] = 0.0;
                        ClmChk[l][m] = 0.0;
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (l-m>=0 && l-m<=mx) {
                C1 = mat4[l-m+1][n1][n2]/mat4[1][n1][n2];

                Clm[l][m] = C1;
                ClmChk[l][m] = Clm[l][m];

                //B
                if (scanopt(cmd->options, "full-sky")) {
                    if (l!=m) {
                        Clm[l][m] = 0.0;
                        ClmChk[l][m] = 0.0;
                    }
                }
                //E
                if (cmd->verbose_log>=3)
                    verb_log_print(cmd->verbose_log, gd->outlog,
                                   "%g ", Clm[l][m]);
                continue;
            }
            if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log, gd->outlog,
                               "%g ", Clm[l][m]);
        }
        if (cmd->verbose_log>=3)
            verb_log_print(cmd->verbose_log, gd->outlog,"\n\n");
    }

    ludcmp(Clm,neqs,indx,&p);
    lubksb(Clm,neqs,indx,bl);

    // vector solutions
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "\nVector solution:\n");
        for (l=1;l<=neqs;l++) {
            verb_log_print(cmd->verbose_log, gd->outlog,"%8s %g\n"," ", bl[l]);
        }
    }

    // check A x = b
    real v;
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\nA x = b:\n");
        for (l=1; l<=neqs; l++) {
            v = 0.0;
            for (m=0; m<neqs; m++) {
                v += ClmChk[l][m]*bl[m];
            }
            t[l] = v;
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%8s %g %g\n"," ",blChk[l],t[l]);
        }
    }

//B correction 2025-04-06
    for (l=1; l<=neqs; l++) {
        if (l>=mx+2) {
            lmx = (l-2)+2-mx;
            mat5[lmx][n1][n2] = bl[l];
            verb_log_print(cmd->verbose_log, gd->outlog,
                           "%g\n", mat5[lmx][n1][n2]);
        }

    }
//E correction 2025-04-06

    free_ivector(indx,1,cmd->mChebyshev+1);
    free_dvector(t,1,cmd->mChebyshev+1);
    free_dvector(blChk,1,cmd->mChebyshev+1);
    free_dvector(bl,1,cmd->mChebyshev+1);
    free_dmatrix(u,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(ClmChk,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);
    free_dmatrix(Clm,1,cmd->mChebyshev+1,1,cmd->mChebyshev+1);

    return SUCCESS;
}

#endif // ! USEGSL

#undef MHISTZETAHEADER
#undef MHISTZETA
#undef MHISTZETAHEADERSTDDEV
#undef MHISTZETASTDDEV

//E

//B socket:
#ifdef ADDONS
#include "cballsutils_include_02.h"
#endif
//E
