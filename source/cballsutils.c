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

#include "globaldefs.h"


global int search_init_sincos_omp(struct  cmdline_data* cmd, 
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
    int n;
    int m;

    if (cmd->computeTPCF) {
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
    }
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNNSubthread = dvector(1,cmd->sizeHistN);
// 2pcf
    hist->histNNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
//B kappa Avg Rmin
    hist->histNNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
//E
//
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    if (cmd->computeTPCF) {
        hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        
        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
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
    }

   for (n = 1; n <= cmd->sizeHistN; n++) {
       hist->histNthread[n] = 0.0;
       hist->histNNSubthread[n] = 0.0;
       hist->histNNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
       hist->histNNSubXi2pcfthreadp[n] = 0.0;
       hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
//E
       hist->histXi2pcfthread[n] = 0.0;
       hist->histXi2pcfthreadsub[n] = 0.0;
   }

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(hist->histZetaMthreadcossin[m], cmd->sizeHistN);
        }
    }

    return SUCCESS;
}

global int search_free_sincos_omp(struct  cmdline_data* cmd,
                                  struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
    if (cmd->computeTPCF) {
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
    }
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
//B kappa Avg Rmin
    free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
//E
    free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
    if (cmd->computeTPCF) {
        free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
    }
    
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
        xi = Kappa(p)/nbody;
        xi_2p = Kappa(p);
        //B kappa Avg Rmin
        if (scanopt(cmd->options, "smooth-pivot")) {
            xi_2p = KappaRmin(p);
            xi = NbRmin(p)*xi_2p/nbody;
        }
        //E
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
    //E Normalization of histograms

    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            //B Normalization of histograms
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            }
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
            MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            MULMS_ext(hist->histZetaMtmpcossin,hist->xiOUTVPcossin,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],
                     hist->histZetaMthreadcos[m],hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],
                     hist->histZetaMthreadsin[m],hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],
                    hist->histZetaMthreadsincos[m],hist->histZetaMtmpsincos,cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            ADDM_ext(hist->histZetaMthreadcossin[m],
                    hist->histZetaMthreadcossin[m],hist->histZetaMtmpcossin,cmd->sizeHistN);
        }
    }

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

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaM[m], cmd->sizeHistN);
        }
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gd->histNNSubXi2pcftotal[n] = 0.0;
//E
        gd->histXi2pcf[n] = 0.0;
        if (cmd->computeTPCF) {
            for (m = 1; m <= cmd->mChebyshev+1; m++)
                gd->histXi[m][n] = 0.0;
        }
    }
    
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaGmRe[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaGmIm[m], cmd->sizeHistN);
        }
    }

    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

global int search_init_gd_hist_sincos(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
            // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
            CLRM_ext(gd->histZetaMcossin[m], cmd->sizeHistN);
        }
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histNN[n] = 0.0;
        gd->histNNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gd->histNNSubXi2pcftotal[n] = 0.0;
//E
        gd->histXi2pcf[n] = 0.0;
        if (cmd->computeTPCF) {
            for (m = 1; m <= cmd->mChebyshev+1; m++)
                // HERE MUST BE gd->histXicos and gd->histXisin
                gd->histXi[m][n] = 0.0;
        }
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return SUCCESS;
}

//B Computation of histogram of all B-B encounters
// The correlation function is estimated as:
//    xi=(V/v(r))*(DD(r)/N^2)
// where v(r)=4*pi*((r+dr/2)^3-(r-dr/2)^3)/3, V=box_size^3 and N is the
// total # particles.
local int search_compute_Xi(struct  cmdline_data* cmd, struct  global_data* gd, int nbody)
{
    int k;
    int n;
    real normFac;
    real Vol;

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
                r0=(real)n*gd->deltaR;
                r1=(real)(n+1)*gd->deltaR;
            }

#if (NDIM==3)
            if (!scanopt(cmd->options, "cute-box")) {
                normFac = Vol/(2.0*PI*rpow(gd->deltaR,3.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5);
// This line gives results in agreement with CB:
                gd->histCF[n] = gd->histNN[n] * normFac / rsqr((int)n-0.5) -1.0;
            } else {
//B This is version CB
                vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
                rho_r=gd->histNN[n]/((real)nbody*vr);
                corr[n]=rho_r/rho_av-1;             // Correlation function
                ercorr[n]=(1+corr[n])*edd[n];       // Poisson errors
                gd->histCF[n] = corr[n];            // Original line
//E
            }
#else
            if (!scanopt(cmd->options, "cute-box")) {
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
// This line gives results for rdf (radial distribution function):
//                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5);
// This line gives results in agreement with CB:
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5) - 1.0;
            } else { // This should be CB version...
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
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

// Check this factor is correct
    normFac = 0.5;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNN[n] *= normFac;

    if (scanopt(cmd->options, "and-CF"))
        search_compute_Xi(cmd, gd, nbody);

    return SUCCESS;
}


// No se usa qsize: update
global bool reject_cell(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q, real qsize)
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

global bool reject_bodycell(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q)
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
#ifdef ADDONS
#include "cballsutils_include_01.h"
#endif

    return SUCCESS;
}
//E

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
global int spherical_to_cartesians(struct  cmdline_data* cmd, 
                                   struct  global_data* gd,
                                   real theta, real phi, vector xyz)
{
    real ra, dec;

    if (scanopt(cmd->options, "no-arfken")) {
        ra = theta;
        dec = PIO2 - phi;
        xyz[0] = rcos(dec)*rcos(ra);
        xyz[1] = rcos(dec)*rsin(ra);
        xyz[2] = rsin(dec);
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

global int spherical_periodic_condition(real *thetaL, real *thetaR, real *phiL, real *phiR)
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

//B section of several routines to do pre/post processing

#define IN 1
#define OUT 0
#define SI 1
#define NO 0

// Column vector input
//  offset as NR
local int inout_InputDataVector(struct cmdline_data*, struct  global_data*,
                                string, real *, int *);

// Column vector input
//  offset as NR
local int inout_InputDataVector(struct cmdline_data* cmd, struct  global_data* gd,
                                string filename, real *vec, int *npts)
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    int col1=1;
    
    instr = stropen(filename, "r");

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
               "\nReading matrix elements from file %s... ",filename);

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                   "\n\tGeneral statistics : ");
        verb_log_print(cmd->verbose_log, gd->outlog,
                   "number of lines, words, and characters : %d %d %d", nl, nw, nc);
    }

    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\tValid numbers statistics : ");
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "nrow, ncol, nvalues : %d %d %d", nrow, ncol, nw);
    }

    rewind(instr);
    
    npoint=nrow;
    row = (realptr) allocate(ncol*sizeof(real));

    *npts = npoint;
    
    ip = 1;
    for (i=0; i<nl; i++) {
        if (lineQ[i]) {
            in_vector_ndim(instr, row, ncol);
            vec[ip] = row[col1-1];
            ++ip;
        } else {
            while ((c = getc(instr)) != EOF)        // Reading dummy line ...
                if (c=='\n') break;
        }
    }
    fclose(instr);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\t... done.\n");

    return SUCCESS;
}

local int inout_InputDataMatrix(struct cmdline_data*, struct  global_data*,
                                string, real **, int *);

local int inout_InputDataMatrix(struct cmdline_data* cmd, struct  global_data* gd,
                                string filename, real **mat, int *npts)
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    int col1=1;
    int col2=2;
    
    instr = stropen(filename, "r");

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,
               "\nReading matrix elements from file %s... ",filename);

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,
                   "\n\tGeneral statistics : ");
        verb_log_print(cmd->verbose_log, gd->outlog,
                   "number of lines, words, and characters : %d %d %d", nl, nw, nc);
    }

    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    if (cmd->verbose_log>=3) {
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\tValid numbers statistics : ");
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "nrow, ncol, nvalues : %d %d %d", nrow, ncol, nw);
    }

    rewind(instr);
    
    npoint=nrow;
    row = (realptr) allocate(ncol*sizeof(real));

    *npts = npoint;
    
    ip = 1;
    for (i=0; i<nl; i++) {
        if (lineQ[i]) {
            in_vector_ndim(instr, row, ncol);
            for (int j=0; j<ncol; j++) {
                mat[ip][j+1] = row[j];
            }
            ++ip;
        } else {
            while ((c = getc(instr)) != EOF)        // Reading dummy line ...
                if (c=='\n') break;
        }
    }
    fclose(instr);

    if (cmd->verbose_log>=3)
        verb_log_print(cmd->verbose_log, gd->outlog,"\n\t... done.\n");

    return SUCCESS;
}

#undef IN
#undef OUT
#undef SI
#undef NO

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
    int m;
    int n1;
    int n2;
    int npts1;
    int npts2;
    int npts3;
    real matElement;
    real matElementIm;
    real **mat1;
    real **mat2;
    real ***mat3;
    real ***matAvg;
    real ***matStdDev;
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

    sizeHistN = cmd->sizeHistN;
    mat1 = dmatrix(1,sizeHistN,1,sizeHistN);
    mat2 = dmatrix(1,sizeHistN,1,sizeHistN);

//B read one file to got needed info to go...
    verb_print(cmd->verbose,
               "\nstatHistograms: first read two file to get some needed info:");
    //B we use rootDir string to construct input/output file names
    sprintf(rootDirPath,"%s/%s%03d/%s/",
            cmd->rootDir,"r",nrealization,cmd->suffixOutFiles);
    //E
    m=1;
    //B file 1
    sprintf(namebuf1,"%s%s_%d%s",
            rootDirPath, gd->infilenames[ifile], m, EXTFILES);
    verb_print(cmd->verbose,
                "\nstatHistograms: opening file %s...",namebuf1);
    if (stat(namebuf1, &buf) == 0)               // no input file exists?
        instr1 = stropen(namebuf1, "r");
    else {
        verb_print(cmd->verbose,
                   "\nstatHistograms: Input file does not exist: %s\n",
                    namebuf1);
        status1 = 0;
    }
    //E file 1
    //B file 2
    sprintf(namebuf2,"%s%s_%d%s",
            rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
    verb_print(cmd->verbose,
               "\nstatHistograms: opening file %s...",namebuf2);
    if (stat(namebuf2, &buf) == 0)               // no input file exists?
        instr2 = stropen(namebuf2, "r");
    else {
        verb_print(cmd->verbose,
                   "\nstatHistograms: Input file does not exist: %s\n",
                    namebuf2);
        status2 = 0;
    }
    //E file 2

    if (status1 == 0 || status2 == 0) {
        if (status1 != 0) fclose(instr1);
        if (status2 != 0) fclose(instr2);
        error("\nstatHistogram: one of the input file does not exist: %s %s",
              namebuf1, namebuf2);
    } else {
        //B processing input files
        inout_InputDataMatrix(cmd, gd, namebuf1, mat1, &npts1);
        inout_InputDataMatrix(cmd, gd, namebuf2, mat2, &npts2);
        if (npts1 != npts2) {
            error("\nstatHistogram: in realization r%03d: %s %d %d",
                  nrealization, "npts are different:", npts1,npts2);
        }
        cmd->sizeHistN = npts1;
        fclose(instr1);
        fclose(instr2);
        //E
    } // ! status
    verb_print(cmd->verbose,
               "\nstatHistograms: mChebyshev and sizeHistN: %d %d\n",
               cmd->mChebyshev, cmd->sizeHistN);
    free_dmatrix(mat2,1,sizeHistN,1,sizeHistN);
    free_dmatrix(mat1,1,sizeHistN,1,sizeHistN);
//E read one file

//B read rBin file
    rBin = dvector(1,sizeHistN);
    //B file 3
    verb_print(cmd->verbose,
               "\nstatHistograms: reading rbins...",
               cmd->mChebyshev, cmd->sizeHistN);
    sprintf(namebuf5,"%s%s%s",
            rootDirPath, "rbins", EXTFILES);
    verb_print(cmd->verbose,
               "\nstatHistograms: opening file %s...",namebuf5);
    if (stat(namebuf5, &buf) == 0)               // no input file exists?
        instr3 = stropen(namebuf5, "r");
    else {
        verb_print(cmd->verbose,
                   "\nstatHistograms: Input file does not exist: %s\n",
                    namebuf5);
        status3 = 0;
    }
    //E file 3
    if (status3 == 0) {
        error("\nstatHistogram: the input file does not exist: %s",
              namebuf5);
    } else {
        //B processing input file
        inout_InputDataVector(cmd, gd, namebuf5, rBin, &npts3);
        if (npts3 != cmd->sizeHistN) {
            error("\nstatHistogram: in realization r%03d: %s %d %d",
                  nrealization, "npts is not equal to sizeHistN:",
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
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        CLRM_ext(mat3[m], cmd->sizeHistN);
        CLRM_ext(matAvg[m], cmd->sizeHistN);
        CLRM_ext(matStdDev[m], cmd->sizeHistN);
    }

//B read realization files and compute mean values
    verb_print(cmd->verbose,
               "\nstatHistograms: mean value computation:");
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
                       "\nstatHistograms: opening file %s...",namebuf1);
            if (stat(namebuf1, &buf) == 0)               // no input file exists?
                instr1 = stropen(namebuf1, "r");
            else {
                verb_print(cmd->verbose,
                           "\nstatHistograms: Input file does not exist: %s\n",
                           namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\nstatHistograms: opening file %s...",namebuf2);
            if (stat(namebuf2, &buf) == 0)               // no input file exists?
                instr2 = stropen(namebuf2, "r");
            else {
                verb_print(cmd->verbose,
                           "\nstatHistograms: Input file does not exist: %s\n",
                           namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
                if (status1 != 0) fclose(instr1);
                if (status2 != 0) fclose(instr2);
                break;
            } else {
                //B output file 1
                sprintf(namebuf3, "%s%s_%d%s", rootDirPath, cmd->outfile, m, EXTFILES);
                verb_print(cmd->verbose,
                           "\nstatHistograms: opening file %s... to save statistics",
                           namebuf3);
                outstr1 = stropen(namebuf3, "w!");
                //E output file 1
                //B output file 2
                sprintf(namebuf4, "%s%s%s_%d%s",
                        rootDirPath, "m", cmd->outfile, m, EXTFILES);
                verb_print(cmd->verbose,
                           "\nstatHistograms: opening file %s... to save statistics",
                           namebuf4);
                outstr2 = stropen(namebuf4, "w!");
                //E output file 2

                //B processing input files
                inout_InputDataMatrix(cmd, gd, namebuf1, mat1, &npts1);
                inout_InputDataMatrix(cmd, gd, namebuf2, mat2, &npts2);
                if (npts1 != npts2) {
                error("\nstatHistogram: in realization r%03d: %s %d %d",
                          nrealization, "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        matElement = mat1[n1][n2]+mat2[n1][n2];
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
                    fprintf(outstr2,MHISTZETA,rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
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
               "\nstatHistograms: number of realization analyzed: %d\n",nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            sprintf(namebuf3, "%s%s_%s_%d%s",
                    rootDirPath, cmd->outfile, "Avg", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\nstatHistograms: opening file %s... to save statistics",
                       namebuf3);
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
                       "\nstatHistograms: opening file %s... to save statistics",
                       namebuf6);
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
               "\n\nstatHistograms: standard deviation computation:");
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
                       "\nstatHistograms: opening file %s...",namebuf1);
            if (stat(namebuf1, &buf) != 0)               // input file exists?
//                instr1 = stropen(namebuf1, "r");
//            else
            {
                verb_print(cmd->verbose,
                           "\nstatHistograms: Input file does not exist: %s\n",
                           namebuf1);
                status1 = 0;
            }
            //E file 1
            //B file 2
            sprintf(namebuf2,"%s%s_%d%s",
                    rootDirPath, gd->infilenames[ifile+1], m, EXTFILES);
            verb_print(cmd->verbose,
                       "\nstatHistograms: opening file %s...",namebuf2);
            if (stat(namebuf2, &buf) != 0)               // input file exists?
//                instr2 = stropen(namebuf2, "r");
//            else
            {
                verb_print(cmd->verbose,
                           "\nstatHistograms: Input file does not exist: %s\n",
                           namebuf2);
                status2 = 0;
            }
            //E file 2
            
            if (status1 == 0 || status2 == 0) {
//                if (status1 != 0) fclose(instr1);
//                if (status2 != 0) fclose(instr2);
                break;
            } else {
                //B processing input files
                inout_InputDataMatrix(cmd, gd, namebuf1, mat1, &npts1);
                inout_InputDataMatrix(cmd, gd, namebuf2, mat2, &npts2);
                if (npts1 != npts2) {
                error("\nstatHistogram: in realization r%03d: %s %d %d",
                          nrealization, "npts are different:", npts1,npts2);
                }

                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    for (n2=1; n2<=cmd->sizeHistN; n2++) {
                        matElement = rabs(
                                           mat1[n1][n2]+mat2[n1][n2]
                                           -matAvg[m][n1][n2]
                                           );
//                        mat3[m][n1][n2] = matElement;
                        matStdDev[m][n1][n2] += matElement;
//                        fprintf(outstr1,"%16.8e ",matElement);
                    }
//                    fprintf(outstr1,"\n");
                }
/*
                fprintf(outstr2,MHISTZETAHEADER);
                for (n1=1; n1<=cmd->sizeHistN; n1++) {
                    Zeta = mat3[m][n1][n1];
                    Zeta2 = mat3[m][n1][(int)(Nbins/4.0)];
                    Zeta3 = mat3[m][n1][(int)(2.0*Nbins/4.0)];
                    Zeta4 = mat3[m][n1][(int)(3.0*Nbins/4.0)];
                    Zeta5 = mat3[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                    fprintf(outstr2,MHISTZETA,rBin[n1],Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
                }
*/
//                fclose(instr1);
//                fclose(instr2);
//                fclose(outstr1);
//                fclose(outstr2);
                //E
            } // ! status
        } // ! end m loop
        nrealization++;
    } while (status1 || status2);

    verb_print(cmd->verbose,
               "\nstatHistograms: number of realization analyzed: %d\n",nrealization-1);

    if (nrealization-1 > 2) {
        //B we use rootDir string to construct input/output file names
        sprintf(rootDirPath,"%s/", cmd->rootDir);
        //E
        for (m=1; m<=cmd->mChebyshev+1; m++) {
            //B output file 1
            sprintf(namebuf3, "%s%s_%s_%d%s",
                    rootDirPath, cmd->outfile, "StdDev", m, EXTFILES);
            verb_print(cmd->verbose,
                       "\nstatHistograms: opening file %s... to save statistics",
                       namebuf3);
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
                       "\nstatHistograms: opening file %s... to save statistics",
                       namebuf6);
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

#undef MHISTZETAHEADER
#undef MHISTZETA
#undef MHISTZETAHEADERSTDDEV
#undef MHISTZETASTDDEV

//E

#ifdef ADDONS
#include "cballsutils_include_02.h"
#endif
