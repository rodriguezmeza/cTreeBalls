/* ==============================================================================
!	MODULE: treeutils.c			[cTreeBalls]									!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:	april 2023                                                  !
!	Purpose: 3-point correlation function computation       					!
!	Language: C																	!
!	Use:				    								                    !
!	Major revisions:															!
!==============================================================================*/
//        1          2          3          4          5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"


global int search_init_sincos_omp(struct  cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
    int n;
    int m;

    if (cmd->computeTPCF) {
        hist->ChebsT = dvector(1,cmd->mchebyshev+1);
        hist->ChebsU = dvector(1,cmd->mchebyshev+1);
    }
    hist->histNthread = dvector(1,cmd->sizeHistN);
    hist->histNSubthread = dvector(1,cmd->sizeHistN);
// 2pcf
    hist->histNSubXi2pcfthread = dvector(1,cmd->sizeHistN);
//B kappa Avg Rmin
    hist->histNSubXi2pcfthreadp = dvector(1,cmd->sizeHistN);
    hist->histNSubXi2pcfthreadtotal = dvector(1,cmd->sizeHistN);
//E
//
    hist->histXi2pcfthread = dvector(1,cmd->sizeHistN);
    hist->histXi2pcfthreadsub = dvector(1,cmd->sizeHistN);
    if (cmd->computeTPCF) {
        hist->histXithreadcos = dmatrix(1,cmd->mchebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mchebyshev+1,1,cmd->sizeHistN);
        
        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos = dmatrix3D(1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        
        hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
    }

   for (n = 1; n <= cmd->sizeHistN; n++) {
       hist->histNthread[n] = 0.0;
       hist->histNSubthread[n] = 0.0;
       hist->histNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
       hist->histNSubXi2pcfthreadp[n] = 0.0;
       hist->histNSubXi2pcfthreadtotal[n] = 0.0;
//E
       hist->histXi2pcfthread[n] = 0.0;
       hist->histXi2pcfthreadsub[n] = 0.0;
   }

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mchebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        }
    }

    return _SUCCESS_;
}

global int search_free_sincos_omp(struct  cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_sincos_omp hist)
{
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsincos,1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsin,1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadcos,1,cmd->mchebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadsin,1,cmd->mchebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadcos,1,cmd->mchebyshev+1,1,cmd->sizeHistN);
    }
    free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
    free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
//B kappa Avg Rmin
    free_dvector(hist->histNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
    free_dvector(hist->histNSubXi2pcfthreadp,1,cmd->sizeHistN);
//E
    free_dvector(hist->histNSubXi2pcfthread,1,cmd->sizeHistN);
    free_dvector(hist->histNSubthread,1,cmd->sizeHistN);
    free_dvector(hist->histNthread,1,cmd->sizeHistN);
    if (cmd->computeTPCF) {
        free_dvector(hist->ChebsU,1,cmd->mchebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mchebyshev+1);
    }
    
    return _SUCCESS_;
}


global int computeBodyProperties_sincos_omp(struct  cmdline_data* cmd, struct  global_data* gd,
                                            bodyptr p, int nbody,
                                            gdhistptr_sincos_omp hist)
{
    int n;
    int m;
    real xi, xi_2p;

// BODY3
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
//

    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mchebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= MAX(hist->histNSubthread[n],1.0);
                hist->histXithreadsin[m][n] /= MAX(hist->histNSubthread[n],1.0);
            }
        
        for (m=1; m<=cmd->mchebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos, hist->histXithreadcos[m], hist->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin, hist->histXithreadsin[m], hist->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos, hist->histXithreadsin[m], hist->histXithreadcos[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],hist->histZetaMthreadcos[m],hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],hist->histZetaMthreadsin[m],hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],hist->histZetaMthreadsincos[m],hist->histZetaMtmpsincos,cmd->sizeHistN);
        }
    }

    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
    }

    return _SUCCESS_;
}


global int search_init_gd_hist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mchebyshev+1; m++) {
            CLRM_ext(gd->histZetaM[m], cmd->sizeHistN);
        }
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histN[n] = 0.0;
        gd->histNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gd->histNSubXi2pcftotal[n] = 0.0;
//E
        gd->histXi2pcf[n] = 0.0;
        if (cmd->computeTPCF) {
            for (m = 1; m <= cmd->mchebyshev+1; m++)
                gd->histXi[m][n] = 0.0;
        }
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return _SUCCESS_;
}

global int search_init_gd_hist_sincos(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;
    int m;

    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mchebyshev+1; m++) {
            CLRM_ext(gd->histZetaMcos[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsin[m], cmd->sizeHistN);
            CLRM_ext(gd->histZetaMsincos[m], cmd->sizeHistN);
        }
    }
    for (n = 1; n <= cmd->sizeHistN; n++) {
        gd->histN[n] = 0.0;
        gd->histNSubXi2pcf[n] = 0.0;
//B kappa Avg Rmin
        gd->histNSubXi2pcftotal[n] = 0.0;
//E
        gd->histXi2pcf[n] = 0.0;
        if (cmd->computeTPCF) {
            for (m = 1; m <= cmd->mchebyshev+1; m++)
                gd->histXi[m][n] = 0.0;
        }
    }
    gd->actmax = gd->nbbcalc = gd->nbccalc = gd->ncccalc = 0;

    return _SUCCESS_;
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
    if (!(scanopt(cmd->options, "Rapaport"))) {
        gd->histN[1]-=nbody;
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
        edd[n] = 1./rsqrt(gd->histN[n]);

    for (n = 1; n <= cmd->sizeHistN; n++) {
        if(gd->histN[n]==0) {
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
            if (scanopt(cmd->options, "Rapaport")) {
                normFac = Vol/(2.0*PI*rpow(gd->deltaR,3.0)*nbody*nbody);
                gd->histCF[n] = gd->histN[n] * normFac / rsqr((int)n-0.5);
            } else {
//B This is version CB
                vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
                rho_r=gd->histN[n]/((real)nbody*vr);
                corr[n]=rho_r/rho_av-1;                 // Correlation function
                ercorr[n]=(1+corr[n])*edd[n];           // Poisson errors
                gd->histCF[n] = corr[n];
//E
            }
#else
            if (scanopt(cmd->options, "Rapaport")) {
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
                gd->histCF[n] = gd->histN[n] * normFac / ((int)n-0.5);
            } else { // This should be CB version...
                normFac = Vol/(PI*rpow(gd->deltaR,2.0)*nbody*nbody);
                gd->histCF[n] = gd->histN[n] * normFac / ((int)n-0.5) - 1.0;
            }
#endif // ! NDIM
        }
    }

    free_dvector(ercorr,1,cmd->sizeHistN);
    free_dvector(corr,1,cmd->sizeHistN);
    free_dvector(edd,1,cmd->sizeHistN);

    return _SUCCESS_;
}


global int search_compute_HistN(struct  cmdline_data* cmd, struct  global_data* gd, int nbody)
{
    int n;
    real normFac;

// Check this factor is correct
    normFac = 0.5;

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histN[n] *= normFac;

    if (scanopt(cmd->options, "and-CF"))
        search_compute_Xi(cmd, gd, nbody);

    return _SUCCESS_;
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
global bool reject_balls(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q, real *drpq, vector dr)
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

global bool reject_cell_balls(struct  cmdline_data* cmd, struct  global_data* gd, nodeptr p, nodeptr q, real *drpq, vector dr)
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

global int compute_cosphi(real dr1, vector dr, real *cosphi, gdhist hist)
{
    real s;

    if (NDIM == 3) {
        DOTVP(s, dr, hist.dr0);
        *cosphi = -s/(dr1*hist.drpq);
    } else
        *cosphi = -dr[1]/dr1;        // x,y

    if (rabs(*cosphi)>1.0) error("\ncossphi must be in (-1,1): %g\n",cosphi);

    return _SUCCESS_;
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
#ifdef BALLS
#include "treeutils_balls_omp.h"
#endif

    return _SUCCESS_;
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
global int spherical_to_cartesians(struct  cmdline_data* cmd, struct  global_data* gd,
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

    return _SUCCESS_;
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

    return _SUCCESS_;
}

#ifdef ADDONS
#include "tree_include.h"
#endif
