/* ==============================================================================
!	MODULE: search_omp.c		[cTreeBalls]                                    !
!	Written by: M.A. Rodriguez-Meza.											!
!    Starting date:    april 2023                                               !
!    Purpose: 3-point correlation function computation                          !
!	Language: C																	!
!	Use: searchcalc_xxx();												!
!	Major revisions:															!
!==============================================================================*/
//        1          2          3          4          5          6          7

//B Search omp parallel routine using tree method:
//
// search=tree-omp-sincos ::
//                  searchcalc_normal_omp_sincos(btab, nbody, ipmin, ipmax,
//                                                              cat1, cat2)
//E

// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"


//B COMIENZA METODO NORMAL DE BUSQUEDA (SINCOS-OMP)

local void normal_walktree_sincos(bodyptr, nodeptr, real,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
#ifdef BODY3ON
local void sumnode_sincos_body3(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
#endif
local void sumnode_sincos_cell(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);

// search=tree-omp-sincos
global void searchcalc_normal_omp_sincos(bodyptr *btable, INTEGER *nbody,
                            INTEGER ipmin, INTEGER *ipmax, int cat1, int cat2)
{
    double cpustart;
    int nn;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "\nevalHistN: Running ... (tree omp-sincos) \n");

    if (scanopt(cmd.options, "behavior-ball"))
        verb_print(cmd.verbose, "with option behavior-ball... \n");
    if (scanopt(cmd.options, "no-one-ball"))
        verb_print(cmd.verbose, "with option no-one-ball... \n");
    if (scanopt(cmd.options, "smooth-pivot"))
        verb_print(cmd.verbose,
                   "with option smooth-pivot... rsmooth=%g\n",gd.rsmooth[0]);
#ifndef TPCF
    verb_print(cmd.verbose, "computing only 2pcf... \n");
#endif

    ThreadCount(nbody[cat1], cat1);

    search_init_gd_hist_sincos();

    INTEGER ipfalse;
    ipfalse=0;
    
    verb_print(cmd.verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd.bytes_tot/(1024.0*1024.0));


#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,root,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode, cat1, cat2, ipfalse)
 {
     bodyptr p;
     int n;
     INTEGER nbbcalcthread = 0;
     INTEGER nbccalcthread = 0;
     gdhist_sincos_omp hist;

     search_init_sincos_omp(&hist);

    INTEGER ipfalsethreads;
    ipfalsethreads = 0;

#pragma omp for nowait schedule(dynamic)
     for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
// Correct it because p and q are in differents node structures!!!... cat1!=cat2
        if (scanopt(cmd.options, "smooth-pivot")) {
//            NbRmin(p) = 0.;
//            NbRminOverlap(p) = 0.;
            if (Update(p) == FALSE) {
                ipfalsethreads++;
                continue;
            }
        }
//
        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;           // Affects only 3pcf
            hist.histXi2pcfthreadsub[n] = 0.0;      // Affects only 2pcf
        }
#ifdef TPCF
        CLRM_ext_ext(hist.histXithreadcos, cmd.mchebyshev+1, cmd.sizeHistN);
        CLRM_ext_ext(hist.histXithreadsin, cmd.mchebyshev+1, cmd.sizeHistN);
         
#if NDIM == 3
          dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
#ifdef SINGLEP
         hist.drpq = sqrt(hist.drpq2);
#else
          hist.drpq = rsqrt(hist.drpq2);
#endif
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION
          real rtheta;
          vector dr0rot;
          rtheta = xrandom(0.0, TWOPI);
          RotationVecAWRtoVecB(dr0rot, hist.dr0, Pos(p), rtheta);
          SETV(hist.dr0, dr0rot);
#endif
//E
#endif // ! NDIM
#endif // ! TPCF

         normal_walktree_sincos(p, ((nodeptr) roottable[cat2]), gd.rSizeTable[cat2], &nbbcalcthread, &nbccalcthread, &hist);
         computeBodyProperties_sincos_omp(p, nbody[cat1], &hist);

         INTEGER ip;
         ip = p - btable[cat1] + 1;
         if (ip%cmd.stepState == 0) {
             verb_log_print(cmd.verbose_log, gd.outlog, " - Completed pivot: %ld\n", ip);
         }
     } // end do body p // end pragma omp DO_BODY p

#pragma omp critical
 {
    for (n = 1; n <= cmd.sizeHistN; n++) {
        gd.histN[n] += hist.histNthread[n];
        gd.histNSub[n] += hist.histNSubthread[n];
        gd.histNSubXi2pcf[n] += hist.histNSubXi2pcfthread[n];
        gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
    }
#ifdef TPCF
     int m;
     for (m=1; m<=cmd.mchebyshev+1; m++) {
         ADDM_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],hist.histZetaMthreadcos[m],cmd.sizeHistN);
         ADDM_ext(gd.histZetaMsin[m],gd.histZetaMsin[m],hist.histZetaMthreadsin[m],cmd.sizeHistN);
         ADDM_ext(gd.histZetaMsincos[m],gd.histZetaMsincos[m],hist.histZetaMthreadsincos[m],cmd.sizeHistN);
     }
#endif
    gd.nbbcalc += nbbcalcthread;
    gd.nbccalc += nbccalcthread;
 }
    search_free_sincos_omp(&hist);

    ipfalse += ipfalsethreads;
 } // end pragma omp parallel

    real xi, den, num;
    int m;
    if (scanopt(cmd.options, "smooth-pivot")) {
        num = nbody[cat1];
        den = nbody[cat1]-ipfalse;
        xi = num/den;
        verb_print(cmd.verbose, 
                   "tree-omp-sincos: p falses found = %ld and %ld %ld\n",
                   ipfalse, num, den);
#ifdef TPCF
        for (m=1; m<=cmd.mchebyshev+1; m++) {
            MULMS_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],xi,cmd.sizeHistN);
            MULMS_ext(gd.histZetaMsin[m],gd.histZetaMsin[m],xi,cmd.sizeHistN);
            MULMS_ext(gd.histZetaMsincos[m],gd.histZetaMsincos[m],xi,cmd.sizeHistN);
        }
#endif
    }

    if (!scanopt(cmd.options, "asymmetric")) {
        for (nn = 1; nn <= cmd.sizeHistN; nn++) {
            gd.histXi2pcf[nn] /= 2.0;
            gd.histNSubXi2pcf[nn] /= 2.0;
            gd.histXi2pcf[nn] /= MAX(gd.histNSubXi2pcf[nn],1.0);
        }
    }

    if (scanopt(cmd.options, "compute-HistN")) {
        if (scanopt(cmd.options, "smooth-pivot")) {
            search_compute_HistN(nbody[cat1]-ipfalse);
        } else {
            search_compute_HistN(nbody[cat1]);
            //        search_compute_HistN(nbody[cat2]);
        }
    }

    verb_print(cmd.verbose, "tree-omp-sincos: p falses found = %ld\n",ipfalse);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void normal_walktree_sincos(bodyptr p, nodeptr q, real qsize,
                        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    nodeptr l;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    real dr1;
    vector dr;
#endif
    int n;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell((nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd.options, "no-one-ball")) {
                        accept_body(p, (nodeptr)q, &dr1, dr);
#ifdef LOGHIST
//B Segment as original
                        if (cmd.rminHist==0)
                            n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                                      + cmd.sizeHistN);
                        else
                            n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR);
                        if (n<=cmd.sizeHistN-1 && n>=1) {
                            if ( gd.deltaRV[n] < dr1 - Radius(q)
                                && dr1 + Radius(q) < gd.deltaRV[n+1]) {
                                sumnode_sincos_cell(p, ((cellptr) q),
                                        ( (cellptr) q+1),
                                        nbbcalcthread, nbccalcthread, hist);
                            } else {
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(p,l,qsize/2,
                                        nbbcalcthread, nbccalcthread, hist);
                            }
                        } else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(p,l,qsize/2,
                                        nbbcalcthread, nbccalcthread, hist);
//E Segment as original

#else // ! LOGHIST
                        real rBin;
                        n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                        rBin = cmd.rminHist + ((real)n)*gd.deltaR;
                        if ( rBin-gd.deltaR < dr1 - Radius(q)
                            && dr1 + Radius(q) < rBin ) {
                            sumnode_sincos_cell(p, ((cellptr) q),
                            ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(p,l,qsize/2,
                                        nbbcalcthread, nbccalcthread, hist);
                        }
#endif
                    } else {    // ! no-one-ball
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(p,l,qsize/2,
                                    nbbcalcthread, nbccalcthread, hist);
                    }           // ! no-one-ball
                } // ! reject_cell
            } else { // ! Type(q) == CELL
                if (Type(q) == BODY)
                    sumnode_sincos(p, ((cellptr) q),
                        ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                else
                    if (Type(p) == BODY3) {
#ifdef BODY3ON
                        sumnode_sincos_body3(p, ((cellptr) q),
                                             ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
#endif
                    }
            } // ! Type(q) == CELL
        } // ! p != q
    } // ! Update
}

local void sumnode_sincos(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
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

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {

#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                              + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNSubXi2pcfthread[n] =
                            hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(q);
#ifdef TPCF
                    REAL cosphi,sinphi;
#if NDIM == 3
                    REAL s, sy;
#ifdef SINGLEP
                    float pr0[NDIM];
#else
                    vector pr0;
#endif
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif // ! PTOPIVOTROTATION2
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
#ifdef SINGLEP
                    if (rabs(cosphi)>1.0)
                        sinphi = 0.0;
                    else
                        sinphi = sqrt(1.0 - cosphi*cosphi);
#else
                    sinphi = rsqrt(1.0 - rsqr(cosphi));
#endif
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
#ifdef SINGLEP
                    if (cosphi>1.0) cosphi = 1.0;
                    if (cosphi<-1.0) cosphi = -1.0;
#else
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
#endif
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNSubXi2pcfthread[n] =
                                hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(q);
#ifdef TPCF
                    real cosphi,sinphi;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}

#ifdef BODY3ON
local void sumnode_sincos_body3(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
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

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nbb(q);
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + Nbb(q);
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                    xi = Kappa(q);
#ifdef TPCF
                    real cosphi,sinphi;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                       "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                    *nbbcalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nbb(q);
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + Nbb(q);
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                    xi = Kappa(q);
#ifdef TPCF
                    real cosphi, sinphi;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                    *nbbcalcthread += 1;
                }
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}
#endif

local void sumnode_sincos_cell(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
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

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                              + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                    hist->histNSubXi2pcfthread[n] =
                            hist->histNSubXi2pcfthread[n] + 1.0;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.0;

#ifdef KappaAvgON
                    xi = KappaAvg(q)/Nb(q);
#else
                    xi = Kappa(q);
#endif

#ifdef TPCF
                    REAL cosphi,sinphi;
#if NDIM == 3
                    REAL s, sy;
#ifdef SINGLEP
                    float pr0[NDIM];
#else
                    vector pr0;
#endif
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
#ifdef SINGLEP
                    if (rabs(cosphi)>1.0)
                        sinphi = 0.0;
                    else
                        sinphi = sqrt(1.0 - cosphi*cosphi);
#else
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
#endif
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
#ifdef SINGLEP
                    if (cosphi>1.0) cosphi = 1.0;
                    if (cosphi<-1.0) cosphi = -1.0;
#else
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
#endif
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                    hist->histNSubXi2pcfthread[n] =
                            hist->histNSubXi2pcfthread[n] + 1.0;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.0;
                    xi = Kappa(q);
#ifdef TPCF
                    real cosphi, sinphi;
#if NDIM == 3
                    real s, sy, phi;
                    vector pr0;
//B Random rotation of dr0:
#ifdef PTOPIVOTROTATION2
                    real rtheta;
                    vector dr0rot;
                    rtheta = xrandom(0.0, TWOPI);
                    RotationVecAWRtoVecB(dr0rot, hist->dr0, Pos(p), rtheta);
                    DOTVP(s, dr, dr0rot);
#else
                    DOTVP(s, dr, hist->dr0);
#endif
//E
                    cosphi = s/(dr1*hist->drpq);
                    CROSSVP(pr0,hist->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVTUOMPSINCOS;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                }
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}

//E TERMINA METODO NORMAL DE BUSQUEDA (SINCOS-OMP)



#ifdef ADDONS
#include "search_include.h"
#endif

