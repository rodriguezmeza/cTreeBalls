/* ==============================================================================
!	MODULE: search_omp.c		[tpcf]                                          !
!	Written by: M.A. Rodriguez-Meza.											!
!    Starting date:    april 2023                                               !
!    Purpose: 3-point correlation function computation                          !
!	Language: C																	!
!	Use: searchcalc_normal_omp();												!
!	Major revisions:															!
!==============================================================================*/

/*
NOTA: El orden de los terminos en la busqueda del árbol de cada particula
	es importante. Diferentes ordenamientos (diferentes metodos de calculo)
	dan diferentes resultados despues de un numero grande de iteraciones.
	Por ejemplo, un calculo con 
		tpcf searchMethod=direct nbody=512
	da un resultado diferente a
		tpcf searchMethod=normal nbody=512
	debido a que el orden de aparicion de los vecinos cercanos a cada particula
	es diferente. Despues de unas 3600 iteraciones el momento angular comienza
	a mostrarse diferente en los resultados del momento angular de los dos casos.
*/

//B Set of search omp parallel routines using brute force and tree methods:
//
// search=tree-omp :: searchcalc_normal_omp(void)
// search=direct-3pcf-omp :: search_direct_3pcf_omp(void)
// search=tree-3pcf-direct-omp :: searchcalc_normal_3pcf_direct_omp(btab, nbody, ipmin, ipmax)
//
//E



// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"


// Bloque para normal (OMP)
local void normal_walktree(bodyptr, nodeptr, real, INTEGER *, INTEGER *, gdhistptr_omp);
local void sumnode(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);
// BODY3
local void sumnode_body3(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);
local void sumnode_cell(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);



//B COMIENZA METODO NORMAL DE BUSQUEDA (OMP)
// search=tree-omp
global void searchcalc_normal_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "searchcalc_normal_omp: Running ... (tree-omp) \n");

    ThreadCount();

    search_init_gd_hist();

#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,root,ipmin,ipmax)
  {
      bodyptr p;
      int m, n;
      INTEGER nbbcalcthread = 0;
      INTEGER nbccalcthread = 0;
      gdhist_omp hist;
      int ip;

      search_init_omp(&hist);

#pragma omp for nowait schedule(dynamic)
      for (p = btab + ipmin -1; p < btab + ipmax; p++) {
//B Set histograms to zero for the pivot
          for (n = 1; n <= cmd.sizeHistN; n++) {
              hist.histNSubthread[n] = 0.0;             // Affect only 3pcf evaluation
              hist.histXi2pcfthreadsub[n] = 0.0;
          }
#ifdef TPCF
          CLRM_ext_ext(hist.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
#endif
//E

#ifdef TPCF
#if NDIM == 3
              dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
              DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
              hist.drpq = rsqrt(hist.drpq2);
#endif
#endif
//E
          normal_walktree(p, ((nodeptr) root), gd.rSize, &nbbcalcthread, &nbccalcthread, &hist);
          computeBodyProperties_omp(p, nbody, &hist);

          ip = p - btab + 1;
          if (ip%cmd.stepState == 0) {
              verb_log_print(cmd.verbose_log, gd.outlog, " - Completed pivot: %ld\n", ip);
          }
    } // end do body p // end pragma omp DO_BODY p

#pragma omp critical
    {
      for (n = 1; n <= cmd.sizeHistN; n++) {
          gd.histN[n] += hist.histNthread[n];
          gd.histNSub[n] += hist.histNSubthread[n];
          gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
      }

#ifdef TPCF
      for (m=1; m<=cmd.mchebyshev+1; m++)
          ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],hist.histZetaMthread[m],cmd.sizeHistN);
#endif
        gd.nbbcalc += nbbcalcthread;
        gd.nbccalc += nbccalcthread;
    }

      search_free_omp(&hist);
  } // end pragma omp parallel

    int nn;
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
        gd.histNSub[nn] = gd.histN[nn];
        gd.histNSub[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSub[nn],1.0);
    }

//B Computation of histogram of all B-B encounters
    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN(nbody);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void normal_walktree(bodyptr p, nodeptr q, real qsize,
                           INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_omp hist)
{
    nodeptr l;
    real dr1;
    vector dr;
    int n;

    if (Update(p)) {
		if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell((nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd.options, "no-check-bin-cell")) {
                        accept_body(p, (nodeptr)q, &dr1, dr);
#ifdef LOGHIST
                        if (cmd.rminHist==0)
                            n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN);
                        else
                            n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR);
                        if (n<=cmd.sizeHistN-1 && n>=1) {
                if ( gd.deltaRV[n] < dr1 - Radius(q) && dr1 + Radius(q) < gd.deltaRV[n+1]) {
                        sumnode_cell(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                        }
                        } else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
#else // ! LOGHIST
                        real rBin;
                        n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                        rBin = cmd.rminHist + ((real)n)*gd.deltaR;
                    if ( rBin-gd.deltaR < dr1 - Radius(q) && dr1 + Radius(q) < rBin ) {
                            sumnode_cell(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                        }
#endif
                    } else {    // ! no-check-bin-cell
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                    }           // ! no-check-bin-cell
                }
            } else {
// BODY3
                if (Type(q) == BODY)
                sumnode(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                else
                    if (Type(p) == BODY3)
                sumnode_body3(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
            }
		}
	}
}

local void sumnode(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real cosphi;
    real xicosmphi;

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;

                if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

              CHEBYSHEVOMP;
//                NOCHEBYSHEVOMP;
#endif // ! TPCF
                hist->histXi2pcfthreadsub[n] += xi;
                *nbbcalcthread += 1;
                }
            }
#else // LOGHIST
            if(dr1>cmd.rminHist) {
            n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
            hist->histNthread[n] = hist->histNthread[n] + 1.;
            hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
            xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
            real s;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
#else
                cosphi = -dr[1]/dr1;
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd.verbose, gd.outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

            CHEBYSHEVOMP;
//            NOCHEBYSHEVOMP;
#endif // ! TPCF
            hist->histXi2pcfthreadsub[n] += xi;
            *nbbcalcthread += 1;
            }
#endif // LOGHIST
		}
    }
}

// BODY3
local void sumnode_body3(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real cosphi;
    real xicosmphi;

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
                hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

              CHEBYSHEVOMP;
//                NOCHEBYSHEVOMP;
#endif // ! TPCF
                hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                *nbbcalcthread += 1;
                }
            }
#else // LOGHIST
            if(dr1>cmd.rminHist) {
            n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                hist->histNthread[n] = hist->histNthread[n] + Nbb(q);
                hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
            xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
            real s;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
#else
                cosphi = -dr[1]/dr1;
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd.verbose, gd.outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

            CHEBYSHEVOMP;
//            NOCHEBYSHEVOMP;
#endif // ! TPCF
            hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
            *nbbcalcthread += 1;
            }
#endif // LOGHIST
        }
    }
}

local void sumnode_cell(bodyptr p, cellptr start, cellptr finish,
                    INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real cosphi;
    real xicosmphi;

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>0) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(q);
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nb(q);
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    CHEBYSHEVOMP;
//                  NOCHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nb(q);
                    *nbccalcthread += 1;
                }
            }
#else // LOGHIST
                n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
            hist->histNthread[n] = hist->histNthread[n] + Nb(q);
            hist->histNSubthread[n] = hist->histNSubthread[n] + Nb(q);
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                CHEBYSHEVOMP;
//            NOCHEBYSHEVOMP;
#endif // ! TPCF
            hist->histXi2pcfthreadsub[n] += xi*Nb(q);
            *nbccalcthread += 1;

#endif // LOGHIST
            }
    }
}


// search=direct-3pcf-omp       // Brute force (no harmonic approximation)
global int search_direct_3pcf_omp(bodyptr btab, int nbody,
                                  INTEGER ipmin, INTEGER ipmax)
{
    int nn, mm, ll;
    double cpustart;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "evalHistN: Running ... (direct simple, 3pcf-omp) \n");

    ThreadCount();

//B Init:
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histN[nn] = 0.0;
        gd.histXi2pcf[nn] = 0.0;
    }
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histNNN[nn] = 0.0;
        for (mm = 1; mm <= cmd.sizeHistN; mm++)
            for (ll = 1; ll <= cmd.sizeHistTheta; ll++) {
                gd.histXi3pcf[nn][mm][ll] = 0.0;
                gd.histNNNSub[nn][mm][ll] = 0.0;
            }
    }
    gd.nbbcalc = 0;
//E

#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,ipmin,ipmax)
  {
    bodyptr j1, j2, j3;
    real rr12;
    real rr13;
    int n, m, l;
    vector dr12, dr13;
    real xi1, xi2, xi3;
    real theta1, theta2, theta;
    real dr12_2, xi;
    gdhist_omp_3pcfbf hist;
    int nbbcalcthread = 0;

    search_init_omp_3pcfbf(&hist);

#pragma omp for nowait schedule(dynamic)
    DO_BODY(j1, btab, btab+nbody) {

        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;
            hist.histXi2pcfthreadsub[n] = 0.0;
        }

        for (n = 1; n <= cmd.sizeHistN; n++)
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++)
                    hist.histNNNSubthread[n][m][l] = 0.0;

        DO_BODY(j2, btab, btab+nbody) {
            if (j1 != j2) {
//B 2pcf:
                if (accept_body(j1, (nodeptr)j2, &dr12_2, dr12)) {
                    n = (int) ( (dr12_2-cmd.rminHist) * gd.i_deltaR) + 1;
                    if (n<=cmd.sizeHistN && n>=1) {
                    hist.histNthread[n] = hist.histNthread[n] + 1.;
                    hist.histNSubthread[n] = hist.histNSubthread[n] + 1.;
                    xi = Kappa(j2);
                    hist.histXi2pcfthreadsub[n] += xi;
                    }
                }
//E

                DO_BODY(j3, btab, btab+nbody) {
                    if (j1 != j3) {
                        DOTPSUBV(rr13, dr13, Pos(j1), Pos(j3));
#ifdef PERIODIC
                        VWrapAll(dr13);
                        DOTVP(rr13, dr13, dr13);
#endif
                        if (dr12_2 < gd.Rcut && rr13 < gd.RcutSq) {
                            theta1 = angle_dxdy(dr12[0], dr12[1]);
                            theta2 = angle_dxdy(dr13[0], dr13[1]);
                            theta = rabs(theta2-theta1);
                            n = (int) (rsqrt(rr12) / gd.deltaR) + 1;
                            m = (int) (rsqrt(rr13) / gd.deltaR) + 1;
                            l = (int) (theta / gd.deltaTheta) + 1;
                            if ( (n<=cmd.sizeHistN && n>=1)
                                && (m<=cmd.sizeHistN && m>=1) && (l<=cmd.sizeHistTheta && l>=1)) {
                            hist.histNNNthread[n] = hist.histNNNthread[n] + 1.;
                            hist.histNNNSubthread[n][m][l] = hist.histNNNSubthread[n][m][l] + 1.;
                            xi3 = Kappa(j3);
                            xi2 = Kappa(j2);
                            xi1 = Kappa(j1);
                            hist.histXi3pcfthread[n][m][l] += xi3 * xi2 * xi1;
                            nbbcalcthread += 1;
                            }
                        }
                    }
                } // end do body j3
            }
        } // end do body j2
        computeBodyProperties_omp_3pcfbf(j1, nbody, hist);
        for (n = 1; n <= cmd.sizeHistN; n++)
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++)
                    hist.histXi3pcfthread[n][m][l] /= MAX(hist.histNNNSubthread[n][m][l],1.0);
    } // end do body j1 // end pragma omp DO_BODY j1
#pragma omp critical
    {
        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histN[n] += hist.histNthread[n];
            gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
        }

        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histNNN[n] += hist.histNNNthread[n];
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++) {
                    gd.histNNNSub[n][m][l] += hist.histNNNSubthread[n][m][l];
                    gd.histXi3pcf[n][m][l] += hist.histXi3pcfthread[n][m][l];
                }
        }
        gd.nbbcalc += nbbcalcthread;
    }
    search_free_omp_3pcfbf(&hist);
  } // end pragma omp parallel

    int i;
    double data[cmd.sizeHistTheta];

    gsl_fft_real_wavetable * real;
    gsl_fft_real_workspace * work;

    work = gsl_fft_real_workspace_alloc (cmd.sizeHistTheta);
    real = gsl_fft_real_wavetable_alloc (cmd.sizeHistTheta);

    for (nn = 1; nn <= cmd.sizeHistN; nn++)
        for (i = 1; i <= cmd.sizeHistN; i++) {
            for (ll = 0; ll < cmd.sizeHistTheta; ll++)
                data[ll] = gd.histXi3pcf[nn][i][ll+1];
            gsl_fft_real_transform (data, 1, cmd.sizeHistTheta, real, work);
            for (mm = 1; mm <= cmd.mchebyshev+1; mm++) {
                gd.histZetaM[mm][nn][i] = data[mm-1];
            }
        }

    gsl_fft_real_wavetable_free (real);

//B Computation of histogram of all B-B encounters
    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN_3pcfbf(cmd.nbody);
    for (nn = 1; nn <= cmd.sizeHistN; nn++)
        gd.histXi2pcf[nn] /= 2.0;

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);

    return _SUCCESS_;
}


local void normal_walktree_nblist_omp(bodyptr, nodeptr, real);
local void find_nblist_omp(bodyptr, cellptr, cellptr);
local void sumnode_nblist_omp(bodyptr p, INTEGER *, gdhistptr_omp_3pcfbf hist);

#if !defined(FACTIVE)
#  define FACTIVE  0.75
#endif
 
local int actlen;
local int *activenb;
local int nblist;


//B COMIENZA METODO NORMAL DE BUSQUEDA
// search=tree-3pcf-direct-omp
global void searchcalc_normal_3pcf_direct_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;
    int nn, mm, ll;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "searchcalc_normal: Running 3pcf-direct ... (tree-omp) \n");

    ThreadCount();

//B Init:
    search_init_gd_hist();
#ifdef TPCF
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histNNN[nn] = 0.0;
        for (mm = 1; mm <= cmd.sizeHistN; mm++)
            for (ll = 1; ll <= cmd.sizeHistTheta; ll++) {
                gd.histXi3pcf[nn][mm][ll] = 0.0;
                gd.histNNNSub[nn][mm][ll] = 0.0;
            }
    }
#endif
//E
    
    actlen = FACTIVE * 216 * 512 * gd.tdepth;
    verb_log_print(cmd.verbose,gd.outlog,"searchcalc: actlen = %ld",actlen);
    activenb = (int *) allocate(actlen * sizeof(int));

#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,root,nblist,actlen,activenb,ipmin,ipmax)
  {
    bodyptr p;
    int n, m, l, ip;
    gdhist_omp_3pcfbf hist;

    INTEGER nbbcalcthread = 0;
    search_init_omp_3pcfbf(&hist);

#pragma omp for nowait schedule(dynamic)
    DO_BODY(p, btab+ipmin-1, btab+ipmax) {
//B Set histograms to zero for the pivot
        hist.ipcount = 0;
        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;               // Affect only 3pcf evaluation
            hist.histXi2pcfthreadsub[n] = 0.0;
        }
#ifdef TPCF
        for (n = 1; n <= cmd.sizeHistN; n++)
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++)
                    hist.histNNNSubthread[n][m][l] = 0.0;
#endif
//E
        nblist=0;

/*
#ifdef TPCF
#if NDIM == 3
// IS NEEDED THIS SEGMENT? CHECK!!!
        dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
        DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
        hist.drpq = rsqrt(hist.drpq2);
#endif
#endif
*/
        normal_walktree_nblist_omp(p, ((nodeptr) root), gd.rSize);
        int_piksrt(nblist, activenb);
        verb_log_print(cmd.verbose_log, gd.outlog, " - Summing nblist: %ld\n", nblist);
        sumnode_nblist_omp(p,&nbbcalcthread,  &hist);
        computeBodyProperties_omp_3pcfbf(p, ipmax-ipmin+1, hist);
#ifdef TPCF
        for (n = 1; n <= cmd.sizeHistN; n++)
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++)
                    hist.histXi3pcfthread[n][m][l] = Kappa(p)*MAX(hist.histNNNSubthread[n][m][l],1.0);
#endif
        ip = p - btab + 1;
        if (ip%cmd.stepState == 0) {
            verb_log_print(cmd.verbose_log, gd.outlog, " - Completed pivot: %ld\n", ip);
        }
    } // end do body p

#pragma omp critical
    {
        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histN[n] += hist.histNthread[n];
            gd.histNSub[n] += hist.histNSubthread[n];
            gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
        }
#ifdef TPCF
        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histNNN[n] += hist.histNNNthread[n];
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++) {
                    gd.histNNNSub[n][m][l] += hist.histNNNSubthread[n][m][l];
                    gd.histXi3pcf[n][m][l] += hist.histXi3pcfthread[n][m][l];
                }
        }
#endif
        gd.nbbcalc += nbbcalcthread;
    }

    search_free_omp_3pcfbf(&hist);

  } // end pragma omp parallel

    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
        gd.histNSub[nn] = gd.histN[nn];
        gd.histNSub[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSub[nn],1.0);
    }

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN_3pcfbf(ipmax-ipmin+1);

//B Computation of tpcf:
#ifdef TPCF
    int i;
    double data[cmd.sizeHistTheta];

    gsl_fft_real_wavetable * real;
    gsl_fft_real_workspace * work;

    work = gsl_fft_real_workspace_alloc (cmd.sizeHistTheta);
    real = gsl_fft_real_wavetable_alloc (cmd.sizeHistTheta);

    for (nn = 1; nn <= cmd.sizeHistN; nn++)
        for (i = 1; i <= cmd.sizeHistN; i++) {
            for (ll = 0; ll < cmd.sizeHistTheta; ll++)
                data[ll] = gd.histXi3pcf[nn][i][ll+1];
            gsl_fft_real_transform (data, 1, cmd.sizeHistTheta, real, work);
            for (mm = 1; mm <= cmd.mchebyshev+1; mm++) {
                gd.histZetaM[mm][nn][i] = data[mm-1];
            }
        }

    gsl_fft_real_wavetable_free (real);
#endif
//E

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void normal_walktree_nblist_omp(bodyptr p, nodeptr q, real qsize)
{
    nodeptr l;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (reject_cell((nodeptr)p, q, qsize)) {
                } else {
                    for (l = More(q); l != Next(q); l = Next(l)) {
                        normal_walktree_nblist_omp(p,l,qsize/2);
                    }
                }
            } else
                find_nblist_omp(p, ((cellptr) q),( (cellptr) q+1));
        }
    }
}

local void find_nblist_omp(bodyptr p, cellptr start, cellptr finish)
{
    cellptr q;
    real dr1;
    vector dr;
    int iq;

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
            gd.nbbcalc += 1;
            iq = (bodyptr)q-bodytab;
            activenb[nblist]=iq;
            nblist +=1;
            if (nblist > actlen)
                error("find_nblist: too many neighbors\n");
        }
    }
}

local void sumnode_nblist_omp(bodyptr p, INTEGER *nbbcalcthread, gdhistptr_omp_3pcfbf hist)
{
    bodyptr q, h;
    real dr1, dr1_h;
    vector dr, dr_h;
    int i, j;
    int n, m, l;
    real theta, theta1, theta2;

    for (i = 0; i < nblist-1; i++) {
        q = bodytab + activenb[i];
        accept_body(p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd.rminHist) {
            theta1 = angle_dxdy(dr[0], dr[1]);
            n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
            if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                hist->histXi2pcfthreadsub[n] += Kappa(q);
#ifdef TPCF
                hist->histNNNthread[n] = hist->histNNNthread[n] + 1.;
#endif
            }
        }

/*

// IS NEEDED THIS SEGMENT? CHECK!!!

#ifdef TPCF
#if NDIM == 3
//            if (NDIM == 3) {
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
//            } else
#else
                    cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#endif
*/

#ifdef TPCF
//        for (j = 0; j < nblist; j++) {
        for (j = i+1; j < nblist; j++) {
//            if (i != j) {
                h = bodytab + activenb[j];
                accept_body(p, (nodeptr)h, &dr1_h, dr_h);
            if(dr1_h>cmd.rminHist) {
//                theta1 = angle_dxdy(dr[0], dr[1]);
                theta2 = angle_dxdy(dr_h[0], dr_h[1]);
                theta = rabs(theta2-theta1);
//                n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                m = (int) ((dr1_h-cmd.rminHist) * gd.i_deltaR) + 1;
                l = (int) (theta / gd.deltaTheta) + 1;
                if ( (n<=cmd.sizeHistN && n>=1)
                    && (m<=cmd.sizeHistN && m>=1) && (l<=cmd.sizeHistTheta && l>=1)) {
//                    hist->histNNNthread[n] = hist->histNNNthread[n] + 1.;
                    hist->histNNNSubthread[n][m][l] = hist->histNNNSubthread[n][m][l] + 1.;
//                    hist->histXi3pcfthread[n][m][l] += Kappa(h) * Kappa(q) * Kappa(p);
                    hist->histXi3pcfthread[n][m][l] += Kappa(h) * Kappa(q);
                    *nbbcalcthread += 1;
                }
            }
//            }
        }
#endif
    }

//
    i = nblist;
        q = bodytab + activenb[i];
        accept_body(p, (nodeptr)q, &dr1, dr);
        if(dr1>cmd.rminHist) {
            n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
            if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                hist->histXi2pcfthreadsub[n] += Kappa(q);
            }
        }

}

#ifdef BALLS

//B COMIENZA METODO BALLS (OMP) DE BUSQUEDA

local void walktree_balls_omp(nodeptr, nodeptr, gdhistptr_omp_balls, INTEGER *);
local void walktree_balls_omp_nodes(nodeptr, nodeptr, gdhistptr_omp_balls, INTEGER *);
local void sumnodes_bb_omp(nodeptr p, nodeptr q, gdhistptr_omp_balls hist);
local void sumnodes_bc_omp(nodeptr, nodeptr, gdhistptr_omp_balls);
local void sumnodes_cc_omp(nodeptr, nodeptr, gdhistptr_omp_balls);

// search=balls-omp
global void searchcalc_balls_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;
    int nn;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "searchcalc_balls6: Running ... (balls-omp) \n");

    ThreadCount();

    search_init_gd_hist();
    gd.nsmoothcount = 0;

    verb_print(cmd.verbose,
               "\nsearchcalc_balls6: Total allocated %g MByte storage so far.\n",
               gd.bytes_tot/(1024.0*1024.0));

    if (cmd.scanLevel<3)
        error("-- too small scanLevel: %d. Use a value > 2... stopping \n",cmd.scanLevel);

#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,root,ipmin,ipmax,nodetabscanlev,rootnode)
  {
    gdhist_omp_balls hist;
    nodeptr p,q;
    int i, j;
    int n,m;
      INTEGER nsmoothcountthread = 0;

    search_init_balls_omp(&hist);

#pragma omp for nowait schedule(dynamic)
        for (i=0; i< gd.nnodescanlev; i++) {
            p = nodetabscanlev[i];
//B Set histograms to zero for the pivot
            hist.ipcount = 0;
            for (n = 1; n <= cmd.sizeHistN; n++) {
                hist.histNSubthread[n] = 0.0;             // Affect only 3pcf evaluation
                hist.histXi2pcfthreadsub[n] = 0.0;
            }
#ifdef TPCF
            CLRM_ext_ext(hist.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
#endif
//E

// Change p with a valid point on the sphere. Use a greater scanLevel
#ifdef TPCF
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
            DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
            hist.drpq = rsqrt(hist.drpq2);
#endif
#endif

#ifdef TREENODE
            walktree_balls_omp_nodes(p, (nodeptr) rootnode, &hist, &nsmoothcountthread);
#else
            if (scanopt(cmd.options, "compute-j-no-eq-i")) {
                for (j=0; j< gd.nnodescanlev; j++) {
                    q = nodetabscanlev[j];
                    walktree_balls_omp(p, q, &hist, &nsmoothcountthread);
                }
            } else {
                for (j=i; j< gd.nnodescanlev; j++) {
                    q = nodetabscanlev[j];
                    walktree_balls_omp(p, q, &hist, &nsmoothcountthread);
                }
            }
#endif

            computeBodyProperties_balls_omp((bodyptr)p, nbody, &hist);

            if (i%cmd.stepState == 0) {
                verb_log_print(cmd.verbose_log, gd.outlog,
                    " - Completed pivot node: %d\n", i);
            }
        } // end loop i
//E

#pragma omp critical
    {
      for (n = 1; n <= cmd.sizeHistN; n++) {
          gd.histN[n] += hist.histNthread[n];
          gd.histNSub[n] += hist.histNSubthread[n];
// 2pcf
          gd.histNSubXi2pcf[n] += hist.histNSubXi2pcfthread[n];
//
          gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
      }

#ifdef TPCF
      for (m=1; m<=cmd.mchebyshev+1; m++)
          ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],hist.histZetaMthread[m],cmd.sizeHistN);
#endif
        gd.nsmoothcount += nsmoothcountthread;
    }

    search_free_balls_omp(&hist);
  } // end pragma omp parallel

    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
// 2pcf
        gd.histNSubXi2pcf[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSubXi2pcf[nn],1.0);
//
    }

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN_balls(nbody);

    verb_print(cmd.verbose, "balls: nsmoothcount = %ld\n",gd.nsmoothcount);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void walktree_balls_omp(nodeptr p, nodeptr q, gdhistptr_omp_balls hist, INTEGER *nsmoothcountthread)
{
    nodeptr h,l;
    real qsize=1.0;     // Dummy
    int n;
    real dr1;
    vector dr;

    if (Type(p) == CELL && Type(q) == CELL) {
        if (!reject_cell(p, q, qsize)) {

#ifdef BUCKET
            if (Nb(p)<=gd.nsmooth[0] && Nb(q)<=gd.nsmooth[0]) {
                *nsmoothcountthread += 1;
                
                if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                    if (n==-1)
                        error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                    sumnodes_cc_omp(p, q, hist);
                } else {
                    for (h = More(p); h != Next(p); h = Next(h))
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
                }

            } else {
#endif
//B As original::
            if (!scanopt(cmd.options, "no-two-balls")) {
                if (nodes_condition(p, q)) {
                    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_cc_omp(p, q, hist);
                    } else {
                        for (h = More(p); h != Next(p); h = Next(h))
                            for (l = More(q); l != Next(q); l = Next(l))
                                walktree_balls_omp(h,l,hist, nsmoothcountthread);
                    }

//                    }

                } else {
                    for (h = More(p); h != Next(p); h = Next(h))
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
                }

            } else {
                for (h = More(p); h != Next(p); h = Next(h))
                    for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
            }
//E

#ifdef BUCKET
            }
#endif
        }
    } else {
        if (Type(p) == BODY && Type(q) == CELL) {
            if (!reject_cell(p, q, qsize)) {

#ifdef BUCKET
                if (Nb(q)<=gd.nsmooth[0]) {
                    *nsmoothcountthread += 1;
                    
                    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_bc_omp(p, q, hist);
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, nsmoothcountthread);
                    }

                } else {
#endif
                if (!scanopt(cmd.options, "no-one-balls")) {
                    if (nodes_condition(p, q)) {
                        if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                            if (n==-1)
                                error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                            sumnodes_bc_omp(p, q, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                walktree_balls_omp(p,l,hist, nsmoothcountthread);
                        }
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, nsmoothcountthread);
                    }
                } else {
                    for (l = More(q); l != Next(q); l = Next(l))
                        walktree_balls_omp(p,l,hist, nsmoothcountthread);
                }

#ifdef BUCKET
                }
#endif
            }
        }
        if (Type(p) == BODY && Type(q) == BODY) {
#ifdef DEBUG
                HIT(p) = TRUE;
                HIT(q) = TRUE;
#endif
            sumnodes_bb_omp(p, q, hist);
        }
        if (Type(p) == CELL && Type(q) == BODY) {
            if (!reject_cell(q, p, qsize)) {

#ifdef BUCKET
                if (Nb(p)<=gd.nsmooth[0]) {
                    *nsmoothcountthread += 1;

                    if (nodes_set_bin(q, p, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_bc_omp(q, p, hist);
                    } else {
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(q,l,hist, nsmoothcountthread);
                    }

                } else {
#endif
                if (!scanopt(cmd.options, "no-one-balls")) {
                    if (nodes_condition(q, p)) {
                        if (nodes_set_bin(q, p, &n, &dr1, dr)) {
                            if (n==-1)
                                error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                            sumnodes_bc_omp(q, p, hist);
                        } else {
                            for (l = More(p); l != Next(p); l = Next(l))
                                walktree_balls_omp(q,l,hist, nsmoothcountthread);
                        }
                    } else {
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(q,l,hist, nsmoothcountthread);
                    }
                } else {
                    for (l = More(p); l != Next(p); l = Next(l))
                        walktree_balls_omp(q,l,hist, nsmoothcountthread);
                }

#ifdef BUCKET
            }
#endif
            }
        }
    }
}

local void walktree_balls_omp_nodes(nodeptr p, nodeptr q, gdhistptr_omp_balls hist, INTEGER *nsmoothcountthread)
{
    nodeptr h,l;
    real qsize=1.0;     // Dummy
    int n;
    real dr1;
    vector dr;

    if (Type(p) == CELL && Type(q) == CELL) {
        if (!reject_cell(p, q, qsize)) {

#ifdef BUCKET
            if (Nb(p)<=gd.nsmooth[0] && Nb(q)<=gd.nsmooth[0]) {
                *nsmoothcountthread += 1;
                
                if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                    if (n==-1)
                        error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                    sumnodes_cc_omp(p, q, hist);
                } else {
                    for (h = More(p); h != Next(p); h = Next(h))
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
                }

            } else {
#endif
//B As original::
            if (!scanopt(cmd.options, "no-two-balls")) {
                if (nodes_condition(p, q)) {
                    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_cc_omp(p, q, hist);
                    } else {
                        for (h = More(p); h != Next(p); h = Next(h))
                            for (l = More(q); l != Next(q); l = Next(l))
                                walktree_balls_omp(h,l,hist, nsmoothcountthread);
                    }

//                    }

                } else {
                    for (h = More(p); h != Next(p); h = Next(h))
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
                }

            } else {
                for (h = More(p); h != Next(p); h = Next(h))
                    for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(h,l,hist, nsmoothcountthread);
            }
//E

#ifdef BUCKET
            }
#endif
        }
    } else {
        if (Type(p) == BODY && Type(q) == CELL) {
            if (!reject_cell(p, q, qsize)) {

#ifdef BUCKET
                if (Nb(q)<=gd.nsmooth[0]) {
                    *nsmoothcountthread += 1;
                    
                    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_bc_omp(p, q, hist);
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, nsmoothcountthread);
                    }

                } else {
#endif
                if (!scanopt(cmd.options, "no-one-balls")) {
                    if (nodes_condition(p, q)) {
                        if (nodes_set_bin(p, q, &n, &dr1, dr)) {
                            if (n==-1)
                                error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                            sumnodes_bc_omp(p, q, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                walktree_balls_omp(p,l,hist, nsmoothcountthread);
                        }
                    } else {
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, nsmoothcountthread);
                    }
                } else {
                    for (l = More(q); l != Next(q); l = Next(l))
                        walktree_balls_omp(p,l,hist, nsmoothcountthread);
                }

#ifdef BUCKET
                }
#endif
            }
        }
        if (Type(p) == BODY && Type(q) == BODY) {
#ifdef DEBUG
                HIT(p) = TRUE;
                HIT(q) = TRUE;
#endif
            sumnodes_bb_omp(p, q, hist);
        }
        if (Type(p) == CELL && Type(q) == BODY) {
            if (!reject_cell(q, p, qsize)) {

#ifdef BUCKET
                if (Nb(p)<=gd.nsmooth[0]) {
                    *nsmoothcountthread += 1;

                    if (nodes_set_bin(q, p, &n, &dr1, dr)) {
                        if (n==-1)
                            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                        sumnodes_bc_omp(q, p, hist);
                    } else {
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(q,l,hist, nsmoothcountthread);
                    }

                } else {
#endif
                if (!scanopt(cmd.options, "no-one-balls")) {
                    if (nodes_condition(q, p)) {
                        if (nodes_set_bin(q, p, &n, &dr1, dr)) {
                            if (n==-1)
                                error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
                            sumnodes_bc_omp(q, p, hist);
                        } else {
                            for (l = More(p); l != Next(p); l = Next(l))
                                walktree_balls_omp(q,l,hist, nsmoothcountthread);
                        }
                    } else {
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(q,l,hist, nsmoothcountthread);
                    }
                } else {
                    for (l = More(p); l != Next(p); l = Next(l))
                        walktree_balls_omp(q,l,hist, nsmoothcountthread);
                }

#ifdef BUCKET
            }
#endif
            }
        }
    }
}

local void sumnodes_bb_omp(nodeptr p, nodeptr q, gdhistptr_omp_balls hist)
{
    real dr1;
    vector dr;
    int n;
    int m;
    real xi, xj;
    real cosphi;
    real xicosmphi;

    if (accept_body((bodyptr)p, q, &dr1, dr)) {
        if(dr1>cmd.rminHist) {
            if (cmd.rminHist==0)
                n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;

            if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                xj = Kappa(p);
                xi = Kappa(q);
                hist->histXi2pcfthreadsub[n] += xj*xi;

#ifdef TPCF
#if NDIM == 3
                real s;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
#else
                cosphi = -dr[1]/dr1;        // x,y
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                    "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                CHEBYSHEVOMPBALLS;
#endif // ! TPCF
                gd.nbbcalc += 1;
            }
        }
    } // ! accept_body
}

local void sumnodes_bc_omp(nodeptr p, nodeptr q, gdhistptr_omp_balls hist)
{
    real dr1;
    vector dr;
    int n;
    int m;
    real xi, xj;
    real cosphi;
    real xicosmphi;

    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
        if (n==-1)
            error("\nsumnodes_bc: error in setting the bin '%d' \n",n);
        if(dr1>cmd.rminHist) {
            if (cmd.rminHist==0)
                n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
            if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + Nb(q);
                hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + Nb(q);
                hist->histNSubthread[n] = hist->histNSubthread[n] + Nb(q);
                xj = Kappa(p);
                xi = Kappa(q);
                hist->histXi2pcfthreadsub[n] += xj*xi*Nb(q);
#ifdef TPCF
#if NDIM == 3
                real s;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
#else
                cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenodes_bc: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                CHEBYSHEVOMPBALLS;
#endif // ! TPCF
                gd.nbccalc += 1;
            }
        }
    }
}

local void sumnodes_cc_omp(nodeptr p, nodeptr q, gdhistptr_omp_balls hist)
{
    real dr1;
    vector dr;
    int n;
    int m;
    real xi, xj;
    real cosphi;
    real xicosmphi;

    if (nodes_set_bin(p, q, &n, &dr1, dr)) {
        if (n==-1)
            error("\nsumnodes_cc: error in setting the bin '%d' \n",n);
        if(dr1>cmd.rminHist) {
            if (cmd.rminHist==0)
                n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
            else
                n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
            if (n<=cmd.sizeHistN && n>=1) {
                hist->histNthread[n] = hist->histNthread[n] + Nb(p)*Nb(q);
                hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + Nb(p)*Nb(q);
                hist->histNSubthread[n] = hist->histNSubthread[n] + Nb(q);
                xj = Kappa(p);
                xi = Kappa(q);
                hist->histXi2pcfthreadsub[n] += xj*Nb(p)*xi*Nb(q);
#ifdef TPCF
#if NDIM == 3
                real s;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
#else
                cosphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                        "sumenodes_cc: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                CHEBYSHEVOMPBALLS;
                gd.ncccalc += 1;
#endif // ! TPCF
            }
        }
    }
}

//BE TERMINA METODO BALLS (OMP) DE BUSQUEDA

#endif // ! BALLS


//B COMIENZA METODO NORMAL DE BUSQUEDA (SINCOS-OMP)

local void normal_walktree_sincos(bodyptr, nodeptr, real,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos_body3(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);
local void sumnode_sincos_cell(bodyptr, cellptr, cellptr,
                    INTEGER *, INTEGER *, gdhistptr_sincos_omp);

// search=tree-omp-sincos
global void searchcalc_normal_omp_sincos(bodyptr btab, int nbody,
                                         INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;
    int nn;
    int mm;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "evalHistN: Running ... (tree omp-sincos) \n");

    ThreadCount();

//B Init:: gd_hist
#ifdef TPCF
    for (mm = 1; mm <= cmd.mchebyshev+1; mm++) {
        CLRM_ext(gd.histZetaMcos[mm], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsin[mm], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsincos[mm], cmd.sizeHistN);
    }
#endif
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histN[nn] = 0.0;
        gd.histXi2pcf[nn] = 0.0;
#ifdef TPCF
        for (mm = 1; mm <= cmd.mchebyshev+1; mm++)
            gd.histXi[mm][nn] = 0.0;
#endif
    }
    gd.nbbcalc = gd.nbccalc = 0;
//E

#pragma omp parallel default(none)   shared(cmd,gd, btab, nbody, root, ipmin, ipmax)
 {
     bodyptr p;
     int m, n;
     INTEGER nbbcalcthread = 0;
     INTEGER nbccalcthread = 0;
     gdhist_sincos_omp hist;

     search_init_sincos_omp(&hist);

#pragma omp for nowait schedule(dynamic)
     for (p = btab + ipmin -1; p < btab + ipmax; p++) {
//B
        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;               // Affect only 3pcf evaluation
            hist.histXi2pcfthreadsub[n] = 0.0;
        }
#ifdef TPCF
        CLRM_ext_ext(hist.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
        CLRM_ext_ext(hist.histXithreadsin, cmd.mchebyshev+1, cmd.sizeHistN);
#if NDIM == 3
             dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
             DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
             hist.drpq = rsqrt(hist.drpq2);
#endif
#endif
//E
         normal_walktree_sincos(p, ((nodeptr) root), gd.rSize, &nbbcalcthread, &nbccalcthread, &hist);
        computeBodyProperties_sincos_omp(p, nbody, &hist);
    } // end do body p // end pragma omp DO_BODY p

#pragma omp critical
 {
    for (n = 1; n <= cmd.sizeHistN; n++) {
        gd.histN[n] += hist.histNthread[n];
        gd.histNSub[n] += hist.histNSubthread[n];
        gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
    }
#ifdef TPCF
     for (m=1; m<=cmd.mchebyshev+1; m++) {
         ADDM_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],hist.histZetaMthreadcos[m],cmd.sizeHistN);
         ADDM_ext(gd.histZetaMsin[m],gd.histZetaM[m],hist.histZetaMthreadsin[m],cmd.sizeHistN);
         ADDM_ext(gd.histZetaMsincos[m],gd.histZetaM[m],hist.histZetaMthreadsincos[m],cmd.sizeHistN);
     }
#endif
    gd.nbbcalc += nbbcalcthread;
 }

    search_free_sincos_omp(&hist);

 } // end pragma omp parallel

    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
        gd.histNSub[nn] = gd.histN[nn];
        gd.histNSub[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSub[nn],1.0);
    }

//B Computation of histogram of all B-B encounters
    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN(nbody);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void normal_walktree_sincos(bodyptr p, nodeptr q, real qsize,
                        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
/*
    nodeptr l;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell((nodeptr)p, q, qsize)) {
                    for (l = More(q); l != Next(q); l = Next(l)) {
                        normal_walktree_sincos(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                    }
                }
            } else {
                sumnode_sincos(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread,hist);
            }
        }
    }
*/

    nodeptr l;
    real dr1;
    vector dr;
    int n;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell((nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd.options, "no-check-bin-cell")) {
                        accept_body(p, (nodeptr)q, &dr1, dr);
#ifdef LOGHIST
                        if (cmd.rminHist==0)
                            n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN);
                        else
                            n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR);
                        if (n<=cmd.sizeHistN-1 && n>=1) {
                if ( gd.deltaRV[n] < dr1 - Radius(q) && dr1 + Radius(q) < gd.deltaRV[n+1]) {
                        sumnode_sincos_cell(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                        }
                        } else
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
#else // ! LOGHIST
                        real rBin;
                        n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                        rBin = cmd.rminHist + ((real)n)*gd.deltaR;
                    if ( rBin-gd.deltaR < dr1 - Radius(q) && dr1 + Radius(q) < rBin ) {
                            sumnode_sincos_cell(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                        } else {
                            for (l = More(q); l != Next(q); l = Next(l))
                                normal_walktree_sincos(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                        }
#endif
                    } else {    // ! no-check-bin-cell
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree_sincos(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                    }           // ! no-check-bin-cell
                }
            } else {
// BODY3
                if (Type(q) == BODY)
                sumnode_sincos(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                else
                    if (Type(p) == BODY3)
                sumnode_sincos_body3(p, ((cellptr) q),( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
            }
        }
    }

}

local void sumnode_sincos(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real Cheb,ChebU;
    real cosphi,sinphi;
    real xicosmphi,xisinmphi;

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;

                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s, sy;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
                    sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                    sinphi = sy/(dr1*hist->drpq);
#else
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                       "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    for (m=1; m<=cmd.mchebyshev+1; m++){
                        Cheb = rcos((real)(m-1) * racos(cosphi));
                        ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                        xicosmphi = xi * Cheb;
                        xisinmphi = xi * ChebU;
                        hist->histXithread[m][n] += xicosmphi;
                        hist->histXithreadsin[m][n] += xisinmphi;
                    }
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                hist->histNthread[n] = hist->histNthread[n] + 1.;
                hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                real s, sy;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
                sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                sinphi = sy/(dr1*hist->drpq);
#else
                cosphi = -dr[0]/dr1;
                sinphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                for (m=1; m<=cmd.mchebyshev+1; m++){
                    Cheb = rcos((real)(m-1) * racos(cosphi));
                    ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                    xicosmphi = xi * Cheb;
                    xisinmphi = xi * ChebU;
                    hist->histXithread[m][n] += xicosmphi;
                    hist->histXithreadsin[m][n] += xisinmphi;
                }
#endif // ! TPCF
                hist->histXi2pcfthreadsub[n] += xi;
                *nbbcalcthread += 1;
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}

// BODY3
local void sumnode_sincos_body3(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real Cheb,ChebU;
    real cosphi,sinphi;
    real xicosmphi,xisinmphi;

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
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s, sy;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
                    sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                    sinphi = sy/(dr1*hist->drpq);
#else
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                       "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    for (m=1; m<=cmd.mchebyshev+1; m++){
                        Cheb = rcos((real)(m-1) * racos(cosphi));
                        ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                        xicosmphi = xi * Cheb;
                        xisinmphi = xi * ChebU;
                        hist->histXithread[m][n] += xicosmphi;
                        hist->histXithreadsin[m][n] += xisinmphi;
                    }
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                    *nbbcalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                hist->histNthread[n] = hist->histNthread[n] + Nbb(q);
                hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                real s, sy;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
                sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                sinphi = sy/(dr1*hist->drpq);
#else
                cosphi = -dr[0]/dr1;
                sinphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                for (m=1; m<=cmd.mchebyshev+1; m++){
                    Cheb = rcos((real)(m-1) * racos(cosphi));
                    ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                    xicosmphi = xi * Cheb;
                    xisinmphi = xi * ChebU;
                    hist->histXithread[m][n] += xicosmphi;
                    hist->histXithreadsin[m][n] += xisinmphi;
                }
#endif // ! TPCF
                hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                *nbbcalcthread += 1;
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}

local void sumnode_sincos_cell(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
    int n;
    int m;
    real xi;
    real Cheb,ChebU;
    real cosphi,sinphi;
    real xicosmphi,xisinmphi;

    for (q = start; q < finish; q++) {
        if (accept_body(p, (nodeptr)q, &dr1, dr)) {
#ifdef LOGHIST
            if(dr1>cmd.rminHist) {
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN)) + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;

                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                    hist->histNSubthread[n] = hist->histNSubthread[n] +  Nb(q);
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s, sy;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = -s/(dr1*hist->drpq);
                    sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                    sinphi = sy/(dr1*hist->drpq);
#else
                    cosphi = -dr[0]/dr1;
                    sinphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                       "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    for (m=1; m<=cmd.mchebyshev+1; m++){
                        Cheb = rcos((real)(m-1) * racos(cosphi));
                        ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                        xicosmphi = xi * Cheb;
                        xisinmphi = xi * ChebU;
                        hist->histXithread[m][n] += xicosmphi;
                        hist->histXithreadsin[m][n] += xisinmphi;
                    }
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi* Nb(q);
                    *nbbcalcthread += 1;
                }
            }
# else // ! LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                hist->histNthread[n] = hist->histNthread[n] +  Nb(q);
                hist->histNSubthread[n] = hist->histNSubthread[n] +  Nb(q);
                xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                real s, sy;
                DOTVP(s, dr, hist->dr0);
                cosphi = -s/(dr1*hist->drpq);
                sy = rsqrt( rsqr(dr1*hist->drpq) + rsqr(s) );
                sinphi = sy/(dr1*hist->drpq);
#else
                cosphi = -dr[0]/dr1;
                sinphi = -dr[1]/dr1;
#endif
                if (rabs(cosphi)>1.0)
                    verb_log_print(cmd.verbose, gd.outlog,
                    "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                for (m=1; m<=cmd.mchebyshev+1; m++){
                    Cheb = rcos((real)(m-1) * racos(cosphi));
                    ChebU = sinphi * rsin((real)(m-2) * racos(cosphi))/rsin(racos(cosphi));
                    xicosmphi = xi * Cheb;
                    xisinmphi = xi * ChebU;
                    hist->histXithread[m][n] += xicosmphi;
                    hist->histXithreadsin[m][n] += xisinmphi;
                }
#endif // ! TPCF
                hist->histXi2pcfthreadsub[n] += xi* Nb(q);
                *nbbcalcthread += 1;
            }
#endif // ! LOGHIST
        } // ! accept_body
    }
}

//E TERMINA METODO NORMAL DE BUSQUEDA (SINCOS-OMP)






