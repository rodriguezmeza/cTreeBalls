/* ==============================================================================
!	MODULE: search_omp.c		[cTreeBalls]                                    !
!	Written by: M.A. Rodriguez-Meza.											!
!    Starting date:    april 2023                                               !
!    Purpose: 3-point correlation function computation                          !
!	Language: C																	!
!	Use: searchcalc_normal_omp();												!
!	Major revisions:															!
!==============================================================================*/
//        1          2          3          4          5          6          7

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

//B Set of search omp parallel routines using brute force
//                  and tree methods (and balls):
//
// search=tree-omp :: searchcalc_normal_omp(btab, nbody, ipmin, ipmax)
// search=tree-3pcf-direct-omp ::
//                  searchcalc_normal_3pcf_direct_omp(btab, nbody, ipmin, ipmax)
// search=balls-omp :: searchcalc_balls_omp(btab, nbody, ipmin, ipmax)
// search=tree-omp-sincos ::
//                  searchcalc_normal_omp_sincos(btab, nbody, ipmin, ipmax)
//
//E

// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

//B COMIENZA METODO NORMAL DE BUSQUEDA (OMP)

// Bloque para normal (OMP)
local void normal_walktree(bodyptr, nodeptr, real, INTEGER *, INTEGER *, gdhistptr_omp);
local void sumnode(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);
// BODY3
local void sumnode_body3(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);
local void sumnode_cell(bodyptr, cellptr, cellptr, INTEGER *, INTEGER *, gdhistptr_omp);

//ADDONS:
#ifdef ADDONSDEVELOP
#include "search_omp_00a.h"
#endif


// search=tree-omp
global void searchcalc_normal_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "searchcalc_normal_omp: Running ... (tree-omp) \n");

    ThreadCount(nbody, 0);

    search_init_gd_hist();

//B 2023.11.29
//#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,root,ipmin,ipmax)
#pragma omp parallel default(none)   shared(cmd,gd,btab,nbody,roottable,ipmin,ipmax)
  {
      bodyptr p;
      int n;
      INTEGER nbbcalcthread = 0;
      INTEGER nbccalcthread = 0;
      gdhist_omp hist;
      int ip;

      search_init_omp(&hist);

#pragma omp for nowait schedule(dynamic)
      for (p = btab + ipmin -1; p < btab + ipmax; p++) {
//B Set histograms to zero for the pivot
          for (n = 1; n <= cmd.sizeHistN; n++) {
              hist.histNSubthread[n] = 0.0;         // Affect only 3pcf evaluation
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
//E

//B 2023.11.29
//          normal_walktree(p, ((nodeptr) root), gd.rSize, &nbbcalcthread, &nbccalcthread, &hist);
          normal_walktree(p, ((nodeptr) roottable[0]), gd.rSizeTable[0], &nbbcalcthread, &nbccalcthread, &hist);
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
          gd.histNSubXi2pcf[n] += hist.histNSubXi2pcfthread[n];
          gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
      }

#ifdef TPCF
        int m;
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
        gd.histNSubXi2pcf[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSubXi2pcf[nn],1.0);
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
                    if (!scanopt(cmd.options, "no-one-ball")) {
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
                    } else {    // ! no-one-ball
                        for (l = More(q); l != Next(q); l = Next(l))
                            normal_walktree(p,l,qsize/2, nbbcalcthread, nbccalcthread, hist);
                    }           // ! no-one-ball
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
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(q);

#ifndef PTOPIVOTROTATION3
#ifdef TPCF
                    real cosphi;
#if NDIM == 3
                    real s;
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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    CHEBYSHEVOMP;
#endif // ! TPCF

#else // ! PTOPIVOTROTATION3
                    cosphi = -dr[1]/dr1;
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVOMP;
#endif // ! PTOPIVOTROTATION3
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
            }
#else // LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + 1.;
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);

                    CHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbbcalcthread += 1;
                }
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
                    real cosphi;
#if NDIM == 3
                    real s;

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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                        "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                    *nbbcalcthread += 1;
                }
            }
#else // LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ( (dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nbb(q);
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + Nbb(q);
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nbb(q);
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*Nbb(q);
                    *nbbcalcthread += 1;
                }
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
                    hist->histNthread[n] = hist->histNthread[n] + Nb(q); // This changes wrt tpcf
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.0;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.0;
                    xi = Kappa(q);
#ifdef TPCF
                    real cosphi;
#if NDIM == 3
                    real s;
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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                }
            }
#else // LOGHIST
            if(dr1>cmd.rminHist) {
                n = (int) ((dr1-cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(q); // This changes wrt tpcf
                    hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.0;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.0;
                    xi = Kappa(q);
#ifdef TPCF
#if NDIM == 3
                    real s;
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
#else // ! NDIM
                    cosphi = -dr[1]/dr1;
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
                    CHEBYSHEVOMP;
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                    *nbccalcthread += 1;
                }
            }
#endif // LOGHIST
        }
    }
}

//E TERMINA METODO NORMAL DE BUSQUEDA (OMP)


//B COMIENZA METODO NORMAL (3pcf - brute force) DE BUSQUEDA

local void normal_walktree_nblist_omp(bodyptr, nodeptr, real);
local void find_nblist_omp(bodyptr, cellptr, cellptr);
local void sumnode_nblist_omp(bodyptr p, INTEGER *, gdhistptr_omp_3pcfbf hist);

#if !defined(FACTIVE)
#  define FACTIVE  0.75
#endif
 
local int actlen;
local int *activenb;
local int nblist;


// search=tree-3pcf-direct-omp
global void searchcalc_normal_3pcf_direct_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax)
{
    double cpustart;
    int nn;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "searchcalc_normal: Running 3pcf-direct ... (tree-omp) \n");

    ThreadCount(nbody, 0);

//B Init:
    search_init_gd_hist();
#ifdef TPCF
    int mm, ll;
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
    int n, ip;
#ifdef TPCF
        int m, l;
#endif
    gdhist_omp_3pcfbf hist;

    INTEGER nbbcalcthread = 0;
    search_init_omp_3pcfbf(&hist);

#pragma omp for nowait schedule(dynamic)
    DO_BODY(p, btab+ipmin-1, btab+ipmax) {
//B Set histograms to zero for the pivot
        hist.ipcount = 0;
        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;           // Affect only 3pcf evaluation
            hist.histXi2pcfthreadsub[n] = 0.0;
        }
#ifdef TPCF
//        int m, l;
        for (n = 1; n <= cmd.sizeHistN; n++)
            for (m = 1; m <= cmd.sizeHistN; m++)
                for (l = 1; l <= cmd.sizeHistTheta; l++)
                    hist.histNNNSubthread[n][m][l] = 0.0;
#endif
//E
        nblist=0;

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
//        int m, l;
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
    bodyptr q;
    real dr1;
    vector dr;
    int i;
    int n;
    real theta1;

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

#ifdef TPCF
        bodyptr h;
        int m, j, l;
        real theta, theta2, dr1_h;
        vector dr_h;
        for (j = i+1; j < nblist; j++) {
                h = bodytab + activenb[j];
                accept_body(p, (nodeptr)h, &dr1_h, dr_h);
            if(dr1_h>cmd.rminHist) {
                theta2 = angle_dxdy(dr_h[0], dr_h[1]);
                theta = rabs(theta2-theta1);
                m = (int) ((dr1_h-cmd.rminHist) * gd.i_deltaR) + 1;
                l = (int) (theta / gd.deltaTheta) + 1;
                if ( (n<=cmd.sizeHistN && n>=1)
                    && (m<=cmd.sizeHistN && m>=1) && (l<=cmd.sizeHistTheta && l>=1)) {
                    hist->histNNNSubthread[n][m][l] = hist->histNNNSubthread[n][m][l] + 1.;
                    hist->histXi3pcfthread[n][m][l] += Kappa(h) * Kappa(q);
                    *nbbcalcthread += 1;
                }
            }
        } // ! end loop j
#endif
    } // ! end loop i

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

//E TERMINA METODO NORMAL (3pcf - brute force) DE BUSQUEDA



#ifdef BALLS

//B COMIENZA METODO BALLS (OMP) DE BUSQUEDA
local void walktree_balls_omp(nodeptr, nodeptr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *, INTEGER *);
//B 2023.11.22
local void sumnodes_bb_omp(nodeptr, nodeptr, real *dr1, vector dr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_bc_omp(nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cb_omp(nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cc_omp(nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);

//B BALLS4: TREENODEBALLS4
local void walktree_balls6_omp(nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local void walksub6_omp(nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *, INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local void sum_balls6_omp(cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6, INTEGER);
local void sumnode_balls6_omp(cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
local void sumcell_balls6_omp(cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
local void sumcellcell_balls6_omp(cellptr, cellptr, nodeptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
//E


//B BALLS4
local INTEGER ihit;
//E
local int imiss;

// search=balls-omp
//global void searchcalc_balls_omp(bodyptr btab, int nbody, INTEGER ipmin, INTEGER ipmax, int cat1, int cat2)
global void searchcalc_balls_omp(bodyptr *btable, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax, int cat1, int cat2)
{
    double cpustart;
//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
    INTEGER ibfcount=0;
#endif
#endif
    int nn;

//B BALLS4
    ihit=0;

    imiss = 0;

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "\nsearchcalc_balls: Running ... (balls-omp) \n");
#ifdef TREENODEALLBODIES
    verb_print(cmd.verbose, "treenode with all bodies... \n");
#else
#ifdef TREENODEBALLS4
    verb_print(cmd.verbose, "treenode with 2 balls using balls4 method...\n");
    verb_print(cmd.verbose,
               "finding at the same time lists of neighbour cells and bodies...\n");
#else
    verb_print(cmd.verbose, "treenode with 2 balls method...\n");
    verb_print(cmd.verbose, "using tree of nodes at two levels: root to search and lower for cell pivots...\n");
#endif // ! TREENODEBALLS4
#endif // ! TREENODEALLBODIES

#ifdef SINCOS
    verb_print(cmd.verbose, "sincos base... \n");
#endif
    if (scanopt(cmd.options, "no-two-balls"))
        verb_print(cmd.verbose, "with option no-two-balls... \n");
    if (scanopt(cmd.options, "no-one-ball"))
        verb_print(cmd.verbose, "with option no-one-ball... \n");

    if (cmd.verbose==3)
    verb_print(cmd.verbose,
        "\ncat1, cat2, nbody[cat1], nbody[cat2], ipmin, imax[cat1], ipmax[cat1]:\
        %d %d %ld %ld %ld %ld %ld\n\n",
        cat1,cat2,nbody[cat1],nbody[cat2],ipmin,ipmax[cat1],ipmax[cat2]);
    ThreadCount(nbody[cat1], cat1);

// Here we clear: histZetaM, histXi, histN, histNSubXi2pcf, histXi2pcf...
    search_init_gd_hist();
#ifdef SINCOS
    search_init_gd_hist_sincos();
#endif

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
    bodytabbf = (bodyptr) allocate(nbody[cat1] * sizeof(body));
    gd.bytes_tot += nbody[cat1]*sizeof(body);
    fprintf(stdout,"\n\nAllocated %g MByte for particle (found) storage.\n",
            nbody[cat1]*sizeof(body)/(1024.0*1024.0));

    verb_print(cmd.verbose,
               "\nsearchcalc_balls_omp: Total allocated %g MByte storage so far.\n",
               gd.bytes_tot/(1024.0*1024.0));
#endif
#endif

    gd.nsmoothcount = 0;

    verb_print(cmd.verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd.bytes_tot/(1024.0*1024.0));

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
#pragma omp parallel default(none)  \
    shared(cmd,gd,btable,nbody,root,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, \
    bodytabbf,ibfcount,ihit,cat1, cat2)
#else
#pragma omp parallel default(none) \
    shared(cmd,gd,btable,nbody,root,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, cat1, cat2)
#endif
#else // ! TREENODEBALLS4
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,root,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, cat1, cat2)
#endif // ! TREENODEBALLS4
  {
    int i;
    int n;
    INTEGER nbbcalcthread = 0;
    INTEGER nbccalcthread = 0;
    INTEGER ncccalcthread = 0;
//#ifdef DEBUG
    INTEGER nsmoothcountthread = 0;
//#endif

//B BALLS4
#ifdef TREENODEBALLS4
      INTEGER ibfcountthread = 0;
      gdhist_omp_balls6 histb;
      gdhist_omp_balls6 histccb;
      gdhist_sincos_omp_balls6 histbsincos;

      search_init_omp_balls6(&histb, cat1);
      search_init_omp_balls6_cc(&histccb);
      search_init_sincos_omp_balls6(&histbsincos);
#endif
//E

    gdhist_omp_balls hist;
    gdhist_omp_balls histcc;
//#ifdef SINCOS
      gdhist_sincos_omp histsincos;
      search_init_sincos_omp(&histsincos);
//#endif

// Here we clear threads of: histZetaM, histXi, histNSub, histN, histNSubXi2pcf, histXi2pcf, histXi2pcf()sub...
    search_init_balls_omp(&hist);
    search_init_balls_omp_cc(&histcc);

#ifdef TREENODEALLBODIES
    bodyptr p;
#else
    nodeptr p;
#endif // ! TREENODEALLBODIES


#ifdef PIVOTEXTERNAL
//B (1)
// Setting reference axis here instead in the loop i,
//  does affect the 3pcf results for m=2 or higher
#ifdef TPCF
#if NDIM == 3
#ifdef TREENODEALLBODIES
    p = (bodyptr)nodetablescanlev[cat1][0];
#else
    p = nodetablescanlev[cat1][0];
#endif
#ifdef SINCOS
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histsincos.q0);
    DOTPSUBV(histsincos.drpq2, histsincos.dr0, Pos(p), histsincos.q0);
    histsincos.drpq = rsqrt(histsincos.drpq2);
#ifdef PTOPIVOTROTATION
    real rtheta;
    vector dr0rot;
    rtheta = xrandom(0.0, TWOPI);
    RotationVecAWRtoVecB(dr0rot, histsincos.dr0, Pos(p), rtheta);
    SETV(histsincos.dr0, dr0rot);
#endif
#else // ! SINCOS
    dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histcc.q0);
    DOTPSUBV(histcc.drpq2, histcc.dr0, Pos(p), histcc.q0);
    histcc.drpq = rsqrt(histcc.drpq2);
#ifdef PTOPIVOTROTATION
    real rtheta;
    vector dr0rot;
    rtheta = xrandom(0.0, TWOPI);
    RotationVecAWRtoVecB(dr0rot, histcc.dr0, Pos(p), rtheta);
    SETV(histcc.dr0, dr0rot);
#endif
#endif // // ! SINCOS
#endif // ! NDIM
#endif // ! TPCF
//E
#endif // ! PIVOTEXTERNAL

#pragma omp for nowait schedule(dynamic)
#ifdef TREENODEALLBODIES
        for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
#else // TREENODEALLBODIES exclude TREENODEBALLS4. Then this is for TREENODEBALLS4
        for (i=0; i< gd.nnodescanlevTable[cat1]; i++) {
            if (cmd.scanLevel==0)
                p = (nodeptr) roottable[cat1];
            else
                p = nodetablescanlev[cat1][i];
#endif
            hist.ipcount = 0;

#ifndef PIVOTEXTERNAL
//B (2)
// Set histograms to zero for the pivot
#ifdef TPCF
#ifdef SINCOS
            for (n = 1; n <= cmd.sizeHistN; n++) {
                histsincos.histNSubthread[n] = 0.0;
            }
            CLRM_ext_ext(histsincos.histXithreadcos, cmd.mchebyshev+1,
                         cmd.sizeHistN);
            CLRM_ext_ext(histsincos.histXithreadsin, cmd.mchebyshev+1,
                         cmd.sizeHistN);

#ifdef TREENODEBALLS4
            for (n = 1; n <= cmd.sizeHistN; n++) {
                histbsincos.histNSubthread[n] = 0.0;
            }
            CLRM_ext_ext(histbsincos.histXithreadcos, cmd.mchebyshev+1,
                         cmd.sizeHistN);
            CLRM_ext_ext(histbsincos.histXithreadsin, cmd.mchebyshev+1,
                         cmd.sizeHistN);
#endif

#else // ! SINCOS
            for (n = 1; n <= cmd.sizeHistN; n++) {
                histcc.histNSubthread[n] = 0.0;
//B BALLS4
#ifdef TREENODEBALLS4
                histb.histNSubthread[n] = 0.0;      // Affect only 3pcf evaluation
                histb.histXi2pcfthreadsub[n] = 0.0;
//B two-balls
                histccb.histNSubthread[n] = 0.0;    // Affect only 3pcf evaluation
                histccb.histXi2pcfthreadsub[n] = 0.0;
//E
#endif // ! TREENODEBALLS4
            }
            CLRM_ext_ext(histcc.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
#ifdef TREENODEBALLS4
            CLRM_ext_ext(histb.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
            CLRM_ext_ext(histccb.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
//E
#endif
#endif // ! SINCOS
#endif // ! TPCF
//E
#endif // ! PIVOTEXTERNAL

#ifndef PIVOTEXTERNAL
//B (3)
// Change p with a valid point on the sphere. Use a greater scanLevel
#ifdef TPCF
#if NDIM == 3
#ifdef SINCOS
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histsincos.q0);
            DOTPSUBV(histsincos.drpq2, histsincos.dr0, Pos(p), histsincos.q0);
            histsincos.drpq = rsqrt(histsincos.drpq2);
#ifdef PTOPIVOTROTATION
            real rtheta;
            vector dr0rot;
            rtheta = xrandom(0.0, TWOPI);
            RotationVecAWRtoVecB(dr0rot, histsincos.dr0, Pos(p), rtheta);
            SETV(histsincos.dr0, dr0rot);
#endif
#else // ! SINCOS
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histcc.q0);
            DOTPSUBV(histcc.drpq2, histcc.dr0, Pos(p), histcc.q0);
            histcc.drpq = rsqrt(histcc.drpq2);
#ifdef PTOPIVOTROTATION
          real rtheta;
          vector dr0rot;
          rtheta = xrandom(0.0, TWOPI);
          RotationVecAWRtoVecB(dr0rot, histcc.dr0, Pos(p), rtheta);
          SETV(histcc.dr0, dr0rot);
#endif
            
#ifdef TREENODEBALLS4
#if NDIM == 3
            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histb.q0);
            DOTPSUBV(histb.drpq2, histb.dr0, Pos(p), histb.q0);
            histb.drpq = rsqrt(histb.drpq2);

            dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, histccb.q0);
            DOTPSUBV(histccb.drpq2, histccb.dr0, Pos(p), histccb.q0);
            histccb.drpq = rsqrt(histccb.drpq2);
#endif // ! NDIM
#endif // ! TREENODEBALLS4

#endif // // ! SINCOS
#endif // ! NDIM
#endif // ! TPCF
//E
#endif // ! PIVOTEXTERNAL

#ifdef TREENODEALLBODIES
            nodeptr q;
            int j;
            if (cmd.scanLevelRoot==0) {
                walktree_balls_omp((nodeptr)p, (nodeptr) roottable[cat2],
                                   &hist, &histcc, &histsincos,
                                   &nsmoothcountthread, &nbbcalcthread,
                                   &nbccalcthread, &ncccalcthread);
            } else {
                for (j=0; j< gd.nnodescanlev_rootTable[cat2]; j++) {
                    q = nodetablescanlev_root[cat2][j];
                    walktree_balls_omp((nodeptr)p, q, &hist, &histcc, &histsincos,
                            &nsmoothcountthread, &nbbcalcthread,
                                       &nbccalcthread, &ncccalcthread);
                }
            }
//E
#else // ! TREENODEALLBODIES
#ifdef TREENODEBALLS4
            nodeptr q;
            int j;
// This method is faster when scanLevel is higher...
            histb.active[0] = (nodeptr) (roottable[cat2]);
//            histb.active[0] = (nodeptr) root;
//            for (j=0; j< gd.nnodescanlev_rootTable[cat2]; j++) {
//               q = nodetablescanlev_root[cat2][j];
//                histb.active[0] = q;

//                verb_print_debug(1, "\nAqui voy (0): %g\n", Kappa(histb.active[0]));
//            q = (nodeptr) (roottable[cat2]);
//                histb.active[0] = q;
//                verb_print_debug(1, "\nAqui voy (1): %g %g %g\n",Pos(p)[0],Pos(p)[1],Pos(p)[2]);


            walktree_balls6_omp(histb.active, histb.active + 1,
                    histb.interact, histb.interact + histb.actlen,
                    p, Size(p), Pos(p),
                    &nbbcalcthread, &nbccalcthread, &ncccalcthread,
                    &histb, &histccb, &histbsincos,
                    &ibfcountthread, &nsmoothcountthread, nbody[cat1]);
//            }

            if (!scanopt(cmd.options, "no-two-balls")) {
#ifndef SINCOS
                computeBodyProperties_omp_balls6_cc((bodyptr)p,
                                                    nbody[cat1], &histccb);
#else
                computeBodyProperties_sincos_omp_balls6_cc((bodyptr)p,
                                                    nbody[cat1], &histbsincos);
#endif
            }
#else // ! TREENODEBALLS4
            nodeptr q;
            int j;
//B Faster if first option is used:
            if (cmd.scanLevelRoot==0)
                walktree_balls_omp(p, (nodeptr) roottable[cat2], &hist, &histcc,
                                   &histsincos, &nsmoothcountthread,
                                   &nbbcalcthread, &nbccalcthread, &ncccalcthread);
//E Above one.
            else {
// Must be the same node's trees
                if (gd.nnodescanlevTable[cat1] == gd.nnodescanlev_rootTable[cat2])
                    for (j=i; j< gd.nnodescanlevTable[cat2]; j++) {
                        q = nodetablescanlev_root[cat2][j];
                        walktree_balls_omp(p, q, &hist, &histcc, &histsincos,
                                &nsmoothcountthread, &nbbcalcthread, &nbccalcthread, &ncccalcthread);
                    }
                else
                for (j=0; j< gd.nnodescanlev_rootTable[cat2]; j++) {
                q = nodetablescanlev_root[cat2][j];
                walktree_balls_omp(p, q, &hist, &histcc, &histsincos,
                    &nsmoothcountthread, &nbbcalcthread, &nbccalcthread, &ncccalcthread);
                }
            }
#endif // ! TREENODEBALLS4
//E
#endif // ! TREENODEALLBODIES

#ifdef TPCF
#ifndef PIVOTEXTERNAL
//B (4)
#ifdef SINCOS
            computeBodyProperties_balls_omp_cc_sincos((bodyptr)p, nbody[cat1], &histsincos);
#else
            computeBodyProperties_balls_omp_cc((bodyptr)p, nbody[cat1], &histcc);
#endif
//E
#endif // ! PIVOTEXTERNAL
#endif

#ifndef TREENODEBALLS4
#ifdef TREENODEALLBODIES
            i = p - btable[cat1] + 1;
            if (i%cmd.stepState == 0) {
                verb_log_print(cmd.verbose_log, gd.outlog,
                    " - Completed pivot node: %d\n", i);
            }
#else
            if (i%cmd.stepState == 0) {
                verb_log_print(cmd.verbose_log, gd.outlog,
                    " - Completed pivot node: %d\n", i);
            }
#endif
#endif
        } // end loop i

#ifndef TREENODEBALLS4
      computeBodyProperties_balls_omp((bodyptr)p, nbody[cat1], &hist);
#endif

#ifdef TPCF
#ifdef PIVOTEXTERNAL
//B (5)
#ifdef SINCOS
            computeBodyProperties_balls_omp_cc_sincos((bodyptr)p, nbody[cat1], &histsincos);
#else
            computeBodyProperties_balls_omp_cc((bodyptr)p, nbody[cat1], &histcc);
#endif
//E
#endif
#endif

#pragma omp critical
    {
#ifdef TREENODEBALLS4
        int n;
        for (n = 1; n <= cmd.sizeHistN; n++) {
            gd.histN[n] += histb.histNthread[n];
            gd.histNSub[n] += histb.histNSubthread[n];
            gd.histNSubXi2pcf[n] += histb.histNSubXi2pcfthread[n];
            gd.histXi2pcf[n] += histb.histXi2pcfthread[n];
            if (!scanopt(cmd.options, "no-two-balls")) {
                gd.histN[n] += histccb.histNthread[n];
                gd.histNSub[n] += histccb.histNSubthread[n];
                gd.histNSubXi2pcf[n] += histccb.histNSubXi2pcfthread[n];
                gd.histXi2pcf[n] += histccb.histXi2pcfthread[n];
            }
        }
#else
      for (n = 1; n <= cmd.sizeHistN; n++) {
          gd.histN[n] += hist.histNthread[n];
          gd.histNSubXi2pcf[n] += hist.histNSubXi2pcfthread[n];
          gd.histXi2pcf[n] += hist.histXi2pcfthread[n];
      }
#endif // ! TREENODEBALLS4

#ifdef TPCF
#ifdef TREENODEBALLS4
#ifdef SINCOS
        int m;
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        ADDM_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],
                 histbsincos.histZetaMthreadcos[m],cmd.sizeHistN);
        ADDM_ext(gd.histZetaMsin[m],gd.histZetaMsin[m],
                 histbsincos.histZetaMthreadsin[m],cmd.sizeHistN);
        ADDM_ext(gd.histZetaMsincos[m],gd.histZetaMsincos[m],
                 histbsincos.histZetaMthreadsincos[m],cmd.sizeHistN);
//        if (!scanopt(cmd.options, "no-two-balls")) {
//            ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],histccb.histZetaMthread[m],
//                     cmd.sizeHistN);
//        }
    }
#else // ! SINCOS
        int m;
        for (m=1; m<=cmd.mchebyshev+1; m++) {
            ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],histb.histZetaMthread[m],
                     cmd.sizeHistN);
            if (!scanopt(cmd.options, "no-two-balls")) {
                ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],histccb.histZetaMthread[m],
                         cmd.sizeHistN);
            }
        }
#endif // ! SINCOS
#else // ! TREENODEBALLS4
#ifdef SINCOS
        int m;
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        ADDM_ext(gd.histZetaMcos[m],gd.histZetaMcos[m],
                 histsincos.histZetaMthreadcos[m],cmd.sizeHistN);
        ADDM_ext(gd.histZetaMsin[m],gd.histZetaMsin[m],
                 histsincos.histZetaMthreadsin[m],cmd.sizeHistN);
        ADDM_ext(gd.histZetaMsincos[m],gd.histZetaMsincos[m],
                 histsincos.histZetaMthreadsincos[m],cmd.sizeHistN);
    }
#else
        int m;
    for (m=1; m<=cmd.mchebyshev+1; m++)
        ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],histcc.histZetaMthread[m],
                 cmd.sizeHistN);
#endif // ! SINCOS
#endif // ! TREENODEBALLS4
#endif // ! TPCF
        gd.nsmoothcount += nsmoothcountthread;
        gd.nbbcalc += nbbcalcthread;
        gd.nbccalc += nbccalcthread;
        gd.ncccalc += ncccalcthread;
#ifdef TREENODEBALLS4
#ifdef DEBUG
        ibfcount += ibfcountthread;
#endif
#endif
    }

//B BALLS4
//#ifdef SINCOS
    search_free_sincos_omp(&histsincos);
//#endif
    search_free_balls_omp_cc(&histcc);
    search_free_balls_omp(&hist);
#ifdef TREENODEBALLS4
    search_free_sincos_omp_balls6(&histbsincos);
    search_free_omp_balls6_cc(&histccb);
    search_free_omp_balls6(&histb);
#endif
  } // end pragma omp parallel

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
    gd.nbodybf = cmd.ntosave;
    verb_print(cmd.verbose, "\nsearchcalc_balls_omp: Total bodies found: %ld %ld\n",ibfcount, gd.nbodybf);
#endif
#endif

    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
        gd.histNSubXi2pcf[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSubXi2pcf[nn],1.0);
    }

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN_balls(nbody[cat1]);

    verb_print(cmd.verbose, "balls: nsmoothcount = %ld\n",gd.nsmoothcount);

    verb_print(cmd.verbose, "balls: imiss = %d\n",imiss);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

#ifdef TREENODEALLBODIES

local void walktree_balls_omp(nodeptr p, nodeptr q,
        gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc,
        gdhistptr_sincos_omp histccsincos, INTEGER *nsmoothcountthread,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
 nodeptr l;
 real dr1;
 vector dr;

    if (!reject_balls(p, q, &dr1, dr)) {
        if (Type(q) == CELL) {
            if (!scanopt(cmd.options, "no-one-ball")) {

                if ( ((Nb(q)<=gd.nsmooth[0]) || (Size(q)<=gd.rminCell[0]))
                    && (dr1 > gd.rminCell[1]) ) {
                    *nsmoothcountthread += 1;
                    sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                 nbbcalcthread, nbccalcthread, ncccalcthread);
                } else {
                    
                    if ( Radius(q)/(dr1) < gd.deltaR)
                        sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc,
                                        histccsincos, nbbcalcthread,
                                        nbccalcthread, ncccalcthread);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, histcc, histccsincos,
                                               nsmoothcountthread, nbbcalcthread,
                                               nbccalcthread, ncccalcthread);
                }

            } else {    // ! no-one-ball
                for (l = More(q); l != Next(q); l = Next(l))
                    walktree_balls_omp(p,l,hist, histcc, histccsincos,
                            nsmoothcountthread, nbbcalcthread,
                            nbccalcthread, ncccalcthread);
                 }           // ! no-one-ball
        } else { // ! Type(q) == CELL
            sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                    nbbcalcthread, nbccalcthread, ncccalcthread);
             } // ! Type(q) == CELL
    } // ! reject_cell
}

#else // ! TREENODEALLBODIES

#define CELLCELL 1
#define CELLBODY 2
#define BODYCELL 3
#define BODYBODY 4

local void walktree_balls_omp(nodeptr p, nodeptr q,
        gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc,
        gdhistptr_sincos_omp histccsincos, INTEGER *nsmoothcountthread,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
    nodeptr l;
    nodeptr h;
    real dr1;
    vector dr;

    int SWITCHVALUE=0;
    if (Type(p) == CELL && Type(q) == CELL) SWITCHVALUE = 1;
    if (Type(p) == CELL && Type(q) == BODY) SWITCHVALUE = 2;
    if (Type(p) == BODY && Type(q) == CELL) SWITCHVALUE = 3;
    if (Type(p) == BODY && Type(q) == BODY) SWITCHVALUE = 4;

//B 2023.11.22
    if (!reject_balls(p, q, &dr1, dr)) {
//        if (dr1 > cmd.rminHist) {                   // 2024.12.04
        switch (SWITCHVALUE){
            case CELLCELL:
                 if (!scanopt(cmd.options, "no-two-ball")) {

                     if ( ( (Size(p)<=gd.rminCell[0] && Size(q)<=gd.rminCell[0])
                           || (Nb(p)<=gd.nsmooth[0] && Nb(q)<=gd.nsmooth[0]) )
                         && (dr1 > gd.rminCell[1]) ) {
//                         && (dr1 > gd.scanLevelMin[1]*cmd.rangeN) ) {
                         *nsmoothcountthread += 1;
                         sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                  nbbcalcthread, nbccalcthread, ncccalcthread);
                     } else {

                     if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if ( (Size(p)<=gd.rminCell[0] && Size(q)<=gd.rminCell[0])
                                 && (Nb(p)<=gd.nsmooth[0]
                                     && Nb(q)<=gd.nsmooth[0]) ) {
                             *nsmoothcountthread += 1;
                             sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc,
                                             histccsincos, nbbcalcthread,
                                             nbccalcthread, ncccalcthread);
                        } else
                             sumnodes_cc_omp(p, q, &dr1, dr, hist, histcc,
                                             histccsincos, nbbcalcthread,
                                             nbccalcthread, ncccalcthread);
                     } else // ! nodes_condition
                         for (h = More(p); h != Next(p); h = Next(h))
                         for (l = More(q); l != Next(q); l = Next(l))
                             walktree_balls_omp(h, l, hist, histcc, histccsincos,
                              nsmoothcountthread, nbbcalcthread, nbccalcthread, ncccalcthread);

                     }

                 } else // ! no-two-ball
                     for (h = More(p); h != Next(p); h = Next(h))
                     for (l = More(q); l != Next(q); l = Next(l))
                             walktree_balls_omp(h, l, hist, histcc, histccsincos, nsmoothcountthread,
                                     nbbcalcthread, nbccalcthread, ncccalcthread);
                break;
            case CELLBODY:
                if (!scanopt(cmd.options, "no-one-ball")) {

                    if ( ((Nb(p)<=gd.nsmooth[0]) || (Size(p)<=gd.rminCell[0]))
                        && (dr1 > gd.rminCell[1]) ) {
//                        && (dr1 > gd.scanLevelMin[1]*cmd.rangeN) ) {
                        *nsmoothcountthread += 1;
                        sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                 nbbcalcthread, nbccalcthread, ncccalcthread);
                    } else {

                    if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if (Size(p)<=gd.rminCell[0] || Nb(p)<=gd.nsmooth[0]) {
                            *nsmoothcountthread += 1;
                            sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                     nbbcalcthread, nbccalcthread, ncccalcthread);
                        } else
                            sumnodes_cb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                         nbbcalcthread, nbccalcthread, ncccalcthread);
                    } else // ! nodes_condition
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(l, q, hist, histcc, histccsincos,
                             nsmoothcountthread, nbbcalcthread, nbccalcthread, ncccalcthread);
                    }

                } else // ! no-one-ball
                    for (l = More(p); l != Next(p); l = Next(l)) {
//                        if (Type(l) == CELL) {
                            walktree_balls_omp(l, q, hist, histcc, histccsincos, nsmoothcountthread,
                                    nbbcalcthread, nbccalcthread, ncccalcthread);
//                        } else // Found BODYBODY
//                            sumnodes_bb_omp(l, q, &dr1, dr, hist, histcc, histccsincos,
//                                            nbbcalcthread, nbccalcthread, ncccalcthread);
                    }
                break;
            case BODYCELL:
                if (!scanopt(cmd.options, "no-one-ball")) {

                    if ( ((Nb(q)<=gd.nsmooth[0]) || (Size(q)<=gd.rminCell[0]))
                        && (dr1 > gd.rminCell[1]) ) {
//                        && (dr1 > gd.scanLevelMin[1]) ) {
//                        && (dr1 > gd.scanLevelMin[1]*cmd.rangeN) ) {
                        *nsmoothcountthread += 1;
//                      sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
//                                 nbbcalcthread, nbccalcthread, ncccalcthread);
                        sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                     nbbcalcthread, nbccalcthread, ncccalcthread);
                    } else {

                    if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if (Size(q)<=gd.rminCell[0] || Nb(q)<=gd.nsmooth[0]) {
                            *nsmoothcountthread += 1;
//                       sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
//                                     nbbcalcthread, nbccalcthread, ncccalcthread);
                            sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc,
                                            histccsincos, nbbcalcthread,
                                            nbccalcthread, ncccalcthread);
                        } else
                            sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc,
                                            histccsincos, nbbcalcthread,
                                            nbccalcthread, ncccalcthread);
                    } else // ! nodes_condition
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(p,l,hist, histcc, histccsincos,
                                               nsmoothcountthread, nbbcalcthread,
                                               nbccalcthread, ncccalcthread);

                    }

                } else // ! no-one-ball
//B Gives consistence results with only walktree-balls_omp
                    for (l = More(q); l != Next(q); l = Next(l)) {
//                        if (Type(l) == CELL) {
                            walktree_balls_omp(p,l,hist, histcc, histccsincos, nsmoothcountthread,
                                    nbbcalcthread, nbccalcthread, ncccalcthread);
//                        } else // Found BODYBODY
//                            sumnodes_bb_omp(p, l, &dr1, dr, hist, histcc, histccsincos,
//                                            nbbcalcthread, nbccalcthread, ncccalcthread);
                    }
                break;
            case BODYBODY: // Found BODYBODY
#ifdef DEBUG
                HIT(p) = TRUE;
                HIT(q) = TRUE;
#endif
                sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                    nbbcalcthread, nbccalcthread, ncccalcthread);
                break;
        } // ! switch
//        } // ! dr1 > rminHist
    } // ! reject_cell
}

#undef CELLCELL
#undef CELLBODY
#undef BODYCELL
#undef BODYBODY

#endif // ! TREENODEALLBODIES

local void sumnodes_bb_omp(nodeptr p, nodeptr q, real *dr1, vector dr,
            gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc, gdhistptr_sincos_omp histccsincos,
            INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
    int n;
    real xi, xj;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif

    if (nodes_set_bin(p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(p)*Nb(q);
        hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
        xj = Kappa(p);
        xi = Kappa(q);
        hist->histXi2pcfthreadsub[n] += xj*xi;
#ifndef PTOPIVOTROTATION3
#ifdef TPCF
        real cosphi;
#ifdef SINCOS
        real sinphi;
        histccsincos->histNSubthread[n] = histccsincos->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
//B
#if NDIM == 3
        real s, sy;
        vector pr0;
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        CROSSVP(pr0,histccsincos->dr0,Pos(p));
        DOTVP(sy, dr, pr0);
        sinphi = rsqrt(1.0 - rsqr(cosphi));;
        if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#endif // ! NDIM
//E

#else // ! NDIMMAKE

//B
        real s, sy;
        vector pr0;
        if (gd.dimension == 3) {
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        CROSSVP_3D(pr0,histccsincos->dr0,Pos(p));
        DOTVP(sy, dr, pr0);
        sinphi = rsqrt(1.0 - rsqr(cosphi));;
        if (sy < 0) sinphi *= -1.0;
        } else {
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
        }
//E
#endif // ! NDIMMAKE

#else // ! SINCOS
        histcc->histNSubthread[n] = histcc->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
#if NDIM == 3
        real s;
        DOTVP(s, dr, histcc->dr0);
        cosphi = s/((*dr1)*histcc->drpq);
#else
        cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#else // ! NDIMMAKE
        real s;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
        } else {
            cosphi = -dr[1]/(*dr1);        // x,y
        }
#endif // ! NDIMMAKE


#endif // ! SINCOS
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                    "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
//#ifdef TREENODEALLBODIES
#ifdef SINCOS
        CHEBYSHEVTUOMP
#else
        CHEBYSHEVOMPBALLSCC;
#endif
//#else // ! TREENODEALLBODIES
//#ifdef SINCOS
//        CHEBYSHEVTUOMP
//#else
//        CHEBYSHEVOMPBALLSCC;
//#endif
//#endif // ! TREENODEALLBODIES
#endif // ! TPCF
#else // ! PTOPIVOTROTATION3
#ifdef TPCF
#if SINCOS
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#else
        cosphi = -dr[1]/(*dr1);
#endif
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
        CHEBYSHEVOMPBALLSCC;
#endif // ! TPCF
#endif // ! PTOPIVOTROTATION3
        *nbbcalcthread += 1;
    } // ! accept_body
}

local void sumnodes_bc_omp(nodeptr p, nodeptr q, real *dr1, vector dr,
            gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc, gdhistptr_sincos_omp histccsincos,
            INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
    int n;
    real xi, xj;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif

    if (nodes_set_bin(p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(q);
        hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
        xj = Kappa(p);
        xi = Kappa(q);
        hist->histXi2pcfthreadsub[n] += xj*xi;
#ifndef PTOPIVOTROTATION3
#ifdef TPCF
        real cosphi;
#ifdef SINCOS
        real sinphi;
        histccsincos->histNSubthread[n] = histccsincos->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
//B
#if NDIM == 3
        real s, sy;
        vector pr0;
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        CROSSVP_3D(pr0,histccsincos->dr0,Pos(p));
        DOTVP(sy, dr, pr0);
        sinphi = rsqrt(1.0 - rsqr(cosphi));;
        if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#endif // ! NDIM
//E

#else // ! NDIMMAKE
//B
        real s, sy;
        vector pr0;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histccsincos->dr0);
            cosphi = s/((*dr1)*histccsincos->drpq);
            CROSSVP_3D(pr0,histccsincos->dr0,Pos(p));
            DOTVP(sy, dr, pr0);
            sinphi = rsqrt(1.0 - rsqr(cosphi));;
            if (sy < 0) sinphi *= -1.0;
            } else {
                cosphi = -dr[0]/(*dr1);
                sinphi = -dr[1]/(*dr1);
            }
//E
#endif // ! NDIMMAKE

    
#else // ! SINCOS
        histcc->histNSubthread[n] = histcc->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
#if NDIM == 3
        real s;
        DOTVP(s, dr, histcc->dr0);
        cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
        cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
//E
#else // ! NDIMMAKE
        real s;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
        } else {
            cosphi = -dr[1]/(*dr1);        // x,y
        }
#endif // ! NDIMMAKE

#endif // ! SINCOS
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
//#ifdef TREENODEALLBODIES
#ifdef SINCOS
        CHEBYSHEVTUOMP
#else
        CHEBYSHEVOMPBALLSCC;
#endif
//#else // ! TREENODEALLBODIES
//#ifdef SINCOS
//        CHEBYSHEVTUOMP
//#else
//        CHEBYSHEVOMPBALLSCC;
//#endif
//#endif // ! TREENODEALLBODIES
#endif // ! TPCF
#else // ! PTOPIVOTROTATION3
#ifdef TPCF
#if SINCOS
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#else
        cosphi = -dr[1]/(*dr1);
#endif
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
        CHEBYSHEVOMPBALLSCC;
#endif // ! TPCF
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cb_omp(nodeptr q, nodeptr p, real *dr1, vector dr,
            gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc, gdhistptr_sincos_omp histccsincos,
            INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
    int n;
    real xi, xj;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif

    if (nodes_set_bin(p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(q);
        hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
        xj = Kappa(q);
        xi = Kappa(p);
        hist->histXi2pcfthreadsub[n] += xj*xi;
#ifndef PTOPIVOTROTATION3
#ifdef TPCF
        real cosphi;
#ifdef SINCOS
        real sinphi;
        histccsincos->histNSubthread[n] = histccsincos->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
//B
#if NDIM == 3
        real s, sy;
        vector pr0;
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        CROSSVP(pr0,histccsincos->dr0,Pos(q));
        DOTVP(sy, dr, pr0);
        sinphi = rsqrt(1.0 - rsqr(cosphi));;
        if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#endif // ! NDIM
//E
#else // ! NDIMMAKE
//B
        real s, sy;
        vector pr0;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histccsincos->dr0);
            cosphi = s/((*dr1)*histccsincos->drpq);
            CROSSVP_3D(pr0,histccsincos->dr0,Pos(q));
            DOTVP(sy, dr, pr0);
            sinphi = rsqrt(1.0 - rsqr(cosphi));;
            if (sy < 0) sinphi *= -1.0;
        } else {
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
        }
//E
#endif // ! NDIMMAKE


#else // ! SINCOS
        histcc->histNSubthread[n] = histcc->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
#if NDIM == 3
        real s;
        DOTVP(s, dr, histcc->dr0);
        cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
        cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#else // ! NDIMMAKE
        real s;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
        } else {
            cosphi = -dr[1]/(*dr1);        // x,y
        }
#endif // ! NDIMMAKE

#endif // ! SINCOS
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
//#ifdef TREENODEALLBODIES
#ifdef SINCOS
        CHEBYSHEVTUOMP
#else
        CHEBYSHEVOMPBALLSCC;
#endif
//#else // ! TREENODEALLBODIES
//#ifdef SINCOS
//        CHEBYSHEVTUOMP
//#else
//        CHEBYSHEVOMPBALLSCC;
//#endif
//#endif // ! TREENODEALLBODIES
#endif // ! TPCF
#else // ! PTOPIVOTROTATION3
#ifdef TPCF
#if SINCOS
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#else
        cosphi = -dr[1]/(*dr1);
#endif
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
        CHEBYSHEVOMPBALLSCC;
#endif // ! TPCF
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cc_omp(nodeptr p, nodeptr q, real *dr1, vector dr,
            gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc, gdhistptr_sincos_omp histccsincos,
            INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
    int n;
    real xi, xj;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif

    if (nodes_set_bin(p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(p)*Nb(q);
        hist->histNSubXi2pcfthread[n] = hist->histNSubXi2pcfthread[n] + 1.;
        xj = Kappa(p);
        xi = Kappa(q);
        hist->histXi2pcfthreadsub[n] += xj*xi;
#ifndef PTOPIVOTROTATION3
#ifdef TPCF
        real cosphi;
#ifdef SINCOS
        real sinphi;
        histccsincos->histNSubthread[n] = histccsincos->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
//B
#if NDIM == 3
/*
//B
        real s, sy, phi;
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        phi = racos(cosphi);
        sy = rsqrt( rsqr(*dr1) - rsqr(s) );
        if ( phi > PI || phi < TWOPI) sy *= -1;
        sinphi = sy/(*dr1);
//E
*/
//B
        real s, sy;
        vector pr0;
        DOTVP(s, dr, histccsincos->dr0);
        cosphi = s/((*dr1)*histccsincos->drpq);
        CROSSVP(pr0,histccsincos->dr0,Pos(q));
        DOTVP(sy, dr, pr0);
        sinphi = rsqrt(1.0 - rsqr(cosphi));;
        if (sy < 0) sinphi *= -1.0;
//E
#else // ! NDIM
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#endif // ! NDIM
//E
#else // ! NDIMMAKE
        real s, sy;
        vector pr0;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histccsincos->dr0);
            cosphi = s/((*dr1)*histccsincos->drpq);
            CROSSVP_3D(pr0,histccsincos->dr0,Pos(q));
            DOTVP(sy, dr, pr0);
            sinphi = rsqrt(1.0 - rsqr(cosphi));;
            if (sy < 0) sinphi *= -1.0;
        } else {
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
        }
#endif // ! NDIMMAKE


#else // ! SINCOS
        histcc->histNSubthread[n] = histcc->histNSubthread[n] + 1.;

#ifdef NDIMMAKE
#if NDIM == 3
        real s;
        DOTVP(s, dr, histcc->dr0);
        cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
        cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#else // ! NDIMMAKE
        real s;
        if (gd.dimension == 3) {
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
        } else {
            cosphi = -dr[1]/(*dr1);        // x,y
        }
#endif // ! NDIMMAKE


#endif // ! SINCOS
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
//#ifdef TREENODEALLBODIES
#ifdef SINCOS
        CHEBYSHEVTUOMP
#else
        CHEBYSHEVOMPBALLSCC;
#endif
//#else // ! TREENODEALLBODIES
//#ifdef SINCOS
//        CHEBYSHEVTUOMP
//#else
//        CHEBYSHEVOMPBALLSCC;
//#endif
//#endif // ! TREENODEALLBODIES
#endif // ! TPCF
#else // ! PTOPIVOTROTATION3
#ifdef TPCF
#if SINCOS
        cosphi = -dr[0]/(*dr1);
        sinphi = -dr[1]/(*dr1);
#else
        cosphi = -dr[1]/(*dr1);
#endif
        if (rabs(cosphi)>1.0)
            verb_log_print(cmd.verbose, gd.outlog,
                "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
        CHEBYSHEVOMPBALLSCC;
#endif // ! TPCF
#endif // ! PTOPIVOTROTATION3
        *ncccalcthread += 1;
    } // ! accept_body
}


//B BALLS4 :: METODO BALLS4 DE BUSQUEDA

local void walktree_balls6_omp(nodeptr *aptr, nodeptr *nptr,
            cellptr cptr, cellptr bptr,
            nodeptr p, real psize, vector pmid,
            INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
            gdhistptr_omp_balls6 hist, gdhistptr_omp_balls6 histcc,
            gdhistptr_sincos_omp_balls6 histsincos,
            INTEGER *ibfcountthread, INTEGER *nsmoothcountthread, INTEGER nbody)
{
    nodeptr *np, *ap, q;
    int actsafe;
    bodyptr pbf;
    real dr1;
    vector dr;
    int n;

    if (Update(p)) {
        np = nptr;
        actsafe = hist->actlen - NSUB;
 //B for ap
         for (ap = aptr; ap < nptr; ap++) {
             if (Type(*ap) == CELL) {
                 if (!reject_cell_balls(p, *ap, &dr1, dr)) {
                     if ( ((Nb(*ap)<=gd.nsmooth[0]) || (Size(*ap)<=gd.rminCell[0]))
                         && (dr1 > gd.rminCell[1]) ) {
                             if (np - hist->active >= actsafe)
                                 error("walktree (1): active list overflow\n");
                             *nsmoothcountthread += 1;
                             --bptr;
                             Weight(bptr) = Weight(*ap);
                             Kappa(bptr) = Kappa(*ap);
                             SETV(Pos(bptr), Pos(*ap));
                             Id(bptr) = Id(*ap);
                             Type(bptr) = Type(*ap);
                             Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(p, *ap, &dr1, dr)) {
                             if (!scanopt(cmd.options, "no-two-balls")
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp((cellptr)(*ap),
                                    (cellptr)*ap+1, p,
                                    nbbcalcthread, nbccalcthread, ncccalcthread,
                                    histcc, histsincos);
                             } else {
                                 if (np - hist->active >= actsafe)
                                     error("walktree (2): active list overflow\n");
                                 if (!scanopt(cmd.options, "no-one-ball")) {
                                     Weight(cptr) = Weight(*ap);
                                     Kappa(cptr) = Kappa(*ap);
                                     SETV(Pos(cptr), Pos(*ap));
                                     Id(cptr) = Id(*ap);
                                     Type(cptr) = Type(*ap);
                                     Nb(cptr) = Nb(*ap);
                                     cptr++;
                                 } else // options : ! no-one-ball
                                     for (q = More(*ap); q != Next(*ap);
                                          q = Next(q))
                                         *np++= q;
                             } // meet condition :: no-wo-balls
                         } else // First meet condition
                             for (q = More(*ap); q != Next(*ap); q = Next(q))
                                 *np++= q;
                     } // ! bucket condition
                 } // ! reject_cell
             } else  // ! == CELL
                 if (*ap != p) {
                     --bptr;
                     Weight(bptr) = Weight(*ap);
                     Kappa(bptr) = Kappa(*ap);
                     SETV(Pos(bptr), Pos(*ap));
                     Id(bptr) = Id(*ap);
                     Type(bptr) = Type(*ap);
                     Nb(bptr) = 1;
                 }
         }
//E End loop for ap

        gd.actmax = MAX(gd.actmax, np - hist->active);
        if (np != nptr)
            walksub6_omp(nptr, np, cptr, bptr, p, psize, pmid,
                        nbbcalcthread, nbccalcthread, ncccalcthread,
                            hist, histcc, histsincos, ibfcountthread, nsmoothcountthread, nbody);
        else {
            if (Type(p) != BODY)
                error("walktree: recursion terminated with cell\n");

            sum_balls6_omp(cptr, bptr, (bodyptr) p,
                           nbbcalcthread, nbccalcthread, ncccalcthread,
                           hist, histsincos, nbody);
            Update(p) = FALSE;

#ifdef DEBUG
            pbf = bodytabbf + *ibfcountthread;
            Weight(pbf) = Weight(p);
            Kappa(pbf) = Kappa(p);
            SETV(Pos(pbf), Pos(p));
            *ibfcountthread += 1;
            Id(pbf) = *ibfcountthread;
#endif
        }
    }   // ! update
}

local void walksub6_omp(nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
        nodeptr p, real psize, vector pmid,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_omp_balls6 histcc,
        gdhistptr_sincos_omp_balls6 histsincos,
        INTEGER *ibfcountthread, INTEGER *nsmoothcountthread, INTEGER nbody)
{
    real poff;
    nodeptr q;
    int k;
    vector nmid;

    poff = psize / 4;
    if (Type(p) == CELL) {
        for (q = More(p); q != Next(p); q = Next(q)) {
            for (k = 0; k < NDIM; k++)
                nmid[k] = pmid[k] + (Pos(q)[k] < pmid[k] ? - poff : poff);
            walktree_balls6_omp(nptr, np, cptr, bptr, q, psize / 2, nmid,
                    nbbcalcthread, nbccalcthread, ncccalcthread, hist, histcc,
                    histsincos, ibfcountthread, nsmoothcountthread, nbody);
        }
    } else {
        for (k = 0; k < NDIM; k++)
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        walktree_balls6_omp(nptr, np, cptr, bptr, p, psize / 2, nmid,
                    nbbcalcthread, nbccalcthread, ncccalcthread, hist, histcc,
                    histsincos, ibfcountthread, nsmoothcountthread, nbody);
    }
}

local void sum_balls6_omp(cellptr cptr, cellptr bptr, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histsincos,
                          INTEGER nbody)
{
    int n;
    local INTEGER ip;
    gdhist_omp_balls6 hist1;
    gdhist_sincos_omp_balls6 hist1sincos;

    search_init_omp_balls6_cc(&hist1);
    search_init_sincos_omp_balls6(&hist1sincos);
//B
    for (n = 1; n <= cmd.sizeHistN; n++) {
        hist1.histNSubthread[n] = 0.0;
        hist1.histXi2pcfthreadsub[n] = 0.0;
    }
#ifdef TPCF
#ifdef SINCOS
            for (n = 1; n <= cmd.sizeHistN; n++) {
                hist1sincos.histNSubthread[n] = 0.0;
            }
            CLRM_ext_ext(hist1sincos.histXithreadcos, cmd.mchebyshev+1,
                         cmd.sizeHistN);
            CLRM_ext_ext(hist1sincos.histXithreadsin, cmd.mchebyshev+1,
                         cmd.sizeHistN);
#else // ! SINCOS
    CLRM_ext_ext(hist1.histXithread, cmd.mchebyshev+1, cmd.sizeHistN);
#endif
#endif // ! TPCF

#ifdef TPCF
#if NDIM == 3
#ifdef SINCOS
            dRotation3D(Pos(p0), ROTANGLE, ROTANGLE, ROTANGLE, hist1sincos.q0);
            DOTPSUBV(hist1sincos.drpq2, hist1sincos.dr0, Pos(p0), hist1sincos.q0);
            hist1sincos.drpq = rsqrt(hist1sincos.drpq2);
#ifdef PTOPIVOTROTATION
            real rtheta;
            vector dr0rot;
            rtheta = xrandom(0.0, TWOPI);
            RotationVecAWRtoVecB(dr0rot, hist1sincos.dr0, Pos(p0), rtheta);
            SETV(hist1sincos.dr0, dr0rot);
#endif
#else // ! SINCOS
    dRotation3D(Pos(p0), ROTANGLE, ROTANGLE, ROTANGLE, hist1.q0);
    DOTPSUBV(hist1.drpq2, hist1.dr0, Pos(p0), hist1.q0);
    hist1.drpq = rsqrt(hist1.drpq2);
#ifdef PTOPIVOTROTATION
          real rtheta;
          vector dr0rot;
          rtheta = xrandom(0.0, TWOPI);
          RotationVecAWRtoVecB(dr0rot, hist1.dr0, Pos(p0), rtheta);
          SETV(hist1.dr0, dr0rot);
#endif
#endif // ! SINCOS
#endif // ! NDIM
#endif // ! TPCF
//E
    if (!scanopt(cmd.options, "no-one-ball"))
        sumcell_balls6_omp(hist->interact, cptr, (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);
    sumnode_balls6_omp(bptr, hist->interact + hist->actlen, (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);

//B Section of type: computeBodyProperties_omp_balls6(p0, cmd.nbody, hist)
    real xi, xi_2p;

// BODY3
    if (Type(p0) == BODY) {
        xi = Kappa(p0)/nbody;
        xi_2p = Kappa(p0);
    } else if (Type(p0) == BODY3) {
        xi = Nbb(p0)*Kappa(p0)/nbody;
        xi_2p = Nbb(p0)*Kappa(p0);
    }
//
#ifdef TPCF
#ifdef SINCOS
    int m;
    for (m=1; m<=cmd.mchebyshev+1; m++)
        for (n=1; n<=cmd.sizeHistN; n++) {
            hist1sincos.histXithreadcos[m][n] /=
                    MAX(hist1sincos.histNSubthread[n],1.0);
            hist1sincos.histXithreadsin[m][n] /=
                MAX(hist1sincos.histNSubthread[n],1.0);
        }

    for (m=1; m<=cmd.mchebyshev+1; m++){
        OUTVP_ext(hist1sincos.xiOUTVPcos, hist1sincos.histXithreadcos[m],
                  hist1sincos.histXithreadcos[m], cmd.sizeHistN);
        OUTVP_ext(hist1sincos.xiOUTVPsin, hist1sincos.histXithreadsin[m],
                  hist1sincos.histXithreadsin[m],cmd.sizeHistN);
        OUTVP_ext(hist1sincos.xiOUTVPsincos, hist1sincos.histXithreadsin[m],
                  hist1sincos.histXithreadcos[m],cmd.sizeHistN);
        CLRM_ext(hist1sincos.histZetaMtmpcos,cmd.sizeHistN);
        CLRM_ext(hist1sincos.histZetaMtmpsin,cmd.sizeHistN);
        CLRM_ext(hist1sincos.histZetaMtmpsincos,cmd.sizeHistN);
        MULMS_ext(hist1sincos.histZetaMtmpcos,hist1sincos.xiOUTVPcos,
                  xi,cmd.sizeHistN);
        MULMS_ext(hist1sincos.histZetaMtmpsin,hist1sincos.xiOUTVPsin,
                  xi,cmd.sizeHistN);
        MULMS_ext(hist1sincos.histZetaMtmpsincos,hist1sincos.xiOUTVPsincos,
                  xi,cmd.sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadcos[m],
                 histsincos->histZetaMthreadcos[m],
                 hist1sincos.histZetaMtmpcos,cmd.sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadsin[m],
                 histsincos->histZetaMthreadsin[m],
                 hist1sincos.histZetaMtmpsin,cmd.sizeHistN);
        ADDM_ext(histsincos->histZetaMthreadsincos[m],
                 histsincos->histZetaMthreadsincos[m],
                 hist1sincos.histZetaMtmpsincos,cmd.sizeHistN);
    }
#else // ! SINCOS
    int m;
    for (m=1; m<=cmd.mchebyshev+1; m++)
        for (n=1; n<=cmd.sizeHistN; n++)
            hist1.histXithread[m][n] /= MAX(hist1.histNSubthread[n],1.0);
    for (m=1; m<=cmd.mchebyshev+1; m++){
        OUTVP_ext(hist1.xiOUTVP, hist1.histXithread[m],
                  hist1.histXithread[m],cmd.sizeHistN);
        CLRM_ext(hist1.histZetaMtmp,cmd.sizeHistN);
        MULMS_ext(hist1.histZetaMtmp,hist1.xiOUTVP,xi,cmd.sizeHistN);
        ADDM_ext(hist->histZetaMthread[m],hist->histZetaMthread[m],
                 hist1.histZetaMtmp,cmd.sizeHistN);
    }
#endif // ! SINCOS
#endif // ! TPCF
    for (n=1; n<=cmd.sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist1.histXi2pcfthreadsub[n];
    }
//E
//B
    for (n = 1; n <= cmd.sizeHistN; n++) {
        hist->histNthread[n] += hist1.histNthread[n];
        hist->histNSubthread[n] += hist1.histNSubthread[n];
        hist->histNSubXi2pcfthread[n] += hist1.histNSubXi2pcfthread[n];
    }
//E
    *nbbcalcthread += hist->interact + hist->actlen - bptr;
    *nbccalcthread += cptr - hist->interact;

    ip = p0 - bodytab + 1;
    if (ip%cmd.stepState == 0) {
        verb_log_print(cmd.verbose_log, gd.outlog, " - Completed pivot: %ld\n", ip);
    }

    search_free_omp_balls6_cc(&hist1);
    search_free_sincos_omp_balls6(&hist1sincos);
}

local void sumnode_balls6_omp(cellptr start, cellptr finish, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
    vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif
//
    int n;
    real xi;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(p0, pb, &dr1, dr)) {
            if(dr1>cmd.rminHist) {
                ibodycount++;
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                              + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(p);
                    hist->histNSubXi2pcfthread[n] =
                                hist->histNSubXi2pcfthread[n]+1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(pb);
#ifdef TPCF
#ifndef SINCOS
                    real cosphi;
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                        "sumnode: Warning!... cossphi must be in (-1,1): %g\n",
                                       cosphi);
                    CHEBYSHEVOMP;
#else // ! SINCOS
                    real cosphi;
                    real sinphi;
                    histccsincos->histNSubthread[n] =
                                    histccsincos->histNSubthread[n] + 1.;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, histccsincos->dr0);
                    cosphi = s/((dr1)*histccsincos->drpq);
                    CROSSVP(pr0,histccsincos->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/(dr1);
                    sinphi = -dr[1]/(dr1);
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumnode: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVTUOMP;
#endif // ! SINCOS
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcell_balls6_omp(cellptr start, cellptr finish, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
    vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif
//
    int n;
    real xi;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body(p0, pb, &dr1, dr)) {
            if(dr1>cmd.rminHist) {
                ibodycount++;
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                              + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(pb);
                    hist->histNSubXi2pcfthread[n] =
                            hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + 1.;
                    xi = Kappa(pb);
#ifdef TPCF
#ifndef SINCOS
                    real cosphi;
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumcell: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVOMP;
#else // ! SINCOS
                    real cosphi;
                    real sinphi;
                    histccsincos->histNSubthread[n] =
                                    histccsincos->histNSubthread[n] + 1.;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, histccsincos->dr0);
                    cosphi = s/((dr1)*histccsincos->drpq);
                    CROSSVP(pr0,histccsincos->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/(dr1);
                    sinphi = -dr[1]/(dr1);
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                            "sumcell: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVTUOMP;
#endif // ! SINCOS
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcellcell_balls6_omp(cellptr start, cellptr finish, nodeptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
    vector dr;
    real dr1;
    nodeptr pb;
    INTEGER ibodycount=0;
#ifdef SINCOS
    real Cheb,ChebU;
    real xicosmphi;
    int m;
    real xisinmphi;
#endif
//
    int n;
    real xi;
    real xj;

    for (p = start; p < finish; p++) {
        pb = ((nodeptr) p);
        if (accept_body((bodyptr)p0, pb, &dr1, dr)) {
            if(dr1>cmd.rminHist) {
                ibodycount++;
                if (cmd.rminHist==0)
                    n = (int)(NLOGBINPD*(rlog10(dr1) - rlog10(cmd.rangeN))
                              + cmd.sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd.rminHist) * gd.i_deltaR) + 1;
                if (n<=cmd.sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(pb);
                    hist->histNSubXi2pcfthread[n] =
                            hist->histNSubXi2pcfthread[n] + 1.;
                    hist->histNSubthread[n] = hist->histNSubthread[n] + Nb(pb);
                    xj = Kappa(p0);
                    xi = Kappa(pb);
#ifdef TPCF
#ifndef SINCOS
                    real cosphi;
#if NDIM == 3
                    real s;
                    DOTVP(s, dr, hist->dr0);
                    cosphi = s/(dr1*hist->drpq);
#else
                    cosphi = -dr[1]/dr1;
#endif
                    if (rabs(cosphi)>1.0)
                       verb_log_print(cmd.verbose, gd.outlog,
                        "sumcellcell: Warning!... cossphi must be in (-1,1): %g\n",
                        cosphi);
                    CHEBYSHEVOMPCC;
#else // ! SINCOS
                    real cosphi;
                    real sinphi;
                    histccsincos->histNSubthread[n] =
                                    histccsincos->histNSubthread[n] + 1.;
#if NDIM == 3
                    real s, sy;
                    vector pr0;
                    DOTVP(s, dr, histccsincos->dr0);
                    cosphi = s/((dr1)*histccsincos->drpq);
                    CROSSVP(pr0,histccsincos->dr0,Pos(p));
                    DOTVP(sy, dr, pr0);
                    sinphi = rsqrt(1.0 - rsqr(cosphi));;
                    if (sy < 0) sinphi *= -1.0;
#else // ! NDIM
                    cosphi = -dr[0]/(dr1);
                    sinphi = -dr[1]/(dr1);
#endif // ! NDIM
                    if (rabs(cosphi)>1.0)
                        verb_log_print(cmd.verbose, gd.outlog,
                        "sumcellcell: Warning!... cossphi must be in (-1,1): %g\n",
                            cosphi);
                    CHEBYSHEVTUOMP;
#endif // ! SINCOS
#endif // ! TPCF
                    hist->histXi2pcfthreadsub[n] += xi*xj;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    *ncccalcthread += 1;
}

//E BALLS4 :: DE METODO BALLS4 DE BUSQUEDA


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

    cpustart = CPUTIME;
    verb_print(cmd.verbose, "evalHistN: Running ... (tree omp-sincos) \n");

    ThreadCount(nbody, 0);

//B Init:: gd_hist
#ifdef TPCF
    int mm;
    for (mm = 1; mm <= cmd.mchebyshev+1; mm++) {
        CLRM_ext(gd.histZetaMcos[mm], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsin[mm], cmd.sizeHistN);
        CLRM_ext(gd.histZetaMsincos[mm], cmd.sizeHistN);
    }
#endif
    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histN[nn] = 0.0;
        gd.histNSubXi2pcf[nn] = 0.0;
        gd.histXi2pcf[nn] = 0.0;
#ifdef TPCF
        for (mm = 1; mm <= cmd.mchebyshev+1; mm++)
            gd.histXi[mm][nn] = 0.0;
#endif
    }
    gd.actmax = gd.nbbcalc = gd.nbccalc = gd.ncccalc = 0;
//E


//B 2023.11.29
//#pragma omp parallel default(none)   shared(cmd,gd, btab, nbody, root, ipmin, ipmax)
#pragma omp parallel default(none)   shared(cmd,gd, btab, nbody, roottable, ipmin, ipmax)
 {
     bodyptr p;
     int n;
     INTEGER nbbcalcthread = 0;
     INTEGER nbccalcthread = 0;
     gdhist_sincos_omp hist;

     search_init_sincos_omp(&hist);

#pragma omp for nowait schedule(dynamic)
     for (p = btab + ipmin -1; p < btab + ipmax; p++) {
//B
        for (n = 1; n <= cmd.sizeHistN; n++) {
            hist.histNSubthread[n] = 0.0;               // Affect only 3pcf
                                                        //  evaluation
            hist.histXi2pcfthreadsub[n] = 0.0;
        }
#ifdef TPCF
        CLRM_ext_ext(hist.histXithreadcos, cmd.mchebyshev+1, cmd.sizeHistN);
        CLRM_ext_ext(hist.histXithreadsin, cmd.mchebyshev+1, cmd.sizeHistN);
         
#if NDIM == 3
          dRotation3D(Pos(p), ROTANGLE, ROTANGLE, ROTANGLE, hist.q0);
          DOTPSUBV(hist.drpq2, hist.dr0, Pos(p), hist.q0);
          hist.drpq = rsqrt(hist.drpq2);
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
//E
#endif // ! TPCF
//E

//B 2023.11.29
//         normal_walktree_sincos(p, ((nodeptr) root), gd.rSize, &nbbcalcthread, &nbccalcthread, &hist);
         normal_walktree_sincos(p, ((nodeptr) roottable[0]), gd.rSizeTable[0], &nbbcalcthread, &nbccalcthread, &hist);
        computeBodyProperties_sincos_omp(p, nbody, &hist);

         INTEGER ip;
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

 } // end pragma omp parallel

    for (nn = 1; nn <= cmd.sizeHistN; nn++) {
        gd.histXi2pcf[nn] /= 2.0;
        gd.histNSubXi2pcf[nn] /= 2.0;
        gd.histXi2pcf[nn] /= MAX(gd.histNSubXi2pcf[nn],1.0);
    }

    if (scanopt(cmd.options, "compute-HistN"))
        search_compute_HistN(nbody);

    gd.cpusearch = CPUTIME - cpustart;
    verb_print(cmd.verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}

local void normal_walktree_sincos(bodyptr p, nodeptr q, real qsize,
                        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    nodeptr l;
    real dr1;
    vector dr;
    int n;

    if (Update(p)) {
        if ( ((nodeptr) p) != q ) {
            if (Type(q) == CELL) {
                if (!reject_cell((nodeptr)p, q, qsize)) {
                    if (!scanopt(cmd.options, "no-one-ball")) {
                        accept_body(p, (nodeptr)q, &dr1, dr);
#ifdef LOGHIST
//B 2023.12.14
                        if (scanopt(cmd.options, "behavior-ball")) {
                            if ( (Radius(p)+Radius(q))/(dr1) < gd.deltaR)
                                sumnode_sincos_cell(p, ((cellptr) q),
                                        ( (cellptr) q+1),
                                        nbbcalcthread, nbccalcthread, hist);
                            else
                                for (l = More(q); l != Next(q); l = Next(l))
                                    normal_walktree_sincos(p,l,qsize/2,
                                        nbbcalcthread, nbccalcthread, hist);
                        } else {
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
                        }
//E 2023.12.14

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
// BODY3
                if (Type(q) == BODY)
                    sumnode_sincos(p, ((cellptr) q),
                        ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
                else
                    if (Type(p) == BODY3)
                        sumnode_sincos_body3(p, ((cellptr) q),
                            ( (cellptr) q+1), nbbcalcthread, nbccalcthread, hist);
            } // ! Type(q) == CELL
        } // ! p != q
    } // ! Update
}

local void sumnode_sincos(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
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
#endif // ! PTOPIVOTROTATION2
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

// BODY3
local void sumnode_sincos_body3(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
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

local void sumnode_sincos_cell(bodyptr p, cellptr start, cellptr finish,
                   INTEGER *nbbcalcthread, INTEGER *nbccalcthread, gdhistptr_sincos_omp hist)
{
    cellptr q;
    real dr1;
    vector dr;
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



//ADDONS:
#ifdef ADDONSDEVELOP
#include "search_omp_00b.h"
#endif

//ADDONS:
#ifdef ADDONSDEVELOP
#include "search_omp_01.h"
#endif

//ADDONS:
#ifdef PATCHES
#include "direct_simple.c"
#endif

