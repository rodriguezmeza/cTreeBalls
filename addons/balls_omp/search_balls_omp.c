/* ==============================================================================
 MODULE: search_balls_omp.c		[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:    april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: searchcalc_balls_omp();
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

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

// search=balls-omp :: searchcalc_balls_omp(btab, nbody, ipmin, ipmax, cat1, cat2)

// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

//#ifdef BALLS
#if defined(OPENMPCODE) || defined(BALLS)

local void walktree_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr, nodeptr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *, INTEGER *);
//B 2023.11.22
local void sumnodes_bb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *dr1, vector dr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_bc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);

#ifdef TREENODEBALLS4
//B BALLS4: TREENODEBALLS4
local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *, INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local void sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                          cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6, INTEGER);
local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr, cellptr, bodyptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
local void sumcellcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                  cellptr, cellptr, nodeptr, INTEGER *,
            INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_sincos_omp_balls6);
//E
#endif

//B BALLS4
local INTEGER ihit;
//E
local int imiss;

// search=balls-omp
global void searchcalc_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                 bodyptr *btable,
                                 INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
                                 int cat1, int cat2)
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
    verb_print(cmd->verbose, "\nsearchcalc_balls: Running ... (balls-omp) \n");
    if (cmd->useLogHist==FALSE)
    error("CheckParameters: can´t have normal hist and BALLS definition (useLogHist=%d)\nSet useLogHist = true\n",
              cmd->useLogHist);
#ifdef TREENODEALLBODIES
    verb_print(cmd->verbose, "treenode with all bodies... \n");
#else
#ifdef TREENODEBALLS4
    verb_print(cmd->verbose, "treenode with 2 balls using balls4 method...\n");
    verb_print(cmd->verbose,
               "finding at the same time lists of neighbour cells and bodies...\n");
#else
    verb_print(cmd->verbose, "treenode with 2 balls method...\n");
    verb_print(cmd->verbose, 
"using tree of nodes at two levels: root to search and lower for cell pivots...\n");
#endif // ! TREENODEBALLS4
#endif // ! TREENODEALLBODIES

#ifdef SINCOS
//#ifdef TPCF
    if (!cmd->computeTPCF)
    verb_print(cmd->verbose, "sincos base... \n");
//#endif
#endif
    if (scanopt(cmd->options, "no-two-balls"))
        verb_print(cmd->verbose, "with option no-two-balls... \n");
    if (scanopt(cmd->options, "no-one-ball"))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (scanopt(cmd->options, "behavior-tree-omp"))
        verb_print(cmd->verbose, "with option behavior-tree-omp... \n");
    if (scanopt(cmd->options, "center-of-mass"))
    verb_print(cmd->verbose, "with option center-of-mass... \n");

#ifdef TREENODEALLBODIES
    if (scanopt(cmd->options, "smooth-pivot"))
    verb_print(cmd->verbose, 
               "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
#endif

//#ifndef TPCF
    if (!cmd->computeTPCF)
    verb_print(cmd->verbose, "computing only 2pcf... \n");
//#endif

    if (cmd->verbose==3)
    verb_print(cmd->verbose,
        "\ncat1, cat2, nbody[cat1], nbody[cat2], ipmin, imax[cat1], ipmax[cat1]:\
        %d %d %ld %ld %ld %ld %ld\n\n",
        cat1,cat2,nbody[cat1],nbody[cat2],ipmin,ipmax[cat1],ipmax[cat2]);
    ThreadCount(cmd, gd, nbody[cat1], cat1);

// Here we clear: histZetaM, histXi, histN, histNNSubXi2pcf, histXi2pcf...
//    search_init_gd_hist();
#ifdef SINCOS
    search_init_gd_hist_sincos(cmd, gd);
#else
    search_init_gd_hist(cmd, gd);
#endif

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
    bodytabbf = (bodyptr) allocate(nbody[cat1] * sizeof(body));
    gd->bytes_tot += nbody[cat1]*sizeof(body);
    fprintf(stdout,"\n\nAllocated %g MByte for particle (found) storage.\n",
            nbody[cat1]*sizeof(body)/(1024.0*1024.0));

    verb_print(cmd->verbose,
               "\nsearchcalc_balls_omp: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));
#endif
#endif

    gd->nsmoothcount = 0;

    verb_print(cmd->verbose,
               "\nsearchcalc_balls: Total allocated %g MByte storage so far.\n",
               gd->bytes_tot/(1024.0*1024.0));


//B 2024-02-01
#ifdef TREENODEALLBODIES
    INTEGER ipfalse;
    ipfalse=0;
//B kappa Avg Rmin
    INTEGER icountNbRmin;
    icountNbRmin=0;
    INTEGER icountNbRminOverlap;
    icountNbRminOverlap=0;
//E
#endif

//    shared(cmd,gd,btable,nbody,root,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, \

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
#pragma omp parallel default(none)  \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, \
    bodytabbf,ibfcount,ihit,cat1, cat2)
#else
#pragma omp parallel default(none) \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, cat1, cat2)
#endif
#else // ! TREENODEBALLS4
#ifdef TREENODEALLBODIES
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode,imiss, cat1, cat2, ipfalse, \
           icountNbRmin, icountNbRminOverlap)
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, cat1, cat2)
#endif
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

      search_init_omp_balls6(cmd, gd, &histb, cat1);
      search_init_omp_balls6_cc(cmd, gd, &histccb);
      search_init_sincos_omp_balls6(cmd, gd, &histbsincos);
#else
      gdhist_omp_balls hist;
      gdhist_omp_balls histcc;
  //#ifdef SINCOS
        gdhist_sincos_omp histsincos;
//      verb_print_debug(1, "\nAqui voy (10): %g %d %s\n",
//                       cmd->theta, cmd->sizeHistN, cmd->rootDir);
        search_init_sincos_omp(cmd, gd, &histsincos);
//      verb_print_debug(1, "\nAqui voy (11)\n");
  //#endif

  // Here we clear threads of: histZetaM, histXi, histNNSub,
  //  histN, histNNSubXi2pcf, histXi2pcf, histXi2pcf()sub...
      search_init_balls_omp(cmd, gd, &hist, cat1);
      search_init_balls_omp_cc(cmd, gd, &histcc);
#endif
//E
/*
    gdhist_omp_balls hist;
    gdhist_omp_balls histcc;
//#ifdef SINCOS
      gdhist_sincos_omp histsincos;
      search_init_sincos_omp(&histsincos);
//#endif

// Here we clear threads of: histZetaM, histXi, histNNSub, 
//  histN, histNNSubXi2pcf, histXi2pcf, histXi2pcf()sub...
    search_init_balls_omp(&hist);
    search_init_balls_omp_cc(&histcc);
*/

#ifdef TREENODEALLBODIES
    bodyptr p;
#else
    nodeptr p;
#endif // ! TREENODEALLBODIES


#ifdef PIVOTEXTERNAL
//B (1)
// Setting reference axis here instead in loop i
//  does affect the 3pcf results for m=2 or higher
//#ifdef TPCF
      if (cmd->computeTPCF) {
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
//#endif // ! TPCF
      }
//E
#endif // ! PIVOTEXTERNAL

      //B 2024-02-01
#ifdef TREENODEALLBODIES
            INTEGER ipfalsethreads;
            ipfalsethreads = 0;
//B kappa Avg Rmin
            INTEGER icountNbRminthread;
            icountNbRminthread=0;
            INTEGER icountNbRminOverlapthread;
            icountNbRminOverlapthread=0;
//E
#endif
      //E

#pragma omp for nowait schedule(dynamic)
#ifdef TREENODEALLBODIES
        for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
// Correct it because p and q are in differents node structures!!!... cat1!=cat2
//B kappa Avg Rmin
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            KappaRmin(p) = Kappa(p);
            if (scanopt(cmd->options, "smooth-pivot")) {
                if (Update(p) == FALSE) {
                    ipfalsethreads++;
                    continue;
                }
            }
//E
#else // TREENODEALLBODIES exclude TREENODEBALLS4. Then this is for TREENODEBALLS4
        for (i=0; i< gd->nnodescanlevTable[cat1]; i++) {
            if (cmd->scanLevel==0)
                p = (nodeptr) roottable[cat1];
            else
                p = nodetablescanlev[cat1][i];
#endif
            hist.ipcount = 0;

#ifndef PIVOTEXTERNAL
//B (2)
// Set histograms to zero for the pivot
//B kappa Avg Rmin
            for (n = 1; n <= cmd->sizeHistN; n++) {
#ifdef TREENODEBALLS4
                histbsincos.histNNSubXi2pcfthreadp[n] = 0.0; // Affects only 2pcf
                histb.histNNSubXi2pcfthreadp[n] = 0.0;       // Affects only 2pcf
                histccb.histNNSubXi2pcfthreadp[n] = 0.0;     // Affects only 2pcf
#else
                hist.histNNSubXi2pcfthreadp[n] = 0.0;        // Affects only 2pcf
                histsincos.histNNSubXi2pcfthreadp[n] = 0.0;  // Affects only 2pcf
                histcc.histNNSubXi2pcfthreadp[n] = 0.0;      // Affects only 2pcf
#ifdef TREENODEALLBODIES
                hist.histXi2pcfthreadsub[n] = 0.0;      // Affects only 2pcf
//B kappa Avg Rmin
                hist.histNNSubXi2pcfthreadp[n] = 0.0;    // Affects only 2pcf
//E
#endif
#endif
            }
//E
//#ifdef TPCF
            if (cmd->computeTPCF) {
#ifdef SINCOS
                /*
                 for (n = 1; n <= cmd->sizeHistN; n++) {
                 histsincos.histNNSubthread[n] = 0.0;
                 }
                 CLRM_ext_ext(histsincos.histXithreadcos, cmd->mChebyshev+1,
                 cmd->sizeHistN);
                 CLRM_ext_ext(histsincos.histXithreadsin, cmd->mChebyshev+1,
                 cmd->sizeHistN);
                 */
#ifdef TREENODEBALLS4
                for (n = 1; n <= cmd->sizeHistN; n++) {
                    histbsincos.histNNSubthread[n] = 0.0;
                }
                CLRM_ext_ext(histbsincos.histXithreadcos, cmd->mChebyshev+1,
                             cmd->sizeHistN);
                CLRM_ext_ext(histbsincos.histXithreadsin, cmd->mChebyshev+1,
                             cmd->sizeHistN);
#else
                for (n = 1; n <= cmd->sizeHistN; n++) {
                    histsincos.histNNSubthread[n] = 0.0;
                }
                CLRM_ext_ext(histsincos.histXithreadcos, cmd->mChebyshev+1,
                             cmd->sizeHistN);
                CLRM_ext_ext(histsincos.histXithreadsin, cmd->mChebyshev+1,
                             cmd->sizeHistN);
#endif
                
#else // ! SINCOS
                for (n = 1; n <= cmd->sizeHistN; n++) {
                    //                histcc.histNNSubthread[n] = 0.0;
                    //B BALLS4
#ifdef TREENODEBALLS4
                    histb.histNNSubthread[n] = 0.0;      // Affect only 3pcf evaluation
                    histb.histXi2pcfthreadsub[n] = 0.0;
                    //B two-balls
                    histccb.histNNSubthread[n] = 0.0;    // Affect only 3pcf evaluation
                    histccb.histXi2pcfthreadsub[n] = 0.0;
                    //E
#else
                    histcc.histNNSubthread[n] = 0.0;
#endif // ! TREENODEBALLS4
                }
                //            CLRM_ext_ext(histcc.histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
#ifdef TREENODEBALLS4
                CLRM_ext_ext(histb.histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
                CLRM_ext_ext(histccb.histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
                //E
#else
                CLRM_ext_ext(histcc.histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
#endif
#endif // ! SINCOS
//#endif // ! TPCF
            }
//E
#endif // ! PIVOTEXTERNAL

#ifndef PIVOTEXTERNAL
//B (3)
// Change p with a valid point on the sphere. Use a greater scanLevel
//#ifdef TPCF
            if (cmd->computeTPCF) {
#if NDIM == 3
                
#ifdef SINCOS
                
#ifdef TREENODEALLBODIES
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
#endif // ! TREENODEALLBODIES
                
#else // ! SINCOS
                
#ifdef TREENODEALLBODIES
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
#endif // ! TREENODEALLBODIES
                
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
//#endif // ! TPCF
            }
//E
#endif // ! PIVOTEXTERNAL

#ifdef TREENODEALLBODIES
            nodeptr q;
            int j;
            if (cmd->scanLevelRoot==0) {
                walktree_balls_omp(cmd, gd, (nodeptr)p, (nodeptr) roottable[cat2],
                                   &hist, &histcc, &histsincos,
                                   &nsmoothcountthread, &nbbcalcthread,
                                   &nbccalcthread, &ncccalcthread);
//B kappa Avg Rmin
                for (n = 1; n <= cmd->sizeHistN; n++) {
                    hist.histNNSubXi2pcfthreadp[n] =
                                    ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
                    hist.histNNSubXi2pcfthreadtotal[n] +=
                                    hist.histNNSubXi2pcfthreadp[n];
                    histsincos.histNNSubXi2pcfthreadp[n] =
                            ((real)NbRmin(p))*histsincos.histNNSubXi2pcfthreadp[n];
                    histsincos.histNNSubXi2pcfthreadtotal[n] +=
                                    histsincos.histNNSubXi2pcfthreadp[n];
                    if (scanopt(cmd->options, "smooth-pivot")) {
                        hist.histNNSubthread[n] =
                        ((real)NbRmin(p))*hist.histNNSubthread[n];
                        histcc.histNNSubthread[n] =
                        ((real)NbRmin(p))*histcc.histNNSubthread[n];
                        histsincos.histNNSubthread[n] =
                        ((real)NbRmin(p))*histsincos.histNNSubthread[n];
                    }
                }
//E
            } else {
                for (j=0; j< gd->nnodescanlev_rootTable[cat2]; j++) {
                    q = nodetablescanlev_root[cat2][j];
                    walktree_balls_omp(cmd, gd,
                                       (nodeptr)p, q, &hist, &histcc, &histsincos,
                            &nsmoothcountthread, &nbbcalcthread,
                                       &nbccalcthread, &ncccalcthread);
                }
//B kappa Avg Rmin
                for (n = 1; n <= cmd->sizeHistN; n++) {
                    hist.histNNSubXi2pcfthreadp[n] =
                                    ((real)NbRmin(p))*hist.histNNSubXi2pcfthreadp[n];
                    hist.histNNSubXi2pcfthreadtotal[n] +=
                                    hist.histNNSubXi2pcfthreadp[n];
                    histsincos.histNNSubXi2pcfthreadp[n] =
                                ((real)NbRmin(p))*histsincos.histNNSubXi2pcfthreadp[n];
                    histsincos.histNNSubXi2pcfthreadtotal[n] +=
                                histsincos.histNNSubXi2pcfthreadp[n];
                    if (scanopt(cmd->options, "smooth-pivot")) {
                        hist.histNNSubthread[n] =
                        ((real)NbRmin(p))*hist.histNNSubthread[n];
                        histcc.histNNSubthread[n] =
                        ((real)NbRmin(p))*histcc.histNNSubthread[n];
                        histsincos.histNNSubthread[n] =
                        ((real)NbRmin(p))*histsincos.histNNSubthread[n];
                    }
                }
//E
            }
//#ifndef TREENODEBALLS4
      computeBodyProperties_balls_omp(cmd, gd, (bodyptr)p, nbody[cat1], &hist);
//#endif
//E
#else // ! TREENODEALLBODIES
#ifdef TREENODEBALLS4
            nodeptr q;
            int j;
// This method is faster when scanLevel is higher...
            histb.active[0] = (nodeptr) (roottable[cat2]);
//            histb.active[0] = (nodeptr) root;
//            for (j=0; j< gd->nnodescanlev_rootTable[cat2]; j++) {
//               q = nodetablescanlev_root[cat2][j];
//                histb.active[0] = q;

//                verb_print_debug(1, "\nAqui voy (0): %g\n", Kappa(histb.active[0]));
//            q = (nodeptr) (roottable[cat2]);
//                histb.active[0] = q;
//                verb_print_debug(1, "\nAqui voy (1): %g %g %g\n",Pos(p)[0],Pos(p)[1],Pos(p)[2]);


            walktree_balls6_omp(cmd, gd, histb.active, histb.active + 1,
                    histb.interact, histb.interact + histb.actlen,
                    p, Size(p), Pos(p),
                    &nbbcalcthread, &nbccalcthread, &ncccalcthread,
                    &histb, &histccb, &histbsincos,
                    &ibfcountthread, &nsmoothcountthread, nbody[cat1]);
//            }

            if (!scanopt(cmd->options, "no-two-balls")) {
#ifndef SINCOS
                computeBodyProperties_omp_balls6_cc(cmd, gd, (bodyptr)p,
                                                    nbody[cat1], &histccb);
#else
                computeBodyProperties_sincos_omp_balls6_cc(cmd, gd, (bodyptr)p,
                                                    nbody[cat1], &histbsincos);
#endif
            }
#else // ! TREENODEBALLS4
            nodeptr q;
            int j;
//B Faster if first option is used:
            if (cmd->scanLevelRoot==0)
                walktree_balls_omp(cmd, gd,
                                   p, (nodeptr) roottable[cat2], &hist, &histcc,
                                   &histsincos, &nsmoothcountthread,
                                   &nbbcalcthread, &nbccalcthread, &ncccalcthread);
//E Above one.
            else {
// Must be the same node's trees
                if (gd->nnodescanlevTable[cat1] == gd->nnodescanlev_rootTable[cat2])
                    for (j=i; j< gd->nnodescanlevTable[cat2]; j++) {
                        q = nodetablescanlev_root[cat2][j];
                        walktree_balls_omp(cmd, gd,
                                           p, q, &hist, &histcc, &histsincos,
                                &nsmoothcountthread, &nbbcalcthread, &nbccalcthread, &ncccalcthread);
                    }
                else
                for (j=0; j< gd->nnodescanlev_rootTable[cat2]; j++) {
                q = nodetablescanlev_root[cat2][j];
                walktree_balls_omp(cmd, gd, p, q, &hist, &histcc, &histsincos,
                    &nsmoothcountthread, &nbbcalcthread, &nbccalcthread, &ncccalcthread);
                }
            }
#endif // ! TREENODEBALLS4
//E
#endif // ! TREENODEALLBODIES

//#ifdef TPCF
            if (cmd->computeTPCF) {
#ifndef PIVOTEXTERNAL
                //B (4)
#ifdef TREENODEALLBODIES
#ifdef SINCOS
                computeBodyProperties_balls_omp_cc_sincos(cmd, gd, (bodyptr)p,
                                                          nbody[cat1], &histsincos);
#else
                computeBodyProperties_balls_omp_cc(cmd, gd,
                                                   (bodyptr)p, nbody[cat1], &histcc);
#endif // ! SINCOS
#endif // ! TREENODEALLBODIES
                //E
#endif // ! PIVOTEXTERNAL
//#endif // ! TPCF
            }

//B kappa Avg Rmin
#ifdef TREENODEALLBODIES
            icountNbRminthread += NbRmin(p);
            icountNbRminOverlapthread += NbRminOverlap(p);
#endif
//E

#ifndef TREENODEBALLS4
#ifdef TREENODEALLBODIES
            i = p - btable[cat1] + 1;
            if (i%cmd->stepState == 0) {
                verb_log_print(cmd->verbose_log, gd->outlog,
                    " - Completed pivot node: %d\n", i);
/*                if (scanopt(cmd->options, "smooth-pivot"))
                verb_print(cmd->verbose,
                    " - Number of neighbours: %ld %d %d\n",
                           i, NbRmin(p), NbRminOverlap(p)); */
            }

//B 2024-02-01
// Correct it because p and q are in differents node structures!!!... cat1!=cat2
//            if (scanopt(cmd->options, "smooth-pivot")) {
//                verb_print(cmd->verbose,
//                    " - Number of neighbours: %d\n", NbRmin(p));
//            }
//E

#else
            if (i%cmd->stepState == 0) {
                verb_log_print(cmd->verbose_log, gd->outlog,
                    " - Completed pivot node: %d\n", i);
            }
#endif
#endif
        } // end loop i

#ifndef TREENODEBALLS4
#ifndef TREENODEALLBODIES
      computeBodyProperties_balls_omp(cmd, gd, (bodyptr)p, nbody[cat1], &hist);
#endif
#endif

//#ifdef TPCF
            if (cmd->computeTPCF) {
#ifdef PIVOTEXTERNAL
                //B (5)
#ifdef TREENODEALLBODIES
#ifdef SINCOS
                computeBodyProperties_balls_omp_cc_sincos(cmd, gd,
                                                          (bodyptr)p, nbody[cat1], &histsincos);
#else
                computeBodyProperties_balls_omp_cc(cmd, gd,
                                                   (bodyptr)p, nbody[cat1], &histcc);
#endif
#endif // ! TREENODEALLBODIES
                //E
#endif
//#endif
            }

#pragma omp critical
    {
#ifdef TREENODEBALLS4
        int n;
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gd->histNN[n] += histb.histNthread[n];
            gd->histNNSub[n] += histb.histNNSubthread[n];
            gd->histNNSubXi2pcf[n] += histb.histNNSubXi2pcfthread[n];
            gd->histXi2pcf[n] += histb.histXi2pcfthread[n];
            if (!scanopt(cmd->options, "no-two-balls")) {
                gd->histNN[n] += histccb.histNthread[n];
                gd->histNNSub[n] += histccb.histNNSubthread[n];
                gd->histNNSubXi2pcf[n] += histccb.histNNSubXi2pcfthread[n];
                gd->histXi2pcf[n] += histccb.histXi2pcfthread[n];
            }
        }
#else
      for (n = 1; n <= cmd->sizeHistN; n++) {
          gd->histNN[n] += hist.histNthread[n];
          gd->histNNSubXi2pcf[n] += hist.histNNSubXi2pcfthread[n];
//B kappa Avg Rmin
          gd->histNNSubXi2pcftotal[n] += hist.histNNSubXi2pcfthreadtotal[n];
//E
          gd->histXi2pcf[n] += hist.histXi2pcfthread[n];
      }
#endif // ! TREENODEBALLS4

//#ifdef TPCF
        if (cmd->computeTPCF) {
#ifdef TREENODEBALLS4
#ifdef SINCOS
            int m;
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                         histbsincos.histZetaMthreadcos[m],cmd->sizeHistN);
                ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                         histbsincos.histZetaMthreadsin[m],cmd->sizeHistN);
                ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                         histbsincos.histZetaMthreadsincos[m],cmd->sizeHistN);
                //        if (!scanopt(cmd->options, "no-two-balls")) {
                //            ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histccb.histZetaMthread[m],
                //                     cmd->sizeHistN);
                //        }
            }
#else // ! SINCOS
            int m;
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histb.histZetaMthread[m],
                         cmd->sizeHistN);
                if (!scanopt(cmd->options, "no-two-balls")) {
                    ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histccb.histZetaMthread[m],
                             cmd->sizeHistN);
                }
            }
#endif // ! SINCOS
#else // ! TREENODEBALLS4
#ifdef SINCOS
            
            int m;
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                ADDM_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],
                         histsincos.histZetaMthreadcos[m],cmd->sizeHistN);
                ADDM_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],
                         histsincos.histZetaMthreadsin[m],cmd->sizeHistN);
                ADDM_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],
                         histsincos.histZetaMthreadsincos[m],cmd->sizeHistN);
                
            }
            
#else
            int m;
            for (m=1; m<=cmd->mChebyshev+1; m++)
                ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histcc.histZetaMthread[m],
                         cmd->sizeHistN);
#endif // ! SINCOS
#endif // ! TREENODEBALLS4
//#endif // ! TPCF
        }
        gd->nsmoothcount += nsmoothcountthread;
        gd->nbbcalc += nbbcalcthread;
        gd->nbccalc += nbccalcthread;
        gd->ncccalc += ncccalcthread;
#ifdef TREENODEBALLS4
#ifdef DEBUG
        ibfcount += ibfcountthread;
#endif
#endif

#ifdef TREENODEALLBODIES
//        ipfalse += ipfalsethreads;
        //B kappa Avg Rmin
                ipfalse += ipfalsethreads;
                icountNbRmin += icountNbRminthread;
                icountNbRminOverlap += icountNbRminOverlapthread;
        //E
#endif
    } // end pragma omp critical

//B BALLS4
//    search_free_sincos_omp(&histsincos);
//    search_free_balls_omp_cc(&histcc);
//    search_free_balls_omp(&hist);
#ifdef TREENODEBALLS4
    search_free_sincos_omp_balls6(cmd, gd, &histbsincos);
    search_free_omp_balls6_cc(cmd, gd, &histccb);
    search_free_omp_balls6(cmd, gd, &histb);
#else
    search_free_sincos_omp(cmd, gd, &histsincos);
    search_free_balls_omp_cc(cmd, gd, &histcc);
    search_free_balls_omp(cmd, gd, &hist);
#endif
  } // end pragma omp parallel

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
    gd->nbodybf = cmd->ntosave;
    verb_print(cmd->verbose, "\nsearchcalc_balls_omp: Total bodies found: %ld %ld\n",ibfcount, gd->nbodybf);
#endif
#endif
//E
#ifdef TREENODEALLBODIES
      bodyptr p;
      real xi, den, num;
      int m;
      if (scanopt(cmd->options, "smooth-pivot")) {
//          den = 0;
//          for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++)
//              den += NbRminOverlap(p);
          num = (real)nbody[cat1];
          den = (real)(nbody[cat1]-ipfalse);              // Little higher than expec.
//          den = nbody[cat1];                            // Lower than expected
          xi = num/den;
          verb_print(cmd->verbose,
                "balls-omp: p falses found (false, num, den, xi) = %ld %e %e %e\n",
                ipfalse, num, den, xi);
//#ifdef TPCF
          if (cmd->computeTPCF) {
#ifdef SINCOS
              for (m=1; m<=cmd->mChebyshev+1; m++) {
                  MULMS_ext(gd->histZetaMcos[m],gd->histZetaMcos[m],xi,cmd->sizeHistN);
                  MULMS_ext(gd->histZetaMsin[m],gd->histZetaMsin[m],xi,cmd->sizeHistN);
                  MULMS_ext(gd->histZetaMsincos[m],gd->histZetaMsincos[m],xi,cmd->sizeHistN);
              }
#else
              for (m=1; m<=cmd->mChebyshev+1; m++) {
                  MULMS_ext(gd->histZetaM[m],gd->histZetaM[m],xi,cmd->sizeHistN);
              }
#endif
//#endif // ! TPCF
          }

      }
#endif // ! TREENODEALLBODIES

    for (nn = 1; nn <= cmd->sizeHistN; nn++) {
        if (cmd->verbose > 2)
            verb_print(cmd->verbose,"%d %e %e\n",
                       nn, gd->histNNSubXi2pcf[nn], gd->histNNSubXi2pcftotal[nn]);
        gd->histXi2pcf[nn] /= 2.0;
        gd->histNNSubXi2pcf[nn] /= 2.0;
//        gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
//B kappa Avg Rmin
        gd->histNNSubXi2pcftotal[nn] /= 2.0;
        if (scanopt(cmd->options, "smooth-pivot")) {
            gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
        } else {
            gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
        }
//E
    }

#ifdef TREENODEALLBODIES
    if (scanopt(cmd->options, "compute-HistN")) {
        if (scanopt(cmd->options, "smooth-pivot")) {
            search_compute_HistN_balls(cmd, gd, nbody[cat1]-ipfalse);
        } else {
            search_compute_HistN_balls(cmd, gd, nbody[cat1]);
        }
    }
#else
      if (scanopt(cmd->options, "compute-HistN")) {
              search_compute_HistN_balls(cmd, gd, nbody[cat1]);
      }
#endif

    verb_print(cmd->verbose, "balls: nsmoothcount = %ld\n",gd->nsmoothcount);
    verb_print(cmd->verbose, "balls: imiss = %d\n",imiss);

#ifdef TREENODEALLBODIES
      verb_print(cmd->verbose, "balls: p falses found = %ld\n",ipfalse);
      //B kappa Avg Rmin
          verb_print(cmd->verbose,
                     "tree-omp-sincos: count NbRmin found = %ld\n",
                     icountNbRmin);
          verb_print(cmd->verbose,
                     "tree-omp-sincos: count overlap found = %ld\n",
                     icountNbRminOverlap);

//          bodyptr p;
          INTEGER ifalsecount;
          ifalsecount = 0;
          INTEGER itruecount;
          itruecount = 0;
          for (p = btable[cat1] + ipmin -1; p < btable[cat1] + ipmax[cat1]; p++) {
              if (Update(p) == FALSE) {
                  ifalsecount++;
              } else {
                  itruecount++;
              }
          }
          verb_print(cmd->verbose, "tree-omp-sincos: p falses found = %ld\n",
                     ifalsecount);
          verb_print(cmd->verbose, "tree-omp-sincos: p true found = %ld\n",
                     itruecount);
          verb_print(cmd->verbose, "tree-omp-sincos: total = %ld\n",
                     itruecount+ifalsecount);
      //E
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
}



#ifdef TREENODEALLBODIES

local void walktree_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr p, nodeptr q,
        gdhistptr_omp_balls hist, gdhistptr_omp_balls histcc,
        gdhistptr_sincos_omp histccsincos, INTEGER *nsmoothcountthread,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread)
{
 nodeptr l;
 real dr1;
 vector dr;

    if (!reject_balls(cmd, gd, p, q, &dr1, dr)) {
        if (scanopt(cmd->options, "smooth-pivot"))
            if (dr1<=gd->rsmooth[0] && Type(q)==BODY) {
                if (Update(q)=TRUE) {
                    Update(q) = FALSE;
                    NbRmin(p) += 1;
                    KappaRmin(p) += Kappa(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
        if (Type(q) == CELL) {
            if (!scanopt(cmd->options, "no-one-ball")) {
                if ( ((Nb(q)<=gd->nsmooth[0]) || (Size(q)<=gd->rminCell[0]))
                                            && (dr1 > gd->rminCell[1]) ) {
                    *nsmoothcountthread += 1;
                    sumnodes_bc_omp(cmd, gd,
                                    p, q, &dr1, dr, hist, histcc, histccsincos,
                                        nbbcalcthread, nbccalcthread, ncccalcthread);
                } else {
                    if ( Radius(q)/(dr1) < gd->deltaR)
                        sumnodes_bc_omp(cmd, gd, p, q, &dr1, dr, hist, histcc,
                                            histccsincos, nbbcalcthread,
                                            nbccalcthread, ncccalcthread);
                    else
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(cmd, gd,
                                               p,l,hist, histcc, histccsincos,
                                                   nsmoothcountthread, nbbcalcthread,
                                                   nbccalcthread, ncccalcthread);
                }
            } else {    // ! no-one-ball
                for (l = More(q); l != Next(q); l = Next(l))
                    walktree_balls_omp(cmd, gd, p,l,hist, histcc, histccsincos,
                                           nsmoothcountthread, nbbcalcthread,
                                           nbccalcthread, ncccalcthread);
            }           // ! no-one-ball
        } else { // ! Type(q) == CELL
                sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr, hist, histcc, histccsincos,
                                nbbcalcthread, nbccalcthread, ncccalcthread);
        } // ! Type(q) == CELL
    } // ! reject_cell
}

#else // ! TREENODEALLBODIES

#define CELLCELL 1
#define CELLBODY 2
#define BODYCELL 3
#define BODYBODY 4

local void walktree_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr p, nodeptr q,
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

    if (!reject_balls(cmd, gd, p, q, &dr1, dr)) {
        switch (SWITCHVALUE){
            case CELLCELL:
                 if (!scanopt(cmd->options, "no-two-ball")) {
                     if ( ( (Size(p)<=gd->rminCell[0] && Size(q)<=gd->rminCell[0])
                           || (Nb(p)<=gd->nsmooth[0] && Nb(q)<=gd->nsmooth[0]) )
                         && (dr1 > gd->rminCell[1]) ) {
                         *nsmoothcountthread += 1;
                         sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                  nbbcalcthread, nbccalcthread, ncccalcthread);
                     } else {
                     if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if ( (Size(p)<=gd->rminCell[0] && Size(q)<=gd->rminCell[0])
                                 && (Nb(p)<=gd->nsmooth[0]
                                     && Nb(q)<=gd->nsmooth[0]) ) {
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
                if (!scanopt(cmd->options, "no-one-ball")) {
                    if ( ((Nb(p)<=gd->nsmooth[0]) || (Size(p)<=gd->rminCell[0]))
                        && (dr1 > gd->rminCell[1]) ) {
                        *nsmoothcountthread += 1;
                        sumnodes_bb_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                 nbbcalcthread, nbccalcthread, ncccalcthread);
                    } else {

                    if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if (Size(p)<=gd->rminCell[0] || Nb(p)<=gd->nsmooth[0]) {
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
                            walktree_balls_omp(l, q, hist, histcc, histccsincos, nsmoothcountthread,
                                    nbbcalcthread, nbccalcthread, ncccalcthread);
                    }
                break;
            case BODYCELL:
                if (!scanopt(cmd->options, "no-one-ball")) {
                    if ( ((Nb(q)<=gd->nsmooth[0]) || (Size(q)<=gd->rminCell[0]))
                        && (dr1 > gd->rminCell[1]) ) {
                        *nsmoothcountthread += 1;
                        sumnodes_bc_omp(p, q, &dr1, dr, hist, histcc, histccsincos,
                                     nbbcalcthread, nbccalcthread, ncccalcthread);
                    } else {

                    if (nodes_condition_balls(p, q, &dr1, dr)) {
                        if (Size(q)<=gd->rminCell[0] || Nb(q)<=gd->nsmooth[0]) {
                            *nsmoothcountthread += 1;
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
                    for (l = More(q); l != Next(q); l = Next(l)) {
                            walktree_balls_omp(p,l,hist, histcc, histccsincos, nsmoothcountthread,
                                    nbbcalcthread, nbccalcthread, ncccalcthread);
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
    } // ! reject_cell
}

#undef CELLCELL
#undef CELLBODY
#undef BODYCELL
#undef BODYBODY

#endif // ! TREENODEALLBODIES

#include "datastruc_defs_balls_omp.h"


local void sumnodes_bb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, vector dr,
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

    if (nodes_set_bin(cmd, gd, p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(p)*Nb(q);
        hist->histNNSubXi2pcfthread[n] = hist->histNNSubXi2pcfthread[n] + 1.;
//B kappa Avg Rmin
        hist->histNNSubXi2pcfthreadp[n] = hist->histNNSubXi2pcfthreadp[n] + 1.;
//E
        xj = Kappa(p);
        xi = Kappa(q);
#ifdef TREENODEALLBODIES
        hist->histXi2pcfthreadsub[n] += xi;
#else
        hist->histXi2pcfthreadsub[n] += xj*xi;
#endif

#ifndef PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            
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
#else // ! SINCOS
            histcc->histNNSubthread[n] = histcc->histNNSubthread[n] + 1.;
            
#if NDIM == 3
            real s;
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
#else
            cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#endif // ! SINCOS
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbbcalcthread += 1;
    } // ! accept_body
}

local void sumnodes_bc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, vector dr,
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

    if (nodes_set_bin(cmd, gd, p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(q);
        hist->histNNSubXi2pcfthread[n] = hist->histNNSubXi2pcfthread[n] + 1.;
        xj = Kappa(p);
        xi = Kappa(q);
#ifdef TREENODEALLBODIES
        hist->histXi2pcfthreadsub[n] += xi;
#else
        hist->histXi2pcfthreadsub[n] += xj*xi;
#endif

#ifndef PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
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
#else // ! SINCOS
            histcc->histNNSubthread[n] = histcc->histNNSubthread[n] + 1.;
#if NDIM == 3
            real s;
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
            cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
            //E
#endif // ! SINCOS
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr q, nodeptr p, real *dr1, vector dr,
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

    if (nodes_set_bin(cmd, gd, p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(q);
        hist->histNNSubXi2pcfthread[n] = hist->histNNSubXi2pcfthread[n] + 1.;
        xj = Kappa(q);
        xi = Kappa(p);
#ifdef TREENODEALLBODIES
        hist->histXi2pcfthreadsub[n] += xi;
#else
        hist->histXi2pcfthreadsub[n] += xj*xi;
#endif

#ifndef PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
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
#else // ! SINCOS
            histcc->histNNSubthread[n] = histcc->histNNSubthread[n] + 1.;
#if NDIM == 3
            real s;
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
            cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#endif // ! SINCOS
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, vector dr,
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

    if (nodes_set_bin(cmd, gd, p, q, &n, dr1, dr)) {
        hist->histNthread[n] = hist->histNthread[n] + Nb(p)*Nb(q);
        hist->histNNSubXi2pcfthread[n] = hist->histNNSubXi2pcfthread[n] + 1.;
        xj = Kappa(p);
        xi = Kappa(q);
#ifdef TREENODEALLBODIES
        hist->histXi2pcfthreadsub[n] += xi;
#else
        hist->histXi2pcfthreadsub[n] += xj*xi;
#endif

#ifndef PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            //B
#if NDIM == 3
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
#else // ! SINCOS
            histcc->histNNSubthread[n] = histcc->histNNSubthread[n] + 1.;
#if NDIM == 3
            real s;
            DOTVP(s, dr, histcc->dr0);
            cosphi = s/((*dr1)*histcc->drpq);
#else // ! NDIM
            cosphi = -dr[1]/(*dr1);        // x,y
#endif // ! NDIM
#endif // ! SINCOS
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenodes_bb: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (cmd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *ncccalcthread += 1;
    } // ! accept_body
}


//B BALLS4 :: METODO BALLS4 DE BUSQUEDA
#ifdef TREENODEBALLS4

local void walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr *aptr, nodeptr *nptr,
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
                     if ( ((Nb(*ap)<=gd->nsmooth[0]) || (Size(*ap)<=gd->rminCell[0]))
                         && (dr1 > gd->rminCell[1]) ) {
                             if (np - hist->active >= actsafe)
                                 error("walktree (1): active list overflow\n");
                             *nsmoothcountthread += 1;
                             --bptr;
                             Mass(bptr) = Mass(*ap);
                             Kappa(bptr) = Kappa(*ap);
                             SETV(Pos(bptr), Pos(*ap));
                             Id(bptr) = Id(*ap);
                             Type(bptr) = Type(*ap);
                             Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(p, *ap, &dr1, dr)) {
                             if (!scanopt(cmd->options, "no-two-balls")
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp((cellptr)(*ap),
                                    (cellptr)*ap+1, p,
                                    nbbcalcthread, nbccalcthread, ncccalcthread,
                                    histcc, histsincos);
                             } else {
                                 if (np - hist->active >= actsafe)
                                     error("walktree (2): active list overflow\n");
                                 if (!scanopt(cmd->options, "no-one-ball")) {
                                     Mass(cptr) = Mass(*ap);
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
                     Mass(bptr) = Mass(*ap);
                     Kappa(bptr) = Kappa(*ap);
                     SETV(Pos(bptr), Pos(*ap));
                     Id(bptr) = Id(*ap);
                     Type(bptr) = Type(*ap);
                     Nb(bptr) = 1;
                 }
         }
//E End loop for ap

        gd->actmax = MAX(gd->actmax, np - hist->active);
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
            Mass(pbf) = Mass(p);
            Kappa(pbf) = Kappa(p);
            SETV(Pos(pbf), Pos(p));
            *ibfcountthread += 1;
            Id(pbf) = *ibfcountthread;
#endif
        }
    }   // ! update
}

local void walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                        nodeptr *nptr, nodeptr *np, cellptr cptr, cellptr bptr,
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

local void sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                          cellptr cptr, cellptr bptr, bodyptr p0,
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
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1.histNNSubthread[n] = 0.0;
        hist1.histXi2pcfthreadsub[n] = 0.0;
    }
//#ifdef TPCF
    if (cmd->computeTPCF) {
#ifdef SINCOS
        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist1sincos.histNNSubthread[n] = 0.0;
        }
        CLRM_ext_ext(hist1sincos.histXithreadcos, cmd->mChebyshev+1,
                     cmd->sizeHistN);
        CLRM_ext_ext(hist1sincos.histXithreadsin, cmd->mChebyshev+1,
                     cmd->sizeHistN);
#else // ! SINCOS
        CLRM_ext_ext(hist1.histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
#endif
//#endif // ! TPCF
    }

//#ifdef TPCF
    if (cmd->computeTPCF) {
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
//#endif // ! TPCF
    }
//E
    if (!scanopt(cmd->options, "no-one-ball"))
        sumcell_balls6_omp(hist->interact, cptr, (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);
    sumnode_balls6_omp(bptr, hist->interact + hist->actlen, (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);

//B Section of type: computeBodyProperties_omp_balls6(p0, cmd->nbody, hist)
    real xi, xi_2p;

// BODY3
    if (Type(p0) == BODY) {
        xi = Kappa(p0)/nbody;
        xi_2p = Kappa(p0);
    } else if (Type(p0) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p0)*Kappa(p0)/nbody;
        xi_2p = Nbb(p0)*Kappa(p0);
#endif
    }
//
//#ifdef TPCF
    if (cmd->computeTPCF) {
#ifdef SINCOS
        int m;
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist1sincos.histXithreadcos[m][n] /=
                MAX(hist1sincos.histNNSubthread[n],1.0);
                hist1sincos.histXithreadsin[m][n] /=
                MAX(hist1sincos.histNNSubthread[n],1.0);
            }
        
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist1sincos.xiOUTVPcos, hist1sincos.histXithreadcos[m],
                      hist1sincos.histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist1sincos.xiOUTVPsin, hist1sincos.histXithreadsin[m],
                      hist1sincos.histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist1sincos.xiOUTVPsincos, hist1sincos.histXithreadsin[m],
                      hist1sincos.histXithreadcos[m],cmd->sizeHistN);
            CLRM_ext(hist1sincos.histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist1sincos.histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist1sincos.histZetaMtmpsincos,cmd->sizeHistN);
            MULMS_ext(hist1sincos.histZetaMtmpcos,hist1sincos.xiOUTVPcos,
                      xi,cmd->sizeHistN);
            MULMS_ext(hist1sincos.histZetaMtmpsin,hist1sincos.xiOUTVPsin,
                      xi,cmd->sizeHistN);
            MULMS_ext(hist1sincos.histZetaMtmpsincos,hist1sincos.xiOUTVPsincos,
                      xi,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadcos[m],
                     histsincos->histZetaMthreadcos[m],
                     hist1sincos.histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadsin[m],
                     histsincos->histZetaMthreadsin[m],
                     hist1sincos.histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadsincos[m],
                     histsincos->histZetaMthreadsincos[m],
                     hist1sincos.histZetaMtmpsincos,cmd->sizeHistN);
        }
#else // ! SINCOS
        int m;
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++)
                hist1.histXithread[m][n] /= MAX(hist1.histNNSubthread[n],1.0);
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist1.xiOUTVP, hist1.histXithread[m],
                      hist1.histXithread[m],cmd->sizeHistN);
            CLRM_ext(hist1.histZetaMtmp,cmd->sizeHistN);
            MULMS_ext(hist1.histZetaMtmp,hist1.xiOUTVP,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthread[m],hist->histZetaMthread[m],
                     hist1.histZetaMtmp,cmd->sizeHistN);
        }
#endif // ! SINCOS
//#endif // ! TPCF
    }
    for (n=1; n<=cmd->sizeHistN; n++) {
        hist->histXi2pcfthread[n] += xi_2p*hist1.histXi2pcfthreadsub[n];
    }
//E
//B
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist->histNthread[n] += hist1.histNthread[n];
        hist->histNNSubthread[n] += hist1.histNNSubthread[n];
        hist->histNNSubXi2pcfthread[n] += hist1.histNNSubXi2pcfthread[n];
    }
//E
    *nbbcalcthread += hist->interact + hist->actlen - bptr;
    *nbccalcthread += cptr - hist->interact;

    ip = p0 - bodytable[gd->iCatalogs[0]] + 1;
    if (ip%cmd->stepState == 0) {
        verb_log_print(cmd->verbose_log, gd->outlog, " - Completed pivot: %ld\n", ip);
    }

    search_free_omp_balls6_cc(cmd, gd, &hist1);
    search_free_sincos_omp_balls6(cmd, gd, &hist1sincos);
}


#define CHEBYSHEVOMP                                         \
{real xicosmphi; int m;                                      \
    hist->Chebs[1] = 1.0;                                    \
    xicosmphi = xi * hist->Chebs[1];                         \
    hist->histXithread[1][n] += xicosmphi;                   \
    hist->Chebs[2] = cosphi;                                 \
    xicosmphi = xi * hist->Chebs[2];                         \
    hist->histXithread[2][n] += xicosmphi;                   \
    hist->Chebs[3] = 2.0*(cosphi)*(cosphi) - (1.0);          \
    xicosmphi = xi * hist->Chebs[3];                         \
    hist->histXithread[3][n] += xicosmphi;                   \
    for (m=4; m<=cmd->mChebyshev+1; m++){                     \
        hist->Chebs[m] = 2.0*(cosphi)*hist->Chebs[m-1] - hist->Chebs[m-2];  \
        xicosmphi = xi * hist->Chebs[m];                     \
        hist->histXithread[m][n] += xicosmphi;               \
}}


local void sumnode_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
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
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(p);
                    hist->histNNSubXi2pcfthread[n] =
                                hist->histNNSubXi2pcfthread[n]+1.;
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                    xi = Kappa(pb);
//#ifdef TPCF
                    if (cmd->computeTPCF) {
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumnode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVOMP;
#else // ! SINCOS
                        real cosphi;
                        real sinphi;
                        histccsincos->histNNSubthread[n] =
                        histccsincos->histNNSubthread[n] + 1.;
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumnode: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVTUOMP;
#endif // ! SINCOS
//#endif // ! TPCF
                    }
                    hist->histXi2pcfthreadsub[n] += xi;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

local void sumcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              cellptr start, cellptr finish, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
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
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(pb);
                    hist->histNNSubXi2pcfthread[n] =
                            hist->histNNSubXi2pcfthread[n] + 1.;
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + 1.;
                    xi = Kappa(pb);
//#ifdef TPCF
                    if (cmd->computeTPCF) {
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumcell: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVOMP;
#else // ! SINCOS
                        real cosphi;
                        real sinphi;
                        histccsincos->histNNSubthread[n] =
                        histccsincos->histNNSubthread[n] + 1.;
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumcell: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVTUOMP;
#endif // ! SINCOS
//#endif // ! TPCF
                    }
                    hist->histXi2pcfthreadsub[n] += xi;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p
}

#undef CHEBYSHEVOMP

local void sumcellcell_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                  cellptr start, cellptr finish, nodeptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histccsincos)
{
    cellptr p;
#ifdef SINGLEP
    float dr1;
    float dr[NDIM];
#else
    vector dr;
    real dr1;
#endif
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
            if(dr1>cmd->rminHist) {
                ibodycount++;
                if (cmd->rminHist==0)
                    n = (int)(cmd->logHistBinsPD*(rlog10(dr1) - rlog10(cmd->rangeN))
                              + cmd->sizeHistN) + 1;
                else
                    n = (int)(rlog10(dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                if (n<=cmd->sizeHistN && n>=1) {
                    hist->histNthread[n] = hist->histNthread[n] + Nb(pb);
                    hist->histNNSubXi2pcfthread[n] =
                            hist->histNNSubXi2pcfthread[n] + 1.;
                    hist->histNNSubthread[n] = hist->histNNSubthread[n] + Nb(pb);
                    xj = Kappa(p0);
                    xi = Kappa(pb);
//#ifdef TPCF
                    if (cmd->computeTPCF) {
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumcellcell: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVOMPCC;
#else // ! SINCOS
                        real cosphi;
                        real sinphi;
                        histccsincos->histNNSubthread[n] =
                        histccsincos->histNNSubthread[n] + 1.;
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
                            verb_log_print(cmd->verbose, gd->outlog,
                                           "sumcellcell: Warning!... cossphi must be in (-1,1): %g\n",
                                           cosphi);
                        CHEBYSHEVTUOMP;
#endif // ! SINCOS
//#endif // ! TPCF
                    }
                    hist->histXi2pcfthreadsub[n] += xi*xj;
                } // ! 1 < n < HistN
            } // ! dr1 > rmin
        } // ! accept_body
    } // ! loop p

    *ncccalcthread += 1;
}

#endif // ! TREENODEBALLS4
//E BALLS4 :: DE METODO BALLS4 DE BUSQUEDA


//B Routines like in treeutils

//B BALLS4
global int search_init_omp_balls6(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_omp_balls6 hist, int ifile)
{
    int n;
    int m;

#  define FACTIVE  0.75
//#  define FACTOR  1
#  define FACTOR  316
//#  define FACTOR  1024

//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->Chebs = dvector(1,cmd->mChebyshev+1);
        hist->xiOUTVP = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->histXithread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histZetaMthread =
        dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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

//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthread[m], cmd->sizeHistN);
        }
        CLRM_ext_ext(hist->histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
//#endif
    }

    hist->actlen = FACTIVE * 216 * FACTOR * gd->tdepthTable[ifile];
    hist->actlen = hist->actlen * rpow(cmd->theta, -2.5);
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "searchcalc_balls6: actlen=%d\n",hist->actlen);
    hist->active = (nodeptr *) allocate(hist->actlen * sizeof(nodeptr));
    gd->bytes_tot += hist->actlen*sizeof(nodeptr);
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "\nAllocated %g MByte for active list storage.\n",
    hist->actlen*sizeof(nodeptr)/(1024.0*1024.0));
    hist->interact = (cellptr) allocate(hist->actlen * sizeof(cell));
    gd->bytes_tot += hist->actlen*sizeof(cell);
    verb_log_print(cmd->verbose_log, gd->outlog, 
                   "Allocated %g MByte for interact list storage.\n",
                   hist->actlen*sizeof(cell)/(1024.0*1024.0));

#undef FACTOR
#undef FACTIVE

    return SUCCESS;
}

global int search_init_sincos_omp_balls6(struct cmdline_data* cmd, 
                                         struct  global_data* gd,
                                         gdhistptr_sincos_omp_balls6 hist)
{
    int n;
    int m;

//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->ChebsT = dvector(1,cmd->mChebyshev+1);
        hist->ChebsU = dvector(1,cmd->mChebyshev+1);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->histXithreadcos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histXithreadsin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        
        hist->histZetaMthreadcos = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsin = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMthreadsincos = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        
        hist->xiOUTVPcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->xiOUTVPsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpcos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsin = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmpsincos = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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

//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthreadcos[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsin[m], cmd->sizeHistN);
            CLRM_ext(hist->histZetaMthreadsincos[m], cmd->sizeHistN);
        }
//#endif
    }

    return SUCCESS;
}

global int search_init_omp_balls6_cc(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                         gdhistptr_omp_balls6 hist)
{
    int n;
    int m;

//    #ifdef TPCF
    if (cmd->computeTPCF) {
        hist->Chebs = dvector(1,cmd->mChebyshev+1);
        hist->xiOUTVP = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->histXithread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histZetaMthread = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
    }

        for (n = 1; n <= cmd->sizeHistN; n++) {
            hist->histNthread[n] = 0.0;
            hist->histNNSubthread[n] = 0.0;
// 2pcf
            hist->histNNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
            hist->histNNSubXi2pcfthreadp[n] = 0.0;
            hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
//E
//
            hist->histXi2pcfthread[n] = 0.0;
            hist->histXi2pcfthreadsub[n] = 0.0;
        }

//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthread[m], cmd->sizeHistN);
        }
        CLRM_ext_ext(hist->histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
//#endif
    }

        return SUCCESS;
    }

global int computeBodyProperties_omp_balls6(struct cmdline_data* cmd,
                                    struct  global_data* gd,
                                    bodyptr p, int nbody, gdhistptr_omp_balls6 hist)
{
    int n;
    int m;
    real xi, xi_2p;

//B BODY3
    if (Type(p) == BODY) {
        xi = Kappa(p)/nbody;
        xi_2p = Kappa(p);
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    }
//E
//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++)
                hist->histXithread[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVP, hist->histXithread[m], hist->histXithread[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmp,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthread[m],hist->histZetaMthread[m],hist->histZetaMtmp,cmd->sizeHistN);
        }
//#endif
    }
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
        }

        return SUCCESS;
}

global int computeBodyProperties_omp_balls6_cc(struct cmdline_data* cmd, 
                                               struct  global_data* gd,
                                               bodyptr p, INTEGER nbody,
                                               gdhistptr_omp_balls6 hist)
{
    int n;
    int m;
    real xi, xi_2p;

// BODY3
    if (Type(p) == BODY) {
        xi = Kappa(p)/nbody;
        xi_2p = 1.0;
    } else if (Type(p) == BODY3) {
#ifdef BODY3ON
        xi = Nbb(p)*Kappa(p)/nbody;
        xi_2p = Nbb(p)*Kappa(p);
#endif
    } else {
        xi = Kappa(p)/nbody;
        xi_2p = 1.0;
    }
//
//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++)
                hist->histXithread[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (m=1; m<=cmd->mChebyshev+1; m++){
                OUTVP_ext(hist->xiOUTVP, hist->histXithread[m], hist->histXithread[m],cmd->sizeHistN);
                CLRM_ext(hist->histZetaMtmp,cmd->sizeHistN);
                MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd->sizeHistN);
                ADDM_ext(hist->histZetaMthread[m],hist->histZetaMthread[m],hist->histZetaMtmp,cmd->sizeHistN);
            }
//#endif
    }
        for (n=1; n<=cmd->sizeHistN; n++) {
            hist->histXi2pcfthread[n] += xi_2p*hist->histXi2pcfthreadsub[n];
        }

    return SUCCESS;
}

global int computeBodyProperties_sincos_omp_balls6_cc(struct cmdline_data* cmd, 
                                                      struct  global_data* gd,
                                                      bodyptr p, INTEGER nbody,
                                            gdhistptr_sincos_omp_balls6 histsincos)
{
        int n;
        int m;
        real xi, xi_2p;

    // BODY3
            if (Type(p) == BODY) {
                xi = Kappa(p)/nbody;
                xi_2p = 1.0;
            } else if (Type(p) == BODY3) {
#ifdef BODY3ON
                xi = Nbb(p)*Kappa(p)/nbody;
                xi_2p = Nbb(p)*Kappa(p);
#endif
            } else {
                xi = Kappa(p)/nbody;
                xi_2p = 1.0;
            }
//
//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++) {
                histsincos->histXithreadcos[m][n] /=
                MAX(histsincos->histNNSubthread[n],1.0);
                histsincos->histXithreadsin[m][n] /=
                MAX(histsincos->histNNSubthread[n],1.0);
            }
        
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(histsincos->xiOUTVPcos, histsincos->histXithreadcos[m],
                      histsincos->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(histsincos->xiOUTVPsin, histsincos->histXithreadsin[m],
                      histsincos->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(histsincos->xiOUTVPsincos, histsincos->histXithreadsin[m],
                      histsincos->histXithreadcos[m],cmd->sizeHistN);
            CLRM_ext(histsincos->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(histsincos->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(histsincos->histZetaMtmpsincos,cmd->sizeHistN);
            MULMS_ext(histsincos->histZetaMtmpcos,histsincos->xiOUTVPcos,
                      xi,cmd->sizeHistN);
            MULMS_ext(histsincos->histZetaMtmpsin,histsincos->xiOUTVPsin,
                      xi,cmd->sizeHistN);
            MULMS_ext(histsincos->histZetaMtmpsincos,histsincos->xiOUTVPsincos,
                      xi,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadcos[m],
                     histsincos->histZetaMthreadcos[m],
                     histsincos->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadsin[m],
                     histsincos->histZetaMthreadsin[m],
                     histsincos->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(histsincos->histZetaMthreadsincos[m],
                     histsincos->histZetaMthreadsincos[m],
                     histsincos->histZetaMtmpsincos,cmd->sizeHistN);
        }
//#endif
    }
        for (n=1; n<=cmd->sizeHistN; n++) {
            histsincos->histXi2pcfthread[n] += xi_2p*histsincos->histXi2pcfthreadsub[n];
        }

    return SUCCESS;
}


global int search_free_omp_balls6(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_omp_balls6 hist)
{
        free(hist->interact);
        free(hist->active);

//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix3D(hist->histZetaMthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//#endif
    }
        free_dvector(hist->histXi2pcfthreadsub,1,cmd->sizeHistN);
        free_dvector(hist->histXi2pcfthread,1,cmd->sizeHistN);
    // 2pcf
        //B kappa Avg Rmin
            free_dvector(hist->histNNSubXi2pcfthreadtotal,1,cmd->sizeHistN);
            free_dvector(hist->histNNSubXi2pcfthreadp,1,cmd->sizeHistN);
        //E
        free_dvector(hist->histNNSubXi2pcfthread,1,cmd->sizeHistN);
    //
        free_dvector(hist->histNNSubthread,1,cmd->sizeHistN);
        free_dvector(hist->histNthread,1,cmd->sizeHistN);
//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVP,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dvector(hist->Chebs,1,cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_sincos_omp_balls6(struct cmdline_data* cmd, 
                                         struct  global_data* gd,
                                         gdhistptr_sincos_omp_balls6 hist)
{
//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmpsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histZetaMtmpcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsincos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPsin,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVPcos,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsincos,1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(hist->histZetaMthreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadsin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithreadcos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//#endif
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
//    #ifdef TPCF
    if (cmd->computeTPCF) {
        free_dvector(hist->ChebsU,1,cmd->mChebyshev+1);
        free_dvector(hist->ChebsT,1,cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_omp_balls6_cc(struct cmdline_data* cmd, 
                                     struct  global_data* gd,
                                     gdhistptr_omp_balls6 hist)
{
//    #ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix3D(hist->histZetaMthread,1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//#endif
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
//    #ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVP,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dvector(hist->Chebs,1,cmd->mChebyshev+1);
//#endif
    }
        
    return SUCCESS;
}

//E BALLS4


// #ifdef BALLS

global int search_init_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                 gdhistptr_omp_balls hist, int ifile)
{
     int n, m;

 #  define FACTIVE  0.75
 //#  define FACTOR  1
 #  define FACTOR  316
 //#  define FACTOR  1024

// #ifdef TPCF
    if (cmd->computeTPCF) {
        hist->Chebs = dvector(1,cmd->mChebyshev+1);
        hist->xiOUTVP = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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
// #ifdef TPCF
    if (cmd->computeTPCF) {
        hist->histXithread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histZetaMthread = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
    }

     for (n = 1; n <= cmd->sizeHistN; n++) {
         hist->histNthread[n] = 0.0;
         hist->histNNSubthread[n] = 0.0;
// 2pcf
         hist->histNNSubXi2pcfthread[n] = 0.0;
//B kappa Avg Rmin
         hist->histNNSubXi2pcfthreadp[n] = 0.0;
         hist->histNNSubXi2pcfthreadtotal[n] = 0.0;
//E
//
         hist->histXi2pcfthread[n] = 0.0;
         hist->histXi2pcfthreadsub[n] = 0.0;
     }

// #ifdef TPCF
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthread[m], cmd->sizeHistN);
        }
        CLRM_ext_ext(hist->histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
//#endif
    }

//    verb_print_debug(1, "\nAqui voy (12.0): %d %d %d %ld %d\n",
//                     cmd->mChebyshev, cmd->sizeHistN,
//                     cmd->computeTPCF, hist->actlen, gd->tdepthTable[ifile]);
     hist->actlen = FACTIVE * 216 * FACTOR * gd->tdepthTable[ifile];
//    verb_print_debug(1, "\nAqui voy (12): %d %d %d %ld\n",
//                     cmd->mChebyshev, cmd->sizeHistN,
//                     cmd->computeTPCF, hist->actlen);
     hist->actlen = hist->actlen * rpow(cmd->theta, -2.5);
//    verb_print_debug(1, "\nAqui voy (12.1): %d %d %d %ld\n",
//                     cmd->mChebyshev, cmd->sizeHistN,
//                     cmd->computeTPCF, hist->actlen);
     verb_log_print(cmd->verbose_log, gd->outlog,
                    "search_init_balls_omp: actlen=%d\n",hist->actlen);
     hist->active = (nodeptr *) allocate(hist->actlen * sizeof(nodeptr));
//    verb_print_debug(1, "\nAqui voy (13): %d %d\n",
//                     cmd->mChebyshev, cmd->sizeHistN);
     gd->bytes_tot += hist->actlen*sizeof(nodeptr);
     verb_log_print(cmd->verbose_log, gd->outlog, 
                    "\nAllocated %g MByte for active list storage.\n",
                    hist->actlen*sizeof(nodeptr)/(1024.0*1024.0));
     hist->interact = (cellptr) allocate(hist->actlen * sizeof(cell));
     gd->bytes_tot += hist->actlen*sizeof(cell);
     verb_log_print(cmd->verbose_log, gd->outlog, 
                    "Allocated %g MByte for interact list storage.\n",
                    hist->actlen*sizeof(cell)/(1024.0*1024.0));

 #undef FACTOR
 #undef FACTIVE

     return SUCCESS;
}

global int search_init_balls_omp_cc(struct cmdline_data* cmd, 
                                    struct  global_data* gd,
                                    gdhistptr_omp_balls hist)
{
     int n, m;

//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->Chebs = dvector(1,cmd->mChebyshev+1);
        hist->xiOUTVP = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
        hist->histZetaMtmp = dmatrix(1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        hist->histXithread = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        hist->histZetaMthread =
        dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
//#endif
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

//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m = 1; m <= cmd->mChebyshev+1; m++) {
            CLRM_ext(hist->histZetaMthread[m], cmd->sizeHistN);
        }
        CLRM_ext_ext(hist->histXithread, cmd->mChebyshev+1, cmd->sizeHistN);
//#endif
    }

    return SUCCESS;
}

global int computeBodyProperties_balls_omp(struct cmdline_data* cmd, 
                                           struct  global_data* gd,
                                           bodyptr p, int nbody,
                                           gdhistptr_omp_balls hist)
{
    int n;
    real xi_p;

#ifdef TREENODEALLBODIES
//B kappa Avg Rmin
    xi_p = Kappa(p);
    if (scanopt(cmd->options, "smooth-pivot")) {
         xi_p = KappaRmin(p);
//         xi_p = NbRmin(p);
    }
//     printf("%e %e\n", xi_p, Kappa(p));
//E
#else
    xi_p = 1.0;
#endif

    for (n=1; n<=cmd->sizeHistN; n++)
        hist->histXi2pcfthread[n] += xi_p*hist->histXi2pcfthreadsub[n];

    return SUCCESS;
}

global int computeBodyProperties_balls_omp_cc(struct cmdline_data* cmd, 
                                              struct  global_data* gd,
                                              bodyptr p, int nbody,
                                              gdhistptr_omp_balls hist)
{
     int n;
     int m;
     real xi;

     xi = 1.0/((real)nbody);

 #ifdef TREENODEALLBODIES
     xi *= Kappa(p);
 #else
 #ifdef TREENODEBALLS4
     xi *= Kappa(p);
 #else
     if (Kappa(p) == 0) {
         verb_log_print(cmd->verbose, gd->outlog,
         "Warning!... Kappa: %e\n",Kappa(p));
     xi = 0.0;
     } else {
         xi *= ((real)Nb(p))*Kappa(p);
     }
 #endif // ! TREENODEBALLS4
 #endif // ! TREENODEALLBODIES
 //
//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++)
                hist->histXithread[m][n] /= MAX(hist->histNNSubthread[n],1.0);
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVP, hist->histXithread[m],
                      hist->histXithread[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmp,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthread[m],hist->histZetaMthread[m],
                     hist->histZetaMtmp,cmd->sizeHistN);
        }
//#endif
    }

    return SUCCESS;
}

global int computeBodyProperties_balls_omp_cc_sincos(struct cmdline_data* cmd,
                                                     struct  global_data* gd,
                                                     bodyptr p, int nbody,
                                                     gdhistptr_sincos_omp hist)
{
     int n;
     int m;
     real xi;

 #ifdef TREENODEALLBODIES
     xi = Kappa(p)/nbody;
//B kappa Avg Rmin
    if (scanopt(cmd->options, "smooth-pivot")) {
        xi = NbRmin(p)*KappaRmin(p)/nbody;
    }
//E
 #else
 #ifdef TREENODEBALLS4
     xi = Kappa(p)/nbody;
 #else
     if (Kappa(p) == 0) {
         verb_log_print(cmd->verbose, gd->outlog,
         "Warning!... Kappa: %e\n",Kappa(p));
     xi = 0.0;
     } else {
         if (Type(p) == CELL) {
         xi = (1.0/nbody)*Nb(p)/Kappa(p);
         } else {
             xi = (1.0/nbody)/Kappa(p);
         }
     }
 #endif // ! TREENODEBALLS4
 #endif // ! TREENODEALLBODIES

//#ifdef TPCF
    if (cmd->computeTPCF) {
        for (m=1; m<=cmd->mChebyshev+1; m++)
            for (n=1; n<=cmd->sizeHistN; n++) {
                hist->histXithreadcos[m][n] /= MAX(hist->histNNSubthread[n],1.0);
                hist->histXithreadsin[m][n] /= MAX(hist->histNNSubthread[n],1.0);
            }
        for (m=1; m<=cmd->mChebyshev+1; m++){
            OUTVP_ext(hist->xiOUTVPcos, hist->histXithreadcos[m],
                      hist->histXithreadcos[m], cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsin, hist->histXithreadsin[m],
                      hist->histXithreadsin[m],cmd->sizeHistN);
            OUTVP_ext(hist->xiOUTVPsincos, hist->histXithreadsin[m],
                      hist->histXithreadcos[m],cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpcos,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsin,cmd->sizeHistN);
            CLRM_ext(hist->histZetaMtmpsincos,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpcos,hist->xiOUTVPcos,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsin,hist->xiOUTVPsin,xi,cmd->sizeHistN);
            MULMS_ext(hist->histZetaMtmpsincos,hist->xiOUTVPsincos,xi,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadcos[m],hist->histZetaMthreadcos[m],
                     hist->histZetaMtmpcos,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsin[m],hist->histZetaMthreadsin[m],
                     hist->histZetaMtmpsin,cmd->sizeHistN);
            ADDM_ext(hist->histZetaMthreadsincos[m],hist->histZetaMthreadsincos[m],
                     hist->histZetaMtmpsincos,cmd->sizeHistN);
        }
//#endif
    }

    return SUCCESS;
}

global int search_free_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                                 gdhistptr_omp_balls hist)
{
     free(hist->interact);
     free(hist->active);

//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix3D(hist->histZetaMthread,1,cmd->mChebyshev+1,1,
                       cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVP,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dvector(hist->Chebs,1,cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_balls_omp_cc(struct cmdline_data* cmd, 
                                    struct  global_data* gd,
                                    gdhistptr_omp_balls hist)
{
//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix3D(hist->histZetaMthread,1,cmd->mChebyshev+1,1,
                       cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->histXithread,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
//#endif
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
//#ifdef TPCF
    if (cmd->computeTPCF) {
        free_dmatrix(hist->histZetaMtmp,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(hist->xiOUTVP,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dvector(hist->Chebs,1,cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

//B 2023.11.22
global bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                  nodeptr p, nodeptr q, real *dr1, vector dr)
{
//    real drpq, drpq2;
/*
     DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
 #ifdef PERIODIC
     VWrapAll(dr);
     DOTVP(drpq2, dr, dr);
 #endif
     drpq = rsqrt(drpq2);
     *dr1 = drpq;
*/
     int n;

     if ( *dr1 == 0.0)
         return (FALSE);
     else
         if ( (Radius(p)+Radius(q))/(*dr1) < gd->deltaR) {
             if (scanopt(cmd->options, "behavior-tree-omp")) {
//B To behaves as tree-omp
                 if ( (*dr1)<gd->Rcut ) {
                     if((*dr1)>cmd->rminHist) {
                         if (cmd->rminHist==0)
                             n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) -
                                     rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                         else
                             n = (int)(rlog10((*dr1)/cmd->rminHist)
                                       * gd->i_deltaR) + 1;
                         if (n<=cmd->sizeHistN-1 && n>=1) {
                             if ( gd->deltaRV[n] < *dr1 - Radius(q) && *dr1
                                 + Radius(q) < gd->deltaRV[n+1]) {
                                 return (TRUE);
                             } else {
                                 return (FALSE);
                             }
                         } else
                             return (FALSE);
                     } else
                         return (FALSE);
                 } else
                     return (FALSE);
 //E
             } else { // ! behavior-tree-omp
                 return (TRUE);
             }
         } else
             return (FALSE);
}
//E

global bool nodes_condition(struct cmdline_data* cmd, struct  global_data* gd,
                            nodeptr p, nodeptr q, real *dr1, vector dr)
{
     real drpq, drpq2;

     DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
// #ifdef PERIODIC
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
//#endif
    }
     drpq = rsqrt(drpq2);
     *dr1 = drpq;

     int n;

     if ( drpq == 0.0)
         return (FALSE);
     else
         if ( (Radius(p)+Radius(q))/drpq < gd->deltaR) {
             if (scanopt(cmd->options, "behavior-tree-omp")) {
 //B To behaves as tree-omp
             if ( *dr1<gd->Rcut ) {
                 if(*dr1>cmd->rminHist) {
                     if (cmd->rminHist==0)
                         n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) - rlog10(cmd->rangeN)) + cmd->sizeHistN) + 1;
                     else
                         n = (int)(rlog10(*dr1/cmd->rminHist) * gd->i_deltaR) + 1;
                     if (n<=cmd->sizeHistN-1 && n>=1) {
                         if ( gd->deltaRV[n] < *dr1 - Radius(q) && *dr1 + Radius(q) < gd->deltaRV[n+1]) {
                             return (TRUE);
                         } else {
                             return (FALSE);
                         }
                     } else
                         return (FALSE);
                 } else
                     return (FALSE);
             } else
                 return (FALSE);
 //E
             } else {
                 return (TRUE);
             }
         } else
             return (FALSE);
}

global bool nodes_condition5(struct cmdline_data* cmd, struct  global_data* gd,
                             nodeptr p, nodeptr q)
{
     real drpq, drpq2;
     vector dr;

     DOTPSUBV(drpq2, dr, Pos(p), Pos(q));
// #ifdef PERIODIC
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(drpq2, dr, dr);
//#endif
    }
     drpq = rsqrt(drpq2);

     if ( drpq == 0.0)
         return (FALSE);
     else
         if ( (Radius(p)+Radius(q))/drpq < gd->deltaR)
             return (TRUE);
         else
             return (FALSE);
}

global bool nodes_set_bin(struct cmdline_data* cmd, struct  global_data* gd,
                          nodeptr p, nodeptr q, int *n, real *dr1, vector dr)
{
     *n=-1;

     if ( *dr1<gd->Rcut ) {
         if(*dr1>cmd->rminHist) {
             if (cmd->rminHist==0)
                 *n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) - rlog10(cmd->rangeN))
                            + cmd->sizeHistN) + 1;
             else
                 *n = (int)(rlog10(*dr1/cmd->rminHist) * gd->i_deltaR) + 1;
             if (*n<=cmd->sizeHistN && *n>=1)
                     return (TRUE);
             else
                 return (FALSE);
         } else
             return (FALSE);
     } else
         return (FALSE);
}

#ifdef SINGLEP
global bool nodes_set_bin5(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, int *n, float *dr1, float *dr)
#else
    global bool nodes_set_bin5(struct cmdline_data* cmd, struct  global_data* gd,
                               nodeptr p, nodeptr q, int *n, real *dr1, vector dr)
#endif
{
    *n=-1;

    if (accept_body(cmd, gd, (bodyptr)p, q, dr1, dr)) {
         if(*dr1>cmd->rminHist) {
             if (cmd->rminHist==0)
                 *n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) - rlog10(cmd->rangeN)) 
                            + cmd->sizeHistN) + 1;
             else
                 *n = (int)(rlog10(*dr1/cmd->rminHist) * gd->i_deltaR) + 1;
             if (*n<=cmd->sizeHistN && *n>=1) {
                     return (TRUE);
             } else
                 return (FALSE);
         } else
             return (FALSE);
     } else
         return (FALSE);
}

#ifdef SINGLEP
global bool nodes_set_bin_bc(struct cmdline_data* cmd, struct  global_data* gd,
                             nodeptr p, nodeptr q, int *n, float *dr1, float *dr)
#else
global bool nodes_set_bin_bc(struct cmdline_data* cmd, struct  global_data* gd,
                             nodeptr p, nodeptr q, int *n, real *dr1, vector dr)
#endif
 {
     *n=-1;

     if (accept_body(cmd, gd, (bodyptr)p, q, dr1, dr)) {
         if(*dr1>cmd->rminHist) {
             if (cmd->rminHist==0)
                 *n = (int)(cmd->logHistBinsPD*(rlog10(*dr1) - rlog10(cmd->rangeN)) + cmd->sizeHistN);
             else
                 *n = (int)(rlog10(*dr1/cmd->rminHist) * gd->i_deltaR);
             if (*n<=cmd->sizeHistN-1 && *n>=1) {
                     return (TRUE);
             } else
                 return (FALSE);
         } else
             return (FALSE);
     } else
         return (FALSE);
 }

// #endif // ! BALLS

local int search_compute_Xi_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                  int nbody)
{
        int k;
        int n;
        real normFac;
        real Vol;

        Vol = 1.0;
        DO_COORD(k)
            Vol = Vol*gd->Box[k];

    #ifdef RAPHISTNVER
    // Check for log-scale value of deltaR
        if (NDIM == 3)
            normFac = Vol / (2.0 * PI * rpow(gd->deltaR, 3.0) * nbody * nbody);
        else if (NDIM == 2)
            normFac = Vol / (PI * rpow(gd->deltaR, 2.0) * nbody * nbody);
        else error("\n\nWrong NDIM!\n\n");
    #endif

    //#ifdef TREENODE
    //    normFac *= 0.5;
    //#else
    //    if (scanopt(cmd->options, "compute-j-no-eq-i"))
    //        normFac /= 2.0;
    //    else
    //        normFac /= 1.0;
    //#endif

        if (gd->nnodescanlev == gd->nnodescanlev_root)
            normFac *= 0.5;

    #ifdef RAPHISTNVER
    // Check for log-scale
        for (n = 1; n <= cmd->sizeHistN; n++)
            if (NDIM == 3)
            gd->histCF[n] = gd->histNN[n] * normFac / rsqr(n-0.5); // add "- 1.0" to have zeta = rho/rho_avg - 1
            else if (NDIM == 2)
                gd->histCF[n] = gd->histNN[n] * normFac / ((int)n-0.5);
            else error("\n\nWrong NDIM!\n\n");
    #endif

//    #ifndef LOGHIST
    if (cmd->useLogHist) {
        gd->histNN[1]-=nbody;
    }
//    #endif
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
          }
          else {
            double r0,r1,vr,rho_r;
//    #ifdef LOGHIST
              if (cmd->useLogHist) {
                  if (cmd->rminHist==0) {
                      r0 = rpow(10.0, ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD + rlog10(cmd->rangeN) );
                      r1 = rpow(10.0, ((real)(n+1-cmd->sizeHistN))/cmd->logHistBinsPD + rlog10(cmd->rangeN) );
                  } else {
                      r0 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n))*gd->deltaR );
                      r1 = rpow(10.0, rlog10(cmd->rminHist) + ((real)(n+1))*gd->deltaR );
                  }
//#else
              } else {
                  r0=(real)n*gd->deltaR;
                  r1=(real)(n+1)*gd->deltaR;
//    #endif
              }
            vr=4.0*PI*(r1*r1*r1-r0*r0*r0)/3.0;
              rho_r=gd->histNN[n]/((real)nbody*vr);
              corr[n]=rho_r/rho_av-1;
              ercorr[n]=(1+corr[n])*edd[n];
              gd->histCF[n] = corr[n] + 1.0;
          }
        }

    return SUCCESS;
}

global int search_compute_HistN_balls(struct cmdline_data* cmd, 
                                      struct  global_data* gd,
                                      int nbody)
{
    int n;
    real normFac;

#ifdef TREENODEBALLS4
    normFac = 0.5;
#else
    normFac = 0.5;
    if (gd->nnodescanlev == gd->nnodescanlev_root)
        normFac = 0.5;
//        normFac = 1.0;              // Check
#endif

//#ifndef LOGHIST
    if (!cmd->useLogHist) {
        //#ifdef KDLIB
        //    gd->histNN[1]-=nbody; //Substract diagonal
        // Check for no KD routines like normal tree...
        //#endif
//#endif
    }

// CHECK FOR LOGHIST!!!! (Must be correct, but...) 
// Checar también el conteo de las bolas...

    for (n = 1; n <= cmd->sizeHistN; n++)
        gd->histNN[n] *= normFac;

    if (scanopt(cmd->options, "and-CF"))
        search_compute_Xi_balls(cmd, gd, nbody);

    return SUCCESS;
}

    
//E Routines like in treeutils

#endif // ! BALLS
