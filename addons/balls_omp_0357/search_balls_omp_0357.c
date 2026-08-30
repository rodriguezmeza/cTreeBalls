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

local bool nodes_condition_balls(struct cmdline_data* cmd,
                                  struct  global_data* gd,
                                  nodeptr p, nodeptr q, real *dr1, compute_vector dr);

//#ifdef BALLS
#if defined(OPENMPCODE) || defined(BALLS)

local void walktree_balls_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr, nodeptr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *, INTEGER *);
//B 2023.11.22
local void sumnodes_bb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *dr1, compute_vector dr,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_bc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, compute_vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, compute_vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);
local void sumnodes_cc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr, nodeptr, real *, compute_vector,
            gdhistptr_omp_balls, gdhistptr_omp_balls, gdhistptr_sincos_omp,
            INTEGER *, INTEGER *, INTEGER *);

#ifdef TREENODEBALLS4
//B BALLS4: TREENODEBALLS4
local int walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                              nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *,  INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local int walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                       nodeptr *, nodeptr *, cellptr, cellptr,
            nodeptr, real, vector, INTEGER *, INTEGER *, INTEGER *,
            gdhistptr_omp_balls6, gdhistptr_omp_balls6,
            gdhistptr_sincos_omp_balls6, INTEGER *, INTEGER *, INTEGER);
local int sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
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

typedef enum {
    BALLS0357_INIT_BALLS,
    BALLS0357_INIT_BALLS_CC,
    BALLS0357_INIT_BALLS6,
    BALLS0357_INIT_BALLS6_CC,
    BALLS0357_INIT_SINCOS_BALLS6
} balls0357_init_kind;

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    void *hist;
    int ifile;
    balls0357_init_kind kind;
} balls0357_init_context;

local int search_init_balls_0357_callback(void *argument)
{
    balls0357_init_context *context = argument;

    switch (context->kind) {
        case BALLS0357_INIT_BALLS:
            return search_init_balls_omp(context->cmd, context->gd,
                (gdhistptr_omp_balls)context->hist, context->ifile);
        case BALLS0357_INIT_BALLS_CC:
            return search_init_balls_omp_cc(context->cmd, context->gd,
                (gdhistptr_omp_balls)context->hist);
        case BALLS0357_INIT_BALLS6:
            return search_init_omp_balls6(context->cmd, context->gd,
                (gdhistptr_omp_balls6)context->hist, context->ifile);
        case BALLS0357_INIT_BALLS6_CC:
            return search_init_omp_balls6_cc(context->cmd, context->gd,
                (gdhistptr_omp_balls6)context->hist);
        case BALLS0357_INIT_SINCOS_BALLS6:
            return search_init_sincos_omp_balls6(context->cmd, context->gd,
                (gdhistptr_sincos_omp_balls6)context->hist);
    }

    return FAILURE;
}

local int search_init_balls_0357_guarded(struct cmdline_data *cmd,
                                         struct global_data *gd,
                                         void *hist, int ifile,
                                         balls0357_init_kind kind)
{
    balls0357_init_context context;
    ErrorMsg allocation_error;

    context.cmd = cmd;
    context.gd = gd;
    context.hist = hist;
    context.ifile = ifile;
    context.kind = kind;

    return cballs_allocation_guard(search_init_balls_0357_callback,
                                   &context, allocation_error,
                                   sizeof(allocation_error));
}

// search=balls-omp
global int searchcalc_balls_omp_0357(struct cmdline_data* cmd,
                                    struct global_data* gd,
                                    bodyptr *btable,
                                    INTEGER *nbody, INTEGER ipmin,
                                    INTEGER *ipmax, int cat1, int cat2)
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
    if (cmd->useLogHist == FALSE)
        cBALLS_FAIL(cmd,
            "CheckParameters: can´t have normal hist and BALLS0357 "
            "definition (useLogHist=%d)\nSet useLogHist = true\n",
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
    if (!gd->computeTPCF)
    verb_print(cmd->verbose, "sincos base... \n");
//#endif
#endif
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose, "with option no-two-balls... \n");
    if (cballs_opt_no_one_ball(cmd))
        verb_print(cmd->verbose, "with option no-one-ball... \n");
    if (cballs_opt_behavior_tree_omp(cmd))
        verb_print(cmd->verbose, "with option behavior-tree-omp... \n");
    if (cballs_opt_center_of_mass(cmd))
    verb_print(cmd->verbose, "with option center-of-mass... \n");

#ifdef TREENODEALLBODIES
    if (cballs_opt_smooth_pivot(cmd))
    verb_print(cmd->verbose, 
               "with option smooth-pivot... rsmooth=%g\n",gd->rsmooth[0]);
#endif

//#ifndef TPCF
    if (!gd->computeTPCF)
    verb_print(cmd->verbose, "computing only 2pcf... \n");
//#endif

    if (cmd->verbose==3)
    verb_print(cmd->verbose,
        "\ncat1, cat2, nbody[cat1], nbody[cat2], ipmin, imax[cat1], ipmax[cat1]:\
        %d %d %ld %ld %ld %ld %ld\n\n",
        cat1,cat2,nbody[cat1],nbody[cat2],ipmin,ipmax[cat1],ipmax[cat2]);

#if defined(OPENMPCODE)
    if (ThreadCount(cmd, gd, nbody[cat1], cat1) == FAILURE)
        return FAILURE;
#endif

// Here we clear: histZetaM, histXi, histN, histNNSubXi2pcf, histXi2pcf...
//    search_init_gd_hist();
#ifdef SINCOS
    if (search_init_gd_hist_sincos(cmd, gd) == FAILURE)
        return FAILURE;
#else
    if (search_init_gd_hist(cmd, gd) == FAILURE)
        return FAILURE;
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
    int allocation_failed = FALSE;
    int worker_failed = FALSE;

// The OpenMP shared lists are written explicitly below for each tree profile.

//B BALLS4
#ifdef TREENODEBALLS4
#ifdef DEBUG
#pragma omp parallel default(none)  \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss, \
    bodytabbf,ibfcount,ihit,cat1,cat2,allocation_failed,worker_failed)
#else
#pragma omp parallel default(none) \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss,cat1,cat2,allocation_failed,worker_failed)
#endif
#else // ! TREENODEBALLS4
#ifdef TREENODEALLBODIES
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev, \
           nodetablescanlev_root,rootnode,imiss, cat1, cat2, ipfalse, \
           icountNbRmin,icountNbRminOverlap,allocation_failed,worker_failed)
#else
#pragma omp parallel default(none)   \
    shared(cmd,gd,btable,nbody,roottable,ipmin,ipmax,nodetablescanlev,nodetablescanlev_root,rootnode,imiss,cat1,cat2,allocation_failed,worker_failed)
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

      memset(&histb, 0, sizeof(histb));
      memset(&histccb, 0, sizeof(histccb));
      memset(&histbsincos, 0, sizeof(histbsincos));
      int hist_ready =
          search_init_balls_0357_guarded(cmd, gd, &histb, cat1,
                                        BALLS0357_INIT_BALLS6) == SUCCESS &&
          search_init_balls_0357_guarded(cmd, gd, &histccb, cat1,
                                        BALLS0357_INIT_BALLS6_CC) == SUCCESS &&
          search_init_balls_0357_guarded(cmd, gd, &histbsincos, cat1,
                               BALLS0357_INIT_SINCOS_BALLS6) == SUCCESS;
#else
      gdhist_omp_balls hist;
      gdhist_omp_balls histcc;
      gdhist_sincos_omp histsincos;

      memset(&hist, 0, sizeof(hist));
      memset(&histcc, 0, sizeof(histcc));
      memset(&histsincos, 0, sizeof(histsincos));
      int hist_ready =
          search_init_sincos_omp(cmd, gd, &histsincos) == SUCCESS &&
          search_init_balls_0357_guarded(cmd, gd, &hist, cat1,
                                        BALLS0357_INIT_BALLS) == SUCCESS &&
          search_init_balls_0357_guarded(cmd, gd, &histcc, cat1,
                                        BALLS0357_INIT_BALLS_CC) == SUCCESS;

  // Here we clear threads of: histZetaM, histXi, histNNSub,
  //  histN, histNNSubXi2pcf, histXi2pcf, histXi2pcf()sub...
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

    if (!hist_ready) {
#pragma omp atomic write
        allocation_failed = TRUE;
    }

#pragma omp barrier

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
      if (gd->computeTPCF) {
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
            if (allocation_failed) continue;
// Correct it because p and q are in differents node structures!!!... cat1!=cat2
//B kappa Avg Rmin
            NbRmin(p) = 1;
            NbRminOverlap(p) = 0;
            KappaRmin(p) = Kappa(p);
            if (cballs_opt_smooth_pivot(cmd)) {
                if (Update(p) == FALSE) {
                    ipfalsethreads++;
                    continue;
                }
            }
//E
#else // TREENODEALLBODIES exclude TREENODEBALLS4. Then this is for TREENODEBALLS4
        for (i=0; i< gd->nnodescanlevTable[cat1]; i++) {
            if (allocation_failed) continue;
            if (cmd->scanLevel==0)
                p = (nodeptr) roottable[cat1];
            else
                p = nodetablescanlev[cat1][i];
#endif

#ifdef TREENODEALLBODIES
            hist.ipcount = 0;
#endif

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
            if (gd->computeTPCF) {
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
            if (gd->computeTPCF) {
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
                    if (cballs_opt_smooth_pivot(cmd)) {
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
                    if (cballs_opt_smooth_pivot(cmd)) {
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
      if (computeBodyProperties_balls_omp(cmd, gd, (bodyptr)p,
                                          nbody[cat1], &hist) == FAILURE) {
#pragma omp atomic write
          worker_failed = TRUE;
      }
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


            if (walktree_balls6_omp(cmd, gd,
                    histb.active, histb.active + 1,
                    histb.interact, histb.interact + histb.actlen,
                    p, Size(p), Pos(p),
                    &nbbcalcthread, &nbccalcthread, &ncccalcthread,
                    &histb, &histccb, &histbsincos,
                    &ibfcountthread, &nsmoothcountthread,
                    nbody[cat1]) == FAILURE) {
#pragma omp atomic write
                worker_failed = TRUE;
                continue;
            }
//            }

            if (!cballs_opt_no_two_balls(cmd)) {
#ifndef SINCOS
                if (computeBodyProperties_omp_balls6_cc(cmd, gd, (bodyptr)p,
                                      nbody[cat1], &histccb) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
#else
                if (computeBodyProperties_sincos_omp_balls6_cc(cmd, gd,
                              (bodyptr)p, nbody[cat1], &histbsincos) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
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
            if (gd->computeTPCF) {
#ifndef PIVOTEXTERNAL
                //B (4)
#ifdef TREENODEALLBODIES
#ifdef SINCOS
                if (computeBodyProperties_balls_omp_cc_sincos(cmd, gd,
                             (bodyptr)p, nbody[cat1], &histsincos) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
#else
                if (computeBodyProperties_balls_omp_cc(cmd, gd, (bodyptr)p,
                                      nbody[cat1], &histcc) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
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
/*                if (cballs_opt_smooth_pivot(cmd))
                verb_print(cmd->verbose,
                    " - Number of neighbours: %ld %d %d\n",
                           i, NbRmin(p), NbRminOverlap(p)); */
            }

//B 2024-02-01
// Correct it because p and q are in differents node structures!!!... cat1!=cat2
//            if (cballs_opt_smooth_pivot(cmd)) {
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
      if (computeBodyProperties_balls_omp(cmd, gd, (bodyptr)p,
                                          nbody[cat1], &hist) == FAILURE) {
#pragma omp atomic write
          worker_failed = TRUE;
      }
#endif
#endif

//#ifdef TPCF
            if (gd->computeTPCF) {
#ifdef PIVOTEXTERNAL
                //B (5)
#ifdef TREENODEALLBODIES
#ifdef SINCOS
                if (computeBodyProperties_balls_omp_cc_sincos(cmd, gd,
                             (bodyptr)p, nbody[cat1], &histsincos) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
#else
                if (computeBodyProperties_balls_omp_cc(cmd, gd, (bodyptr)p,
                                      nbody[cat1], &histcc) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
#endif // ! SINCOS
#else // ! TREENODEALLBODIES
                if (computeBodyProperties_balls_omp_cc_sincos(cmd, gd,
                             (bodyptr)p, nbody[cat1], &histsincos) == FAILURE) {
#pragma omp atomic write
                    worker_failed = TRUE;
                }
#endif // ! TREENODEALLBODIES
                //E
#endif // ! PIVOTEXTERNAL
//#endif
            }

#pragma omp critical
    {
      if (hist_ready && !allocation_failed && !worker_failed) {
#ifdef TREENODEBALLS4
        int n;
        for (n = 1; n <= cmd->sizeHistN; n++) {
            gd->histNN[n] += histb.histNthread[n];
            gd->histNNSub[n] += histb.histNNSubthread[n];
            gd->histNNSubXi2pcf[n] += histb.histNNSubXi2pcfthread[n];
            gd->histXi2pcf[n] += histb.histXi2pcfthread[n];
            if (!cballs_opt_no_two_balls(cmd)) {
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
        if (gd->computeTPCF) {
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
                //        if (!cballs_opt_no_two_balls(cmd)) {
                //            ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histccb.histZetaMthread[m],
                //                     cmd->sizeHistN);
                //        }
            }
#else // ! SINCOS
            int m;
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                ADDM_ext(gd->histZetaM[m],gd->histZetaM[m],histb.histZetaMthread[m],
                         cmd->sizeHistN);
                if (!cballs_opt_no_two_balls(cmd)) {
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
      }
    } // end pragma omp critical

//B BALLS4
//    search_free_sincos_omp(&histsincos);
//    search_free_balls_omp_cc(&histcc);
//    search_free_balls_omp(&hist);
    int cleanup_failed = FALSE;
#ifdef TREENODEBALLS4
    if (search_free_sincos_omp_balls6(cmd, gd, &histbsincos) == FAILURE)
        cleanup_failed = TRUE;
    if (search_free_omp_balls6_cc(cmd, gd, &histccb) == FAILURE)
        cleanup_failed = TRUE;
    if (search_free_omp_balls6(cmd, gd, &histb) == FAILURE)
        cleanup_failed = TRUE;
#else
    if (search_free_sincos_omp(cmd, gd, &histsincos) == FAILURE)
        cleanup_failed = TRUE;
    if (search_free_balls_omp_cc(cmd, gd, &histcc) == FAILURE)
        cleanup_failed = TRUE;
    if (search_free_balls_omp(cmd, gd, &hist) == FAILURE)
        cleanup_failed = TRUE;
#endif
    if (cleanup_failed) {
#pragma omp atomic write
        worker_failed = TRUE;
    }
  } // end pragma omp parallel

    if (allocation_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "searchcalc_balls_omp_0357: OpenMP histogram allocation failed");
        return FAILURE;
    }
    if (worker_failed) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "searchcalc_balls_omp_0357: OpenMP worker operation failed");
        return FAILURE;
    }

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
      if (cballs_opt_smooth_pivot(cmd)) {
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
          if (gd->computeTPCF) {
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
        if (cballs_opt_smooth_pivot(cmd)) {
            gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcftotal[nn],1.0);
        } else {
            gd->histXi2pcf[nn] /= MAX(gd->histNNSubXi2pcf[nn],1.0);
        }
//E
    }

#ifdef TREENODEALLBODIES
    if (cballs_opt_compute_histn(cmd)) {
        if (cballs_opt_smooth_pivot(cmd)) {
            if (search_compute_HistN_balls(cmd, gd,
                                           nbody[cat1]-ipfalse) == FAILURE)
                return FAILURE;
        } else {
            if (search_compute_HistN_balls(cmd, gd,
                                           nbody[cat1]) == FAILURE)
                return FAILURE;
        }
    }
#else
      if (cballs_opt_compute_histn(cmd)) {
              if (search_compute_HistN_balls(cmd, gd,
                                             nbody[cat1]) == FAILURE)
                  return FAILURE;
      }
#endif

    verb_print(cmd->verbose, "balls: nsmoothcount = %ld\n",gd->nsmoothcount);
    verb_print(cmd->verbose, "balls: imiss = %d\n",imiss);

#ifdef TREENODEALLBODIES
      if (cballs_opt_smooth_pivot(cmd)) {
          verb_print(cmd->verbose, "balls: p falses found = %ld\n",ipfalse);
          //B kappa Avg Rmin
          verb_print(cmd->verbose,
                     "balls: count NbRmin found = %ld\n",
                     icountNbRmin);
          verb_print(cmd->verbose,
                     "balls: count overlap found = %ld\n",
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
          verb_print(cmd->verbose, "balls: p falses found = %ld\n",
                     ifalsecount);
          verb_print(cmd->verbose, "balls: p true found = %ld\n",
                     itruecount);
          verb_print(cmd->verbose, "balls: total = %ld\n",
                     itruecount+ifalsecount);
          //E
      }
#endif

    gd->cpusearch = CPUTIME - cpustart;
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n",CPUTIME-cpustart);
    return SUCCESS;
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
 compute_vector dr;

    if (!reject_balls(cmd, gd, p, q, &dr1, dr)) {
        if (cballs_opt_smooth_pivot(cmd))
            if (dr1<=gd->rsmooth[0] && Type(q)==BODY) {
                if (Update(q) == TRUE) {
                    Update(q) = FALSE;
                    NbRmin(p) += 1;
                    KappaRmin(p) += Kappa(q);
                } else {
                    NbRminOverlap(p) += 1;
                }
            }
        if (Type(q) == CELL) {
            if (!cballs_opt_no_one_ball(cmd)) {
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
    compute_vector dr;

    int SWITCHVALUE=0;
    if (Type(p) == CELL && Type(q) == CELL) SWITCHVALUE = 1;
    if (Type(p) == CELL && Type(q) == BODY) SWITCHVALUE = 2;
    if (Type(p) == BODY && Type(q) == CELL) SWITCHVALUE = 3;
    if (Type(p) == BODY && Type(q) == BODY) SWITCHVALUE = 4;

    if (!reject_balls(cmd, gd, p, q, &dr1, dr)) {
        switch (SWITCHVALUE){
            case CELLCELL:
                 if (!cballs_opt_no_two_ball(cmd)) {
                     if ( ( (Size(p)<=gd->rminCell[0] && Size(q)<=gd->rminCell[0])
                           || (Nb(p)<=gd->nsmooth[0] && Nb(q)<=gd->nsmooth[0]) )
                         && (dr1 > gd->rminCell[1]) ) {
                         *nsmoothcountthread += 1;
                         sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr,
                                         hist, histcc, histccsincos,
                                         nbbcalcthread, nbccalcthread,
                                         ncccalcthread);
                     } else {
                     if (nodes_condition_balls(cmd, gd, p, q, &dr1, dr)) {
                        if ( (Size(p)<=gd->rminCell[0]
                              && Size(q)<=gd->rminCell[0])
                                 && (Nb(p)<=gd->nsmooth[0]
                                     && Nb(q)<=gd->nsmooth[0]) ) {
                             *nsmoothcountthread += 1;
                             sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr,
                                             hist, histcc,
                                             histccsincos, nbbcalcthread,
                                             nbccalcthread, ncccalcthread);
                        } else
                             sumnodes_cc_omp(cmd, gd, p, q, &dr1, dr,
                                             hist, histcc,
                                             histccsincos, nbbcalcthread,
                                             nbccalcthread, ncccalcthread);
                     } else // ! nodes_condition
                         for (h = More(p); h != Next(p); h = Next(h))
                         for (l = More(q); l != Next(q); l = Next(l))
                             walktree_balls_omp(cmd, gd, h, l,
                                                hist, histcc, histccsincos,
                                                nsmoothcountthread,
                                                nbbcalcthread, nbccalcthread,
                                                ncccalcthread);

                     }
                 } else // ! no-two-ball
                     for (h = More(p); h != Next(p); h = Next(h))
                     for (l = More(q); l != Next(q); l = Next(l))
                             walktree_balls_omp(cmd, gd, h, l,
                                                hist, histcc, histccsincos,
                                                nsmoothcountthread,
                                                nbbcalcthread, nbccalcthread,
                                                ncccalcthread);
                break;
            case CELLBODY:
                if (!cballs_opt_no_one_ball(cmd)) {
                    if ( ((Nb(p)<=gd->nsmooth[0]) || (Size(p)<=gd->rminCell[0]))
                        && (dr1 > gd->rminCell[1]) ) {
                        *nsmoothcountthread += 1;
                        sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr,
                                        hist, histcc, histccsincos,
                                        nbbcalcthread, nbccalcthread,
                                        ncccalcthread);
                    } else {

                    if (nodes_condition_balls(cmd, gd, p, q, &dr1, dr)) {
                        if (Size(p)<=gd->rminCell[0] || Nb(p)<=gd->nsmooth[0]) {
                            *nsmoothcountthread += 1;
                            sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr,
                                            hist, histcc, histccsincos,
                                            nbbcalcthread, nbccalcthread,
                                            ncccalcthread);
                        } else
                            sumnodes_cb_omp(cmd, gd, p, q, &dr1, dr,
                                            hist, histcc, histccsincos,
                                            nbbcalcthread, nbccalcthread,
                                            ncccalcthread);
                    } else // ! nodes_condition
                        for (l = More(p); l != Next(p); l = Next(l))
                            walktree_balls_omp(cmd, gd, l, q,
                                               hist, histcc, histccsincos,
                                               nsmoothcountthread,
                                               nbbcalcthread, nbccalcthread,
                                               ncccalcthread);
                    }
                } else // ! no-one-ball
                    for (l = More(p); l != Next(p); l = Next(l)) {
                            walktree_balls_omp(cmd, gd, l, q,
                                               hist, histcc, histccsincos,
                                               nsmoothcountthread,
                                               nbbcalcthread, nbccalcthread,
                                               ncccalcthread);
                    }
                break;
            case BODYCELL:
                if (!cballs_opt_no_one_ball(cmd)) {
                    if ( ((Nb(q)<=gd->nsmooth[0]) || (Size(q)<=gd->rminCell[0]))
                        && (dr1 > gd->rminCell[1]) ) {
                        *nsmoothcountthread += 1;
                        sumnodes_bc_omp(cmd, gd, p, q, &dr1, dr,
                                        hist, histcc, histccsincos,
                                        nbbcalcthread, nbccalcthread,
                                        ncccalcthread);
                    } else {

                    if (nodes_condition_balls(cmd, gd, p, q, &dr1, dr)) {
                        if (Size(q)<=gd->rminCell[0] || Nb(q)<=gd->nsmooth[0]) {
                            *nsmoothcountthread += 1;
                            sumnodes_bc_omp(cmd, gd, p, q, &dr1, dr,
                                            hist, histcc,
                                            histccsincos, nbbcalcthread,
                                            nbccalcthread, ncccalcthread);
                        } else
                            sumnodes_bc_omp(cmd, gd, p, q, &dr1, dr,
                                            hist, histcc,
                                            histccsincos, nbbcalcthread,
                                            nbccalcthread, ncccalcthread);
                    } else // ! nodes_condition
                        for (l = More(q); l != Next(q); l = Next(l))
                            walktree_balls_omp(cmd, gd, p, l,
                                               hist, histcc, histccsincos,
                                               nsmoothcountthread, nbbcalcthread,
                                               nbccalcthread, ncccalcthread);
                    }
                } else // ! no-one-ball
                    for (l = More(q); l != Next(q); l = Next(l)) {
                            walktree_balls_omp(cmd, gd, p, l,
                                               hist, histcc, histccsincos,
                                               nsmoothcountthread,
                                               nbbcalcthread, nbccalcthread,
                                               ncccalcthread);
                    }
                break;
            case BODYBODY: // Found BODYBODY
#ifdef DEBUG
                HIT(p) = TRUE;
                HIT(q) = TRUE;
#endif
                sumnodes_bb_omp(cmd, gd, p, q, &dr1, dr,
                                hist, histcc, histccsincos,
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

//#include "datastruc_defs_balls_omp.h"
#include "datastruc_defs_balls_omp_0357.h"


local void sumnodes_bb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, compute_vector dr,
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
        if (gd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            
            //B
#if NDIM == 3
            real s, sy;
            compute_vector pr0;
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
                               "sumenodes_bb_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (gd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode_bb_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbbcalcthread += 1;
    } // ! accept_body
}

local void sumnodes_bc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, compute_vector dr,
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
        if (gd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            //B
#if NDIM == 3
            real s, sy;
            compute_vector pr0;
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
                               "sumenodes_bc_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (gd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode_bc_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cb_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr q, nodeptr p, real *dr1, compute_vector dr,
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
        if (gd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            //B
#if NDIM == 3
            real s, sy;
            compute_vector pr0;
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
                               "sumenodes_cb_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (gd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode_cb_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *nbccalcthread += 1;
    } // ! accept_body
}

local void sumnodes_cc_omp(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, real *dr1, compute_vector dr,
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
        if (gd->computeTPCF) {
            real cosphi;
#ifdef SINCOS
            real sinphi;
            histccsincos->histNNSubthread[n] = histccsincos->histNNSubthread[n] + 1.;
            //B
#if NDIM == 3
            //B
            real s, sy;
            compute_vector pr0;
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
                               "sumenodes_cc_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
#ifdef SINCOS
            CHEBYSHEVTUOMP
#else
            CHEBYSHEVOMPBALLSCC;
#endif
//#endif // ! TPCF
        }
#else // ! PTOPIVOTROTATION3
//#ifdef TPCF
        if (gd->computeTPCF) {
#if SINCOS
            cosphi = -dr[0]/(*dr1);
            sinphi = -dr[1]/(*dr1);
#else
            cosphi = -dr[1]/(*dr1);
#endif
            if (rabs(cosphi)>1.0)
                verb_log_print(cmd->verbose, gd->outlog,
                               "sumenode_cc_omp: Warning!... cossphi must be in (-1,1): %g\n",cosphi);
            CHEBYSHEVOMPBALLSCC;
//#endif // ! TPCF
        }
#endif // ! PTOPIVOTROTATION3
        *ncccalcthread += 1;
    } // ! accept_body
}


//B BALLS4 :: METODO BALLS4 DE BUSQUEDA
#ifdef TREENODEBALLS4

local int walktree_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
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
    compute_vector dr;
    int n;

    if (Update(p)) {
        np = nptr;
        actsafe = hist->actlen - NSUB;
 //B for ap
         for (ap = aptr; ap < nptr; ap++) {
             if (Type(*ap) == CELL) {
                 if (!reject_cell_balls(cmd, gd, p, *ap, &dr1, dr)) {
                     if ( ((Nb(*ap)<=gd->nsmooth[0]) || (Size(*ap)<=gd->rminCell[0]))
                         && (dr1 > gd->rminCell[1]) ) {
                             if (np - hist->active >= actsafe)
                                 return FAILURE;
                             *nsmoothcountthread += 1;
                             --bptr;
                             Mass(bptr) = Mass(*ap);
                             Kappa(bptr) = Kappa(*ap);
                             SETV(Pos(bptr), Pos(*ap));
                             Id(bptr) = Id(*ap);
                             Type(bptr) = Type(*ap);
                             Nb(bptr) = Nb(*ap);
                     } else { // ! bucket condition
                         if (nodes_condition_balls(cmd, gd, p, *ap, &dr1, dr)) {
                             if (!cballs_opt_no_two_balls(cmd)
                                 && Type(p) == CELL ) {
                                 sumcellcell_balls6_omp(cmd, gd, (cellptr)(*ap),
                                    (cellptr)*ap+1, p,
                                    nbbcalcthread, nbccalcthread, ncccalcthread,
                                    histcc, histsincos);
                             } else {
                                 if (np - hist->active >= actsafe)
                                     return FAILURE;
                                 if (!cballs_opt_no_one_ball(cmd)) {
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
        if (np != nptr) {
            if (walksub6_omp(cmd, gd, nptr, np, cptr, bptr, p, psize, pmid,
                            nbbcalcthread, nbccalcthread, ncccalcthread,
                            hist, histcc, histsincos, ibfcountthread,
                            nsmoothcountthread, nbody) == FAILURE)
                return FAILURE;
        }
        else {
            if (Type(p) != BODY)
                return FAILURE;

            if (sum_balls6_omp(cmd, gd, cptr, bptr, (bodyptr) p,
                              nbbcalcthread, nbccalcthread, ncccalcthread,
                              hist, histsincos, nbody) == FAILURE)
                return FAILURE;
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

    return SUCCESS;
}

local int walksub6_omp(struct cmdline_data* cmd, struct  global_data* gd,
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
            if (walktree_balls6_omp(cmd, gd, nptr, np, cptr, bptr, q,
                    psize / 2, nmid, nbbcalcthread, nbccalcthread,
                    ncccalcthread, hist, histcc, histsincos, ibfcountthread,
                    nsmoothcountthread, nbody) == FAILURE)
                return FAILURE;
        }
    } else {
        for (k = 0; k < NDIM; k++)
            nmid[k] = pmid[k] + (Pos(p)[k] < pmid[k] ? - poff : poff);
        if (walktree_balls6_omp(cmd, gd, nptr, np, cptr, bptr, p,
                psize / 2, nmid, nbbcalcthread, nbccalcthread,
                ncccalcthread, hist, histcc, histsincos, ibfcountthread,
                nsmoothcountthread, nbody) == FAILURE)
            return FAILURE;
    }

    return SUCCESS;
}

local int sum_balls6_omp(struct cmdline_data* cmd, struct  global_data* gd,
                         cellptr cptr, cellptr bptr, bodyptr p0,
        INTEGER *nbbcalcthread, INTEGER *nbccalcthread, INTEGER *ncccalcthread,
        gdhistptr_omp_balls6 hist, gdhistptr_sincos_omp_balls6 histsincos,
                          INTEGER nbody)
{
    int n;
    INTEGER ip;
    gdhist_omp_balls6 hist1;
    gdhist_sincos_omp_balls6 hist1sincos;

    memset(&hist1, 0, sizeof(hist1));
    memset(&hist1sincos, 0, sizeof(hist1sincos));
    if (search_init_balls_0357_guarded(cmd, gd, &hist1, 0,
                                      BALLS0357_INIT_BALLS6_CC) == FAILURE ||
        search_init_balls_0357_guarded(cmd, gd, &hist1sincos, 0,
                               BALLS0357_INIT_SINCOS_BALLS6) == FAILURE) {
        search_free_omp_balls6_cc(cmd, gd, &hist1);
        search_free_sincos_omp_balls6(cmd, gd, &hist1sincos);
        return FAILURE;
    }
//B
    for (n = 1; n <= cmd->sizeHistN; n++) {
        hist1.histNNSubthread[n] = 0.0;
        hist1.histXi2pcfthreadsub[n] = 0.0;
    }
//#ifdef TPCF
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (!cballs_opt_no_one_ball(cmd))
        sumcell_balls6_omp(cmd, gd, hist->interact, cptr, (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);
    sumnode_balls6_omp(cmd, gd, bptr, hist->interact + hist->actlen,
                       (bodyptr) p0,
                nbbcalcthread, nbccalcthread, ncccalcthread, &hist1, &hist1sincos);

//B Section of type: computeBodyProperties_omp_balls6(p0, cmd->nbody, hist)
    real xi = 0.0, xi_2p = 0.0;

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
    if (gd->computeTPCF) {
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

    int cleanup_failed = FALSE;
    if (search_free_omp_balls6_cc(cmd, gd, &hist1) == FAILURE)
        cleanup_failed = TRUE;
    if (search_free_sincos_omp_balls6(cmd, gd, &hist1sincos) == FAILURE)
        cleanup_failed = TRUE;

    return cleanup_failed ? FAILURE : SUCCESS;
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
    compute_vector dr;
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
        if (accept_body(cmd, gd, p0, pb, &dr1, dr)) {
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
                    if (gd->computeTPCF) {
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
                        compute_vector pr0;
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
    compute_vector dr;
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
        if (accept_body(cmd, gd, p0, pb, &dr1, dr)) {
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
                    if (gd->computeTPCF) {
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
                        compute_vector pr0;
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
    compute_vector dr;
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
        if (accept_body(cmd, gd, (bodyptr)p0, pb, &dr1, dr)) {
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
                    if (gd->computeTPCF) {
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
                        compute_vector pr0;
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    real xi = 0.0, xi_2p = 0.0;

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
    if (gd->computeTPCF) {
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
    real xi = 0.0, xi_2p = 0.0;

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
    if (gd->computeTPCF) {
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
        real xi = 0.0, xi_2p = 0.0;

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
    if (gd->computeTPCF) {
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


#define BALLS0357_FREE_DVECTOR(pointer, low, high)                    \
    do {                                                              \
        if ((pointer) != NULL) {                                      \
            free_dvector((pointer), (low), (high));                   \
            (pointer) = NULL;                                         \
        }                                                             \
    } while (0)
#define BALLS0357_FREE_DMATRIX(pointer, rlow, rhigh, clow, chigh)     \
    do {                                                              \
        if ((pointer) != NULL) {                                      \
            free_dmatrix((pointer), (rlow), (rhigh), (clow), (chigh)); \
            (pointer) = NULL;                                         \
        }                                                             \
    } while (0)
#define BALLS0357_FREE_DMATRIX3D(pointer, rlow, rhigh, clow, chigh,  \
                                 dlow, dhigh)                         \
    do {                                                              \
        if ((pointer) != NULL) {                                      \
            free_dmatrix3D((pointer), (rlow), (rhigh),                \
                           (clow), (chigh), (dlow), (dhigh));          \
            (pointer) = NULL;                                         \
        }                                                             \
    } while (0)

global int search_free_omp_balls6(struct cmdline_data* cmd, struct  global_data* gd,
                                  gdhistptr_omp_balls6 hist)
{
        free(hist->interact);
        hist->interact = NULL;
        free(hist->active);
        hist->active = NULL;

//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthread, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithread, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN);
//#endif
    }
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthreadsub, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthread, 1, cmd->sizeHistN);
    // 2pcf
        //B kappa Avg Rmin
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadtotal,
                                  1, cmd->sizeHistN);
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadp,
                                  1, cmd->sizeHistN);
        //E
        BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthread,
                              1, cmd->sizeHistN);
    //
        BALLS0357_FREE_DVECTOR(hist->histNNSubthread, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histNthread, 1, cmd->sizeHistN);
//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmp, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVP, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->Chebs, 1, cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_sincos_omp_balls6(struct cmdline_data* cmd, 
                                         struct  global_data* gd,
                                         gdhistptr_sincos_omp_balls6 hist)
{
//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmpsincos, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmpsin, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmpcos, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVPsincos, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVPsin, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVPcos, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthreadsincos, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthreadsin, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthreadcos, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithreadsin, 1,
                              cmd->mChebyshev+1, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithreadcos, 1,
                              cmd->mChebyshev+1, 1, cmd->sizeHistN);
//#endif
    }
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthreadsub, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthread, 1, cmd->sizeHistN);
        //B kappa Avg Rmin
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadtotal,
                                  1, cmd->sizeHistN);
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadp,
                                  1, cmd->sizeHistN);
        //E
        BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthread,
                              1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histNNSubthread, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histNthread, 1, cmd->sizeHistN);
//    #ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DVECTOR(hist->ChebsU, 1, cmd->mChebyshev+1);
        BALLS0357_FREE_DVECTOR(hist->ChebsT, 1, cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_omp_balls6_cc(struct cmdline_data* cmd, 
                                     struct  global_data* gd,
                                     gdhistptr_omp_balls6 hist)
{
//    #ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthread, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithread, 1,
                              cmd->mChebyshev+1, 1, cmd->sizeHistN);
//#endif
    }
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthreadsub, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histXi2pcfthread, 1, cmd->sizeHistN);
        //B kappa Avg Rmin
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadtotal,
                                  1, cmd->sizeHistN);
            BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadp,
                                  1, cmd->sizeHistN);
        //E
        BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthread,
                              1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histNNSubthread, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->histNthread, 1, cmd->sizeHistN);
//    #ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmp, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVP, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->Chebs, 1, cmd->mChebyshev+1);
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (gd->computeTPCF) {
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
    if (cballs_opt_smooth_pivot(cmd)) {
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
    if (gd->computeTPCF) {
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
    if (cballs_opt_smooth_pivot(cmd)) {
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
    if (gd->computeTPCF) {
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
     hist->interact = NULL;
     free(hist->active);
     hist->active = NULL;

//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthread, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithread, 1,
                              cmd->mChebyshev+1, 1, cmd->sizeHistN);
//#endif
    }
     BALLS0357_FREE_DVECTOR(hist->histXi2pcfthreadsub, 1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histXi2pcfthread, 1, cmd->sizeHistN);
     //B kappa Avg Rmin
         BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadtotal,
                               1, cmd->sizeHistN);
         BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadp,
                               1, cmd->sizeHistN);
     //E
     BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthread,
                           1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histNNSubthread, 1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histNthread, 1, cmd->sizeHistN);
//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmp, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVP, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->Chebs, 1, cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

global int search_free_balls_omp_cc(struct cmdline_data* cmd, 
                                    struct  global_data* gd,
                                    gdhistptr_omp_balls hist)
{
//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX3D(hist->histZetaMthread, 1,
            cmd->mChebyshev+1, 1, cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->histXithread, 1,
                              cmd->mChebyshev+1, 1, cmd->sizeHistN);
//#endif
    }
     BALLS0357_FREE_DVECTOR(hist->histXi2pcfthreadsub, 1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histXi2pcfthread, 1, cmd->sizeHistN);
     //B kappa Avg Rmin
         BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadtotal,
                               1, cmd->sizeHistN);
         BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthreadp,
                               1, cmd->sizeHistN);
     //E
     BALLS0357_FREE_DVECTOR(hist->histNNSubXi2pcfthread,
                           1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histNNSubthread, 1, cmd->sizeHistN);
     BALLS0357_FREE_DVECTOR(hist->histNthread, 1, cmd->sizeHistN);
//#ifdef TPCF
    if (gd->computeTPCF) {
        BALLS0357_FREE_DMATRIX(hist->histZetaMtmp, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DMATRIX(hist->xiOUTVP, 1,
                              cmd->sizeHistN, 1, cmd->sizeHistN);
        BALLS0357_FREE_DVECTOR(hist->Chebs, 1, cmd->mChebyshev+1);
//#endif
    }

    return SUCCESS;
}

#undef BALLS0357_FREE_DMATRIX3D
#undef BALLS0357_FREE_DMATRIX
#undef BALLS0357_FREE_DVECTOR

//B 2023.11.22
local bool nodes_condition_balls(struct cmdline_data* cmd, struct  global_data* gd,
                                  nodeptr p, nodeptr q, real *dr1, compute_vector dr)
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
             if (cballs_opt_behavior_tree_omp(cmd)) {
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
                            nodeptr p, nodeptr q, real *dr1, compute_vector dr)
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
             if (cballs_opt_behavior_tree_omp(cmd)) {
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
     compute_vector dr;

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
                          nodeptr p, nodeptr q, int *n, real *dr1, compute_vector dr)
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

global bool nodes_set_bin5(struct cmdline_data* cmd, struct  global_data* gd,
                           nodeptr p, nodeptr q, int *n, real *dr1,
                           compute_vector dr)
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

global bool nodes_set_bin_bc(struct cmdline_data* cmd, struct  global_data* gd,
                             nodeptr p, nodeptr q, int *n, real *dr1,
                             compute_vector dr)
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
        else
            cBALLS_FAIL(cmd, "search_compute_Xi_balls: wrong NDIM (%d)", NDIM);
    #endif

    //#ifdef TREENODE
    //    normFac *= 0.5;
    //#else
    //    if (cballs_opt_compute_j_no_eq_i(cmd))
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
            else
                cBALLS_FAIL(cmd,
                            "search_compute_Xi_balls: wrong NDIM (%d)", NDIM);
    #endif

//    #ifndef LOGHIST
    if (cmd->useLogHist) {
        gd->histNN[1]-=nbody;
    }
//    #endif
        real rho_av=(real)nbody/Vol;

        for (n = 1; n <= cmd->sizeHistN; n++) {
            if (gd->histNN[n] != 0) {
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
              gd->histCF[n] = rho_r/rho_av;
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

    if (cballs_opt_and_cf(cmd) &&
        search_compute_Xi_balls(cmd, gd, nbody) == FAILURE)
        return FAILURE;

    return SUCCESS;
}

    
//E Routines like in treeutils

#endif // ! BALLS
