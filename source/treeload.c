/* ==============================================================================
 MODULE: treeload.c			[cTreeBalls]
 Written by: M.A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: 3-point correlation function computation
 Language: C
 Use: maketree(cmd, gd, btab, nbody, ifile)
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "tree_contracts.h"

local int treeInfo(struct  cmdline_data* cmd, struct  global_data* gd,
                   int, INTEGER);
local int newtree(struct  cmdline_data* cmd, struct  global_data* gd, int);
local cellptr makecell(struct  cmdline_data* cmd, struct  global_data* gd, int);
local int loadbody(struct  cmdline_data*, struct  global_data*, bodyptr, int);
local int subindex(bodyptr, cellptr);
local int hackcellprop(struct  cmdline_data* cmd, struct  global_data* gd,
                        cellptr, real, int, int);
local int setradius(struct  cmdline_data* cmd, struct  global_data* gd,
                    cellptr, compute_vector, real, int);
local void threadtree(struct  cmdline_data* cmd, struct  global_data* gd, 
                      nodeptr, nodeptr, int);
local void free_unthreaded_tree(nodeptr, long int *);
local long int free_tree_file(struct global_data *, int);
local int statBodies(struct  cmdline_data* cmd, struct  global_data* gd,
                     bodyptr btab, int nbody);

#define MAXLEVEL  32

local int cellhist[MAXLEVEL];
local int subnhist[MAXLEVEL];
//B Smooth(ing) section
local INTEGER NTOT[1];                              // Two sets of cells to smooth
local INTEGER ip;                                   //  bodies
#define NbMax 33
local INTEGER cellhistNb[NbMax];
local int cellRadius[NbMax];
local real deltaRadius;
//E

local void walktree_selected(nodeptr, real);        // To see the bodies belonging
                                                    //  to a cell

local INTEGER inode;
#ifdef DEBUG
local void walktree_hit(struct  cmdline_data* cmd, struct  global_data* gd,
                        nodeptr, real);
#endif

#ifdef DEBUGTREE
local char treeinfofilePath[MAXLENGTHOFFILES];
local FILE *outtreeinfo;
local void walkTree_printInfo(struct  cmdline_data* cmd, struct  global_data* gd,
                              nodeptr q, real qsize, int lev);
#endif

local int scanLevel(struct  cmdline_data* cmd, struct  global_data* gd, int);

#ifdef CBALLS_NEEDS_BALLS4_SCAN
local int scanLevelB4(struct  cmdline_data* cmd,
                      struct  global_data* gd, int ifile);
local int walktree_scan_lev_balls4(struct cmdline_data* cmd,
                                  struct global_data* gd,
                                  nodeptr q, int lev, int ifile,
                                  int scanLevel);
local INTEGER inodelevB4;
local INTEGER ibodyleftoutB4;
#endif

#ifndef MACONLY
//B celltable
local INTEGER ncell;
//E
#endif

#ifdef PRUNING
local int pruningCells(struct  cmdline_data* cmd,
                        struct  global_data* gd,
                        int ifile, cellptr p, int lev);
#endif

local INTEGER isel, inosel;
//B socket:
#ifdef ADDONS
#include "treeload_include_00.h"
#endif
//E


/*
 MakeTree routine to create octtree structure:

 To be called using: search=octree-omp-sincos

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `btab`: Input: catalog of bodies
    * `nbody`: Input: number of points in table array
    * `ifile`: Input: tag index that identify the catalog
    * Global tructures used: gd, cmd
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int MakeTree(struct  cmdline_data* cmd,
                    struct  global_data* gd,
                    bodyptr btab, INTEGER nbody, int ifile)
{
    string routineName = "MakeTree";
    double cpustart;
    double cpustartMiddle;
    bodyptr p;
    int i;
    bool preserve_catalog_frame = cballs_observer_frame(cmd);

    cpustart = CPUTIME;
    gd->bytes_tot_cells = 0;

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n%s: making the tree...\n", routineName);

#ifdef DEBUG
//B To debug cells:
    //B
    if (format_checked(cellsfilePath, sizeof(cellsfilePath),
                       "cellsfilePath", "%s/cells%s.txt",gd->tmpDir,cmd->suffixOutFiles) != 0)
        return FAILURE;
    //E

    if (!(outcells = fopen(cellsfilePath, "w"))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: error opening file '%s'", routineName, cellsfilePath);
        return FAILURE;
    }
//E
#endif

    if (statBodies(cmd, gd, btab, nbody)==FAILURE)
        return FAILURE;
    
    if (newtree(cmd, gd, ifile) == FAILURE) return FAILURE;
    roottable[ifile] = makecell(cmd, gd, ifile);
    if (roottable[ifile] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: failed to allocate root cell for catalog %d",
                 routineName, ifile);
        return FAILURE;
    }
    gd->tree_allocated = TRUE;
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    cpustartMiddle = CPUTIME;
    if (FindRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]) == FAILURE) return FAILURE;
#ifdef OCTREESHEAROMP
    preserve_catalog_frame |= gd->searchMethod_int == OCTREESHEARMETHOD;
#endif
#ifdef OCTREE2BALLSOMP
    preserve_catalog_frame |= gd->searchMethod_int == OCTREE2BALLSMETHOD;
#endif
    if (!preserve_catalog_frame
        && centerBodies(btab, nbody, ifile, roottable[ifile]) == FAILURE)
        return FAILURE;
    if (FindRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]) == FAILURE) return FAILURE;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n%s: centerBodies-FindRootCenter CPU time: %lf %s\n",
                    routineName,
                    CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);
//E

    cpustartMiddle = CPUTIME;
    CLRV(Pos(roottable[ifile]));
    if (expandbox(cmd, gd, btab, nbody, ifile, roottable[ifile]) == FAILURE) return FAILURE;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n%s: expandbox CPU time: %lf %s\n",
                    routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

    cpustartMiddle = CPUTIME;

    treeInfo(cmd, gd, ifile, nbody);

    DO_BODY(p, btab, btab+nbody) {
#ifdef BODY3ON
        Nbb(p) = 1;                                 // Check consistency with
                                                    //  smoothing... Correction
#endif
        if (loadbody(cmd, gd, p, ifile) == FAILURE) return FAILURE;
        Nb(p) = 1;
        Radius(p) = 0.0;
    }
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n%s: loadbody CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

    gd->tdepthTable[ifile] = 0;

    for (i = 0; i < MAXLEVEL; i++)
        cellhist[i] = subnhist[i] = 0;
    for (i = 0; i < NbMax; i++)
        cellhistNb[i] = cellRadius[i] = 0;
    deltaRadius = gd->rSizeTable[ifile]/NbMax;
    verb_log_print(cmd->verbose_log,gd->outlog,
                   "\ndeltaRadius = %g\n",deltaRadius);

    DO_BODY(p,btab,btab+nbody)                      // See bodies belonging to a
        Selected(p) = FALSE;                        //  cell

    NTOT[0] = 0;                                    // Smooth(ing) section
    cpustartMiddle = CPUTIME;
    if(hackcellprop(cmd, gd, roottable[ifile], gd->rSizeTable[ifile], 0, ifile) == FAILURE)
        return FAILURE;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: hackcellprop CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

//B To bin cell's radius... Check!!! (bins logscale)
    real rBin;
    verb_log_print(cmd->verbose_log,gd->outlog,
                   "\n%s: radius histogram:\n", routineName);
    for (i = 0; i < NbMax; i++) {
        rBin = ((int)i)*gd->rSizeTable[ifile]/deltaRadius;
        verb_log_print(cmd->verbose_log,gd->outlog,"%g %d\n", rBin, cellRadius[i]);
    }
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");
//E

//B socket:
#ifdef ADDONS
#include "treeload_include_01.h"
#endif
//E

#ifdef OCTREESMOOTHING
    tree_is_threaded[ifile] = TRUE;
#endif

#ifndef OCTREESMOOTHING
    ip = 0;

    cpustartMiddle = CPUTIME;
    
#ifndef MACONLY
    //B celltable
    ncell=0;
    celltable[ifile] =
        (cellptr *) allocate(gd->ncellTable[ifile] * sizeof(cellptr));
    //E
#endif

    threadtree(cmd, gd, (nodeptr) roottable[ifile], NULL, ifile);
    tree_is_threaded[ifile] = TRUE;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: threadtree CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "threadtree: number ip of selected cells = %ld\n",ip);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%d real node (range of nodes to search: >nc1 && <nc2)\n",inode);
    gd->rnnode = inode;
#endif

    cpustartMiddle = CPUTIME;
    walktree_selected((nodeptr) roottable[ifile],   // Smooth(ing) section
                      gd->rSizeTable[ifile]);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: walktree_selected CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

    isel=0, inosel=0;
    DO_BODY(p,btab,btab+nbody) {                    // See bodies belonging
        if (Selected(p))                            //  to a cell
            isel++;
        else
            inosel++;
    }
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\nSelected vs NotSelected and total: %ld %ld %ld\n\n",
                isel, inosel, isel + inosel);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "tdepth = %d\n\n",gd->tdepthTable[ifile]);

    cpustartMiddle = CPUTIME;
    // Scan pivot and root trees
    if (scanLevel(cmd, gd, ifile) == FAILURE)
        return FAILURE;
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: scanLevel CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

#ifdef CBALLS_NEEDS_BALLS4_SCAN
    if (cballs_method_needs_balls4_scan(gd->searchMethod_int)) {
        cpustartMiddle = CPUTIME;
        if (scanLevelB4(cmd, gd, ifile) == FAILURE)
            return FAILURE;
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: scanLevelB4 CPU time: %lf %s\n",
                        routineName, CPUTIME - cpustartMiddle,
                        PRNUNITOFTIMEUSED);
    }
#endif

//B socket:
#ifdef ADDONS
#include "treeload_include_01b.h"
#endif
//E

//B Histogram useful to smooth cells
    verb_log_print(cmd->verbose_log,gd->outlog,
        "\n%s: Nb histogram:\n", routineName);
    for (i = 0; i < NbMax; i++)
        verb_log_print(cmd->verbose_log,gd->outlog,
            "%d %ld\n", i, cellhistNb[i]);
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");
//E
    gd->bytes_tot += gd->bytes_tot_cells;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\nAllocated %g MByte for (%d) cells storage.\n",
                        gd->bytes_tot_cells*INMB, gd->ncellTable[ifile]);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: root number of bodies = %ld\n",
                        routineName, Nb(roottable[ifile]));

#ifdef DEBUG
    fclose(outcells);                               // Close file to debug cells
#endif

//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
    DO_BODY(p,btab,btab+nbody)
        HIT(p) = FALSE;
    walktree_hit(cmd, gd, (nodeptr) roottable[ifile], Size(roottable[ifile]));
#endif
//E

#ifdef PRUNING
    cpustartMiddle = CPUTIME;
    for (i = 0; i < NbMax; i++)
        cellhistNb[i] = cellRadius[i] = 0;
    if (pruningCells(cmd, gd, ifile, roottable[ifile], 0) == FAILURE)
        return FAILURE;
    verb_log_print(cmd->verbose_log,gd->outlog,
                   "\n%s: radius histogram (after pruning):\n", routineName);
    for (i = 0; i < NbMax; i++) {
        rBin = ((int)i)*gd->rSizeTable[ifile]/deltaRadius;
        verb_log_print(cmd->verbose_log,gd->outlog,"%g %d\n", rBin, cellRadius[i]);
    }
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");
    verb_log_print(cmd->verbose_log,gd->outlog,
        "\n%s: Nb histogram (after pruning):\n", routineName);
    for (i = 0; i < NbMax; i++)
        verb_log_print(cmd->verbose_log,gd->outlog,
            "%d %ld\n", i, cellhistNb[i]);
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: pruningCells CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);
#endif

#ifdef DEBUGTREE
    //B new
    if (format_checked(treeinfofilePath, sizeof(treeinfofilePath),
                       "treeinfofilePath", "%s/treeinfo%s.txt",
                       gd->tmpDir,cmd->suffixOutFiles) != 0)
        return FAILURE;
    //E

    if (!(outtreeinfo = fopen(treeinfofilePath, "w"))) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: error opening file '%s'", routineName, treeinfofilePath);
        return FAILURE;
    }

    
    fprintf(outtreeinfo, "lev, Nb, Radius, Kappa:\n");
    walkTree_printInfo(cmd, gd, (nodeptr) roottable[ifile],
                       Size(roottable[ifile]), 0);
    fclose(outtreeinfo);
#endif

    gd->cputree = CPUTIME - cpustart;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\tdone with the tree.\nmaking tree CPU time : %lf %s\n\n",
                gd->cputree, PRNUNITOFTIMEUSED);

    gd->tree_allocated = TRUE;

    return SUCCESS;
}

//B socket:
#ifdef ADDONS
#include "treeload_include_02.h"
#endif
//E

//B BALLS :: SCANLEV
local int scanLevel(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    int i;

//B socket:
#ifdef ADDONS
#include "treeload_include_03.h"
#endif
//E

    gd->Rcell[0] = gd->rSizeTable[ifile];
    for (i = 1; i <= gd->tdepthTable[ifile]; i++)
        gd->Rcell[i] = gd->Rcell[i-1]/2;

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\nMaximum and minimum cell size: %e %e\n",
                gd->Rcell[0],gd->Rcell[gd->tdepthTable[ifile]-1]);

//B socket:
#ifdef ADDONS
#include "treeload_include_04.h"
#endif
//E

#ifdef SMOOTHPIVOT
        if (cmd->useLogHist) {
            verb_print(cmd->verbose,
                       "deltaRV min and max: %lg %lg\n",
                       gd->deltaRmin, gd->deltaRmax);
        } else {
            verb_print(cmd->verbose,
                       "deltaR=%lf normal scale):\n",gd->deltaR);
        }
#endif

//B Some useful conversion info:
// rad = (pi/(180*60)) arc minutes
// (pi/(180*60)) = 0.000290888208666
// ((180*60)/pi) = 3437.74677
//
// 1 arcmin = 0.000290888208666 rad
// 1 rad = 3437.74677 arcmin
// 1 degree = 60 arcmin = 60*0.000290888208666 rad = 0.004072434921324 rad
//
// rsmooth must be given in radians.
// plots will be given in arcmin.
//E
// For unit sphere rsmooth are given in arcmin. Here we transform to radians
#ifdef SMOOTHPIVOT
#define ARCMINTORAD   0.000290888208666
        if (!cballs_opt_smooth_pivot(cmd)) {
            gd->rsmooth[0] = 0.0;
        } else if (strnull(cmd->rsmooth)) {
            if (cballs_opt_fix_rsmooth(cmd)) {
                //B Leave it as a refereence: green line
//            gd->rsmooth[0] = 0.00416666665;       // (0.25 arcmin)/60
                //E
                gd->rsmooth[0] = 0.01*cmd->rangeN;  // 1% of rangeN
                                                    // 200*0.01 -> 0.00058178 rad
                verb_print(cmd->verbose,
                           "\n\trsmooth is set to: %g\n", gd->rsmooth[0]);
                if (gd->rsmooth[0]>0.05*cmd->rangeN) {
                    verb_print(cmd->verbose,
            "Warning! rsmooth is greatear than 0.05*rangeN (%g %g)... fixing\n",
                               gd->rsmooth[0],0.05*cmd->rangeN);
                }
            } else {
                if (cballs_opt_default_rsmooth(cmd)) {
                    gd->rsmooth[0] = gd->Rcell[gd->irsmooth-1];
                } else {
                    gd->rsmooth[0] = gd->Rcell[gd->tdepthTable[ifile]-1];
                    verb_print(cmd->verbose,
                               "\tfixing rsmooth to: %g\n", gd->rsmooth[0]);
                    if (gd->rsmooth[0]>0.05*cmd->rangeN) {
                        verb_print(cmd->verbose,
            "Warning! rsmooth is greatear than 0.05*rangeN (%g %g)... fixing\n",
                                   gd->rsmooth[0],0.05*cmd->rangeN);
                        gd->rsmooth[0] = 0.75*gd->Rcell[gd->tdepthTable[ifile]-1];
                        verb_print(cmd->verbose,
                                   "\tfixing rsmooth to: %g\n", gd->rsmooth[0]);
                    }
                }
            }
        } else {
            // For Takahasi Nside 4096, rsmooth = 3 arcmin is still a good value
            gd->rsmooth[0] *= ARCMINTORAD;
            if (gd->rsmooth[0]>0.05*cmd->rangeN) {
                verb_print(cmd->verbose,
            "Warning! rsmooth is greatear than 0.05*rangeN (%g %g)... fixing\n",
                           gd->rsmooth[0],0.05*cmd->rangeN);
                gd->rsmooth[0] = 0.75*gd->Rcell[gd->tdepthTable[ifile]-1];
                verb_print(cmd->verbose,
                           "\tfixing rsmooth to: %g\n", gd->rsmooth[0]);
            } else {
                verb_print(cmd->verbose,
                           "\trsmooth is set to: %g\n", gd->rsmooth[0]);
            }
        }
//B added 2025-12-30
        gd->rsmooth[0] *= THETA;
//E
#undef ARCMINTORAD
#else
        gd->rsmooth[0] = 0;                         // setting a safe value
#endif

#ifdef SMOOTHPIVOT
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "rsmooth and rminHist %% of rangeN: %lg %lg\n",
                100*gd->rsmooth[0]/cmd->rangeN, 100*cmd->rminHist/cmd->rangeN);
#endif

    return SUCCESS;
}
//E BALLS :: SCANLEV

global int FindRootCenter(struct  cmdline_data* cmd,
                          struct  global_data* gd,
                          bodyptr btab, int nbody, int ifile, cellptr root)
{
    string routineName = "FindRootCenter";
    real len;
    int k;
    compute_vector xmin, xmax;

    DO_COORD(k)
        xmin[k] = xmax[k] = Pos(btab)[k];

    bodyptr p;
    int kk;
    DO_BODY(p, btab, btab+nbody) {
        DO_COORD(kk) {
            if (Pos(p)[kk] > xmax[kk])
                xmax[kk] = Pos(p)[kk];
            if (Pos(p)[kk] < xmin[kk])
                xmin[kk] = Pos(p)[kk];
        }
    }

    DO_COORD(k) {
        Pos(root)[k] = (xmax[k]+xmin[k])/2;
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: Pos(root) = %lf\n", routineName, Pos(root)[k]);
    }

    for(k=0, len=xmax[0]-xmin[0]; k<NDIM; k++)
        if((xmax[k]-xmin[k])>len)
            len=xmax[k]-xmin[k];

    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "%s: len = %lf\n", routineName, len);

    return SUCCESS;
}

global int centerBodies(bodyptr btab, int nbody, int ifile, cellptr root)
{
    bodyptr p;
    int k;

    DO_BODY(p, btab, btab+nbody)
        DO_COORD(k)
            Pos(p)[k] = Pos(p)[k] - Pos(root)[k];

    return SUCCESS;
}

local int statBodies(struct  cmdline_data* cmd, struct  global_data* gd,
                     bodyptr btab, int nbody)
{
    string routinename = "statBodies";
    bodyptr p;
    int k;
    real kavg = 0.;
    real kstd = 0.;

    DO_BODY(p, btab, btab+nbody) {
        kavg += Kappa(p);
    }
    
    kavg = kavg /((real)nbody);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\t%s: average of kappa (%ld particles) = %le\n",
                            routinename, nbody, kavg);

    DO_BODY(p, btab, btab+nbody) {
            kstd += (Kappa(p) - kavg)*(Kappa(p) - kavg);
    }
    // corrected sample standard deviation
    kstd = rsqrt( kstd / ((real)nbody - 1.0));
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\t%s: corrected sample standard deviation (%ld particles) = %le\n\n",
            routinename, nbody, kstd);

    if (cballs_opt_remove_mean(cmd)) {
        DO_BODY(p, btab, btab+nbody) {
            Kappa(p) = Kappa(p) - kavg;
        }
        
        kavg = 0.0;
        DO_BODY(p, btab, btab+nbody) {
            kavg += Kappa(p);
        }
        
        kavg = kavg /((real)nbody);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
        "\n\t%s: average of kappa without fluctuations (%ld particles) = %le\n",
        routinename, nbody, kavg);

    }

    return SUCCESS;
}

#define invlog102 3.321928094887362
local int treeInfo(struct  cmdline_data* cmd, struct  global_data* gd,
                   int ifile, INTEGER nbody)
{
    string routineName = "treeInfo";
    real averageLength;
    real tmp;
    real tdepthN;

    averageLength = rsqrt(4.0*PI/(real)nbody);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "%s: average length, root size: %g %g\n",
                           routineName, averageLength, gd->rSizeTable[ifile]);
    tmp = rsqrt((double) nbody)*gd->rSizeTable[ifile]
            / (2.0*SQRTPI);
    tdepthN = invlog102*rlog10(tmp);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "%s: tdepthN: %g\n",
                           routineName, tdepthN);

    return SUCCESS;
}
#undef invlog102

local int newtree(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    roottable[ifile] = NULL;
    gd->ncellTable[ifile] = 0;
    tree_is_threaded[ifile] = FALSE;

    return SUCCESS;
}

local void free_unthreaded_tree(nodeptr p, long int *cellcounter)
{
    int i;

    if (p == NULL || Type(p) != CELL)
        return;

    for (i = 0; i < NSUB; i++)
        free_unthreaded_tree(Subp(p)[i], cellcounter);

    free(p);
    ++(*cellcounter);
}

local long int free_tree_file(struct global_data *gd, int ifile)
{
    long int cellcounter = 0;

    if (roottable[ifile] == NULL)
        goto done;

    if (!tree_is_threaded[ifile]) {
        free_unthreaded_tree((nodeptr) roottable[ifile], &cellcounter);
#ifndef MACONLY
        free(celltable[ifile]);
        celltable[ifile] = NULL;
#endif
        goto done;
    }

#ifndef MACONLY
    if (celltable[ifile] != NULL) {
        for (INTEGER i = gd->ncellTable[ifile] - 1; i >= 0; i--) {
            free(celltable[ifile][i]);
            celltable[ifile][i] = NULL;
            ++cellcounter;
        }
        free(celltable[ifile]);
        celltable[ifile] = NULL;
    } else {
        free_unthreaded_tree((nodeptr) roottable[ifile], &cellcounter);
    }
#else
    nodeptr p = (nodeptr) roottable[ifile];
    nodeptr freecell = NULL;

    while (p != NULL) {
        if (Type(p) == CELL) {
            nodeptr more = More(p);
            Next(p) = freecell;
            freecell = p;
            p = more;
        } else {
            p = Next(p);
        }
    }

    while (freecell != NULL) {
        nodeptr next = Next(freecell);
        free(freecell);
        freecell = next;
        ++cellcounter;
    }
#endif

done:
    roottable[ifile] = NULL;
    gd->ncellTable[ifile] = 0;
    tree_is_threaded[ifile] = FALSE;
    return cellcounter;
}

// deallocate memory tree
global int freeTree(struct  cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "freeTree";
    int ifile;
    long int cellcounter=0;

    if (cballs_opt_read_mask(cmd)) {
        ifile=0;
#ifdef MACONLY
        INTEGER allocated_cells = gd->ncellTable[ifile];
#endif
        cellcounter = free_tree_file(gd, ifile);
#ifdef MACONLY
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                              "%s: allocated cells vs freed cells: %ld %ld\n",
                               routineName, (long int) allocated_cells,
                               cellcounter);
#endif
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
#ifdef MACONLY
            INTEGER allocated_cells = gd->ncellTable[ifile];
#endif
            cellcounter = free_tree_file(gd, ifile);
#ifdef MACONLY
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                  "%s: allocated cells vs freed cells: %ld %ld\n",
                                   routineName, (long int) allocated_cells,
                                   cellcounter);
#endif
        }
    }

    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        roottable[ifile] = NULL;
        gd->ncellTable[ifile] = 0;
        tree_is_threaded[ifile] = FALSE;
    }

    return SUCCESS;
}

local cellptr makecell(struct  cmdline_data* cmd,
                       struct  global_data* gd, int ifile)
{
    cellptr c;
    int i;
 
    c = (cellptr) allocate(sizeof(cell));
    Type(c) = CELL;
    Nb(c) = 0;                                      // To smooth cells
    Update(c) = FALSE;
    Selected(c) = FALSE;
    if (cballs_opt_read_mask(cmd))
        Mask(c) = FALSE;                            // check that FALSE is ok
    for (i = 0; i < NSUB; i++)
        Subp(c)[i] = NULL;
    gd->ncellTable[ifile] = gd->ncellTable[ifile] + 1;
    gd->bytes_tot_cells += sizeof(cell);
    return (c);
}

global int expandbox(struct  cmdline_data* cmd,
                    struct  global_data* gd,
                    bodyptr btab, int nbody, int ifile, cellptr root)
{
    string routineName = "expandbox (treeload)";
    real dmax, d;
    bodyptr p;
    int k;

    dmax = 0.0;
	DO_BODY(p, btab, btab+nbody)
		DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(root)[k]);
            if (d > dmax)
                dmax = d;                       
        }
    while (gd->rSizeTable[ifile] < 2 * dmax)
        gd->rSizeTable[ifile] = 2 * gd->rSizeTable[ifile];

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "%s: rSize = %lf\n",
                           routineName, gd->rSizeTable[ifile]);

    return SUCCESS;
}

//#define EPSILON 1.0E-7                            // Choose well.
#define EPSILON 1.0E-5                              // Choose well.
                                                    //  If not tdepth and
                                                    //  no. of cell will
                                                    //  grow badly...
#define EPSILONFLOAT 1.0E-0                         // this is for float pos

local int loadbody(struct  cmdline_data* cmd,
                   struct  global_data* gd, bodyptr p, int ifile)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;

    int sameposcount = 0;
startagain:
    q = roottable[ifile];
    qind = subindex(p, q);
    if (sameposcount == 0)
        Nb(q) += 1;                                 // Smooth
    qsize = gd->rSizeTable[ifile];
    while (Subp(q)[qind] != NULL) {
// BODY3
        if (Type(Subp(q)[qind]) == BODY || Type(Subp(q)[qind]) == BODY3) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0) { // look for these events in tmp/log file
                verb_log_print(cmd->verbose_log,gd->outlog,
                               "\nIds: %ld and %ld have the same position\n",
                               Id(p),Id(Subp(q)[qind]));
                DO_COORD(k)
                    verb_log_print(cmd->verbose_log,gd->outlog,
                                   "Pos[k]: %le %le\n",
                               Pos(p)[k],Pos(Subp(q)[qind])[k]);
                if (cballs_opt_no_check_equal_positions(cmd)) {
                    DO_COORD(k) {
//                        Pos(p)[k] += EPSILON*grandom(0.0, 0.01*qsize);
                        Pos(p)[k] += EPSILONLOADBODY*grandom(0.0, 1.0);
                    }
                    DO_COORD(k) {
                        verb_log_print(cmd->verbose_log,gd->outlog,
                        "CorrectedPos[k]: %g %g and correction: %le\n",
                        Pos(p)[k],Pos(Subp(q)[qind])[k],
                        EPSILONLOADBODY*grandom(0.0, 1.0));
//                           EPSILON*grandom(0.0, 0.01*qsize));
                    }
                    Update(p) = FALSE;
                    gd->sameposcount++;
                    sameposcount++;
                    goto startagain;
                } else {
                    snprintf(cmd->error_message, _ERRORMSGSIZE_,
                             "loadbody: two bodies have same position; consider option "
                             "'no-check-two-bodies-eq-pos'");
                    return FAILURE;
                }
            }
            c = makecell(cmd, gd, ifile);
            Nb(c) += 1;                             // Smooth
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
            Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
            Nb(c) += 1;                             // Smooth
            Subp(q)[qind] = (nodeptr) c;
        } else {
            Nb(Subp(q)[qind]) += 1;                 // Smooth
        }
        q = (cellptr) Subp(q)[qind];
        qind = subindex(p, q);
        qsize = qsize / 2;
    }
    Subp(q)[qind] = (nodeptr) p;

    return SUCCESS;
}
#undef EPSILON

local int subindex(bodyptr p, cellptr q)
{
    int ind, k;
 
    ind = 0;                                    
	DO_COORD(k)
        if (Pos(q)[k] <= Pos(p)[k])             
            ind += NSUB >> (k + 1);             
    return (ind);
}

local int hackcellprop(struct  cmdline_data* cmd, struct  global_data* gd,
                       cellptr p, real psize, int lev, int ifile)
{
    string routineName = "hackcellprop";
    compute_vector cmpos;
    int i, k;
    nodeptr q;
    short cell_mask = MASK_NODE_EMPTY;
    bool read_mask = cballs_opt_read_mask(cmd);
    static INTEGER mixed_cells;
    static INTEGER contaminated_children;

    if (lev == 0) {
        mixed_cells = 0;
        contaminated_children = 0;
    }

    gd->tdepthTable[ifile] = MAX(gd->tdepthTable[ifile], lev);
    if (lev < 0 || lev >= MAXLEVEL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: tree level out of range: lev=%d MAXLEVEL=%d\n",
                 routineName, lev, MAXLEVEL);
        return FAILURE;
    }
    cellhist[lev]++;
    Level(p) = lev;                                 // To set scanLevel
    Mass(p) = 0.0;
    Weight(p) = 0.0;
    Nb(p) = 0;
    Kappa(p) = 0.0;
    Selected(p) = FALSE;
    Update(p) = FALSE;
    CLRV(cmpos);

    for (i = 0; i < NSUB; i++) {
        if ((q = Subp(p)[i]) != NULL) {
            subnhist[lev]++;
            if (Type(q) == CELL) {
                if (hackcellprop(cmd, gd, (cellptr) q, psize/2, lev+1, ifile)
                    == FAILURE)
                    return FAILURE;
            }
            if (cballs_cell_accumulate_child((nodeptr)p, q, read_mask,
                                             &cell_mask, cmpos)
                && cmd->verbose_log >= 3)
                contaminated_children++;
        } // ! q no NULL
    } // ! loop i: 0 -> NSUB-1
    if (read_mask) {
        Mask(p) = cell_mask;
        if (Mask(p) == MASK_NODE_MIXED)
            mixed_cells++;
    }
//B Smooth(ing) section
    if (Nb(p)==cmd->nsmooth) {                      // Correct to <=
        NTOT[0]=NTOT[0]+1;
    }
//E
    if (cballs_opt_center_of_mass(cmd)) {
        if (Mass(p) > 0.0) {
            DIVVS(cmpos, cmpos, Mass(p));
        } else {
            SETV(cmpos, Pos(p));
        }
    } else
        SETV(cmpos, Pos(p));

//B this is the one working
//#define EPSILON 1.0E-16
//E
//B testing this one...
#define EPSILON 1.0E-12
//E
// Here there appears an error for big numbers of points such 201 millions...
// See line above and uncomment DIVVS(cmpos, cmpos, Mass(p)); line (not working!)
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 || Pos(p)[k] + psize/2 <= cmpos[k]) {
            if (psize/2 > 2.710505e-20 + EPSILONHACKCELL) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "%s: tree structure error: %d %le %le %le %le\n",
                               routineName, k, cmpos[k],
                               Pos(p)[k] - psize/2, Pos(p)[k] + psize/2, psize/2);
                return FAILURE;
            } else {
                if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log,gd->outlog,
                        "%s: tree structure warning! psize/2 to small: %le \n",
                        routineName, psize/2);
            }
        }
#undef EPSILON

#ifdef SINGLEP
    DO_COORD(k)
        cmpos[k] = (real)(cballs_storage_real)cmpos[k];
#endif
    if (setradius(cmd, gd, p, cmpos, psize, ifile) == FAILURE)
        return FAILURE;

    SETV(Pos(p), cmpos);
    if (Nb(p)>0) {
        cballs_cell_finalize_averages((nodeptr)p);
    } else {
#if defined(DEBUGTREE)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: Nb = 0: %ld\n", routineName, Nb(p));
        return FAILURE;
#endif
    }

    if (lev == 0 && cballs_opt_read_mask(cmd)
        && cmd->verbose_log >= 3)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "%s: mixed cells=%ld, masked children with "
                       "nonzero aggregates=%ld\n",
                       routineName, mixed_cells, contaminated_children);

    return SUCCESS;
}

// Parameter theta controls size of the cell.
// theta from 0 to 5:
//  0 always open cells (complexity N^2);
//  1 is the default value.
local int setradius(struct  cmdline_data* cmd, struct  global_data* gd,
                    cellptr p, compute_vector cmpos, real psize, int ifile)
{
    string routineName = "setradius";
    real bmax2, d, radius;
    int k;

    if (cmd->theta == 0.0)
        radius = 2 * gd->rSizeTable[ifile];
    else if (gd->sw94) {
        bmax2 = 0.0;
		DO_COORD(k) {
            d = cmpos[k] - Pos(p)[k] + psize/2; 
            bmax2 += rsqr(MAX(d, psize - d));   
        }
        radius = rsqrt(bmax2) / cmd->theta;
    } else if (gd->bh86)
        radius = psize / cmd->theta;
    else {
        radius = (psize/cmd->theta) * rsqrt((real)(NDIM))/2.0;
//        DISTV(d, cmpos(p), Pos(p));                 // find offset from center
//        Radius(p) = psize / cmd->theta + d;         // use size plus offset
    }

    Radius(p) = cballs_store_search_bound(radius);
    Size(p) = cballs_store_upper_bound(psize);

    int n;

    //B
    if (deltaRadius <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid deltaRadius=%g\n", routineName, deltaRadius);
        return FAILURE;
    }

    if (!isfinite((double)Radius(p)) || Radius(p) < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid cell radius=%g\n", routineName, Radius(p));
        return FAILURE;
    }
    if (cmd->theta == 0.0
        || Radius(p) >= (real)(NbMax - 1) * deltaRadius)
        n = NbMax - 1;                              // radius overflow bin
    else
        n = (int)(Radius(p) / deltaRadius);
    if (n < 0 || n >= NbMax) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: cellRadius index out of range: n=%d NbMax=%d Radius=%g deltaRadius=%g\n",
                 routineName, n, NbMax, Radius(p), deltaRadius);
        return FAILURE;
    }
    cellRadius[n]++;
    //E

    return SUCCESS;
}

local void threadtree(struct  cmdline_data* cmd,
                      struct  global_data* gd, nodeptr p, nodeptr n, int ifile)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    bodyptr q;

    Next(p) = n;
    if (Type(p) == CELL) {
#ifndef MACONLY
        //B celltable
        celltable[ifile][ncell] = (cellptr)p;
        ncell++;
        //E
#endif
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        ndesc = 0;
        for (i = 0; i < NSUB; i++)
            if (Subp(p)[i] != NULL)
                desc[ndesc++] = Subp(p)[i];
        More(p) = desc[0];
        desc[ndesc] = n;
        for (i = 0; i < ndesc; i++)
            threadtree(cmd, gd, desc[i], desc[i+1], ifile);
    } // ! p = CELL
}

// To see the bodies belonging to a cell
local void walktree_selected(nodeptr q, real qsize)
{
    nodeptr l;

    if (Selected(q)) {
        if (Type(q) == CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                Selected(l) = TRUE;
                walktree_selected(l,qsize/2);
            }
        } else {
            Selected(q) = TRUE;
        }
    } else {
        if (Type(q) == CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                walktree_selected(l,qsize/2);
            }
        } else {
            Selected(q) = FALSE;
        }
    }
}

//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
local void walktree_hit(struct  cmdline_data* cmd, struct  global_data* gd,
                        nodeptr q, real qsize)
{
    nodeptr l;

    if (Type(q) == CELL) {
        for (l = More(q); l != Next(q); l = Next(l)) {
            walktree_hit(cmd, gd, l,qsize/2);
        }
    } else {
        HIT(q) = TRUE;
    }
}
#endif
//E


#ifdef CBALLS_NEEDS_BALLS4_SCAN
local int scanLevelB4(struct  cmdline_data* cmd,
                      struct  global_data* gd, int ifile)
{
    string routine_name = "scanLevelB4";
    INTEGER covered_bodies = 0;
    bool read_mask = cballs_opt_read_mask(cmd);
    INTEGER expected_bodies = read_mask ? Nb(roottable[ifile])
                                        : gd->nbodyTable[ifile];

    gd->nnodescanlevTableB4[ifile] =
        gd->ncellTable[ifile]+gd->nbodyTable[ifile];
    nodetablescanlevB4[ifile] =
        (nodeptr *) allocate(gd->nnodescanlevTableB4[ifile] * sizeof(nodeptr));
    if (nodetablescanlevB4[ifile] == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: unable to allocate the B4 pivot table for catalog %d\n",
                 routine_name, ifile);
        gd->nnodescanlevTableB4[ifile] = 0;
        return FAILURE;
    }
    gd->bytes_tot += gd->nnodescanlevTableB4[ifile]*sizeof(nodeptr);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n%s: Allocated %g MByte for (%" INTEGER_FMT
                ") scan nodetab storage.\n",
                routine_name, INMB*gd->nnodescanlevTableB4[ifile]*sizeof(nodeptr),
                gd->nnodescanlevTableB4[ifile]);

    inodelevB4 = 0;
    ibodyleftoutB4 = 0;

//B scan tree up to the smallest cells
    if (walktree_scan_lev_balls4(cmd, gd,
                                (nodeptr)roottable[ifile], 0,
                                ifile, gd->tdepthTable[ifile]) == FAILURE)
        return FAILURE;
//E
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: Found %" INTEGER_FMT
                        " nodes to scan at upper most level %d.\n",
                        routine_name, inodelevB4, gd->tdepthTable[ifile]);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t%" INTEGER_FMT
                        " particles were included to scan at that level.\n",
                        ibodyleftoutB4);

    gd->nnodescanlevTableB4[ifile] = inodelevB4;


    nodeptr p;
    real rmax=0, rmin=cmd->rangeN;
    INTEGER numCells=0, numBodies=0;
    for (INTEGER i = 0; i < gd->nnodescanlevTableB4[ifile]; i++) {
        p = nodetablescanlevB4[ifile][i];
        if (p == NULL || Nb(p) <= 0
            || (read_mask && Mask(p) != MASK_NODE_VALID)) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: invalid B4 pivot at index %" INTEGER_FMT
                     " for catalog %d\n",
                     routine_name, i, ifile);
            return FAILURE;
        }
        covered_bodies += Nb(p);
        if (Type(p)==CELL) {
            if (Radius(p)>rmax) rmax=Radius(p);
            if (Radius(p)<rmin) rmin=Radius(p);
            numCells++;
        } else {
            numBodies++;
        }
    }
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: cell radius min and max: %lg %lg\n",
                           routine_name, rmin, rmax);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: number of cell and of bodies and total nodes: %"
                INTEGER_FMT " %" INTEGER_FMT " %" INTEGER_FMT "\n",
                routine_name, numCells, numBodies, numCells+numBodies);

    if (covered_bodies != expected_bodies) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: B4 pivots cover %" INTEGER_FMT
                 " bodies, expected %" INTEGER_FMT " in catalog %d\n",
                 routine_name, covered_bodies, expected_bodies, ifile);
        return FAILURE;
    }

    return SUCCESS;
}

local int walktree_scan_lev_balls4(struct cmdline_data* cmd,
                                  struct global_data* gd,
                                  nodeptr q, int lev, int ifile,
                                  int scanLevel)
{
    nodeptr child;
    INTEGER capacity = gd->ncellTable[ifile] + gd->nbodyTable[ifile];
    bool read_mask = cballs_opt_read_mask(cmd);

    if (q == NULL)
        return SUCCESS;

    if (read_mask && Mask(q) == MASK_NODE_MASKED)
        return SUCCESS;

    if ((Type(q) != CELL || Radius(q) < gd->deltaRmin*THETA
         || lev >= scanLevel)
        && (!read_mask || Mask(q) == MASK_NODE_VALID)) {
        if (inodelevB4 >= capacity) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "scanLevelB4: pivot table overflow for catalog %d\n",
                     ifile);
            return FAILURE;
        }
        nodetablescanlevB4[ifile][inodelevB4++] = q;
        if (Type(q) != CELL)
            ibodyleftoutB4++;
        return SUCCESS;
    }

    for (child = More(q); child != Next(q); child = Next(child))
        if (walktree_scan_lev_balls4(cmd, gd, child, lev+1, ifile,
                                    scanLevel) == FAILURE)
            return FAILURE;

    return SUCCESS;
}
#endif // ! CBALLS_NEEDS_BALLS4_SCAN


#ifdef PRUNING
// use in combination with options=center-of-mass
local void normal_walktree(struct  cmdline_data* cmd,
                           struct  global_data* gd,
                           cellptr p, real psize,  nodeptr q, real qsize, int lev)
{
    nodeptr l;
    real drpq, drpq2;
    compute_vector dr;
    real d;
    int k;
    bool flag=TRUE;                                 // inside cell

    if (Type(q) == CELL) {                          // is a cell, zoom in
        for (l = More(q); l != Next(q); l = Next(l))
            normal_walktree(cmd, gd, p, psize, l, qsize/2, lev+1);
    } else {                                        // found body, process it
        for (k = 0; k < NDIM; k++) {
            d = rabs(Pos(q)[k] - Pos(p)[k]);
            if (d > psize/2) flag=FALSE;
        }
        if (flag==TRUE) {
            if (cballs_opt_read_mask(cmd)
                && Mask(q) == MASK_NODE_MASKED)
                return;
            Nb(p) += 1;
            Weight(p) += Weight(q);
#ifndef NOWKAvg
            Kappa(p) += Weight(q)*Kappa(q);
#else
            Kappa(p) += Kappa(q);
#endif
            DOTPSUBV(drpq2, dr, Pos(q), Pos(p));
            if (cmd->usePeriodic) {
                VWrapAll(dr);
                DOTVP(drpq2, dr, dr);
            }
            drpq = rsqrt(drpq2);
            if ((real)Radius(p) < drpq)
                Radius(p) = cballs_store_search_bound(drpq);
        }
    }
}

local int pruningCells(struct  cmdline_data* cmd,
                        struct  global_data* gd,
                        int ifile, cellptr p, int lev)
{
    string routineName = "pruningCells";
    int i;
    nodeptr q;
    real qsize;
    real psize;

/*
// correct to consider these properties be computed...
        Mass(p) = 0.0;
        Weight(p) = 0.0;
        CLRV(cmpos);
*/

    Weight(p) = 0.0;
    Nb(p) = 0;
    Kappa(p) = 0.0;
    Radius(p) = 0.0;
    q = (nodeptr) p;
    qsize = Size(q);
    psize = Size(p);
    normal_walktree(cmd, gd, p, psize, q, qsize, 0);

    if (cmd->theta == 0.0)
        Radius(p) = cballs_store_search_bound(2.0*psize);
    else
        Radius(p) = cballs_store_search_bound((real)Radius(p)/cmd->theta);

    int n;
    
    //B
    if (deltaRadius <= 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid deltaRadius=%g\n", routineName, deltaRadius);
        return FAILURE;
    }

    if (!isfinite((double)Radius(p)) || Radius(p) < 0.0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid cell radius=%g\n", routineName, Radius(p));
        return FAILURE;
    }
    if (cmd->theta == 0.0
        || Radius(p) >= (real)(NbMax - 1) * deltaRadius)
        n = NbMax - 1;                              // radius overflow bin
    else
        n = (int)(Radius(p) / deltaRadius);
    if (n < 0 || n >= NbMax) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: cellRadius index out of range: n=%d NbMax=%d Radius=%g deltaRadius=%g\n",
                 routineName, n, NbMax, Radius(p), deltaRadius);
        return FAILURE;
    }
    cellRadius[n]++;
    //E


    
    if (Nb(p)<NbMax)
        cellhistNb[Nb(p)]++;

    if (Nb(p)>0) {
#ifndef NOWKAvg
        if (Weight(p) > 0.0) {
            Kappa(p) /= Weight(p);
        } else {
            Kappa(p) = 0.0;
        }
#else
        Kappa(p) /= Nb(p);
        Weight(p) /= Nb(p);
#endif
    } else {
#if defined(DEBUGTREE)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Nb = 0: %ld\n", Nb(p));
        return FAILURE;
#endif
    }

    for (i = 0; i < NSUB; i++) {                    // loop over existing subnodes
        if ((q = Subp(p)[i]) != NULL) {             // access each one in turn
            if (Type(q) == CELL) {                  // if also a cell, prune it
                if (pruningCells(cmd, gd, ifile, (cellptr) q, lev+1) == FAILURE)
                    return FAILURE;
            }
        }
    }

    return SUCCESS;
}

#endif

#ifdef DEBUGTREE
local void walkTree_printInfo(struct  cmdline_data* cmd, struct  global_data* gd,
                        nodeptr q, real qsize, int lev)
{
    nodeptr l;

    if (Type(q) == CELL) {
        fprintf(outtreeinfo, "%08d \t%08ld \t%16.8e \t%16.8e\n",
               lev, Nb(q), Radius(q), Kappa(q));
        for (l = More(q); l != Next(q); l = Next(l)) {
            walkTree_printInfo(cmd, gd, l,qsize/2, lev+1);
        }
    }
}
#endif

