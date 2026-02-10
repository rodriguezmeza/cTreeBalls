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

local int newtree(struct  cmdline_data* cmd, struct  global_data* gd, int);
local cellptr makecell(struct  cmdline_data* cmd, struct  global_data* gd, int);
local int loadbody(struct  cmdline_data*, struct  global_data*, bodyptr, int);
local int subindex(bodyptr, cellptr);
local void hackcellprop(struct  cmdline_data* cmd, struct  global_data* gd,
                        cellptr, real, int, int);
local int setradius(struct  cmdline_data* cmd, struct  global_data* gd,
                    cellptr, vector, real, int);
local void threadtree(struct  cmdline_data* cmd, struct  global_data* gd, 
                      nodeptr, nodeptr);

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

#ifdef BALLS4SCANLEV
local int scanLevelB4(struct  cmdline_data* cmd,
                      struct  global_data* gd, int ifile);
local void walktree_scan_lev_balls4(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    nodeptr q, int lev, int ifile, int scanLevel);
local INTEGER inodelevB4;
local INTEGER ibodyleftoutB4;
#endif

#ifdef PRUNING
local void pruningCells(struct  cmdline_data* cmd,
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

 To be called using: search=tree-omp-sincos

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

    debug_tracking_s("001", routineName);

    cpustart = CPUTIME;
    gd->bytes_tot_cells = 0;

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n%s: making the tree...\n", routineName);

#ifdef DEBUG
//B To debug cells:
    sprintf(cellsfilePath,"%s/cells%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outcells=fopen(cellsfilePath, "w")))
        error("\n%s: error opening file '%s' \n",routineName, cellsfilePath);
//E
#endif

    newtree(cmd, gd, ifile);
    roottable[ifile] = makecell(cmd, gd, ifile);
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    cpustartMiddle = CPUTIME;
    FindRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]);
    centerBodies(btab, nbody, ifile, roottable[ifile]);
    FindRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n%s: centerBodies-FindRootCenter CPU time: %lf %s\n",
                    routineName,
                    CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);
//E

    cpustartMiddle = CPUTIME;
    CLRV(Pos(roottable[ifile]));
    expandbox(cmd, gd, btab, nbody, ifile, roottable[ifile]);
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n%s: expandbox CPU time: %lf %s\n",
                    routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

    cpustartMiddle = CPUTIME;
    debug_tracking("002");
    DO_BODY(p, btab, btab+nbody) {
#ifdef BODY3ON
        Nbb(p) = 1;                                 // Check consistency with
                                                    //  smoothing... Correction
#endif
        loadbody(cmd, gd, p, ifile);                // Set body in a cell
        Nb(p) = 1;
        Radius(p) = 0.0;
#ifdef KappaAvgON
        KappaAvg(p) = Kappa(p);
#endif
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
    debug_tracking("003");
    hackcellprop(cmd, gd, roottable[ifile], gd->rSizeTable[ifile], 0, ifile);
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

#ifndef OCTREESMOOTHING
    ip = 0;

    cpustartMiddle = CPUTIME;
    debug_tracking("004");
    threadtree(cmd, gd, (nodeptr) roottable[ifile], NULL);
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
    debug_tracking("005");
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
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "tdepth = %d\n\n",gd->tdepthTable[ifile]);

    cpustartMiddle = CPUTIME;
    debug_tracking("006");
    scanLevel(cmd, gd, ifile);                      // Scan pivot and root trees
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%s: scanLevel CPU time: %lf %s\n",
                routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

#ifdef BALLS4SCANLEV
    if (gd->flagBalls4Scanlevel == FALSE) {
        cpustartMiddle = CPUTIME;
        scanLevelB4(cmd, gd, ifile);                // Scan pivot and root trees
        verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%s: scanLevelB4 CPU time: %lf %s\n",
                        routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);
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
    debug_tracking("007");
    pruningCells(cmd, gd, ifile, roottable[ifile], 0);
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
    sprintf(treeinfofilePath,"%s/treeinfo%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outtreeinfo=fopen(treeinfofilePath, "w")))
        error("\n%s: error opening file '%s' \n",routineName, treeinfofilePath);
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

    debug_tracking_s("008... final",routineName);

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
        if (strnull(cmd->rsmooth)) {
            if (scanopt(cmd->options, "fix-rsmooth")) {
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
                if (scanopt(cmd->options, "default-rsmooth")) {
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
    vector xmin, xmax;

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

local int newtree(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    roottable[ifile] = NULL;
    gd->ncellTable[ifile] = 0;

    return SUCCESS;
}

// deallocate memory tree
global int freeTree(struct  cmdline_data* cmd, struct  global_data* gd)
{
    nodeptr p;
    nodeptr freecell = NULL;
    int ifile;

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        p = (nodeptr) roottable[ifile];
        while (p != NULL)
            if (Type(p) == CELL) {
                Next(p) = freecell;
                freecell = p;
                p = More(p);
            } else
                p = Next(p);

        p = freecell;
        while (p != NULL) {
            free(p);
            p = Next(p);
        }
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++) {
            p = (nodeptr) roottable[ifile];
            while (p != NULL)
                if (Type(p) == CELL) {
                    Next(p) = freecell;
                    freecell = p;
                    p = More(p);
                } else
                    p = Next(p);
            
            p = freecell;
            while (p != NULL) {
                free(p);
                p = Next(p);
            }
        }
    }

    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        roottable[ifile] = NULL;
        gd->ncellTable[ifile] = 0;
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
    if (scanopt(cmd->options, "read-mask"))
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

    if (cmd->verbose>=VERBOSENORMALINFO)
        verb_print(cmd->verbose, "treeload expandbox: rSize = %lf\n",
                   gd->rSizeTable[ifile]);

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
                if (scanopt(cmd->options, "no-check-two-bodies-eq-pos")) {
                    DO_COORD(k) {
//                        Pos(p)[k] += EPSILON*grandom(0.0, 0.01*qsize);
                        Pos(p)[k] += EPSILON*grandom(0.0, 1.0);
                        DO_COORD(k) {
                            verb_log_print(cmd->verbose_log,gd->outlog,
                            "CorrectedPos[k]: %g %g and correction: %le\n",
                            Pos(p)[k],Pos(Subp(q)[qind])[k],
                            EPSILON*grandom(0.0, 1.0));
//                            EPSILON*grandom(0.0, 0.01*qsize));
                        }
                    }
                    Update(p) = FALSE;
                    gd->sameposcount++;
                    sameposcount++;
                    goto startagain;
                } else {
                    error("loadbody: two bodies have same position...\n\t%s",
                          "consider the option: 'no-check-two-bodies-eq-pos'");
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

local void hackcellprop(struct  cmdline_data* cmd, struct  global_data* gd,
                       cellptr p, real psize, int lev, int ifile)
{
    string routineName = "hackcellprop";
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;

    gd->tdepthTable[ifile] = MAX(gd->tdepthTable[ifile], lev);
    cellhist[lev]++;
    Level(p) = lev;                                 // To set scanLevel
    Mass(p) = 0.0;
    Weight(p) = 0.0;
    Nb(p) = 0;
    Kappa(p) = 0.0;
#ifdef KappaAvgON
    KappaAvg(p) = 0.0;
#endif
    CLRV(cmpos);

    for (i = 0; i < NSUB; i++) {
        if ((q = Subp(p)[i]) != NULL) {
            subnhist[lev]++;
            if (Type(q) == CELL) {
                hackcellprop(cmd, gd, (cellptr) q, psize/2, lev+1, ifile);
            }
            Selected(p) |= Selected(q);             // bodies belonging to a cell
            Update(p) |= Update(q);
            if (scanopt(cmd->options, "read-mask"))
                Mask(p) |= Mask(q);

            Mass(p) += Mass(q);
            Weight(p) += Weight(q);                 // sum of all weight of q
            if ( Type(q) == CELL) {
                if (scanopt(cmd->options, "read-mask"))
                    Mask(p) |= Mask(q);
                Nb(p) += Nb(q);
                Kappa(p) += Weight(q)*Kappa(q);     // Kappa(q) average at cell q
#ifdef KappaAvgON
                KappaAvg(p) += KappaAvg(q);         // sum of all kappa at cell q
#endif
            } else {
                if (Type(q) == BODY) {
                    Nb(p) += 1;
                    Kappa(p) += Weight(q)*Kappa(q);
#ifdef KappaAvgON
                    KappaAvg(p) += Kappa(q);
#endif
                } else if (Type(q) == BODY3) {      // To set smoothing body
                    Nb(p) += 1;
                    Kappa(p) += Weight(q)*Kappa(q);
                }
            }
            MULVS(tmpv, Pos(q), Mass(q));
            ADDV(cmpos, cmpos, tmpv);
        } // ! q no NULL
    } // ! loop i: 0 -> NSUB-1
//B Smooth(ing) section
    if (Nb(p)==cmd->nsmooth) {                      // Correct to <=
        NTOT[0]=NTOT[0]+1;
    }
//E
    if (scanopt(cmd->options, "center-of-mass")) {
        if (Mass(p) > 0.0) {
            DIVVS(cmpos, cmpos, Mass(p));
        } else {
            SETV(cmpos, Pos(p));
        }
    } else
        SETV(cmpos, Pos(p));

#define EPSILON 1.0E-16
// Here there appears an error for big numbers of points such 201 millions...
// See line above and uncomment DIVVS(cmpos, cmpos, Mass(p)); line (not working!)
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 || Pos(p)[k] + psize/2 <= cmpos[k]) {
            if (psize/2 > 2.710505e-20 + EPSILON)
            error("%s: tree structure error: %d %le %le %le %le\n",
                  routineName, k, cmpos[k],
                  Pos(p)[k] - psize/2, Pos(p)[k] + psize/2, psize/2);
            else {
                if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log,gd->outlog,
                        "%s: tree structure warning! psize/2 to small: %le \n",
                        routineName, psize/2);
            }
        }
#undef EPSILON

    setradius(cmd, gd, p, cmpos, psize, ifile);

    SETV(Pos(p), cmpos);
    if (Nb(p)>0) {
        Kappa(p) /= Nb(p);
    } else {
#if defined(DEBUGTREE)
        error("%s: Nb = 0: %ld\n", routineName, Nb(p));
#endif
    }
}

// Parameter theta controls size of the cell.
// theta from 0 to 5:
//  0 always open cells (complexity N^2);
//  1 is the default value.
local int setradius(struct  cmdline_data* cmd, struct  global_data* gd,
                    cellptr p, vector cmpos, real psize, int ifile)
{
    real bmax2, d;
    int k;

    if (cmd->theta == 0.0)
        Radius(p) = 2 * gd->rSizeTable[ifile];
    else if (gd->sw94) {
        bmax2 = 0.0;
		DO_COORD(k) {
            d = cmpos[k] - Pos(p)[k] + psize/2; 
            bmax2 += rsqr(MAX(d, psize - d));   
        }
        Radius(p) = rsqrt(bmax2) / cmd->theta;
    } else if (gd->bh86)
        Radius(p) = psize / cmd->theta;
    else {
        Radius(p) = (psize/cmd->theta) * rsqrt((real)(NDIM))/2.0;
//        DISTV(d, cmpos(p), Pos(p));                 // find offset from center
//        Radius(p) = psize / cmd->theta + d;         // use size plus offset
    }

    Size(p) = psize;

    int n;
    n = (int) (Radius(p) / deltaRadius);
    (cellRadius[n])++;

    return SUCCESS;
}

local void threadtree(struct  cmdline_data* cmd,
                      struct  global_data* gd, nodeptr p, nodeptr n)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    bodyptr q;

    Next(p) = n;
    if (Type(p) == CELL) {
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        ndesc = 0;
        for (i = 0; i < NSUB; i++)
            if (Subp(p)[i] != NULL)
                desc[ndesc++] = Subp(p)[i];
        More(p) = desc[0];
        desc[ndesc] = n;
        for (i = 0; i < ndesc; i++)
            threadtree(cmd, gd, desc[i], desc[i+1]);
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


#ifdef BALLS4SCANLEV
// need to free allocated here memory
local int scanLevelB4(struct  cmdline_data* cmd,
                      struct  global_data* gd, int ifile)
{
    string routine_name = "scanLevelB4";

    gd->nnodescanlevTableB4[ifile] =
        gd->ncellTable[ifile]+gd->nbodyTable[ifile];
    nodetablescanlevB4[ifile] =
        (nodeptr *) allocate(gd->nnodescanlevTableB4[ifile] * sizeof(nodeptr));
    gd->bytes_tot += gd->nnodescanlevTableB4[ifile]*sizeof(nodeptr);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n%s: Allocated %g MByte for (%d) scan nodetab storage.\n",
                routine_name, INMB*gd->nnodescanlevTableB4[ifile]*sizeof(nodeptr),
                gd->nnodescanlevTableB4[ifile]);

    inodelevB4 = 0;
    ibodyleftoutB4 = 0;

//B scan tree up to the smallest cells
    walktree_scan_lev_balls4(cmd, gd,
                             (nodeptr)roottable[ifile], 0,
                             ifile, gd->tdepthTable[ifile]);
//E
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\n%s: Found %d nodes to scan at upper most level %d.\n",
                        routine_name, inodelevB4, gd->tdepthTable[ifile]);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t%ld particles were included to scan at that level.\n",
                        ibodyleftoutB4);

    gd->nnodescanlevTableB4[ifile] = inodelevB4;


    nodeptr p;
    real rmax=0, rmin=cmd->rangeN;
    INTEGER numCells=0, numBodies=0;
    for (INTEGER i = 0; i < gd->nnodescanlevTableB4[ifile]; i++) {
        p = nodetablescanlevB4[ifile][i];
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
                "%s: number of cell and of bodies and total nodes: %ld %ld %ld\n",
                routine_name, numCells, numBodies, numCells+numBodies);

    return SUCCESS;
}

local void walktree_scan_lev_balls4(struct  cmdline_data* cmd,
                                    struct  global_data* gd,
                                    nodeptr q, int lev, int ifile, int scanLevel)
{
    nodeptr p,g,h,l;

    if ( lev+1 <= scanLevel ) {
        if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Radius(l) < gd->deltaRmin*THETA) {
                    nodetablescanlevB4[ifile][inodelevB4] = l;
                    inodelevB4++;
                    ibodyleftoutB4++;
                } else {
                    if (Type(l)==CELL)
                        walktree_scan_lev_balls4(cmd, gd,
                                                 l, lev+1, ifile, scanLevel);
                    else
                        ibodyleftoutB4++;
                }
            }
        } else {
            nodetablescanlevB4[ifile][inodelevB4] = q;
            inodelevB4++;
            ibodyleftoutB4++;
        }
    }
}
#endif // ! BALLS4SCANLEV


#ifdef PRUNING
// use in combination with options=center-of-mass
local void normal_walktree(struct  cmdline_data* cmd,
                           struct  global_data* gd,
                           cellptr p, real psize,  nodeptr q, real qsize, int lev)
{
    nodeptr l;
    real drpq, drpq2;
    vector dr;
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
            Nb(p) += 1;
            Kappa(p) += Weight(q)*Kappa(q);
            DOTPSUBV(drpq2, dr, Pos(q), Pos(p));
            if (cmd->usePeriodic) {
                VWrapAll(dr);
                DOTVP(drpq2, dr, dr);
            }
            drpq = rsqrt(drpq2);
            if (Radius(p) < drpq)
                Radius(p) = drpq;
        }
    }
}

local void pruningCells(struct  cmdline_data* cmd,
                        struct  global_data* gd,
                        int ifile, cellptr p, int lev)
{
    int i;
    nodeptr q;
    real qsize;
    real psize;

/*
// correct to consider these properties be computed...
        Mass(p) = 0.0;
        Weight(p) = 0.0;
    #ifdef KappaAvgON
        KappaAvg(p) = 0.0;
    #endif
        CLRV(cmpos);
*/

    Nb(p) = 0;
    Kappa(p) = 0.0;
    Radius(p) = 0.0;
    q = (nodeptr) p;
    qsize = Size(q);
    psize = Size(p);
    normal_walktree(cmd, gd, p, psize, q, qsize, 0);

    if (cmd->theta == 0.0)
        Radius(p) = 2.0*psize;                      // always open cell
    else
        Radius(p) /= cmd->theta;

    int n;
    n = (int) (Radius(p) / deltaRadius);
    (cellRadius[n])++;
    if (Nb(p)<NbMax)
        cellhistNb[Nb(p)]++;

    if (Nb(p)>0) {
        Kappa(p) /= Nb(p);
    } else {
#if defined(DEBUGTREE)
        error("Nb = 0: %ld\n", Nb(p));
#endif
    }

    for (i = 0; i < NSUB; i++) {                    // loop over existing subnodes
        if ((q = Subp(p)[i]) != NULL) {             // access each one in turn
            if (Type(q) == CELL)                    // if also a cell, prune it
                pruningCells(cmd, gd, ifile, (cellptr) q, lev+1);
        }
    }
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

