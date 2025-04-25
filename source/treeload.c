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
local void hackCellProp(struct  cmdline_data* cmd, struct  global_data* gd, 
                        cellptr, real, int, int);
local int setRadius(struct  cmdline_data* cmd, struct  global_data* gd, 
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
local int cellhistNb[NbMax];
local int cellRadius[NbMax];
local real deltaRadius;
//E

local void walktree_selected(nodeptr, real);        // To see the bodies belonging
                                                    //  to a cell

//local INTEGER icell;                                // To debug cells
local INTEGER inode;
#ifdef DEBUG
local void walktree_hit(struct  cmdline_data* cmd, struct  global_data* gd,
                        nodeptr, real);
#endif

//B not needed in the public version... they are part of an addon
//local INTEGER Nc1;
//local INTEGER Nc2;
//E

local int scanLevel(struct  cmdline_data* cmd, struct  global_data* gd, int);
local void walktree_index_scan_lev(nodeptr, int, int, int);

local INTEGER inodelev;
local INTEGER ibodyleftout;

//B Root nodes:
local void walktree_index_scan_lev_root(struct cmdline_data* cmd, 
                                        struct  global_data* gd,
                                        nodeptr, int, int);
local INTEGER inodelev_root;
local INTEGER ibodyleftout_root;
//E

local INTEGER isel, inosel;
local int smoothBodies(struct  cmdline_data*, 
                       struct  global_data* gd, bodyptr, INTEGER);

local int save_nodes(struct  cmdline_data*, 
                     struct  global_data* gd, int ifile);
local int save_nodes_root(struct  cmdline_data*, 
                          struct  global_data* gd, int ifile);

local char cellsfilePath[MAXLENGTHOFFILES];
local FILE *outcells;

global int MakeTree(struct  cmdline_data* cmd, 
                    struct  global_data* gd,
                    bodyptr btab, INTEGER nbody, int ifile)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = CPUTIME;
    gd->bytes_tot_cells = 0;

#ifdef DEBUG
//B To debug cells:
    sprintf(cellsfilePath,"%s/cells%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(outcells=fopen(cellsfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",cellsfilePath);
//E
#endif

    newtree(cmd, gd, ifile);
    roottable[ifile] = makecell(cmd, gd, ifile);
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    findRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]);
    centerBodies(btab, nbody, ifile, roottable[ifile]);
    findRootCenter(cmd, gd, btab, nbody, ifile, roottable[ifile]);
//E

    CLRV(Pos(roottable[ifile]));
    expandbox(cmd, gd, btab, nbody, ifile, roottable[ifile]);

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
    hackCellProp(cmd, gd, roottable[ifile], gd->rSizeTable[ifile], 0, ifile);

//B To bin cell's radius... Check!!! (bins logscale)
    real rBin;
    verb_log_print(cmd->verbose_log,gd->outlog,"\nmaketree: radius histogram:\n");
    for (i = 0; i < NbMax; i++) {
        rBin = ((int)i)*gd->rSizeTable[ifile]/deltaRadius;
        verb_log_print(cmd->verbose_log,gd->outlog,"%g %d\n", rBin, cellRadius[i]);
    }
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");
//E

//B Smooth(ing) section
    if (!gd->flagSmooth)                             // Flag to smooth bodies
        if (scanopt(cmd->options, "smooth")) {
            printf("NTOT = %ld \t%ld \t%ld\n",
               gd->nsmooth[0], NTOT[0], gd->nsmooth[0]*NTOT[0]);
            gd->nbodysm = NTOT[0];
            bodytabsm = (bodyptr) allocate((NTOT[0]) * sizeof(body));
            gd->bytes_tot += (NTOT[0])*sizeof(body);
            verb_print(cmd->verbose,
                "Allocated %g MByte for (smooth %ld) particle storage.\n",
                (NTOT[0])*sizeof(body)*INMB, (NTOT[0]));
        }
//E

    ip = 0;
    threadtree(cmd, gd, (nodeptr) roottable[ifile], NULL);
    if (cmd->verbose>=VERBOSEDEBUGINFO) {
        verb_print(cmd->verbose,
                   "threadtree: number ip of selected cells = %ld\n",ip);
        verb_print(cmd->verbose,
                   "%d real node (range of nodes to search: >nc1 && <nc2)\n",inode);
    }
    gd->rnnode = inode;

    walktree_selected((nodeptr) roottable[ifile],   // Smooth(ing) section
                      gd->rSizeTable[ifile]);
    isel=0, inosel=0;
    DO_BODY(p,btab,btab+nbody) {                    // See bodies belonging
        if (Selected(p))                            //  to a cell
            isel++;
        else
            inosel++;
    }
    if (cmd->verbose>=VERBOSEDEBUGINFO) {
        verb_print(cmd->verbose,
                   "\nSelected vs NotSelected and total: %ld %ld %ld\n\n",
                   isel, inosel, isel + inosel);
        verb_print(cmd->verbose,
                   "tdepth = %d\n\n",gd->tdepthTable[ifile]);
    }

    scanLevel(cmd, gd, ifile);                      // Scan pivot and root trees
    smoothBodies(cmd, gd, btab, nbody);             // Smooth cells

//B Histogram useful to smooth cells
    verb_log_print(cmd->verbose_log,gd->outlog,
        "\nmaketree: Nb histogram:\n");
    for (i = 0; i < NbMax; i++)
        verb_log_print(cmd->verbose_log,gd->outlog,
            "%d %ld\n", i, cellhistNb[i]);
    verb_log_print(cmd->verbose_log,gd->outlog,"\n");
//E
    gd->bytes_tot += gd->bytes_tot_cells;
    if (cmd->verbose>=VERBOSENORMALINFO) {
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for (%d) cells storage.\n",
                   gd->bytes_tot_cells*INMB, gd->ncellTable[ifile]);
        verb_print(cmd->verbose,
                   "\nmaketree: root number of bodies = %ld\n",
                   Nb(roottable[ifile]));
    }

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

    gd->cputree = CPUTIME - cpustart;
    verb_print(cmd->verbose, "\nmaking tree CPU time : %lf %s\n\n",
               gd->cputree, PRNUNITOFTIMEUSED);

    return SUCCESS;
}


//B Smooth(ing) section
local int smoothBodies(struct  cmdline_data* cmd, 
                       struct  global_data* gd, bodyptr btab, INTEGER nbody)
{
    bodyptr p;

    if (!gd->flagSmooth || !gd->flagSetNbNoSel) {
        if ( (scanopt(cmd->options, "smooth")
              && scanopt(cmd->options, "set-Nb-noSel")) ) {
            printf("NTOT = %ld \t%ld \t%ld\n\n",
                   gd->nsmooth[0], NTOT[0]+inosel,
                   gd->nsmooth[0]*NTOT[0]+inosel);
            bodytabSel = (bodyptr) allocate((NTOT[0]+inosel) * sizeof(body));
            gd->bytes_tot += (NTOT[0]+inosel)*sizeof(body);
            verb_print(cmd->verbose,
                       "Allocated %g MByte for (smooth) particle (%ld) storage.\n",
                       (NTOT[0]+inosel)*sizeof(body)*INMB, (NTOT[0]+inosel));
            int ipcount=0;
            bodyptr q = bodytabSel;
            DO_BODY(p,bodytabsm,bodytabsm+gd->nbodysm) {
                ipcount++;
                Id(q) = ipcount;
// BODY3
                Type(q) = BODY3;
#ifdef BODY3ON
                Nbb(q) = Nbb(p);
#endif
                Mass(q) = Mass(p);
                Weight(q) = Weight(p);
                SETV(Pos(q), Pos(p));
                Kappa(q) = Kappa(p);
                q++;
            }
            verb_print(cmd->verbose,"Added %ld smoothed cells...\n",ipcount);
            DO_BODY(p,btab,btab+nbody) {
                if (!Selected(p)) {
                    ipcount++;
                    Id(q) = ipcount;
                    Type(q) = BODY;
// BODY3
#ifdef BODY3ON
                    Nbb(q) = 1;
#endif
                    Mass(q) = Mass(p);
                    Weight(q) = Weight(p);
                    SETV(Pos(q), Pos(p));
                    Kappa(q) = Kappa(p);
                    q++;
                }
            }
            verb_print(cmd->verbose,"Added %ld total bodies...\n",ipcount);
            gd->nbodySel = ipcount;
        }
    }

    return SUCCESS;
}
//E

//B BALLS :: SCANLEV
local int scanLevel(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    int i;
    
#ifdef BALLS
    
    if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") ) {
        cmd->scanLevel = gd->tdepthTable[ifile] + 1 + gd->scanLevelMin[0];
        verb_print(cmd->verbose,
                   "\tsmoothCellMin: fixing scanLevel to: %d\n",
                   cmd->scanLevel);
    } else {
        if (scanopt(cmd->options, "set-default-param")) {
            verb_print(cmd->verbose, "\tfixing scanLevel to tdepth-1...\n");
            cmd->scanLevel = MAX(gd->tdepthTable[ifile]-1,3);
            verb_print(cmd->verbose, "\tfinal value is %d.\n", cmd->scanLevel);
        } else {
            if (cmd->scanLevel > gd->tdepthTable[ifile]) {
                verb_print(cmd->verbose,
                           "Warning! tree depth (%d) is less than scanLevel (%d)...\n",
                           gd->tdepthTable[ifile], cmd->scanLevel);
                verb_print(cmd->verbose, "\tfixing to tdepth-1...\n");
                cmd->scanLevel = MAX(gd->tdepthTable[ifile]-1,3);
                verb_print(cmd->verbose, "\tfinal value is %d.\n", cmd->scanLevel);
            }
        }
    }
    
    if (cmd->useLogHist) { // ! rminHist = 0 not allowed when useLogHist is true
        if (cmd->verbose >= VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
            "\n(Only in log-scale) deltaR is %g and root size at scanLevel is %g.\n",
            gd->deltaR, gd->rSizeTable[ifile]/rpow(2.0,cmd->scanLevel));
        }
        i = 0;
        while (gd->deltaR < gd->rSizeTable[ifile]/rpow(2.0,i)) i++;
        verb_print(cmd->verbose,
                   "\t\t\tSuggested scanLevel is %d, where root size is %g.\n",
                   i, gd->rSizeTable[ifile]/rpow(2.0,i));
    }
    inodelev = 0;
    ibodyleftout = 0;
    if (cmd->scanLevel==0) {
        gd->nnodescanlevTable[ifile] = 0;
        gd->bytes_tot += gd->nnodescanlevTable[ifile]*sizeof(nodeptr);
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for (%d) scan nodetab storage.\n",
                   INMB*gd->nnodescanlevTable[ifile]*sizeof(nodeptr),
                   gd->nnodescanlevTable[ifile]);
        gd->nnodescanlev = 1;
    } else {
        gd->nnodescanlevTable[ifile] =
        gd->ncellTable[ifile]+gd->nbodyTable[ifile];
        nodetablescanlev[ifile] =
        (nodeptr *) allocate(gd->nnodescanlevTable[ifile] * sizeof(nodeptr));
        //
        gd->bytes_tot += gd->nnodescanlevTable[ifile]*sizeof(nodeptr);
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for (%d) scan nodetab storage.\n",
                   INMB*gd->nnodescanlevTable[ifile]*sizeof(nodeptr),
                   gd->nnodescanlevTable[ifile]);
        
        if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") ) {
            walktree_index_scan_lev((nodeptr)roottable[ifile], 0,
                                    ifile, cmd->scanLevel);
        } else {
            // Check this line. May be it is correct:
            //            cmd->scanLevel = MAX(gd->tdepthTable[ifile],3);
            //            verb_print(cmd->verbose, "\tfinal value is %d.\n", cmd->scanLevel);
            if (cmd->scanLevel > gd->tdepthTable[ifile]) {
                verb_print(cmd->verbose,
                           "Warning! tree depth (%d) is less than scanLevel (%d)...\n",
                           gd->tdepthTable[ifile], cmd->scanLevel);
                verb_print(cmd->verbose, "\tfixing to tdepth-1...\n");
                cmd->scanLevel = MAX(gd->tdepthTable[ifile]-1,3);
                verb_print(cmd->verbose, "\tfinal value is %d.\n", cmd->scanLevel);
            }
            walktree_index_scan_lev((nodeptr)roottable[ifile], 0,
                                    ifile, cmd->scanLevel);
        }

        if (cmd->verbose > VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "\nFound %d nodes to scan at level %d.\n",
                       inodelev, cmd->scanLevel);
        }
        // Freeing some segment of memory will be necessary
        gd->nnodescanlevTable[ifile] = inodelev;
        save_nodes(cmd, gd, ifile);
    } // ! scanLevel==0

//B socket:
#ifdef ADDONS
#include "tree_include.h"
#endif
//E

    //B Root nodes to scan:
    inodelev_root = 0;
    ibodyleftout_root = 0;
    
    if (cmd->scanLevelRoot==0) {
        gd->nnodescanlev_rootTable[ifile] = 0;
        gd->bytes_tot += gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr);
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for (%d) scan root nodetab storage.\n",
                   INMB*gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr),gd->nnodescanlev_rootTable[ifile]);
    } else {
        gd->nnodescanlev_rootTable[ifile] =
        gd->ncellTable[ifile]+gd->nbodyTable[ifile];
        nodetablescanlev_root[ifile] =
        (nodeptr *) allocate(gd->nnodescanlev_rootTable[ifile]
                             * sizeof(nodeptr));
        gd->bytes_tot += gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr);
        verb_print(cmd->verbose,
                   "\nAllocated %g MByte for (%d) scan root nodetab storage (%ld cells).\n",
                   INMB*gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr),
                   gd->nnodescanlev_rootTable[ifile],gd->ncellTable[ifile]);
        walktree_index_scan_lev_root(cmd, gd, (nodeptr)roottable[ifile], 0, ifile);
        if (cmd->verbose >= VERBOSENORMALINFO) {
            verb_print(cmd->verbose,
                       "\nFound %d root nodes to scan at level %d.\n",
                       inodelev_root, cmd->scanLevelRoot);
            verb_print(cmd->verbose,
                       "%ld particles were left out at root level.\n",
                       ibodyleftout_root);
        }
        gd->nnodescanlev_rootTable[ifile] = inodelev_root;
        if (cmd->verbose==3)
            save_nodes_root(cmd, gd, ifile);
    }
    
#endif // ! BALLS
    
    gd->Rcell[0] = gd->rSizeTable[ifile];
    for (i = 1; i <= gd->tdepthTable[ifile]; i++)
        gd->Rcell[i] = gd->Rcell[i-1]/2;
    
    verb_print(cmd->verbose, "\nMaximum and minimum cell size: %e %e\n",
               gd->Rcell[0],gd->Rcell[gd->tdepthTable[ifile]-1]);
    
#ifdef BALLS

    if (cmd->useLogHist==TRUE) {

    if (gd->flagSmoothCellMin) {
        gd->scanLevelMin[0] = gd->scanLevelMin[1];
        verb_print(cmd->verbose,
                   "\tsmoothCellMin: fixing scanLevelMin[0] to: %d\n",
                   gd->scanLevelMin[0]);
        if (gd->scanLevelMin[0] == 0)
            gd->rminCell[0] = 0.;
        else
            gd->rminCell[0] = gd->Rcell[gd->tdepthTable[ifile]+gd->scanLevelMin[0]];
    } else {
        if (gd->scanLevelMin[0] == 0)
            gd->rminCell[0] = 0.;
        else
            gd->rminCell[0] =
            gd->Rcell[gd->tdepthTable[ifile]+gd->scanLevelMin[0]];
    }

    verb_print(cmd->verbose,
               "Cell size at scanLevelMin (%d) and scale factor: %e %e\n",
               gd->scanLevelMin[0], gd->rminCell[0], gd->rminCell[1]);
    if (gd->rminCell[0] > gd->deltaRmax)
        verb_print(cmd->verbose,
        "Warning! Cell size at scanLevelMin is greatear than deltaRmax: %e %e\n",
        gd->rminCell[0], gd->deltaRmax);
    if (gd->rminCell[0] > gd->deltaRmin)
        verb_print(cmd->verbose,
        "Warning! Cell size at scanLevelMin is greatear than deltaRmin: %e %e\n",
        gd->rminCell[0], gd->deltaRmin);
    for (i=1; i<=cmd->sizeHistN-1; i++)
        if (gd->rminCell[0] < gd->ddeltaRV[i])
            break;
    verb_print(cmd->verbose,
               "rminCell, ddeltaRV and deltaRV (at n = %d): %g %g %g\n",
               i+1, gd->rminCell[0], gd->ddeltaRV[i], gd->deltaRV[i+1]);
    
    if (gd->scanLevelMin[0] == 0) {
        gd->rminCell[1] = cmd->rangeN;
        verb_print(cmd->verbose,
                   "\tfixing rminCell[1] to: %g\n", gd->rminCell[1]);
    } else {
        gd->rminCell[1] = gd->deltaRV[i+1];
        verb_print(cmd->verbose,
                   "\tfixing rminCell[1] to: %g\n", gd->rminCell[1]);
        if (gd->rminCell[1]>cmd->rangeN)
            verb_print(cmd->verbose,
                       "Warning! rminCell[1] is greatear than rangeN\n");
    }
    verb_print(cmd->verbose,
    "Cell size at scanLevelMin (%d) and scale factor (Modified values): %e %e\n",
    gd->scanLevelMin[0], gd->rminCell[0], gd->rminCell[1]);
    //E Root nodes
    
    if (gd->infilefmt_int == INTAKAHASI)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "Unit sphere (Takahasi): (S/N)^(1/2): %g\n\n",
                       rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0));
    
} // ! useLogHist = true
    else {
//        error("CheckParameters: can´t have normal hist and BALLS definition (useLogHist=%d)\nSet useLogHist = true\n",
//              cmd->useLogHist);
        verb_print(cmd->verbose,
                   "Warning! can´t have normal hist and BALLS definition (useLogHist=%d)\nSet useLogHist = true\n",
                   cmd->useLogHist);
    }
#endif // ! BALLS

    if (scanopt(cmd->options, "smooth-pivot")) {
        if (cmd->useLogHist) {
            verb_print(cmd->verbose,
                       "deltaRV min and max: %lg %lg\n",
                       gd->deltaRmin, gd->deltaRmax);
        } else {
            verb_print(cmd->verbose,
                       "deltaR=%lf normal scale):\n",gd->deltaR);
        }
    }

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
    if (scanopt(cmd->options, "smooth-pivot")) {
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
#undef ARCMINTORAD
    } else // ! smooth-pivot
        gd->rsmooth[0] = 0;                         // setting a safe value

    if (scanopt(cmd->options, "smooth-pivot")) {
        verb_print(cmd->verbose,
                   "rsmooth and rminHist %% of rangeN: %lg %lg\n",
                   100*gd->rsmooth[0]/cmd->rangeN, 100*cmd->rminHist/cmd->rangeN);
    }

    return SUCCESS;
}
//E BALLS :: SCANLEV

global int findRootCenter(struct  cmdline_data* cmd, 
                          struct  global_data* gd,
                          bodyptr btab, int nbody, int ifile, cellptr root)
{
    real len;
    bodyptr p;
    int k;
    vector xmin, xmax;

    DO_COORD(k)
        xmin[k] = xmax[k] = Pos(btab)[k];

    DO_BODY(p, btab, btab+nbody)
        DO_COORD(k) {
            if (Pos(p)[k] > xmax[k])
                xmax[k] = Pos(p)[k];
            if (Pos(p)[k] < xmin[k])
                xmin[k] = Pos(p)[k];
        }

    DO_COORD(k) {
        Pos(root)[k] = (xmax[k]+xmin[k])/2;
        if (cmd->verbose>=2)
            verb_print(cmd->verbose,
                       "findRootCenter: Pos(root) = %lf\n", Pos(root)[k]);
    }

    for(k=0, len=xmax[0]-xmin[0]; k<NDIM; k++)
        if((xmax[k]-xmin[k])>len)
            len=xmax[k]-xmin[k];

    if (cmd->verbose>=2)
        verb_print(cmd->verbose, "findRootCenter: len = %lf\n", len);

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

local cellptr makecell(struct  cmdline_data* cmd,
                       struct  global_data* gd, int ifile)
{
    cellptr c;
    int i;
 
    c = (cellptr) allocate(sizeof(cell));
    Type(c) = CELL;
    Nb(c) = 0;                                      // To smooth cells
    Update(c) = FALSE;
    for (i = 0; i < NSUB; i++)                  
        Subp(c)[i] = NULL;
//    (gd->ncellTable[ifile])++;
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

#define EPSILON 1.0E-7                              // Choose well.
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

// Keep it in order to study MPI programming...
//    verb_print_debug(1, "\nloadbody:: Aqui voy (2)\n");
//    MPI_Barrier(MPI_COMM_WORLD);

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
//                        Pos(p)[k] += EPSILON*grandom(0.0, 0.01*gd->Box[k]);
                        Pos(p)[k] += EPSILON*grandom(0.0, 0.01*qsize);
//                        Pos(p)[k] += EPSILONFLOAT*grandom(0.0, 0.01*gd->Box[k]);
//                        Pos(p)[k] += grandom(0.0, 0.1*gd->Box[k]);
                        DO_COORD(k) {
                            verb_log_print(cmd->verbose_log,gd->outlog,
                            "CorrectedPos[k]: %g %g and correction: %le\n",
                            Pos(p)[k],Pos(Subp(q)[qind])[k],
                            EPSILON*grandom(0.0, 0.01*qsize));
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

local void hackCellProp(struct  cmdline_data* cmd, struct  global_data* gd,
                       cellptr p, real psize, int lev, int ifile)
{
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
            if (Type(q) == CELL)
                hackCellProp(cmd, gd, (cellptr) q, psize/2, lev+1, ifile);
            Selected(p) |= Selected(q);             // To see the bodies belonging
                                                    //  to a cell
            Update(p) |= Update(q);
            Mass(p) += Mass(q);
            Weight(p) += Weight(q);
            if ( Type(q) == CELL) {
                Nb(p) += Nb(q);
                Kappa(p) += Mass(q)*Kappa(q);
#ifdef KappaAvgON
                KappaAvg(p) += KappaAvg(q);
#endif
            } else {
                if (Type(q) == BODY) {
                    Nb(p) += 1;
                    Kappa(p) += Mass(q)*Kappa(q);
#ifdef KappaAvgON
                    KappaAvg(p) += Kappa(q);
#endif
                } else if (Type(q) == BODY3) {      // To set smoothing body
                    Nb(p) += 1;
                    Kappa(p) += Mass(q)*Kappa(q);
                }
            }
            MULVS(tmpv, Pos(q), Mass(q));
            ADDV(cmpos, cmpos, tmpv);           
        }
    }
//B Smooth(ing) section
    if (Nb(p)==gd->nsmooth[0]) {                     // Correct to <=
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
            error("hackCellProp: tree structure error: %d %le %le %le %le\n",
                  k, cmpos[k], Pos(p)[k] - psize/2, Pos(p)[k] + psize/2, psize/2);
            else {
                if (cmd->verbose_log>=3)
                verb_log_print(cmd->verbose_log,gd->outlog,
                "hackCellProp: tree structure warning! psize/2 to small: %le \n",
                psize/2);
            }
        }
#undef EPSILON

    setRadius(cmd, gd, p, cmpos, psize, ifile);
    SETV(Pos(p), cmpos);
    if (Nb(p)>0) {
        Kappa(p) /= Nb(p);
    } else
        error("hackCellProp: Nb = 0: %ld\n", Nb(p));
}

// Parameter theta controls size of the cell.
// theta from 0 to 5:
//  0 always open cells (complexity N^2);
//  1 is the default value.
local int setRadius(struct  cmdline_data* cmd, struct  global_data* gd,
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
        //B Smooth(ing) section
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        if (!gd->flagSmooth) {
            if (scanopt(cmd->options, "smooth")) {
                if (Nb(p)==gd->nsmooth[0]) {
                    ip++;
                    q = ip + bodytabsm -1;
                    Id(q) = ip;
                    Type(q) = BODY3;                // To set smoothing body
#ifdef BODY3ON
                    Nbb(q) = Nb(p);
                    Nbb(p) = Nb(p);
#endif
                    Mass(q) = Mass(p);
                    Weight(q) = Weight(p);
                    SETV(Pos(q), Pos(p));
                    Kappa(q) = Kappa(p);
                    Selected(p) = TRUE;             // To see the bodies belonging
                    Selected(q) = TRUE;             //  to a cell
                }
            }
        }
        //E
        ndesc = 0;
        for (i = 0; i < NSUB; i++)
            if (Subp(p)[i] != NULL)
                desc[ndesc++] = Subp(p)[i];
        More(p) = desc[0];
        desc[ndesc] = n;
        for (i = 0; i < ndesc; i++)
            threadtree(cmd, gd, desc[i], desc[i+1]);
    }
}

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

#ifdef BALLS
//
// Unify: walktree_index_scan_lev and walktree_index_scan_lev_root
//
local void walktree_index_scan_lev(nodeptr q, int lev, int ifile, int scanLevel)
{
    nodeptr p,g,h,l;
//    int i;

    if (scanLevel==2) {
        p = (nodeptr)roottable[ifile];
        for (g = More(p); g != Next(p); g = Next(g))
            for (l = More(g); l != Next(g); l = Next(l)) {
                    nodetablescanlev[ifile][inodelev] = l;
                    inodelev++;
            }
        return;
    }

    if (scanLevel==3) {
        p = (nodeptr)roottable[ifile];
        for (g = More(p); g != Next(p); g = Next(g))
            for (h = More(g); h != Next(g); h = Next(h))
                for (l = More(h); l != Next(h); l = Next(l)) {
                        nodetablescanlev[ifile][inodelev] = l;
                        inodelev++;
                }
        return;
    }

    if (lev == scanLevel-1) {
        if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Type(l)==CELL) {
                    nodetablescanlev[ifile][inodelev] = l;
                    inodelev++;
                } else {
                    nodetablescanlev[ifile][inodelev] = l;
                    inodelev++;
                    ibodyleftout++;
                }
            }
        } else {
            ibodyleftout++;
        }
    } else {
        if ( lev+1 <= scanLevel ) {
            if (Type(q)==CELL) {
                for (l = More(q); l != Next(q); l = Next(l)) {
                    if (Type(l)==CELL)
                        walktree_index_scan_lev(l, lev+1, ifile, scanLevel);
                    else
                        ibodyleftout++;
                }
            } else {
                nodetablescanlev[ifile][inodelev] = q;
                inodelev++;
                ibodyleftout++;
            }
        }
    }
}

local void walktree_index_scan_lev_root(struct cmdline_data* cmd, struct  global_data* gd,
                                        nodeptr q, int lev, int ifile)
{
    nodeptr p,g,h,l;
//    int i;

    if (cmd->scanLevelRoot==2) {
        p = (nodeptr)roottable[ifile];
        for (g = More(p); g != Next(p); g = Next(g))
            for (l = More(g); l != Next(g); l = Next(l)) {
                nodetablescanlev_root[ifile][inodelev_root] = l;
                    inodelev_root++;
            }
        return;
    }

    if (cmd->scanLevelRoot==3) {
        p = (nodeptr)roottable[ifile];
        for (g = More(p); g != Next(p); g = Next(g))
            for (h = More(g); h != Next(g); h = Next(h))
                for (l = More(h); l != Next(h); l = Next(l)) {
                    nodetablescanlev_root[ifile][inodelev_root] = l;
                        inodelev_root++;
                }
        return;
    }

    if (lev == cmd->scanLevelRoot-1) {
        if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Type(l)==CELL) {
                    nodetablescanlev_root[ifile][inodelev_root] = l;
                    inodelev_root++;
                } else {
                    nodetablescanlev_root[ifile][inodelev_root] = l;
                    inodelev_root++;
                    ibodyleftout_root++;
                }
            }
        } else {
            ibodyleftout_root++;
        }
    } else {
        if ( lev+1 <= cmd->scanLevelRoot ) {
            if (Type(q)==CELL) {
                for (l = More(q); l != Next(q); l = Next(l)) {
                    if (Type(l)==CELL)
                        walktree_index_scan_lev_root(cmd, gd, l, lev+1, ifile);
                    else
                        ibodyleftout_root++;
                }
            } else {
                nodetablescanlev_root[ifile][inodelev_root] = q;
                inodelev_root++;
                ibodyleftout_root++;
            }
        }
    }
}
#endif // ! BALLS

//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
local void walktree_hit(struct  cmdline_data* cmd, struct  global_data* gd,
                        nodeptr q, real qsize)
{
    nodeptr l;
//    int i;

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

local int save_nodes(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    if (cmd->verbose==VERBOSEDEBUGINFO) {
        sprintf(gd->nodesfilePath,"%s/nodes%s.txt",gd->tmpDir,cmd->suffixOutFiles);
        if(!(gd->outnodelev=fopen(gd->nodesfilePath, "w")))
            error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
    }
//B
    if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") ) {
        if (cmd->verbose==3) {
            sprintf(gd->bodiesfilePath,"%s/bodies%s.txt",
                    gd->tmpDir,cmd->suffixOutFiles);
            if(!(gd->outbodylev=fopen(gd->bodiesfilePath, "w")))
                error("\nstart_Common: error opening file '%s' \n",
                      gd->nodesfilePath);
        }
    }
    int in;
//E
    INTEGER nodescount=0;
    bodyptr pn;
    INTEGER nodescount_smooth=0;
    INTEGER sumbodies=0, sumcells=0;
    INTEGER nodescount_thread=0, nodescount_thread_total=0;

    if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") )
        nodetable = (bodyptr) allocate(gd->nnodescanlevTable[ifile] * sizeof(body));

    for (in=0; in<inodelev; in++) {
//B
        if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") ) {
            pn = nodetable + in;
            Id(pn) = in;
            Type(pn) = BODY;
            Update(pn) = TRUE;
            Mass(pn) = Mass(nodetablescanlev[ifile][in]);
            Weight(pn) = Weight(nodetablescanlev[ifile][in]);
            Kappa(pn) = Kappa(nodetablescanlev[ifile][in]);
            SETV(Pos(pn),Pos(nodetablescanlev[ifile][in]));
            Nb(pn) = Nb(nodetablescanlev[ifile][in]);
            Size(pn) = 0.;
            Radius(pn) = 0.;
            if (cmd->verbose==VERBOSEDEBUGINFO) {
                out_int_mar(gd->outbodylev, Id(pn));
                out_int_mar(gd->outbodylev, Type(pn));
                out_vector_mar(gd->outbodylev, Pos(pn));
                out_real_mar(gd->outbodylev, Kappa(pn));
                out_real_mar(gd->outbodylev, Radius(pn));
                out_real_mar(gd->outbodylev, Mass(pn));
                out_int_long(gd->outbodylev, Nb(pn));
            }
        }
//E
        if (Type(nodetablescanlev[ifile][in]) == CELL) {
            if (cmd->verbose==VERBOSEDEBUGINFO) {
                out_int_mar(gd->outnodelev, IDXSCAN(nodetablescanlev[ifile][in]));
                out_int_mar(gd->outnodelev, Type(nodetablescanlev[ifile][in]));
                out_vector_mar(gd->outnodelev, Pos(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Kappa(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Size(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Mass(nodetablescanlev[ifile][in]));
                out_int_long(gd->outnodelev, Nb(nodetablescanlev[ifile][in]));
            }
            nodescount += Nb(nodetablescanlev[ifile][in]);
//B Smooth(ing) section
            if (Nb(nodetablescanlev[ifile][in]) <= gd->nsmooth[0]) nodescount_smooth++;
//E
            sumbodies += Nb(nodetablescanlev[ifile][in]);
            sumcells += 1;
            nodescount_thread += Nb(nodetablescanlev[ifile][in]);
        } else {
            if (cmd->verbose==VERBOSEDEBUGINFO) {
                out_int_mar(gd->outnodelev, IDXSCAN(nodetablescanlev[ifile][in]));
                out_int_mar(gd->outnodelev, Type(nodetablescanlev[ifile][in]));
                out_vector_mar(gd->outnodelev, Pos(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Kappa(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Radius(nodetablescanlev[ifile][in]));
                out_real_mar(gd->outnodelev, Mass(nodetablescanlev[ifile][in]));
                out_int_long(gd->outnodelev, Nb(nodetablescanlev[ifile][in]));
            }
            sumbodies += Nb(nodetablescanlev[ifile][in]);
            nodescount_thread += Nb(nodetablescanlev[ifile][in]);
        }

#ifdef OPENMPCODE
        if (cmd->verbose_log>=VERBOSEDEBUGINFO)
        if (in%cmd->numthreads == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog, " -(%ld) Chunk: %ld\n", in, nodescount_thread);
            nodescount_thread_total += nodescount_thread;
            nodescount_thread = 0;
        }
#else
        int numthreads=1;
        if (cmd->verbose_log>=VERBOSEDEBUGINFO)
        if (in%numthreads == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog, " -(%ld) Chunk: %ld\n", in, nodescount_thread);
            nodescount_thread_total += nodescount_thread;
            nodescount_thread = 0;
        }
#endif
    } // ! loop in

    if (cmd->verbose_log>=VERBOSEDEBUGINFO)
    verb_log_print(cmd->verbose_log, gd->outlog, " -Total Chunk: %ld\n", nodescount_thread_total);

    if (cmd->verbose >= VERBOSENORMALINFO) {
        verb_print(cmd->verbose,
                   "Found %ld particles in %ld nodes... vs number of total cells %ld\n",
                   nodescount, gd->nnodescanlevTable[ifile], gd->ncellTable[ifile]);
        verb_print(cmd->verbose,
                   "...%ld cells at scan level...\n",
                   sumcells);
        verb_print(cmd->verbose,
                   "...and %ld bodies. Bodies in upper levels: %ld\n",
                   sumbodies, cmd->nbody-sumbodies);
        verb_print(cmd->verbose,
                   "%ld particles were left out of cells at scan level.\n",
                   ibodyleftout);
        //B Smooth(ing) section
        verb_print(cmd->verbose,
                   "%ld cells were with at much %d particles in them.\n",
                   nodescount_smooth,gd->nsmooth[0]);
        //E
        verb_print(cmd->verbose,
                   "Checking sums (bodyleftout+nodescount): %ld.\n",
                   ibodyleftout+nodescount);
    }

    if (cmd->verbose==VERBOSEDEBUGINFO)
        fclose(gd->outnodelev);
//B
    if (!gd->flagSmoothCellMin && scanopt(cmd->options, "smooth-min-cell") )
        if (cmd->verbose==VERBOSEDEBUGINFO)
            fclose(gd->outbodylev);
//E
    return SUCCESS;
}

local int save_nodes_root(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    sprintf(gd->nodesfilePath,"%s/nodes_root%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(gd->outnodelev=fopen(gd->nodesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
    sprintf(gd->bodiesfilePath,"%s/bodies_root%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if(!(gd->outbodylev=fopen(gd->bodiesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
    int in;
    INTEGER nodescount=0;
    bodyptr pn;
    INTEGER nodescount_smooth=0;
    INTEGER sumbodies=0, sumcells=0;
    INTEGER nodescount_thread=0, nodescount_thread_total=0;

    nodetable_root =
        (bodyptr) allocate(gd->nnodescanlev_rootTable[ifile] * sizeof(body));

    for (in=0; in<inodelev_root; in++) {
        pn = nodetable_root + in;
        Id(pn) = in;
        Type(pn) = BODY;
        Update(pn) = TRUE;
        Mass(pn) = Mass(nodetablescanlev_root[ifile][in]);
        Weight(pn) = Weight(nodetablescanlev_root[ifile][in]);
        Kappa(pn) = Kappa(nodetablescanlev_root[ifile][in]);
        SETV(Pos(pn),Pos(nodetablescanlev_root[ifile][in]));

        out_int_mar(gd->outbodylev, Id(pn));
        out_int_mar(gd->outbodylev, Type(pn));
        out_vector_mar(gd->outbodylev, Pos(pn));
        out_real_mar(gd->outbodylev, Kappa(pn));
        out_real_mar(gd->outbodylev, Radius(pn));
        out_real_mar(gd->outbodylev, Mass(pn));
        out_int_long(gd->outbodylev, Nb(pn));

        if (Type(nodetablescanlev_root[ifile][in]) == CELL) {
            out_int_mar(gd->outnodelev, IDXSCAN(nodetablescanlev_root[ifile][in]));
            out_int_mar(gd->outnodelev, Type(nodetablescanlev_root[ifile][in]));
            out_vector_mar(gd->outnodelev, Pos(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Kappa(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Size(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Mass(nodetablescanlev_root[ifile][in]));
            out_int_long(gd->outnodelev, Nb(nodetablescanlev_root[ifile][in]));
            nodescount += Nb(nodetablescanlev_root[ifile][in]);
//B Smooth(ing) section
            if (Nb(nodetablescanlev_root[ifile][in]) <= gd->nsmooth[0]) nodescount_smooth++;
//E
            sumbodies += Nb(nodetablescanlev_root[ifile][in]);
            sumcells += 1;
            nodescount_thread += Nb(nodetablescanlev_root[ifile][in]);
        } else {
            out_int_mar(gd->outnodelev, IDXSCAN(nodetablescanlev_root[ifile][in]));
            out_int_mar(gd->outnodelev, Type(nodetablescanlev_root[ifile][in]));
            out_vector_mar(gd->outnodelev, Pos(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Kappa(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Radius(nodetablescanlev_root[ifile][in]));
            out_real_mar(gd->outnodelev, Mass(nodetablescanlev_root[ifile][in]));
            out_int_long(gd->outnodelev, Nb(nodetablescanlev_root[ifile][in]));
            sumbodies += Nb(nodetablescanlev_root[ifile][in]);
            nodescount_thread += Nb(nodetablescanlev_root[ifile][in]);
        }

#ifdef OPENMPCODE
        if (cmd->verbose_log>=3)
        if (in%cmd->numthreads == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog, " -(%ld) Chunk: %ld\n", in, nodescount_thread);
            nodescount_thread_total += nodescount_thread;
            nodescount_thread = 0;
        }
#else
        int numthreads=1;
        if (cmd->verbose_log>=3)
        if (in%numthreads == 0) {
            verb_log_print(cmd->verbose_log, gd->outlog, " -(%ld) Chunk: %ld\n", in, nodescount_thread);
            nodescount_thread_total += nodescount_thread;
            nodescount_thread = 0;
        }

#endif
    } // ! loop in

    if (cmd->verbose_log==3)
    verb_log_print(cmd->verbose_log, gd->outlog, " -Total Chunk: %ld\n", nodescount_thread_total);

    verb_print(cmd->verbose,
        "Found %ld particles in %ld nodes... vs number of total cells %ld\n",
        nodescount, gd->nnodescanlevTable[ifile], gd->ncellTable[ifile]);
    verb_print(cmd->verbose,
        "...%ld cells at scan level...\n",
        sumcells);
    verb_print(cmd->verbose,
        "...and %ld bodies. Bodies in upper levels: %ld\n",
        sumbodies, cmd->nbody-sumbodies);
    verb_print(cmd->verbose, "%ld particles were left out of cells at scan level.\n",ibodyleftout);
//B Smooth(ing) section
    verb_print(cmd->verbose, "%ld cells were with at much %d particles in them.\n",
               nodescount_smooth,gd->nsmooth[0]);
//E
    verb_print(cmd->verbose, "Checking sums (bodyleftout+nodescount): %ld.\n",
               ibodyleftout+nodescount);
    fclose(gd->outnodelev);
    fclose(gd->outbodylev);
    free(nodetable_root);

    return SUCCESS;
}

//E
