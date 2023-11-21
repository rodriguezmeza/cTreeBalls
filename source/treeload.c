/* ==============================================================================
!	MODULE: treeload.c			[cTreeBalls]									!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:	april 2023                                                  !
!	Purpose: 3-point correlation function computation       					!
!	Language: C																	!
!	Use: maketree(btab, nbody)				    								!
!	Major revisions:															!
!==============================================================================*/
//        1          2          3          4          5          6          7

// Work to do in order to use with boxes not centered at (0,0,...)

#include "globaldefs.h"

local void newtree(void);
local cellptr makecell(void);
local void expandbox(bodyptr, int);
local void loadbody(bodyptr);
local int subindex(bodyptr, cellptr);
local void hackCellProp(cellptr, real, int);
local void setRadius(cellptr, vector, real);
local void threadtree(nodeptr, nodeptr);

local void findRootCenter(bodyptr, int);
local void centerBodies(bodyptr btab, int nbody);

#define MAXLEVEL  32

local int cellhist[MAXLEVEL];
local int subnhist[MAXLEVEL];
//B Smooth(ing) section
//B To smooth bodies
local INTEGER NTOT[1];                      // Two sets of cells to smooth
local INTEGER ip;
#define NbMax 33
local int cellhistNb[NbMax];
local int cellRadius[NbMax];
local real deltaRadius;
//E
//E

//B To see the bodies belonging to a cell:
local void walktree(nodeptr q, real qsize);
//E


//#ifdef BALLS
//B To debug cells:
local INTEGER icell;
local INTEGER inode;
local INTEGER Nc1;
local INTEGER Nc2;
local void walktree_index_scan_lev(nodeptr q, int lev);
#ifdef DEBUG
local void walktree_hit(nodeptr q, real qsize);
#endif
local INTEGER inodelev;
local INTEGER ibodyleftout;
// Root nodes:
local void walktree_index_scan_lev_root(nodeptr q, int lev);
local INTEGER inodelev_root;
local INTEGER ibodyleftout_root;
//local int scanLevel_root;
//E

//B Balls-correction. 2023-11-10
local int save_nodes(void);
//E
//#endif

global void maketree(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = CPUTIME;
    gd.bytes_tot_cells = 0;

//#ifdef BALLS
#ifdef DEBUG
//B To debug cells:
    sprintf(gd.cellsfilePath,"%s/cells%s.txt",gd.tmpDir,cmd.suffixOutFiles);
    if(!(gd.outcells=fopen(gd.cellsfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd.cellsfilePath);
//E
//#endif
#endif

    newtree();
    root = makecell();
//B Set (0,0,...) as the center of the box
// By now it is only working with boxes centered at (0,0,...)
    findRootCenter(btab, nbody);
    centerBodies(btab, nbody);
    findRootCenter(btab, nbody);
//E

    CLRV(Pos(root));
    expandbox(btab, nbody);
    DO_BODY(p, btab, btab+nbody) {
        Nbb(p) = 1;                         // Check consistency with smoothing... Correction
        loadbody(p);
//B Balls-correction. 2023-11-10
        Nb(p) = 1;
        Radius(p) = 0.0;
//E
    }
    gd.tdepth = 0;

    for (i = 0; i < MAXLEVEL; i++)
        cellhist[i] = subnhist[i] = 0;
    for (i = 0; i < NbMax; i++)
        cellhistNb[i] = cellRadius[i] = 0;
    deltaRadius = gd.rSize/NbMax;
    verb_log_print(cmd.verbose_log,gd.outlog,"\ndeltaRadius = %g\n",deltaRadius);

//B To see the bodies belonging to a cell:
    DO_BODY(p,btab,btab+nbody)
        Selected(p) = FALSE;
//E

//B Smooth(ing) section
    NTOT[0] = 0;
//E
    ip = 0;
    hackCellProp(root, gd.rSize, 0);

//B To bin cell's radius... Check!!! (bins logscale)
    real rBin;
    verb_log_print(cmd.verbose_log,gd.outlog,"\nmaketree: radius histogram:\n");
    for (i = 0; i < NbMax; i++) {
        rBin = ((int)i)*gd.rSize/deltaRadius;
        verb_log_print(cmd.verbose_log,gd.outlog,"%g %d\n", rBin, cellRadius[i]);
    }
    verb_log_print(cmd.verbose_log,gd.outlog,"\n");
//E

//B Smooth(ing) section
//B To smooth bodies
    if (!gd.flagSmooth)
        if (scanopt(cmd.options, "smooth")) {
            printf("NTOT = %ld \t%ld \t%ld\n",
               gd.nsmooth[0], NTOT[0], gd.nsmooth[0]*NTOT[0]);
            gd.nbodysm = NTOT[0];
            bodytabsm = (bodyptr) allocate((NTOT[0]) * sizeof(body));
            gd.bytes_tot += (NTOT[0])*sizeof(body);
            verb_print(cmd.verbose,
                "Allocated %g MByte for (smooth %ld) particle storage.\n",
                (NTOT[0])*sizeof(body)*INMB, (NTOT[0]));
        }
//E
//E

//#ifdef BALLS
//B To debug cells:
    icell = 0;
    inode = 0;
    gd.nnode = gd.ncell/cmd.stepNodes;
    nodetab = (nodeptr *) allocate(gd.nnode * sizeof(nodeptr));
    gd.bytes_tot += gd.nnode*sizeof(nodeptr);
    verb_print(cmd.verbose,
               "\nAllocated %g MByte for (%d cells) nodetab storage.\n",
               gd.nnode*sizeof(nodeptr)*INMB,gd.nnode);
    Nc1 = gd.ncritical[0];
    Nc2 = gd.ncritical[1];
//E
//#endif

    threadtree((nodeptr) root, NULL);
    verb_print(cmd.verbose, "threadtree: number ip of selected cells = %ld\n",ip);
//#ifdef BALLS
    verb_print(cmd.verbose, "%d real node (range of nodes to search: >nc1 && <nc2)\n",inode);
    gd.rnnode = inode;
//#endif

//B Smooth(ing) section
//B To see the bodies belonging to a cell:
    walktree((nodeptr) root, gd.rSize);
    int isel=0, inosel=0;
    DO_BODY(p,btab,btab+nbody) {
        if (Selected(p))
            isel++;
        else
            inosel++;
    }
    verb_print(cmd.verbose,
               "\nSelected vs NotSelected and total: %ld %ld %ld\n\n",
               isel, inosel, isel + inosel);
//E

#ifdef BALLS
#ifndef TREENODEALLBODIES
//B BALLS :: SCANLEV
    if (scanopt(cmd.options, "set-default-param")) {
        verb_print(cmd.verbose, "\tfixing scanLevel to tdepth-1...\n");
        cmd.scanLevel = MAX(gd.tdepth-1,3);
        verb_print(cmd.verbose, "\tfinal value is %d.\n", cmd.scanLevel);
    } else {
    if (cmd.scanLevel > gd.tdepth) {
        verb_print(cmd.verbose,
                   "Warning! tree depth (%d) is less than scanLevel (%d)...\n",
                   gd.tdepth, cmd.scanLevel);
        verb_print(cmd.verbose, "\tfixing to tdepth-1...\n");
        cmd.scanLevel = MAX(gd.tdepth-1,3);
        verb_print(cmd.verbose, "\tfinal value is %d.\n", cmd.scanLevel);
    }
    }
#ifdef LOGHIST
    verb_print(cmd.verbose,
               "(Only in log-scale) deltaR is %g and root size at scanLevel is %g.\n",
               gd.deltaR, gd.rSize/rpow(2.0,cmd.scanLevel));
    i = 0;
    while (gd.deltaR < gd.rSize/rpow(2.0,i)) i++;
    verb_print(cmd.verbose,
               "\t\t\tSuggested scanLevel is %d, where root size is %g.\n\n",
               i, gd.rSize/rpow(2.0,i));
#endif
    inodelev = 0;
    ibodyleftout = 0;
    if (cmd.scanLevel==0) {
        gd.nnodescanlev = 0;
        gd.bytes_tot += gd.nnodescanlev*sizeof(nodeptr);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for (%d) scan nodetab storage.\n",
                   INMB*gd.nnodescanlev*sizeof(nodeptr),gd.nnodescanlev);
        gd.nnodescanlev = 1;
    } else {
/*
#if NDIM == 3
    gd.nnodescanlev = rpow(8, cmd.scanLevel);
#else
    gd.nnodescanlev = rpow(4, cmd.scanLevel);
#endif
*/
// Use the node storage for cells. Bodies are left out.
//    gd.nnodescanlev = gd.ncell;
    gd.nnodescanlev = 2.0*gd.ncell;
    nodetabscanlev = (nodeptr *) allocate(gd.nnodescanlev * sizeof(nodeptr));
//
    gd.bytes_tot += gd.nnodescanlev*sizeof(nodeptr);
//    verb_print(cmd.verbose,
//               "Allocated %g MByte for (%d) scan nodetab storage.\n",
//               gd.bytes_tot*INMB,gd.nnodescanlev);
        verb_print(cmd.verbose,
               "Allocated %g MByte for (%d) scan nodetab storage.\n",
               INMB*gd.nnodescanlev*sizeof(nodeptr),gd.nnodescanlev);

//    if (!cmd.scanLevel==0)                  // 0 <= lev <= cmd.scanLevel
        walktree_index_scan_lev((nodeptr)root, 0);
    verb_print(cmd.verbose,
               "Found %d nodes to scan at level %d.\n",
               inodelev, cmd.scanLevel);
    gd.nnodescanlev = inodelev;
    nodetable = (bodyptr) allocate(gd.nnodescanlev * sizeof(body));
    save_nodes();
    }

#endif // ! TREENODEALLBODIES

//B Root nodes to scan:
    inodelev_root = 0;
    ibodyleftout_root = 0;
    if (cmd.scanLevelRoot==0) {
        gd.nnodescanlev_root = 0;
        gd.bytes_tot += gd.nnodescanlev_root*sizeof(nodeptr);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for (%d) scan root nodetab storage.\n",
                   INMB*gd.nnodescanlev_root*sizeof(nodeptr),gd.nnodescanlev_root);
    } else {
        gd.nnodescanlev_root = 2.0*gd.ncell;
        nodetabscanlev_root = (nodeptr *) allocate(gd.nnodescanlev_root * sizeof(nodeptr));
        gd.bytes_tot += gd.nnodescanlev_root*sizeof(nodeptr);
        verb_print(cmd.verbose,
               "Allocated %g MByte for (%d) scan root nodetab storage.\n",
               INMB*gd.nnodescanlev_root*sizeof(nodeptr),gd.nnodescanlev_root);
        walktree_index_scan_lev_root((nodeptr)root, 0); // 0 <= lev <= cmd.scanLevelRoot
        verb_print(cmd.verbose,
           "Found %d root nodes to scan at level %d.\n",
            inodelev_root, cmd.scanLevelRoot);
        verb_print(cmd.verbose, "%ld particles were left out at root level.\n",ibodyleftout_root);
        gd.nnodescanlev_root = inodelev_root;
    }
    
#ifdef BUCKET
    gd.Rcell[0] = gd.rSize;
    for (i = 1; i < gd.tdepth; i++)
        gd.Rcell[i] = gd.Rcell[i-1]/2;
    verb_print(cmd.verbose, "Maximum and minimum cell size: %e %e\n",
               gd.Rcell[0],gd.Rcell[gd.tdepth-1]);
    if (cmd.scanLevelMin == 0)
        gd.rminCell = 0.;
    else
        gd.rminCell = gd.Rcell[gd.tdepth-1+cmd.scanLevelMin+1];
    verb_print(cmd.verbose, "Cell size at scanLevelMin (%d): %e\n",
               cmd.scanLevelMin, gd.rminCell);
#endif
//E Root nodes

//E BALLS :: SCANLEV
#endif // ! BALLS


//B Smooth(ing) section

    if (!gd.flagSmooth || !gd.flagSetNbNoSel) {
    if ( (scanopt(cmd.options, "smooth") && scanopt(cmd.options, "set-Nb-noSel")) ) {
        printf("NTOT = %ld \t%ld \t%ld\n\n",
               gd.nsmooth[0], NTOT[0]+inosel,
               gd.nsmooth[0]*NTOT[0]+ inosel);
        bodytabSel = (bodyptr) allocate((NTOT[0]+inosel) * sizeof(body));
        gd.bytes_tot += (NTOT[0]+inosel)*sizeof(body);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for (smooth) particle (%ld) storage.\n",
                   (NTOT[0]+inosel)*sizeof(body)*INMB, (NTOT[0]+inosel));
        int ipcount=0;
        bodyptr q = bodytabSel;
        DO_BODY(p,bodytabsm,bodytabsm+gd.nbodysm) {
            ipcount++;
            Id(q) = ipcount;
// BODY3
            Type(q) = BODY3;
            Nbb(q) = Nbb(p);
            Weight(q) = Weight(p);
            SETV(Pos(q), Pos(p));
            Kappa(q) = Kappa(p);
            q++;
        }
        verb_print(cmd.verbose,"Added %ld smoothed cells...\n",ipcount);
        DO_BODY(p,btab,btab+nbody) {
            if (!Selected(p)) {
            ipcount++;
            Id(q) = ipcount;
            Type(q) = BODY;
// BODY3
            Nbb(q) = 1;
            Weight(q) = Weight(p);
            SETV(Pos(q), Pos(p));
            Kappa(q) = Kappa(p);
            q++;
            }
        }
        verb_print(cmd.verbose,"Added %ld total bodies...\n",ipcount);
        gd.nbodySel = ipcount;
    }
    }
//E

//B To smooth bodies
    verb_log_print(cmd.verbose_log,gd.outlog,"\nmaketree: Nb histogram:\n");
    for (i = 0; i < NbMax; i++)
        verb_log_print(cmd.verbose_log,gd.outlog,"%d %ld\n", i, cellhistNb[i]);
    verb_log_print(cmd.verbose_log,gd.outlog,"\n");
//E
    gd.bytes_tot += gd.bytes_tot_cells;
    verb_print(cmd.verbose,
               "Allocated %g MByte for (%d) cells storage.\n",
                       gd.bytes_tot_cells*INMB, gd.ncell);
    verb_print(cmd.verbose, "\nmaketree : root number of bodies = %ld\n", Nb(root));
    
#ifdef DEBUG
//B To debug cells:
    fclose(gd.outcells);
//E
#endif

//E Smooth(ing) section


//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
    DO_BODY(p,btab,btab+nbody)
        HIT(p) = FALSE;
//    walktree_hit((nodeptr) root, Size(root));
#endif
//E

    gd.cputree = CPUTIME - cpustart;
    verb_print(cmd.verbose, "\nCPU tree time : %lf\n", gd.cputree);
}

local void findRootCenter(bodyptr btab, int nbody)
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
        verb_print(cmd.verbose, "findRootCenter: Pos(root) = %lf\n", Pos(root)[k]);
    }

    for(k=0, len=xmax[0]-xmin[0]; k<NDIM; k++)
        if((xmax[k]-xmin[k])>len)
            len=xmax[k]-xmin[k];

    verb_print(cmd.verbose, "findRootCenter: len = %lf\n", len);
}

local void centerBodies(bodyptr btab, int nbody)
{
    bodyptr p;
    int k;

    DO_BODY(p, btab, btab+nbody)
        DO_COORD(k)
            Pos(p)[k] = Pos(p)[k] - Pos(root)[k];
}

local void newtree(void)
{
    root = NULL;
    gd.ncell = 0;
}

local cellptr makecell(void)
{
    cellptr c;
    int i;
 
    c = (cellptr) allocate(sizeof(cell));
    Type(c) = CELL;
//B To smooth bodies
    Nb(c) = 0;
//E
    Update(c) = FALSE;
    for (i = 0; i < NSUB; i++)                  
        Subp(c)[i] = NULL;
    gd.ncell++;
    gd.bytes_tot_cells += sizeof(cell);
    return (c);
}

local void expandbox(bodyptr btab, int nbody)
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
    while (gd.rSize < 2 * dmax)
        gd.rSize = 2 * gd.rSize;

    verb_print(cmd.verbose, "treeload expandbox: rSize = %lf\n", gd.rSize);
}

local void loadbody(bodyptr p)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;

// Keep it in order to study MPI programming...
//    verb_print_debug(1, "\nloadbody:: Aqui voy (2)\n");
//    MPI_Barrier(MPI_COMM_WORLD);

    q = root;
    qind = subindex(p, q);
    Nb(q) += 1;                 // Smooth
    qsize = gd.rSize;
    while (Subp(q)[qind] != NULL) {
// BODY3
        if (Type(Subp(q)[qind]) == BODY || Type(Subp(q)[qind]) == BODY3) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0) {
                verb_log_print(cmd.verbose_log,gd.outlog,
                               "\nIds: %ld and %ld have the same position\n",
                               Id(p),Id(Subp(q)[qind]));
                DO_COORD(k)
                    verb_log_print(cmd.verbose_log,gd.outlog,"\nPos[k]: %le %le\n",
                               Pos(p)[k],Pos(Subp(q)[qind])[k]);
                error("loadbody: two bodies have same position\n");
            }
            c = makecell();                     
            Nb(c) += 1;                     // Smooth
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
            Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
            Nb(c) += 1;                     // Smooth
            Subp(q)[qind] = (nodeptr) c;
        } else {
            Nb(Subp(q)[qind]) += 1;         // Smooth
        }
        q = (cellptr) Subp(q)[qind];
        qind = subindex(p, q);
        qsize = qsize / 2;
    }
    Subp(q)[qind] = (nodeptr) p;
}

local int subindex(bodyptr p, cellptr q)
{
    int ind, k;
 
    ind = 0;                                    
	DO_COORD(k)
        if (Pos(q)[k] <= Pos(p)[k])             
            ind += NSUB >> (k + 1);             
    return (ind);
}

local void hackCellProp(cellptr p, real psize, int lev)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;
 
    gd.tdepth = MAX(gd.tdepth, lev);
    cellhist[lev]++;
    Level(p) = lev;                         // To set scanLevel
    Weight(p) = 0.0;
    Nb(p) = 0;
    Kappa(p) = 0.0;
    CLRV(cmpos);

    for (i = 0; i < NSUB; i++) {
        if ((q = Subp(p)[i]) != NULL) {
            subnhist[lev]++;                    
            if (Type(q) == CELL)
                hackCellProp((cellptr) q, psize/2, lev+1);
            Selected(p) |= Selected(q);     // To see the bodies belonging to a cell
            Update(p) |= Update(q);
            Weight(p) += Weight(q);
            if ( Type(q) == CELL) {
                Nb(p) += Nb(q);
                Kappa(p) += Weight(q)*Kappa(q);
            } else {
                if (Type(q) == BODY) {
                    Nb(p) += 1;
                    Kappa(p) += Weight(q)*Kappa(q);
                } else if (Type(q) == BODY3) {  // To set smoothing body
                    Nb(p) += 1;
                    Kappa(p) += Weight(q)*Kappa(q);
                }
            }
            MULVS(tmpv, Pos(q), Weight(q));
            ADDV(cmpos, cmpos, tmpv);           
        }
    }
//B Smooth(ing) section
    if (Nb(p)==gd.nsmooth[0]) {                 // Correct to <=
        NTOT[0]=NTOT[0]+1;
    }
//B Balls-correction. 2023-11-10
    if (scanopt(cmd.options, "center-of-mass")) {
        if (Weight(p) > 0.0) {
            DIVVS(cmpos, cmpos, Weight(p));
        } else {
            SETV(cmpos, Pos(p));
        }
    } else
        SETV(cmpos, Pos(p));
//E

#define EPSILON 1.0E-16
// Here there appears an error for big numbers of points such 201 millions...
// See line above and uncomment DIVVS(cmpos, cmpos, Weight(p)); line (not working!)
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 || Pos(p)[k] + psize/2 <= cmpos[k]) {
            if (psize/2 > 2.710505e-20 + EPSILON)
            error("hackCellProp: tree structure error: %d %le %le %le %le\n",
                  k, cmpos[k], Pos(p)[k] - psize/2, Pos(p)[k] + psize/2, psize/2);
            else {
                if (cmd.verbose_log>=3)
                verb_log_print(cmd.verbose_log,gd.outlog,
                    "hackCellProp: tree structure warning! psize/2 to small: %le \n", psize/2);
            }
        }
#undef EPSILON

    setRadius(p, cmpos, psize);
    SETV(Pos(p), cmpos);
    if (Nb(p)>0) {
        Kappa(p) /= Nb(p);
    } else
        error("hackCellProp: Nb = 0: %ld\n", Nb(p));
}

// Parameter theta controls size of the cell.
// theta from 0 to 1: 0 always open cells (BF, complexity N^2); 1 is the default value.
local void setRadius(cellptr p, vector cmpos, real psize)
{
    real bmax2, d;
    int k;

    if (cmd.theta == 0.0)
        Radius(p) = 2 * gd.rSize;
    else if (gd.sw94) {
        bmax2 = 0.0;
		DO_COORD(k) {
            d = cmpos[k] - Pos(p)[k] + psize/2; 
            bmax2 += rsqr(MAX(d, psize - d));   
        }
        Radius(p) = rsqrt(bmax2) / cmd.theta;
    } else if (gd.bh86)
        Radius(p) = psize / cmd.theta;
    else {
        Radius(p) = (psize/cmd.theta) * rsqrt((real)(NDIM))/2.0;
    }

    Size(p) = psize;

    int n;
    n = (int) (Radius(p) / deltaRadius);
    (cellRadius[n])++;
}

local void threadtree(nodeptr p, nodeptr n)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    bodyptr q;

    Next(p) = n;
    if (Type(p) == CELL) {

#ifdef DEBUG
//B To debug cells:
        icell++;
        if (icell%cmd.stepNodes == 0 && Nb(p)>= Nc1 && Nb(p)<=Nc2) {
            out_vector_mar(gd.outcells, Pos(p));
            out_real_mar(gd.outcells, Kappa(p));
            out_real_mar(gd.outcells, Size(p));
            out_int_long(gd.outcells, Nb(p));
            nodetab[inode] = p;
            inode++;
        }
//E
#else
        icell++;
        if (icell%cmd.stepNodes == 0 && Nb(p)>= Nc1 && Nb(p)<=Nc2) {
            nodetab[inode] = p;
            inode++;
        }
#endif

//B Smooth(ing) section
//B To smooth bodies
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        if (!gd.flagSmooth) {
            if (scanopt(cmd.options, "smooth")) {
                if (Nb(p)==gd.nsmooth[0]) {
                    ip++;
                    q = ip + bodytabsm -1;
                    Id(q) = ip;
                    Type(q) = BODY3;        // To set smoothing body
                    Nbb(q) = Nb(p);
                    Nbb(p) = Nb(p);
                    Weight(q) = Weight(p);
                    SETV(Pos(q), Pos(p));
                    Kappa(q) = Kappa(p);
//B To see the bodies belonging to a cell:
                    Selected(p) = TRUE;
                    Selected(q) = TRUE;
//E
                }
            }
        }
//E
//E
        ndesc = 0;
        for (i = 0; i < NSUB; i++)
            if (Subp(p)[i] != NULL)
                desc[ndesc++] = Subp(p)[i];
        More(p) = desc[0];
        desc[ndesc] = n;
        for (i = 0; i < ndesc; i++)
            threadtree(desc[i], desc[i+1]);
    }
}

local void walktree(nodeptr q, real qsize)
{
    nodeptr l;

    if (Selected(q)) {
        if (Type(q) == CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                Selected(l) = TRUE;
                walktree(l,qsize/2);
            }
        } else {
            Selected(q) = TRUE;
        }
    } else {
        if (Type(q) == CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                walktree(l,qsize/2);
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
local void walktree_index_scan_lev(nodeptr q, int lev)
{
    nodeptr p,g,h,l;
    int i;

    if (cmd.scanLevel==2) {
        p = (nodeptr)root;
        for (g = More(p); g != Next(p); g = Next(g))
            for (l = More(g); l != Next(g); l = Next(l)) {
                    nodetabscanlev[inodelev] = l;
                    inodelev++;
            }
        return;
    }

    if (cmd.scanLevel==3) {
        p = (nodeptr)root;
        for (g = More(p); g != Next(p); g = Next(g))
            for (h = More(g); h != Next(g); h = Next(h))
                for (l = More(h); l != Next(h); l = Next(l)) {
                        nodetabscanlev[inodelev] = l;
                        inodelev++;
                }
        return;
    }

// lev must be <= cmd.scanLevel
    if (lev == cmd.scanLevel-1) {
        if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Type(l)==CELL) {
                    nodetabscanlev[inodelev] = l;
                    inodelev++;
                } else {
                    nodetabscanlev[inodelev] = l;   // Correction
                    inodelev++;                     // Correction
                    ibodyleftout++;
                }
            }
        } else {
            ibodyleftout++;
        }
    } else {
//        if ( lev+1 != cmd.scanLevel )
        if ( lev+1 <= cmd.scanLevel )       // Correction
            if (Type(q)==CELL) {
                for (l = More(q); l != Next(q); l = Next(l)) {
                    if (Type(l)==CELL)
                        walktree_index_scan_lev(l, lev+1);
                    else
                        ibodyleftout++;
                }
            } else {
                nodetabscanlev[inodelev] = q;   // Correction
                inodelev++;                     // Correction
                ibodyleftout++;
            }
    }
}

local void walktree_index_scan_lev_root(nodeptr q, int lev)
{
    nodeptr p,g,h,l;
    int i;

    if (cmd.scanLevelRoot==2) {
        p = (nodeptr)root;
        for (g = More(p); g != Next(p); g = Next(g))
            for (l = More(g); l != Next(g); l = Next(l)) {
                    nodetabscanlev_root[inodelev_root] = l;
                    inodelev_root++;
            }
        return;
    }

    if (cmd.scanLevelRoot==3) {
        p = (nodeptr)root;
        for (g = More(p); g != Next(p); g = Next(g))
            for (h = More(g); h != Next(g); h = Next(h))
                for (l = More(h); l != Next(h); l = Next(l)) {
                        nodetabscanlev_root[inodelev_root] = l;
                        inodelev_root++;
                }
        return;
    }

// lev must be <= cmd.scanLevelRoot
    if (lev == cmd.scanLevelRoot-1) {
        if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Type(l)==CELL) {
                    nodetabscanlev_root[inodelev_root] = l;
                    inodelev_root++;
                } else {
                    nodetabscanlev_root[inodelev_root] = l;   // Correction
                    inodelev_root++;                     // Correction
                    ibodyleftout_root++;
                }
            }
        } else {
            ibodyleftout_root++;
        }
    } else {
        if ( lev+1 <= cmd.scanLevelRoot )       // Correction
            if (Type(q)==CELL) {
                for (l = More(q); l != Next(q); l = Next(l)) {
                    if (Type(l)==CELL)
                        walktree_index_scan_lev_root(l, lev+1);
                    else
                        ibodyleftout_root++;
                }
            } else {
                nodetabscanlev_root[inodelev_root] = q;   // Correction
                inodelev_root++;                     // Correction
                ibodyleftout_root++;
            }
    }
}
#endif // ! BALLS

//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
local void walktree_hit(nodeptr q, real qsize)
{
    nodeptr l;
    int i;

    if (Type(q) == CELL) {
        for (l = More(q); l != Next(q); l = Next(l)) {
            walktree_hit(l,qsize/2);
        }
    } else {
        HIT(q) = TRUE;
    }
}
#endif
//E

local int save_nodes(void)
{
    sprintf(gd.nodesfilePath,"%s/nodes%s.txt",gd.tmpDir,cmd.suffixOutFiles);
    if(!(gd.outnodelev=fopen(gd.nodesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd.nodesfilePath);
    sprintf(gd.bodiesfilePath,"%s/bodies%s.txt",gd.tmpDir,cmd.suffixOutFiles);
    if(!(gd.outbodylev=fopen(gd.bodiesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd.nodesfilePath);
    int in;
    INTEGER nodescount=0;
    bodyptr pn;
    INTEGER nodescount_smooth=0;
    INTEGER sumbodies=0, sumcells=0;

    for (in=0; in<inodelev; in++) {

        pn = nodetable + in;
        Id(pn) = in;
        Type(pn) = BODY;
        Update(pn) = TRUE;
        Weight(pn) = Weight(nodetabscanlev[in]);
        Kappa(pn) = Kappa(nodetabscanlev[in]);
        SETV(Pos(pn),Pos(nodetabscanlev[in]));

        out_int_mar(gd.outbodylev, Id(pn));
        out_int_mar(gd.outbodylev, Type(pn));
        out_vector_mar(gd.outbodylev, Pos(pn));
        out_real_mar(gd.outbodylev, Kappa(pn));
        out_real_mar(gd.outbodylev, Radius(pn));
        out_real_mar(gd.outbodylev, Weight(pn));
        out_int_long(gd.outbodylev, Nb(pn));
//E
        
        if (Type(nodetabscanlev[in]) == CELL) {
            out_int_mar(gd.outnodelev, IDXSCAN(nodetabscanlev[in]));
            out_int_mar(gd.outnodelev, Type(nodetabscanlev[in]));
            out_vector_mar(gd.outnodelev, Pos(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Kappa(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Size(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Weight(nodetabscanlev[in]));
            out_int_long(gd.outnodelev, Nb(nodetabscanlev[in]));
            nodescount += Nb(nodetabscanlev[in]);
//B Smooth(ing) section
            if (Nb(nodetabscanlev[in]) <= gd.nsmooth[0]) nodescount_smooth++;
//E
            sumbodies += Nb(nodetabscanlev[in]);
            sumcells += 1;
        } else {
            out_int_mar(gd.outnodelev, IDXSCAN(nodetabscanlev[in]));
            out_int_mar(gd.outnodelev, Type(nodetabscanlev[in]));
            out_vector_mar(gd.outnodelev, Pos(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Kappa(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Radius(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Weight(nodetabscanlev[in]));
            out_int_long(gd.outnodelev, Nb(nodetabscanlev[in]));
//            nodescount += Nb(nodetabscanlev[in]);          // Correction
            sumbodies += Nb(nodetabscanlev[in]);
        }
    }

    verb_print(cmd.verbose,
        "Found %ld particles in %ld nodes... vs number of total cells %ld\n",
        nodescount, gd.nnodescanlev, gd.ncell);
    verb_print(cmd.verbose,
        "...%ld cells at scan level...\n",
        sumcells);
    verb_print(cmd.verbose,
        "...and %ld bodies. Bodies in upper levels: %ld\n",
        sumbodies, cmd.nbody-sumbodies);
    verb_print(cmd.verbose, "%ld particles were left out of cells at scan level.\n",ibodyleftout);
//B Smooth(ing) section
    verb_print(cmd.verbose, "%ld cells were with at much %d particles in them.\n",
               nodescount_smooth,gd.nsmooth[0]);
//E
    verb_print(cmd.verbose, "Checking sums (bodyleftout+nodescount): %ld.\n",
               ibodyleftout+nodescount);
    fclose(gd.outnodelev);
    fclose(gd.outbodylev);

    return _SUCCESS_;
}
//E
