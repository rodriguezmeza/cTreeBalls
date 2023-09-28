/* ==============================================================================
!	MODULE: treeload.c			[tpcf]											!
!	Written by: M.A. Rodriguez-Meza.											!
!	Starting date:	april 2023                                                  !
!	Purpose: 3-point correlation function computation       					!
!	Language: C																	!
!	Use: maketree(btab, nbody)				    								!
!	Major revisions:															!
!==============================================================================*/

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
//B To smooth bodies
local INTEGER NTOT[3];          // Two sets of cells to smooth
local INTEGER ip;
#define NbMax 33
local int cellhistNb[NbMax];
local int cellRadius[NbMax];
local real deltaRadius;
//E

//B To see the bodies belonging to a cell:
//local cellptr celltabSel;
//local INTEGER icellSel;
//local bodyptr bodytabSel;
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

//    verb_print_debug(1, "\nAqui voy (2)\n");

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
	DO_BODY(p, btab, btab+nbody)
        loadbody(p);
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

//    verb_print_debug(1, "\nAqui voy (3)\n");

    NTOT[0] = 0; NTOT[1] = 0; NTOT[2] = 0;
    ip = 0;
    hackCellProp(root, gd.rSize, 0);

//B To bin cell's radius... Check!!! (bins logscale)
    real rBin;
    verb_log_print(cmd.verbose_log,gd.outlog,"\nmaketree: radius histogram:\n");
    for (i = 0; i < NbMax; i++) {
//        rBin = ((int)i)*gd.rSize/NbMax;
        rBin = ((int)i)*gd.rSize/deltaRadius;
        verb_log_print(cmd.verbose_log,gd.outlog,"%g %d\n", rBin, cellRadius[i]);
    }
    verb_log_print(cmd.verbose_log,gd.outlog,"\n");
//E

//B To smooth bodies
    if (!gd.flagSmooth) {
    if (scanopt(cmd.options, "smooth")) {
//        printf("NTOT = %ld \t%ld\n", NTOT[0]+NTOT[1]+NTOT[2],
//               gd.nsmooth[0]*NTOT[0]+gd.nsmooth[1]*NTOT[1]+gd.nsmooth[2]*NTOT[2]);
        printf("NTOT = %ld \t%ld \t%ld\n",
               gd.nsmooth[0], NTOT[0], gd.nsmooth[0]*NTOT[0]);
        gd.nbodysm = NTOT[0]+NTOT[1]+NTOT[2];
        bodytabsm = (bodyptr) allocate((NTOT[0]+NTOT[1]+NTOT[2]) * sizeof(body));
        gd.bytes_tot += (NTOT[0]+NTOT[1]+NTOT[2])*sizeof(body);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for (smooth %ld) particle storage.\n",
                   (NTOT[0]+NTOT[1]+NTOT[2])*sizeof(body)*INMB, (NTOT[0]+NTOT[1]+NTOT[2]));
    }
    }
//E
//#ifdef BALLS
//B To debug cells:
    icell = 0;
    inode = 0;
    gd.nnode = gd.ncell/cmd.stepNodes;
    nodetab = (nodeptr *) allocate(gd.nnode * sizeof(nodeptr));
    gd.bytes_tot += gd.nnode*sizeof(nodeptr);
    verb_print(cmd.verbose,
               "Allocated %g MByte for (%d) nodetab storage.\n",
               gd.bytes_tot*INMB,gd.nnode);
    Nc1 = gd.ncritical[0];
    Nc2 = gd.ncritical[1];
//E
//#endif

    threadtree((nodeptr) root, NULL);
    verb_print(cmd.verbose, "threadtree: number ip of selected cells = %ld\n",ip);
//#ifdef BALLS
    verb_print(cmd.verbose, "%d real node\n",inode);
    gd.rnnode = inode;
//#endif

//B To see the bodies belonging to a cell:
    walktree((nodeptr) root, gd.rSize);
    int isel=0, inosel=0;
    DO_BODY(p,btab,btab+nbody) {
//        verb_print(cmd.verbose, "%ld ", Selected(p));
        if (Selected(p))
            isel++;
        else
            inosel++;
    }
    verb_print(cmd.verbose, "\nSelected vs NotSelected and total: %ld %ld %ld\n\n",
               isel, inosel, isel + inosel);

//B BALLS :: SCANLEV
// Check segmentation fault at scanLevel>=9!!
// Suggest at cmdline (StartRun: CheckParameters) that scanLevel<9
    if (cmd.scanLevel > gd.tdepth) {
        verb_print(cmd.verbose, "Warning! tree depth (%d) is less than scanLevel (%d)...\n",
                   gd.tdepth, cmd.scanLevel);
        verb_print(cmd.verbose, "\tfixing to tdepth-2...\n");
        cmd.scanLevel = MAX(gd.tdepth-2,3);
        verb_print(cmd.verbose, "\tfinal value is %d.\n", cmd.scanLevel);
    }
    inodelev = 0;
    ibodyleftout = 0;
#if NDIM == 3
    gd.nnodescanlev = rpow(8, cmd.scanLevel);
#else
    gd.nnodescanlev = rpow(4, cmd.scanLevel);
#endif
// Use the node storage use for cells. Bodies are left out. 
//    nodetabscanlev = (nodeptr *) allocate(gd.nnodescanlev * sizeof(nodeptr));
    gd.nnodescanlev = gd.ncell;
    nodetabscanlev = (nodeptr *) allocate(gd.nnodescanlev * sizeof(nodeptr));
//
    gd.bytes_tot += gd.nnodescanlev*sizeof(nodeptr);
    verb_print(cmd.verbose,
               "Allocated %g MByte for (%d) scan nodetab storage.\n",
               gd.bytes_tot*INMB,gd.nnodescanlev);

    if (!cmd.scanLevel==0)
        walktree_index_scan_lev((nodeptr)root, 0); // 0 <= lev <= cmd.scanLevel
    verb_print(cmd.verbose, "Found %d nodes to scan at level %d.\n",inodelev, cmd.scanLevel);

    sprintf(gd.nodesfilePath,"%s/nodes%s.txt",gd.tmpDir,cmd.suffixOutFiles);
    if(!(gd.outnodelev=fopen(gd.nodesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd.nodesfilePath);
    int in;
    INTEGER nodescount=0;
    gd.nnodescanlev = inodelev;

    nodetable = (nodeptr) allocate(gd.nnodescanlev * sizeof(node));
    nodeptr pn;

    INTEGER nodescount_smooth=0;

    for (in=0; in<inodelev; in++) {

        pn = nodetable + in;
        Type(pn) = Type(nodetabscanlev[in]);
        Update(pn) = Update(nodetabscanlev[in]);
        Weight(pn) = Weight(nodetabscanlev[in]);
        Kappa(pn) = Kappa(nodetabscanlev[in]);
        SETV(Pos(pn),Pos(nodetabscanlev[in]));
        Next(pn) = Next(nodetabscanlev[in]);
        More(pn) = More(nodetabscanlev[in]);
        Nb(pn) = Nb(nodetabscanlev[in]);
        Radius(pn) = Radius(nodetabscanlev[in]);
        Size(pn) = Size(nodetabscanlev[in]);

        if (Type(nodetabscanlev[in]) == CELL) {
            out_int_mar(gd.outnodelev, IDXSCAN(nodetabscanlev[in]));
            out_int_mar(gd.outnodelev, Type(nodetabscanlev[in]));
            out_vector_mar(gd.outnodelev, Pos(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Kappa(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Size(nodetabscanlev[in]));
            out_int_long(gd.outnodelev, Nb(nodetabscanlev[in]));
            nodescount += Nb(nodetabscanlev[in]);
            if (Nb(nodetabscanlev[in]) <= gd.nsmooth[0]) nodescount_smooth++;
        } else {
            out_int_mar(gd.outnodelev, IDXSCAN(nodetabscanlev[in]));
            out_int_mar(gd.outnodelev, Type(nodetabscanlev[in]));
            out_vector_mar(gd.outnodelev, Pos(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Kappa(nodetabscanlev[in]));
            out_real_mar(gd.outnodelev, Smooth(nodetabscanlev[in]));
            out_int_long(gd.outnodelev, Nbb(nodetabscanlev[in]));
            nodescount += Nb(nodetabscanlev[in]);
        }
    }

    verb_print(cmd.verbose, "Found %ld particles in %d nodes... vs number of total cells %ld\n",
               nodescount, gd.nnodescanlev, gd.ncell);
    verb_print(cmd.verbose, "%ld particles were left out.\n",ibodyleftout);
    verb_print(cmd.verbose, "%ld cells were with at much %d particles in them.\n",
               nodescount_smooth,gd.nsmooth[0]);
    fclose(gd.outnodelev);
    
//B Creating tree for nodes
    maketreenodes(nodetable, gd.nnodescanlev);
//E

//E

    if (!gd.flagSmooth || !gd.flagSetNbNoSel) {
    if ( (scanopt(cmd.options, "smooth") && scanopt(cmd.options, "set-Nb-noSel")) ) {
//    if ( scanopt(cmd.options, "smooth") ) {
//        printf("NTOT = %ld \t%ld\n\n", NTOT[0]+NTOT[1]+NTOT[2]+inosel,
//               gd.nsmooth[0]*NTOT[0]+gd.nsmooth[1]*NTOT[1]+gd.nsmooth[2]*NTOT[2] + inosel);
        printf("NTOT = %ld \t%ld \t%ld\n\n",
               gd.nsmooth[0], NTOT[0]+NTOT[1]+NTOT[2]+inosel,
               gd.nsmooth[0]*NTOT[0] + inosel);
//        gd.nbodysm = NTOT[0]+NTOT[1]+NTOT[2];
        bodytabSel = (bodyptr) allocate((NTOT[0]+NTOT[1]+NTOT[2]+inosel) * sizeof(body));
        gd.bytes_tot += (NTOT[0]+NTOT[1]+NTOT[2]+inosel)*sizeof(body);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for (smooth) particle (%ld) storage.\n",
                   (NTOT[0]+NTOT[1]+NTOT[2]+inosel)*sizeof(body)*INMB, (NTOT[0]+NTOT[1]+NTOT[2]+inosel));
        int ipcount=0;
        bodyptr q = bodytabSel;
        DO_BODY(p,bodytabsm,bodytabsm+gd.nbodysm) {
            ipcount++;
//            q = ipcount + bodytabsm -1;
            Id(q) = ipcount;
// BODY3
//            Type(q) = BODY;
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
//#ifdef BALLS
//B To debug cells:
    fclose(gd.outcells);
//E
//#endif
#endif

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
//    IdCell(c) += 1;
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

//    if (gd.flagSmooth)
//    verb_print_debug(1, "\nAqui voy (4)\n");

//    verb_print_debug(1, "\nloadbody:: Aqui voy (2)\n");
//    MPI_Barrier(MPI_COMM_WORLD);

    q = root;
    qind = subindex(p, q);
    Nb(q) += 1;                 // Smooth
    qsize = gd.rSize;
    while (Subp(q)[qind] != NULL) {
// BODY3
//        if (Type(Subp(q)[qind]) == BODY) {
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
            Nb(c) += 1;         // Smooth
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;

//#ifdef BALLS
//B To debug cells:
//            out_vector(gd.outcells, Pos(c));
//E
//#endif

            Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
            Nb(c) += 1;         // Smooth
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
//    int sumNb;
 
//    if (gd.flagSmooth)
//    verb_print_debug(1, "\nAqui voy (4)\n");

    gd.tdepth = MAX(gd.tdepth, lev);
    cellhist[lev]++;
// BALLS
    Level(p) = lev;
//
    Weight(p) = 0.0;
//B To smooth bodies :: Always default: smooth and set-Nb
//    if (scanopt(cmd.options, "smooth") || scanopt(cmd.options, "set-Nb")) {
        Nb(p) = 0;
        Kappa(p) = 0.0;
//    }
//E
    CLRV(cmpos);
    for (i = 0; i < NSUB; i++) {
        if ((q = Subp(p)[i]) != NULL) {
            subnhist[lev]++;                    
            if (Type(q) == CELL)
                hackCellProp((cellptr) q, psize/2, lev+1);

//B To see the bodies belonging to a cell:
            Selected(p) |= Selected(q);
//
            Update(p) |= Update(q);
            Weight(p) += Weight(q);
//B To smooth bodies :: Always default: smooth and set-Nb
//  if ( (scanopt(cmd.options, "smooth") || scanopt(cmd.options, "set-Nb")) && Type(q) == CELL) {
                if ( Type(q) == CELL) {
                Nb(p) += Nb(q);
//                    Kappa(p) += Kappa(q);
                    Kappa(p) += Weight(q)*Kappa(q);
                } else {
//  if ( (scanopt(cmd.options, "smooth") || scanopt(cmd.options, "set-Nb")) && Type(q) == BODY) {
// BODY3
                if (Type(q) == BODY) {
                Nb(p) += 1;
//                Kappa(p) += Kappa(q);
                    Kappa(p) += Weight(q)*Kappa(q);
                } else if (Type(q) == BODY3) {
                    Nb(p) += 1;
//                    Kappa(p) += Kappa(q);
                    Kappa(p) += Weight(q)*Kappa(q);
                }
                }
//E
            MULVS(tmpv, Pos(q), Weight(q));
            ADDV(cmpos, cmpos, tmpv);           
        }
    }
    if (Nb(p)==gd.nsmooth[0]) {
        NTOT[0]=NTOT[0]+1;
//        Selected(p) = TRUE; // Propagate Selected value to the fathers...
    }
//    if (Nb(p)==gd.nsmooth[1]) NTOT[1]++;
//    if (Nb(p)==gd.nsmooth[2]) NTOT[2]++;
    if (Weight(p) > 0.0) {                 // We should use geometric center always
//        DIVVS(cmpos, cmpos, Weight(p));   // Activate if you want to use com center
        SETV(cmpos, Pos(p));
    } else {
        SETV(cmpos, Pos(p));
    }

#define EPSILON 1.0E-16
// Here there appears an error for big numbers of points such 201 millions...
// See line above and uncomment DIVVS(cmpos, cmpos, Weight(p)); line (not working!)
	DO_COORD(k)
        if (cmpos[k] < Pos(p)[k] - psize/2 ||   Pos(p)[k] + psize/2 <= cmpos[k]) {
            if (psize/2 > 2.710505e-20 + EPSILON)
            error("hackCellProp: tree structure error: %d %le %le %le %le\n",
                  k, cmpos[k], Pos(p)[k] - psize/2, Pos(p)[k] + psize/2, psize/2);
            else {
                if (cmd.verbose_log>=3)
                verb_log_print(cmd.verbose_log,gd.outlog,
                               "hackCellProp: tree structure warning! psize/2 to small: %le \n",
                               psize/2);
//                verb_print(cmd.verbose,
//                           "hackCellProp: tree structure warning! psize/2 to small: %le \n",
//                           psize/2);
            }
        }
#undef EPSILON

    setRadius(p, cmpos, psize);
    SETV(Pos(p), cmpos);
//B To smooth bodies :: Always default: smooth and set-Nb
//    if (scanopt(cmd.options, "smooth") || scanopt(cmd.options, "set-Nb"))
    if (Nb(p)>0) {
        Kappa(p) /= Nb(p); // times nb(p)?
//        printf("%d \t%le \t%ld\n", Type(p),Kappa(p),Nb(p));
    } else
        error("hackCellProp: Nb = 0: %ld\n", Nb(p));
//E
}

// Parameter theta controls size of the cell.
// theta from 0 to 1: 0 always open cells (BF); 1 default value.
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

//#ifdef BALLS
//B To debug cells:
    Size(p) = psize;
//E
//#endif

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

//#ifdef BALLS
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
//#endif

//B To smooth bodies
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        if (!gd.flagSmooth) {
        if (scanopt(cmd.options, "smooth")) {
//            if (Nb(p)<NbMax)
//                cellhistNb[Nb(p)]++;
            if (Nb(p)==gd.nsmooth[0]) {
                ip++;
                q = ip + bodytabsm -1;
                Id(q) = ip;
// BODY3
//                Type(q) = BODY;
            Type(q) = BODY3;
                Nbb(q) = Nb(p);
                Nbb(p) = Nb(p);
//                verb_print_debug(1, "\nAqui voy (0):: %ld %d\n", Type(q), Nbb(q));
                Weight(q) = Weight(p);
                SETV(Pos(q), Pos(p));
                Kappa(q) = Kappa(p);
                
//B To see the bodies belonging to a cell:
//                celltabSel[icellSel] = p;
//                icellSel++;
                Selected(p) = TRUE;
                Selected(q) = TRUE;
//E

            }
/*
            if (Nb(p)==gd.nsmooth[1]) {
                ip++;
                q = ip + bodytabsm -1;
                Id(q) = ip;
// BODY3
                Type(q) = BODY;
//            Type(q) = BODY3;
//            Nbb(q) = Nb(p);
                Weight(q) = Weight(p);
                SETV(Pos(q), Pos(p));
                Kappa(q) = Kappa(p);

//B To see the bodies belonging to a cell:
//                celltabSel[icellSel] = p;
//                icellSel++;
                Selected(p) = TRUE;
                Selected(q) = TRUE;
//E
            }
            if (Nb(p)==gd.nsmooth[2]) {
                ip++;
                q = ip + bodytabsm -1;
                Id(q) = ip;
// BODY3
                Type(q) = BODY;
//            Type(q) = BODY3;
//            Nbb(q) = Nb(p);
                Weight(q) = Weight(p);
                SETV(Pos(q), Pos(p));
                Kappa(q) = Kappa(p);

//B To see the bodies belonging to a cell:
//                celltabSel[icellSel] = p;
//                icellSel++;
                Selected(p) = TRUE;
                Selected(q) = TRUE;
//E
            }
*/
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
//            verb_print_debug(1, "\nAqui voy (2): %d %d %d\n", lev, inodelev, Type(q));
            for (l = More(q); l != Next(q); l = Next(l)) {
                if (Type(l)==CELL) {
//                verb_print_debug(1, "\nAqui voy (22): %d %d %d\n", lev, inodelev, Type(l));
                nodetabscanlev[inodelev] = l;
                inodelev++;
                } else {
                    ibodyleftout++;
                }
            }
        } else {
            ibodyleftout++;
//        verb_print_debug(1, "\nAqui voy (2): %d %d %d\n", lev, inodelev, Type(q));
//            nodetabscanlev[inodelev] = q;
//            inodelev++;
//            verb_print_debug(1, "\nAqui voy (2): %d %d %d\n", lev, inodelev,
//                             Type(nodetabscanlev[inodelev]));
        }
    } else {
        if ( lev+1 != cmd.scanLevel )
            if (Type(q)==CELL) {
            for (l = More(q); l != Next(q); l = Next(l)) {
                walktree_index_scan_lev(l, lev+1);
            }
            } else {
                ibodyleftout++;
            }
    }
}

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
