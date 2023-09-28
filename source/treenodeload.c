/* ==============================================================================
!	MODULE: treenodeload.c			[cTreeBalls]								!
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
local void expandbox(nodeptr, int);
local void loadnode(nodeptr);
local int subindex(nodeptr, cellptr);
local void hackCellProp(cellptr, real, int);
local void setRadius(cellptr, vector, real);
local void threadtree(nodeptr, nodeptr);
local void findRootCenter(nodeptr, int);
local void centerNodes(nodeptr, int);

#define MAXLEVEL  32  

local int cellhist[MAXLEVEL];
local int subnhist[MAXLEVEL];
local int tdepth;
local int ncell;
local real rSize;
local real cputree;
local INTEGER bytes_tot_cells;

global void maketreenodes(nodeptr ntab, int nnode)
{
    double cpustart;
    nodeptr p;
    int i;

    cpustart = CPUTIME;
    DO_BODY(p, ntab, ntab+nnode) {
        if (Type(p) == CELL)
            Type(p) = NODECELL;
        else
            Type(p) = NODEBODY;
    }
    newtree();
    rootnode = makecell();
    findRootCenter(ntab, nnode);
    centerNodes(ntab, nnode);
    findRootCenter(ntab, nnode);
    CLRV(Pos(rootnode));
    expandbox(ntab, nnode);
	DO_BODY(p, ntab, ntab+nnode)
        loadnode(p);
    tdepth = ncell = bytes_tot_cells = 0;
    for (i = 0; i < MAXLEVEL; i++)
        cellhist[i] = subnhist[i] = 0;
    hackCellProp(rootnode, rSize, 0);
    threadtree((nodeptr) rootnode, NULL);
    DO_BODY(p, ntab, ntab+nnode) {
        if (Type(p) == NODECELL)
            Type(p) = CELL;
        else
            Type(p) = BODY;
    }
    cputree = CPUTIME - cpustart;
    verb_print(cmd.verbose, "\nCPU tree node time : %lf\n", cputree);
}

local void findRootCenter(nodeptr ntab, int nnode)
{
    real len;
    nodeptr p;
    int k;
    vector xmin, xmax;

    DO_COORD(k)
        xmin[k] = xmax[k] = Pos(ntab)[k];

    DO_BODY(p, ntab, ntab+nnode)
        DO_COORD(k) {
            if (Pos(p)[k] > xmax[k])
                xmax[k] = Pos(p)[k];
            if (Pos(p)[k] < xmin[k])
                xmin[k] = Pos(p)[k];
        }

    DO_COORD(k) {
        Pos(rootnode)[k] = (xmax[k]+xmin[k])/2;
        verb_print(cmd.verbose, "findRootCenter: Pos(rootnode) = %lf\n", Pos(rootnode)[k]);
    }

    for(k=0, len=xmax[0]-xmin[0]; k<NDIM; k++)
        if((xmax[k]-xmin[k])>len)
            len=xmax[k]-xmin[k];

    verb_print(cmd.verbose, "findRootCenter: len = %lf\n", len);
}

local void centerNodes(nodeptr ntab, int nnode)
{
    nodeptr p;
    int k;

    DO_BODY(p, ntab, ntab+nnode)
        DO_COORD(k)
            Pos(p)[k] = Pos(p)[k] - Pos(rootnode)[k];
}

local void newtree(void)
{
    rootnode = NULL;
    ncell = 0;
}

local cellptr makecell(void)
{
    cellptr c;
    int i;
 
    c = (cellptr) allocate(sizeof(cell));
    Type(c) = CELL;
    Update(c) = FALSE;
    for (i = 0; i < NSUB; i++)                  
        Subp(c)[i] = NULL;
    ncell++;
    bytes_tot_cells += sizeof(cell);
    return (c);
}

local void expandbox(nodeptr ntab, int nnode)
{
    real dmax, d;
    nodeptr p;
    int k;
 
    rSize = 1.0;
    dmax = 0.0;                                 
	DO_BODY(p, ntab, ntab+nnode)
		DO_COORD(k) {
            d = rabs(Pos(p)[k] - Pos(rootnode)[k]);
            if (d > dmax)
                dmax = d;                       
        }
    while (rSize < 2 * dmax)
        rSize = 2 * rSize;

    verb_print(cmd.verbose, "treenodeload expandbox: rSize = %lf\n", rSize);
}

local void loadnode(nodeptr p)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;

    q = rootnode;
    qind = subindex(p, q);
    qsize = rSize;
    while (Subp(q)[qind] != NULL) {
        if (Type(Subp(q)[qind]) == NODEBODY || Type(Subp(q)[qind]) == NODECELL) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0) {
                verb_log_print(cmd.verbose_log,gd.outlog,
                               "\nIds: %ld and %ld have the same position\n",
                               Id(p),Id(Subp(q)[qind]));
                DO_COORD(k)
                    verb_log_print(cmd.verbose_log,gd.outlog,"\nPos[k]: %le %le\n",
                               Pos(p)[k],Pos(Subp(q)[qind])[k]);
                error("loadnode: two nodes have same position\n");
            }
            c = makecell();
			DO_COORD(k)
                Pos(c)[k] = Pos(q)[k] +
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;

            Subp(c)[subindex((nodeptr) Subp(q)[qind], c)] = Subp(q)[qind];
            Subp(q)[qind] = (nodeptr) c;
        }
        q = (cellptr) Subp(q)[qind];
        qind = subindex(p, q);
        qsize = qsize / 2;
    }
    Subp(q)[qind] = (nodeptr) p;
}

local int subindex(nodeptr p, cellptr q)
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
 
    tdepth = MAX(tdepth, lev);
    cellhist[lev]++;
    Weight(p) = 0.0;
    Kappa(p) = 0.0;
    CLRV(cmpos);
    for (i = 0; i < NSUB; i++) {
        if ((q = Subp(p)[i]) != NULL) {
            subnhist[lev]++;                    
            if (Type(q) == CELL)
                hackCellProp((cellptr) q, psize/2, lev+1);
            Update(p) |= Update(q);
            Weight(p) += Weight(q);
            Kappa(p) += Weight(q)*Kappa(q);
            MULVS(tmpv, Pos(q), Weight(q));
            ADDV(cmpos, cmpos, tmpv);           
        }
    }
    if (Weight(p) > 0.0)
        SETV(cmpos, Pos(p));

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
            }
        }
#undef EPSILON

    setRadius(p, cmpos, psize);
    SETV(Pos(p), cmpos);
}

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
}

local void threadtree(nodeptr p, nodeptr n)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;
    if (Type(p) == CELL) {
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


