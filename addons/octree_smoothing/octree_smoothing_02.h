//=============================================================================
//        1          2          3          4        ^ 5          6          7
// Use:
//#include "octree_smoothing_02.h"

//
// included in file: addons/addons_include/source/tree/treeload_include02.h
//

#ifndef _octree_smoothing_02_h
#define _octree_smoothing_02_h


local void threadtree_smooth(struct  cmdline_data* cmd,
                      struct  global_data* gd, nodeptr p, nodeptr n,
                             int ifile)
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
//        celltable[ifile][gd->ncellt[ifile]] = (cellptr)p;
//        gd->ncellt[ifile] = gd->ncellt[ifile]+1;
        //E
#endif
        //B Smooth(ing) section
        if (Nb(p)<NbMax)
            cellhistNb[Nb(p)]++;
        if (!gd->flagSmooth) {
            if (cballs_opt_smooth(cmd)) {
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
            threadtree_smooth(cmd, gd, desc[i], desc[i+1], ifile);
    } // ! p = CELL
}


//B Smooth(ing) section
local int smoothBodies(struct  cmdline_data* cmd,
                       struct  global_data* gd, bodyptr btab, INTEGER nbody)
{
    bodyptr p;

    if (!gd->flagSmooth || !gd->flagSetNbNoSel) {
        if ( (cballs_opt_smooth(cmd)
              && cballs_opt_set_nb_nosel(cmd)) ) {
            printf("NTOT = %d \t%" INTEGER_FMT " \t%" INTEGER_FMT "\n\n",
                   gd->nsmooth[0], NTOT[0]+inosel,
                   gd->nsmooth[0]*NTOT[0]+inosel);
            bodytabSel = (bodyptr) allocate((NTOT[0]+inosel) * sizeof(body));
            gd->bytes_tot += (NTOT[0]+inosel)*sizeof(body);
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "Allocated %g MByte for (smooth) particle (%" INTEGER_FMT
                    ") storage.\n",
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
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "Added %ld smoothed cells...\n",ipcount);
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
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "Added %ld total bodies...\n",ipcount);
            gd->nbodySel = ipcount;
        }
    }

    return SUCCESS;
}
//E

// Check thoroughly if this flag is not longer needed...
//#ifdef BALLS
//
// Unify: walktree_index_scan_lev and walktree_index_scan_lev_root
//
local void walktree_index_scan_lev(nodeptr q, int lev, int ifile, int scanLevel)
{
    nodeptr l;

    if (q == NULL)
        return;
    if (lev >= scanLevel || Type(q) != CELL) {
        nodetablescanlev[ifile][inodelev++] = q;
        if (Type(q) != CELL)
            ibodyleftout++;
        return;
    }

    for (l = More(q); l != Next(q); l = Next(l))
        walktree_index_scan_lev(l, lev+1, ifile, scanLevel);
}

local void walktree_index_scan_lev_root(struct cmdline_data* cmd, struct  global_data* gd,
                                        nodeptr q, int lev, int ifile)
{
    nodeptr l;

    if (q == NULL)
        return;
    if (lev >= cmd->scanLevelRoot || Type(q) != CELL) {
        nodetablescanlev_root[ifile][inodelev_root++] = q;
        if (Type(q) != CELL)
            ibodyleftout_root++;
        return;
    }

    for (l = More(q); l != Next(q); l = Next(l))
        walktree_index_scan_lev_root(cmd, gd, l, lev+1, ifile);
}
// Check thoroughly if this flag is not longer needed...
//#endif // ! BALLS


local int save_nodes(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
    gd->nodesfilePath[0] = '\0';
    gd->bodiesfilePath[0] = '\0';

    if (cmd->verbose==VERBOSEDEBUGINFO) {
//        sprintf(gd->nodesfilePath,"%s/nodes%s.txt",gd->tmpDir,cmd->suffixOutFiles);
        if (format_checked(gd->nodesfilePath, sizeof(gd->nodesfilePath),
            "nodesfilePath", "%s/nodes%s.txt",gd->tmpDir,cmd->suffixOutFiles) != 0)
            return FAILURE;

        if(!(gd->outnodelev=fopen(gd->nodesfilePath, "w")))
            error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
    }
//B
    if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) ) {
        if (cmd->verbose==3) {
//            sprintf(gd->bodiesfilePath,"%s/bodies%s.txt",
//                    gd->tmpDir,cmd->suffixOutFiles);
            if (format_checked(gd->bodiesfilePath, sizeof(gd->bodiesfilePath),
                "bodiesfilePath", "%s/bodies%s.txt",
                gd->tmpDir,cmd->suffixOutFiles) != 0)
                return FAILURE;

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

    if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) )
        nodetable = (bodyptr) allocate(gd->nnodescanlevTable[ifile] * sizeof(body));

    for (in=0; in<inodelev; in++) {
//B
        if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) ) {
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

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "Found %ld particles in %ld nodes... vs number of total cells %ld\n",
            nodescount, gd->nnodescanlevTable[ifile], gd->ncellTable[ifile]);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "...%ld cells at scan level...\n",
                           sumcells);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "...and %ld bodies. Bodies in upper levels: %ld\n",
                           sumbodies, cmd->nbody-sumbodies);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%ld particles were left out of cells at scan level.\n",
                           ibodyleftout);
    //B Smooth(ing) section
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%ld cells were with at much %d particles in them.\n",
                           nodescount_smooth,gd->nsmooth[0]);
    //E
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "Checking sums (bodyleftout+nodescount): %ld.\n",
                           ibodyleftout+nodescount);

    if (cmd->verbose==VERBOSEDEBUGINFO)
        fclose(gd->outnodelev);
//B
    if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) )
        if (cmd->verbose==VERBOSEDEBUGINFO)
            fclose(gd->outbodylev);
//E
    return SUCCESS;
}

local int save_nodes_root(struct  cmdline_data* cmd, struct  global_data* gd, int ifile)
{
//    sprintf(gd->nodesfilePath,"%s/nodes_root%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if (format_checked(gd->nodesfilePath, sizeof(gd->nodesfilePath),
        "nodesfilePath", "%s/nodes_root%s.txt",gd->tmpDir,cmd->suffixOutFiles) != 0)
        return FAILURE;

    if(!(gd->outnodelev=fopen(gd->nodesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
//    sprintf(gd->bodiesfilePath,"%s/bodies_root%s.txt",gd->tmpDir,cmd->suffixOutFiles);
    if (format_checked(gd->bodiesfilePath, sizeof(gd->bodiesfilePath),
        "bodiesfilePath", "%s/bodies_root%s.txt",gd->tmpDir,cmd->suffixOutFiles) != 0)
        return FAILURE;

    if(!(gd->outbodylev=fopen(gd->bodiesfilePath, "w")))
        error("\nstart_Common: error opening file '%s' \n",gd->nodesfilePath);
    int in;
    INTEGER nodescount=0;
    bodyptr pn;
    INTEGER nodescount_smooth=0;
    INTEGER sumbodies=0, sumcells=0;
    INTEGER nodescount_thread=0, nodescount_thread_total=0;

    nodetable_root =            // check memory allocation freed at the end...
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
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "Checking sums (bodyleftout+nodescount): %ld.\n",
                            ibodyleftout+nodescount);

    fclose(gd->outnodelev);
    fclose(gd->outbodylev);
    free(nodetable_root);

    return SUCCESS;
}

//E

#endif	// ! _octree_smoothing_02_h
