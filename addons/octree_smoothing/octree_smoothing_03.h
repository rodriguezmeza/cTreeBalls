//=============================================================================
//        1          2          3          4        ^ 5          6          7
// Use:
//#include "octree_smoothing_03.h"

//
// included in file: addons/addons_include/source/tree/treeload_include03.h
//

#ifndef _octree_smoothing_03_h
#define _octree_smoothing_03_h

//B update use of this compilation flag
// Check thoroughly if this flag is not longer needed...
//#ifdef BALLS

if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) ) {
    cmd->scanLevel = gd->tdepthTable[ifile] + 1 + gd->scanLevelMin[0];
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\tsmoothCellMin: fixing scanLevel to: %d\n",
                           cmd->scanLevel);
} else {
    if (cballs_opt_set_default_param(cmd)) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\tfixing scanLevel to tdepth-1...\n");
        cmd->scanLevel = MAX(gd->tdepthTable[ifile]-1,3);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\tfinal value is %d.\n", cmd->scanLevel);
    } else {
        if (cmd->scanLevel > gd->tdepthTable[ifile]) {
            verb_print_warning(cmd->verbose,
                "Warning! tree depth (%d) is less than scanLevel (%d)...\n",
                gd->tdepthTable[ifile], cmd->scanLevel);
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "\tfixing to tdepth-1...\n");
            cmd->scanLevel = MAX(gd->tdepthTable[ifile]-1,3);
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "\tfinal value is %d.\n", cmd->scanLevel);
        }
    }
}

if (cmd->useLogHist) { // ! rminHist = 0 not allowed when useLogHist is true
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "\n(Only in log-scale) deltaR is %g and root size at scanLevel is %g.\n",
    gd->deltaR, gd->rSizeTable[ifile]/rpow(2.0,cmd->scanLevel));
    i = 0;
    while (gd->deltaR < gd->rSizeTable[ifile]/rpow(2.0,i)) i++;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\t\t\tSuggested scanLevel is %d, where root size is %g.\n",
                i, gd->rSizeTable[ifile]/rpow(2.0,i));
}

inodelev = 0;
ibodyleftout = 0;
if (cmd->scanLevel==0) {
    gd->nnodescanlevTable[ifile] = 0;
    gd->bytes_tot += gd->nnodescanlevTable[ifile]*sizeof(nodeptr);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
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
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\nAllocated %g MByte for (%d) scan nodetab storage.\n",
                INMB*gd->nnodescanlevTable[ifile]*sizeof(nodeptr),
                gd->nnodescanlevTable[ifile]);

    if (!gd->flagSmoothCellMin && cballs_opt_smooth_min_cell(cmd) ) {
        walktree_index_scan_lev((nodeptr)roottable[ifile], 0,
                                ifile, cmd->scanLevel);
    } else {
        // Check this line. May be it is correct:
        //      cmd->scanLevel = MAX(gd->tdepthTable[ifile],3);
        // verb_print(cmd->verbose, "\tfinal value is %d.\n", cmd->scanLevel);
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

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
               "\nFound %d nodes to scan at level %d.\n",
               inodelev, cmd->scanLevel);

    // Freeing some segment of memory will be necessary
    gd->nnodescanlevTable[ifile] = inodelev;
    save_nodes(cmd, gd, ifile);
} // ! scanLevel==0

//B socket:
//#ifdef ADDONS
//#include "tree_include.h"
//#endif
//E

//B Root nodes to scan:
inodelev_root = 0;
ibodyleftout_root = 0;

if (cmd->scanLevelRoot==0) {
    gd->nnodescanlev_rootTable[ifile] = 0;
    gd->bytes_tot += gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\nAllocated %g MByte for (%d) scan root nodetab storage.\n",
            INMB*gd->nnodescanlev_rootTable[ifile]*sizeof(nodeptr),
            gd->nnodescanlev_rootTable[ifile]);
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

// Check thoroughly if this flag is not longer needed...
//#endif // ! BALLS
//E update use of this compilation flag

#endif	// ! _octree_smoothing_03_h
