//=============================================================================
//        1          2          3          4        ^ 5          6          7
// Use:
//#include "octree_smoothing_01.h"

//
// included in file: addons/addons_include/source/tree/treeload_include01.h
//

#ifndef _octree_smoothing_01_h
#define _octree_smoothing_01_h

//B Smooth(ing) section
    if (!gd->flagSmooth)                             // Flag to smooth bodies
        if (cballs_opt_smooth(cmd)) {
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "NTOT = %ld \t%ld \t%ld\n",
                    gd->nsmooth[0], NTOT[0], gd->nsmooth[0]*NTOT[0]);
            gd->nbodysm = NTOT[0];
            bodytabsm = (bodyptr) allocate((NTOT[0]) * sizeof(body));
            gd->bytes_tot += (NTOT[0])*sizeof(body);
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "Allocated %g MByte for (smooth %ld) particle storage.\n",
                (NTOT[0])*sizeof(body)*INMB, (NTOT[0]));
        }
//E

ip = 0;

cpustartMiddle = CPUTIME;

#ifndef MACONLY
    //B celltable
    ncell=0;
//    gd->ncellt[ifile]=0;
    celltable[ifile] =
        (cellptr *) allocate(gd->ncellTable[ifile] * sizeof(cellptr));
    //E
#endif

threadtree_smooth(cmd, gd, (nodeptr) roottable[ifile], NULL, ifile);
verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "%s: threadtree CPU time: %lf %s\n",
            routineName, CPUTIME - cpustartMiddle, PRNUNITOFTIMEUSED);

verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "threadtree: number ip of selected cells = %ld\n",ip);
verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "%d real node (range of nodes to search: >nc1 && <nc2)\n",inode);
gd->rnnode = inode;

#endif	// ! _octree_smoothing_01_h
