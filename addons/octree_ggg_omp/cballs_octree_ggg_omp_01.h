// Use:
//#include "cballs_octree_ggg_omp_01.h"

//
// included in file: addons/addons_include/source/cballs/cballs_include01.h
//

#ifndef _cballs_octree_ggg_omp_01_h
#define _cballs_octree_ggg_omp_01_h

#ifdef BALLS4SCANLEV
    int k;
    gd->flagBalls4Scanlevel = FALSE;

#ifdef BALLS4SCANLEV2NDROUND

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n\t%s: BALLS4SCANLEVEL: making tree first round...\n\n", routineName);

    DO_BODY(p,bodytable[ifile], bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
    MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);

    verb_print_debug(1, "\nAqui voy (0): %d\n", ifile);

    INTEGER nbodylocal;
    
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "\n\t%s: BALLS4SCANLEVEL: found %ld nodes at upper most level (%ld)...\n\n", routineName, gd->nnodescanlevTableB4[ifile], gd->nbodyTable[ifile]);

    nbodylocal = gd->nnodescanlevTableB4[ifile];

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        free(bodytable[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++)
            free(bodytable[ifile]);
    }

//    gd->nbodyTable[ifile] = gd->nnodescanlevTableB4[ifile];
    ifile=0;
    gd->nbodyTable[ifile] = nbodylocal;
    gd->nnodescanlevTableB4[ifile] = nbodylocal;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "\n\t%s: BALLS4SCANLEVEL: found %ld nodes at upper most level (%ld)...\n\n",
    routineName, gd->nnodescanlevTableB4[ifile], gd->nbodyTable[ifile]);


    verb_print_debug(1, "\nAqui voy (1): %d\n", ifile);

    bodytable[ifile] = (bodyptr) allocate(gd->nbodyTable[ifile] * sizeof(body));
    verb_print_debug(1, "\nAqui voy (2)\n");

    gd->bytes_tot += gd->nbodyTable[ifile] * sizeof(body);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"Allocated %g MByte for final (BALLS4SCANLEVEL) particle (%ld) storage (%ld).\n",
                           gd->nbodyTable[ifile] * sizeof(body)*INMB,
                           gd->nbodyTable[ifile],
                           gd->nbodyTable[ifile]);

    for (INTEGER i=0; i< gd->nnodescanlevTableB4[ifile]; i++) {
        q = (bodyptr)nodetablescanlevB4[ifile][i];
        p = bodytable[ifile] + i;
        DO_COORD(k)
            Pos(p)[k] = Pos(q)[k];
        Kappa(p) = Kappa(q);
        Mass(p) = Mass(q);
        Weight(p) = Weight(q);
        Id(p) = i+1;
        Type(p) = BODY;
        Mask(p) = Mask(q);
    }

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        free(nodetablescanlevB4[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++)
            free(nodetablescanlevB4[ifile]);
    }

    if (!scanopt(cmd->searchMethod, "kdtree-omp")
        && !scanopt(cmd->searchMethod, "kdtree-box-omp") ) {
        freeTree(cmd, gd);
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n\t%s: BALLS4SCANLEVEL: going tree second round...\n\n", routineName);
    gd->flagBalls4Scanlevel = TRUE;
#endif // ! BALLS4SCANLEV2NDROUND
#endif


#endif	// ! _cballs_octree_ggg_omp_01_h
