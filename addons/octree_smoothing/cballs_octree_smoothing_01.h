// Use:
//#include "cballs_octree_smoothing_01.h"

//
// included in file: addons/addons_include/source/cballs/cballs_include01.h
//

#ifndef _cballs_octree_smoothing_01_h
#define _cballs_octree_smoothing_01_h

gd->flagSmoothCellMin = FALSE;
gd->flagSmooth = FALSE;
gd->flagSetNbNoSel = FALSE;

if (cballs_opt_smooth_min_cell(cmd)) {
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\n\t%s: smooth cell min: try making tree...\n\n", routineName);
    DO_BODY(p,bodytable[gd->iCatalogs[0]],
            bodytable[gd->iCatalogs[0]]+gd->nbodyTable[gd->iCatalogs[0]])
        Update(p) = TRUE;
    if (MakeTree(cmd, gd, bodytable[gd->iCatalogs[0]],
                 gd->nbodyTable[gd->iCatalogs[0]], gd->iCatalogs[0])
        == FAILURE)
        return FAILURE;
//B
    free(bodytable[gd->iCatalogs[0]]);
    gd->bytes_tot -= gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body);
    gd->nbodyTable[gd->iCatalogs[0]] = gd->nnodescanlevTable[gd->iCatalogs[0]];
    bodytable[ifile] =
                (bodyptr) allocate(gd->nnodescanlevTable[gd->iCatalogs[0]]
                * sizeof(body));
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"Allocated %g MByte for final (smoothCellMin) particle (%ld) storage (%ld).\n",
                    gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body)*INMB,
                    gd->nnodescanlevTable[gd->iCatalogs[0]],
                    gd->nnodescanlevTable[gd->iCatalogs[0]]);
    kavg = 0;
    q = nodetable;
    bodytable[gd->iCatalogs[0]] = nodetable;
    in = 0;
    DO_BODY(p,bodytable[gd->iCatalogs[0]],
        bodytable[gd->iCatalogs[0]]+gd->nnodescanlevTable[gd->iCatalogs[0]]) {
        kavg += Kappa(p);
        in++;
    }
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "%ld %d %ld %lg %ld %lg\n",
                in,Type(p-1),Id(p-1),Mass(p-1),Nb(p-1),Kappa(p-1));
//E
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "smoothCellMin: %ld particles in nodetable\n", in);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "smoothCellMin: Average of kappa (%ld particles) = %le\n",
                gd->nnodescanlevTable[gd->iCatalogs[0]],
                kavg/((real)gd->nnodescanlevTable[gd->iCatalogs[0]]));
    gd->flagSmoothCellMin = TRUE;
} // ! smooth-min-cell

if (cballs_opt_smooth(cmd) && !cballs_opt_set_nb_nosel(cmd)) {
    verb_print(cmd->verbose, "\n\tMainLoop: smooth: try making tree...\n\n");
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        Update(p) = TRUE;
    if (MakeTree(cmd, gd, bodytable[ifile], cmd->nbody, 0) == FAILURE)
        return FAILURE;
//B
    //B correct these lines following cballsio.c
    //      routine InputData_all_in_one
    //      lines: 270--287
    free(bodytable[ifile]);
    gd->bytes_tot -= cmd->nbody*sizeof(body);
    cmd->nbody = gd->nbodysm;
        bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += cmd->nbody*sizeof(body);
    verb_print(cmd->verbose,
        "Allocated %g MByte for final (smooth) particle (%ld) storage.\n",
        cmd->nbody*sizeof(body)*INMB, cmd->nbody);
    //E
    
    
    kavg = 0;
    q = bodytabsm;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        SETV(Pos(p),Pos(q));
// BODY3
        Type(p) = Type(q);
#ifdef BODY3ON
        Nbb(p) = Nbb(q);
#endif
        Mass(p) = Mass(q);
        Kappa(p) = Kappa(q);
        Id(p) = p-bodytable[ifile]+1;
        Update(p) = TRUE;
        q++;
        kavg += Kappa(p);
    }
    free(bodytabsm);
    gd->bytes_tot -= gd->nbodysm*sizeof(body);
//E
    verb_print(cmd->verbose,
               "smooth: Average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );
    gd->flagSmooth = TRUE;
} // ! smooth && set-Nb-noSel

if ( cballs_opt_smooth(cmd) && cballs_opt_set_nb_nosel(cmd) ) {
    verb_print(cmd->verbose,
           "\n\tMainLoop: smooth & set-Nb-noSel: try making tree...\n\n");
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        Update(p) = TRUE;
    if (MakeTree(cmd, gd, bodytable[ifile], cmd->nbody, 0) == FAILURE)
        return FAILURE;
//B
    free(bodytable[ifile]);
    gd->bytes_tot -= cmd->nbody*sizeof(body);
    cmd->nbody = gd->nbodySel;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
    gd->bytes_tot += cmd->nbody*sizeof(body);
    verb_print(cmd->verbose,
    "Allocated %g MByte for final (smooth-noSel) particle (%ld) storage.\n",
    cmd->nbody*sizeof(body)*INMB, cmd->nbody);
    kavg = 0;
    q = bodytabSel;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        SETV(Pos(p),Pos(q));
// BODY3
        Type(p) = Type(q);
#ifdef BODY3ON
        Nbb(p) = Nbb(q);
#endif
        Mass(p) = Mass(q);
        Kappa(p) = Kappa(q);
        Id(p) = p-bodytable[ifile]+1;
        Update(p) = TRUE;
        q++;
        kavg += Kappa(p);
    }
    free(bodytabSel);
    gd->bytes_tot -= gd->nbodySel*sizeof(body);
//E
    verb_print(cmd->verbose,
        "smooth-set-Nb-noSel: Average of kappa (%ld particles) = %le\n",
        cmd->nbody, kavg/((real)cmd->nbody) );
    gd->flagSmooth = TRUE;
    gd->flagSetNbNoSel = TRUE;
} // ! smooth && set-Nb-noSel


#endif	// ! _cballs_octree_smoothing_01_h
