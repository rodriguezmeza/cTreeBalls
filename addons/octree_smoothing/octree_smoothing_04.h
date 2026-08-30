// Use:
//#include "octree_smoothing_04.h"

//
// included in file: addons/addons_include/source/tree/treeload_include04.h
//

#ifndef _octree_smoothing_04_h
#define _octree_smoothing_04_h

// Check thoroughly if this flag is not longer needed...
//#ifdef BALLS

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

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "Cell size at scanLevelMin (%d) and scale factor: %e %e\n",
                gd->scanLevelMin[0], gd->rminCell[0], gd->rminCell[1]);
    if (gd->rminCell[0] > gd->deltaRmax)
        verb_print_warning(cmd->verbose,
        "Warning! Cell size at scanLevelMin is greatear than deltaRmax: %e %e\n",
        gd->rminCell[0], gd->deltaRmax);
    if (gd->rminCell[0] > gd->deltaRmin)
        verb_print_warning(cmd->verbose,
        "Warning! Cell size at scanLevelMin is greatear than deltaRmin: %e %e\n",
        gd->rminCell[0], gd->deltaRmin);
    for (i=1; i<=cmd->sizeHistN-1; i++)
        if (gd->rminCell[0] < gd->ddeltaRV[i])
            break;
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "rminCell, ddeltaRV and deltaRV (at n = %d): %g %g %g\n",
                i+1, gd->rminCell[0], gd->ddeltaRV[i], gd->deltaRV[i+1]);

    if (gd->scanLevelMin[0] == 0) {
        gd->rminCell[1] = cmd->rangeN;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\tfixing rminCell[1] to: %g\n", gd->rminCell[1]);
    } else {
        gd->rminCell[1] = gd->deltaRV[i+1];
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\tfixing rminCell[1] to: %g\n", gd->rminCell[1]);
        if (gd->rminCell[1]>cmd->rangeN)
            verb_print_warning(cmd->verbose,
                       "Warning! rminCell[1] is greatear than rangeN\n");
    }
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "Cell size at scanLevelMin (%d) and scale factor (Modified values): %e %e\n",
    gd->scanLevelMin[0], gd->rminCell[0], gd->rminCell[1]);
    //E Root nodes

        //B version 1.0.1
    if (gd->infilefmt_int == INTAKAHASHI)
        verb_log_print(cmd->verbose_log, gd->outlog,
                       "Unit sphere (Takahashi): (S/N)^(1/2): %g\n\n",
                       rpow(2.0*TWOPI/gd->nbodyTable[gd->iCatalogs[0]],1.0/2.0));
        //E
    
} // ! useLogHist = true
    else {
        verb_print(cmd->verbose,
                   "Warning! can´t have normal hist and BALLS definition (useLogHist=%d)\nSet useLogHist = true if using balls-omp searching method\n",
                   cmd->useLogHist);
    }
// Check thoroughly if this flag is not longer needed...
//#endif // ! BALLS

#endif	// ! _octree_smoothing_04_h
