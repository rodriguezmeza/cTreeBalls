// Use:
//#include "cballs_boxes_omp_omp.h"

// included in addons_include/source/cballs/cballs_include_02.h

#ifndef _cballs_neighbor_boxes_omp_h
#define _cballs_neighbor_boxes_omp_h

//#define NEIGHBORBOXESOMP         73

    case 73:
//        verb_print(cmd->verbose,
//                   "\n\tevalHist: with neighbor-boxes-omp method\n\n");
        verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n\tevalHist: with neighbor-boxes-omp method\n\n");
        searchcalc_neighbor_boxes_omp(cmd, gd,
                              bodytable, gd->nbodyTable,
                              1, gd->nbodyTable,
                              gd->iCatalogs[0], gd->iCatalogs[1]);
        break;

#endif	// ! _cballs_neighbor_boxes_omp_h
