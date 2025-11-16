// Use:
//#include "cballs_kdtree_box_omp.h"

// included in addons_include/source/cballs/cballs_include_03.h

#ifndef _cballs_print_kdtree_box_omp_h
#define _cballs_print_kdtree_box_omp_h

//#define KDTREEBOXOMP         71

        case 71:
            verb_print(cmd->verbose,
                    "\n\tevalHist: printing kdtree-box-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN"))
                PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            break;

#endif	// ! _cballs_print_kdtree_omp_h
