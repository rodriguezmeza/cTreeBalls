// Use:
//#include "cballs_octree_box_omp.h"

// included in addons_include/source/cballs/cballs_include_03.h

#ifndef _cballs_print_octree_box_omp_h
#define _cballs_print_octree_box_omp_h

//#define OCTREEBOXOMP         72

        case 72:
            verb_print(cmd->verbose,
                    "\n\tevalHist: printing octree-box-omp method\n\n");
            if (scanopt(cmd->options, "compute-HistN"))
                PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            break;

#endif	// ! _cballs_print_octree_omp_h
