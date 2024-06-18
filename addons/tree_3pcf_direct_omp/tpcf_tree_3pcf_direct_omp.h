// Use:
//NOLSST:
//#include "tpcf_tree_3pcf_direct_omp.h"

#ifndef _tpcf_tree_3pcf_direct_omp_h
#define _tpcf_tree_3pcf_direct_omp_h

//#define TREE3PCFBFOMPMETHOD     45

case 45:                   // search=tree-3pcf-direct
    verb_print(cmd->verbose,
    "\n\tevalHist: with normal tree method (tree-3pcf-direct-omp)\n\n");
for (ifile=0; ifile<gd->ninfiles; ifile++) {
    DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
    Update(p) = TRUE;
    maketree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
}
searchcalc_normal_3pcf_direct_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
            gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);

    break;

#endif	// ! _tpcf_tree_3pcf_direct_omp_h
