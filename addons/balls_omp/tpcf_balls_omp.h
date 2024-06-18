// Use:
//ADDONS:
//#include "tpcf_balls_omp.h"

#ifndef _tpcf_balls_omp_h
#define _tpcf_balls_omp_h

//#define BALLSOMPMETHOD         46

#ifdef BALLS
        case 46:
            verb_print(cmd->verbose,
                       "\n\tevalHist: with balls tree-omp method\n\n");
            for (ifile=0; ifile<gd->ninfiles; ifile++) {
                DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
                Update(p) = TRUE;
                maketree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            }
            searchcalc_balls_omp(cmd, gd, bodytable, gd->nbodyTable,
                                 1, gd->nbodyTable,
                                 gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
#endif

#endif	// ! _tpcf_balls_omp_h
