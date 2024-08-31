// Use:
//#include "cballs_balls_omp.h"

#ifndef _cballs_balls_omp_h
#define _cballs_balls_omp_h

//#define BALLSOMPMETHOD         46

#ifdef BALLS
        case 46:
            verb_print(cmd->verbose,
                       "\n\tevalHist: with balls tree-omp method\n\n");
            for (ifile=0; ifile<gd->ninfiles; ifile++) {
                DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
                Update(p) = TRUE;
                MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            }
            searchcalc_balls_omp(cmd, gd, bodytable, gd->nbodyTable,
                                 1, gd->nbodyTable,
                                 gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
#endif

#endif	// ! _cballs_balls_omp_h
