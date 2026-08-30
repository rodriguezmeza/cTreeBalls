// Use:
//#include "cballs_balls_omp_0357.h"

#ifndef _cballs_balls_omp_0357_h
#define _cballs_balls_omp_0357_h

//#define BALLSOMPMETHOD         79

//#ifdef BALLS
case 79:
    verb_print(cmd->verbose,
               "\n\tevalHist: with balls tree-omp method\n\n");
    if (MainLoop_0357(cmd, gd) == FAILURE)
        return FAILURE;
/*
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        MakeTree_0357(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    searchcalc_balls_omp_0357(cmd, gd, bodytable, gd->nbodyTable,
                        1, gd->nbodyTable,
                        gd->iCatalogs[0], gd->iCatalogs[1]);
*/
    break;

/*
case SEARCHNULL:
    verb_print(cmd->verbose, "\n\tEvalHist: null search method.\n");
    verb_print(cmd->verbose,
               "\n\tevalHist: with balls tree-omp method\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    searchcalc_balls_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                    gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;
default:
    verb_print(cmd->verbose, "\n\tEvalHist: dafault search method.\n");
    verb_print(cmd->verbose,
               "\n\tevalHist: with balls tree-omp method\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    searchcalc_balls_omp(cmd, gd, bodytable, gd->nbodyTable, 1,
                gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;
#else
case SEARCHNULL:
    verb_print(cmd->verbose, "\n\tEvalHist: null search method.\n");
    verb_print(cmd->verbose,
            "\n\tevalHist: with normal tree method (sincos-omp)\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    searchcalc_normal_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
                gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;
default:
    verb_print(cmd->verbose, "\n\tEvalHist: dafault search method.\n");
    verb_print(cmd->verbose,
            "\n\tevalHist: with normal tree method (sincos-omp)\n\n");
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
        Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
    }
    searchcalc_normal_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
                gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
    break;
#endif
*/

#endif	// ! _cballs_balls_omp_0357_h
