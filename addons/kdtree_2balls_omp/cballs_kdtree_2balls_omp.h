#ifndef _cballs_kdtree_2balls_omp_h
#define _cballs_kdtree_2balls_omp_h

    case KDTREE2BALLSOMPMETHOD:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with KD two-ball/LogMultipole method\n\n");
        if (searchcalc_kdtree_2balls_omp(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;

#endif /* !_cballs_kdtree_2balls_omp_h */
