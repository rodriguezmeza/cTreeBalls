#ifndef _cballs_balltree_2balls_omp_h
#define _cballs_balltree_2balls_omp_h

    case BALLTREE2BALLSMETHOD:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with TreeCorr-style dual/triple-node method\n\n");
        if (searchcalc_balltree_2balls_omp(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;

#endif /* !_cballs_balltree_2balls_omp_h */
