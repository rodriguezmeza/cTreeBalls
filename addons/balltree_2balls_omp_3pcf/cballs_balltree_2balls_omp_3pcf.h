#ifndef _cballs_balltree_2balls_omp_3pcf_h
#define _cballs_balltree_2balls_omp_3pcf_h

    case BALLTREE2BALLSOMP3PCFMETHOD:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with TreeCorr-style LogMultipole method\n\n");
        if (searchcalc_balltree_2balls_omp_3pcf(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;

#endif /* !_cballs_balltree_2balls_omp_3pcf_h */
