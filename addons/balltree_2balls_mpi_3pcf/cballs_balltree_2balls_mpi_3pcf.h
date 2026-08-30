#ifndef _cballs_balltree_2balls_mpi_3pcf_h
#define _cballs_balltree_2balls_mpi_3pcf_h

case BALLTREE2BALLSMPI3PCFMETHOD:
    verb_print(cmd->verbose,
               "\n\tevalHist: with distributed FCFC ball-tree "
               "LogMultipole 3PCF\n\n");
    if (searchcalc_balltree_2balls_mpi_3pcf(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
        return FAILURE;
    break;

#endif /* !_cballs_balltree_2balls_mpi_3pcf_h */
