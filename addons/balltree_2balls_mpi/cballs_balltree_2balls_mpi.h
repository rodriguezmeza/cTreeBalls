#ifndef _cballs_balltree_2balls_mpi_h
#define _cballs_balltree_2balls_mpi_h

case BALLTREE2BALLSMPIMETHOD:
    verb_print(cmd->verbose,
               "\n\tevalHist: with distributed FCFC ball-tree "
               "dual/triple-node method\n\n");
    if (searchcalc_balltree_2balls_mpi(
            cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
            gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
        return FAILURE;
    break;

#endif /* !_cballs_balltree_2balls_mpi_h */
