#ifndef _cballs_balltree_omp_h
#define _cballs_balltree_omp_h

    case 167:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with FCFC balltree-omp method\n\n");
        if (searchcalc_balltree_omp(cmd, gd, bodytable, gd->nbodyTable,
                                   1, gd->nbodyTable,
                                   gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;

#ifdef BALLTREEMPI
    case 168:
        verb_print(cmd->verbose,
                   "\n\tevalHist: with FCFC balltree-mpi method\n\n");
        if (searchcalc_balltree_mpi(cmd, gd, bodytable, gd->nbodyTable,
                                   1, gd->nbodyTable,
                                   gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;
#endif

#endif /* !_cballs_balltree_omp_h */
