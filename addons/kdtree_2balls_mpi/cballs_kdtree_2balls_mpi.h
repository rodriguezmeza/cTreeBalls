#ifndef _cballs_kdtree_2balls_mpi_h
#define _cballs_kdtree_2balls_mpi_h

    case KDTREE2BALLSMPIMETHOD:
        if (searchcalc_kdtree_2balls_mpi(
                cmd, gd, bodytable, gd->nbodyTable, 1, gd->nbodyTable,
                gd->iCatalogs[0], gd->iCatalogs[1]) == FAILURE)
            return FAILURE;
        break;

#endif /* !_cballs_kdtree_2balls_mpi_h */
