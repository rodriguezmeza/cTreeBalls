#ifndef _protodefs_balltree_2balls_mpi_h
#define _protodefs_balltree_2balls_mpi_h

#define BALLTREE2BALLSMPIMETHOD 179

#include "fcfc_balltree_2balls_mpi.h"

global int searchcalc_balltree_2balls_mpi(
    struct cmdline_data *, struct global_data *, bodyptr *,
    INTEGER *, INTEGER, INTEGER *, int, int);

#endif /* !_protodefs_balltree_2balls_mpi_h */
