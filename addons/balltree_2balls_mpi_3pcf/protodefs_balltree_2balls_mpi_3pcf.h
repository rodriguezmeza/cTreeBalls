#ifndef _protodefs_balltree_2balls_mpi_3pcf_h
#define _protodefs_balltree_2balls_mpi_3pcf_h

#define BALLTREE2BALLSMPI3PCFMETHOD 178

#include "fcfc_balltree_2balls_mpi_3pcf.h"

global int searchcalc_balltree_2balls_mpi_3pcf(
    struct cmdline_data *, struct global_data *, bodyptr *,
    INTEGER *, INTEGER, INTEGER *, int, int);

#endif /* !_protodefs_balltree_2balls_mpi_3pcf_h */
