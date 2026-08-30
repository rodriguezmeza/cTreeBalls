#ifndef _cballs_print_balltree_2balls_mpi_h
#define _cballs_print_balltree_2balls_mpi_h

case BALLTREE2BALLSMPIMETHOD:
    verb_print(cmd->verbose,
               "\n\tprintEvalHist: printing distributed FCFC ball-tree "
               "dual/triple-node method\n\n");
#ifdef TWOPCF
    if (!scanopt(cmd->options, "only-3pcf")) {
        if (cballs_opt_compute_histn(cmd))
            PRINT_OR_FAIL(PrintHistNN(cmd, gd));
    }
#endif
    PRINT_OR_FAIL(PrintHistrBins(cmd, gd));
#ifdef TWOPCF
    if (!scanopt(cmd->options, "only-3pcf"))
        PRINT_OR_FAIL(PrintHistXi2pcf(cmd, gd));
#endif
#ifdef THREEPCFCONVERGENCE
    if (!scanopt(cmd->options, "only-2pcf")) {
        PRINT_OR_FAIL(PrintHistZetaM_sincos(cmd, gd));
        if (cballs_opt_out_m_histzeta(cmd))
            PRINT_OR_FAIL(PrintHistZetaMm_sincos(cmd, gd));
    }
#endif
    break;

#endif /* !_cballs_print_balltree_2balls_mpi_h */
