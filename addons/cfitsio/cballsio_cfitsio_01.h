
#ifndef _cballsio_cfitsio_01_h
#define _cballsio_cfitsio_01_h

case INCFITSIO:
//    verb_print(cmd->verbose,
//               "\n\tInput in fits format...\n");
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits format...\n");
    inputdata_cfitsio(cmd, gd, filename, ifile);
    break;
case INFITSRADECTC:
//    verb_print(cmd->verbose,
//               "\n\tInput in fits-healpix-radec-tc format...\n");
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits-healpix-radec-tc format...\n");
    inputdata_cfitsio_radec_tc(cmd, gd, filename, ifile);
    break;
case INFITSHEALPIX:
//    verb_print(cmd->verbose,
//               "\n\tInput in fits-healpix format...\n");
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits-healpix format...\n");
    inputdata_cfitsio_healpix(cmd, gd, filename, ifile);
    break;

case INNUMPYHEALPIX:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in numpy-healpix format...\n");
    inputdata_numpy_healpix(cmd, gd, filename, ifile);
    break;

#endif	// ! _cballsio_cfitsio_01_h
