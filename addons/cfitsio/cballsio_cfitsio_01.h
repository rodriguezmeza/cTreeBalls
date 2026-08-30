
// included in: addons/addons_include/source/cballsio/cballsio_indlude_01.h

#ifndef _cballsio_cfitsio_01_h
#define _cballsio_cfitsio_01_h
case INCFITSIO:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits format...\n");
    class_call_cballs(inputdata_cfitsio(cmd, gd, filename, ifile), errmsg, errmsg);
    break;
case INFITSRADECFIELD:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits-radec-field format...\n");
    class_call_cballs(inputdata_cfitsio_radec_field(cmd, gd, filename, ifile), errmsg, errmsg);
    break;
case INFITSRADECRFIELD:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits-radecr-field format...\n");
    class_call_cballs(inputdata_cfitsio_radecr_field(cmd, gd, filename, ifile), errmsg, errmsg);
    break;
case INFITSHEALPIX:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in fits-healpix format...\n");
    class_call_cballs(inputdata_cfitsio_healpix(cmd, gd, filename, ifile), errmsg, errmsg);
    break;

case INNUMPYHEALPIX:
    verb_print_debug_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                          "\n\tInput in numpy-healpix format...\n");
    class_call_cballs(inputdata_numpy_healpix(cmd, gd, filename, ifile), errmsg, errmsg);
    break;

#endif	// ! _cballsio_cfitsio_01_h
