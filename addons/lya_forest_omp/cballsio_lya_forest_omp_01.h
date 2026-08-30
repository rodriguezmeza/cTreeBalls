#ifndef _cballsio_lya_forest_omp_01_h
#define _cballsio_lya_forest_omp_01_h

case INLYAASCII:
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tInput in Lyman-alpha six-column ASCII format...\n");
    class_call_cballs(inputdata_lya_ascii(cmd, gd, filename, ifile),
                      errmsg, errmsg);
    break;

#endif
