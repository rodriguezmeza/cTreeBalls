
// included in: addons/addons_include/source/cballsio/cballsio_include_03.h

#ifndef _cballsio_cfitsio_03_h
#define _cballsio_cfitsio_03_h

case OUTCFITSIO:
    verb_print(cmd->verbose,
               "\n\tOutput in fits format...\n");
    class_call_cballs(outputdata_cfitsio(cmd, gd, btable, nbody),
                      cmd->error_message, cmd->error_message);
    break;
case OUTNUMPYHEALPIX:
    verb_print(cmd->verbose,
               "\n\tOutput in numpy-healpix format...\n");
    class_call_cballs(outputdata_numpy_healpix(cmd, gd, btable, nbody),
                      cmd->error_message, cmd->error_message);
    break;

#endif	// ! _cballsio_cfitsio_03_h
