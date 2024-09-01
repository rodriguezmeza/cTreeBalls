
#ifndef _cballsio_cfitsio_03_h
#define _cballsio_cfitsio_03_h

case OUTCFITSIO:
    verb_print(cmd->verbose,
               "\n\tOutput in fits format...\n");
    outputdata_cfitsio(cmd, gd, btable, nbody);
    break;

#endif	// ! _cballsio_cfitsio_03_h
