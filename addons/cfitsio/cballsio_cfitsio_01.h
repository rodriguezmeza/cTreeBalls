
#ifndef _cballsio_cfitsio_01_h
#define _cballsio_cfitsio_01_h

case INCFITSIO:
    verb_print(cmd->verbose,
               "\n\tInput in fits format...\n");
    inputdata_cfitsio(cmd, gd, filename, ifile); 
    break;

#endif	// ! _cballsio_cfitsio_01_h
