
// included in: addons/addons_include/source/cballsio/cballsio_include_03.h

#ifndef _cballsio_iolib_03_h
#define _cballsio_iolib_03_h

case OUTCOLUMNSPOS:
    verb_print(cmd->verbose,
               "\n\tOutput in columns-ascii-pos format...\n");
    class_call_cballs(outputdata_ascii_pos(cmd, gd, btable, nbody),
                      errmsg, errmsg);
    break;

#endif	// ! _cballsio_iolib_03_h
