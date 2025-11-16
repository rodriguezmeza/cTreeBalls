
// included in: addons/addons_include/source/cballsio/cballsio_indlude_01.h

#ifndef _cballsio_iolib_01_h
#define _cballsio_iolib_01_h

case INCOLUMNSPOS:
    verb_print(cmd->verbose,
               "\n\tInput in columns only positions (ascii) format...\n");
    inputdata_ascii_pos(cmd, gd, filename, ifile); break;

#if NDIM == 3
        case INCOLUMNS2DTO3D:
            verb_print(cmd->verbose,
                       "\n\tInput in columns (ascii-2d-to-3d) format...\n");
            inputdata_ascii_2d_to_3d(cmd, gd, filename, ifile); break;
#endif

case INMCOLUMNS:
    printf("\n\tInput in multiple ascii columns format...\n");
    inputdata_ascii_mcolumns(cmd, gd, filename, ifile); break;
case INRADECCOLUMNS:
    printf("\n\tInput in ra-dec ascii columns format...\n");
inputdata_ascii_ra_dec(cmd, gd, filename, ifile); break;

#endif	// ! _cballsio_iolib_01_h
