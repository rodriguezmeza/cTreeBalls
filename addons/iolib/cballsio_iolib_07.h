
//=============================================================================
//        1          2          3          4        ^ 5          6          7

// included in: addons/addons_include/source/cballsio/cballsio_indlude_11b.h

#ifndef _cballsio_iolib_07_h
#define _cballsio_iolib_07_h

// outfileformat: columns-ascii-pos
//  may be used in combination with options=cute-box-fmt
//  no header... the cute-box format input
//  if set pos: 0 < pos < lbox
local int outputdata_ascii_pos(struct cmdline_data* cmd,
                               struct  global_data* gd,
                               bodyptr bodytab, INTEGER nbody)
{
    string routineName = "outputdata_ascii_pos";
    char namebuf[256];
    stream outstr;
    bodyptr p;
    vector vecp;
    int k;

    if (scanopt(cmd->options, "cute-box-fmt")) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\t%s: saving data (box with sides ", routineName);
        DO_COORD(k)
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   " %g", gd->Box[k]);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                ") to file... cute-box-fmt\n");
    } else {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\t%s: saving data (box with sides ", routineName);
        DO_COORD(k)
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   " %g", gd->Box[k]);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                ") to file...\n");
    }

    if (format_checked(namebuf, sizeof(namebuf),
                       "namebuf", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;    //E
    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (scanopt(cmd->options, "cute-box-fmt")) {
            DO_COORD(k)
                vecp[k] = Pos(p)[k] + 0.5*gd->Box[k];

            if (out_vector_checked(outstr, vecp, routineName, namebuf,
                                   cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                if (outstr != NULL) {
                    fclose(outstr);
                    outstr = NULL;
                }
                return FAILURE;
            }
        } else {
            if (out_vector_checked(outstr, Pos(p), routineName, namebuf,
                                   cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                if (outstr != NULL) {
                    fclose(outstr);
                    outstr = NULL;
                }
                return FAILURE;
            }
        }        
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: data output to file %s\n",
                           routineName, namebuf);

    return SUCCESS;
}

#endif	// ! _cballsio_iolib_07_h
