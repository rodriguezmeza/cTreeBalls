
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
    string routineName = "searchcalc_octree_box_omp";
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

    sprintf(namebuf, gd->fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (scanopt(cmd->options, "cute-box-fmt")) {
            DO_COORD(k)
                vecp[k] = Pos(p)[k] + 0.5*gd->Box[k];
            out_vector(outstr, vecp);
        } else {
            out_vector(outstr, Pos(p));
        }
    }
    fclose(outstr);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: data output to file %s\n",
                           routineName, namebuf);

    return SUCCESS;
}

#endif	// ! _cballsio_iolib_07_h
