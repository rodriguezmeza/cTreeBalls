// Use:
//#include "startrun_iolib_03.h"

#ifndef _startrun_iolib_03_h
#define _startrun_iolib_03_h

if (!strnull(cmd->columns)) {
    if ((int)gd->columns[0] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[0] (%s, %d)\n",
                        cmd->columns, gd->columns[0]);
        gd->columns[0] = 1;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "CheckParameters: set it to: %d\n",
                               gd->columns[0]);
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[1] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[1] (%s, %d)\n",
                        cmd->columns, gd->columns[1]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[2] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[2] (%s, %d)\n",
                        cmd->columns, gd->columns[2]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[3] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[3] (%s, %d)\n",
                        cmd->columns, gd->columns[3]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[4] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[4] (%s, %d)\n",
                        cmd->columns, gd->columns[4]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[5] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[5] (%s, %d)\n",
                        cmd->columns, gd->columns[5]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
    if ((int)gd->columns[6] < 1) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "CheckParameters: absurd value for columns[6] (%s, %d)\n",
                        cmd->columns, gd->columns[6]);
        gd->columns[0] = 1;
        gd->columns[1] = 2;
        gd->columns[2] = 3;
        gd->columns[3] = 4;
        gd->columns[4] = 5;
        gd->columns[5] = 6;
        gd->columns[6] = 7;
        goto END;
    }
} else {
    gd->columns[0] = 1;
    gd->columns[1] = 2;
    gd->columns[2] = 3;
    gd->columns[3] = 4;
    gd->columns[4] = 5;
    gd->columns[5] = 6;
    gd->columns[6] = 7;
}
END:

#endif	// ! _startrun_iolib_03_h
