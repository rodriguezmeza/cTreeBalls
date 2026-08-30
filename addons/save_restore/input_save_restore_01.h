// Use:
//#include "input_save_restore_01.h"

#ifndef _input_save_restore_01_h
#define _input_save_restore_01_h

gd->statefileFlag = FALSE;
gd->restorefileFlag = FALSE;

PARSER_READ(parser_read_string(pfc, "statefile", &string1, &flag1, errmsg));

if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"statefile") == 0){
          /* statefile */
          cmd->statefile = copy_param_string(pfc->value[index]);
          if (cmd->statefile == NULL) {
              BASE_FREE_STRINGS_ON_FAILURE();
              return FAILURE;
          }
          gd->statefileFlag = TRUE;

        break;
      }
    }
}

PARSER_READ(parser_read_string(pfc, "restorefile", &string1, &flag1, errmsg));

if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"restorefile") == 0){
          /* restorefile */
          cmd->restorefile = copy_param_string(pfc->value[index]);
          if (cmd->restorefile == NULL) {
              BASE_FREE_STRINGS_ON_FAILURE();
              return FAILURE;
          }
          gd->restorefileFlag = TRUE;

        break;
      }
    }
}


//E

#endif	// ! _input_save_restore_01_h
