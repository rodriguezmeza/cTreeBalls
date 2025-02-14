// Use:
//#include "input_save_restore_01.h"

#ifndef _input_save_restore_01_h
#define _input_save_restore_01_h

/*
class_call(parser_read_string(pfc,"statefile",&string1,&flag1,errmsg),
           errmsg,errmsg);
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"statefile") == 0){
          cmd->statefile = strdup(pfc->value[index]);
        break;
      }
    }
}
*/

class_call(parser_read_string(pfc,"statefile",&string1,&flag1,errmsg),
           errmsg,errmsg);
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"statefile") == 0){
          slen = strlen(pfc->value[index]);
          cmd->statefile = (char*) malloc(slen*sizeof(char));
        memcpy(cmd->statefile,pfc->value[index],slen+1);
        break;
      }
    }
}

class_call(parser_read_string(pfc,"restorefile",&string1,&flag1,errmsg),
           errmsg,errmsg);
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"restorefile") == 0){
          cmd->restorefile = strdup(pfc->value[index]);
        break;
      }
    }
}


#endif	// ! _input_save_restore_01_h
