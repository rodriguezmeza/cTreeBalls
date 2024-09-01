// Use:
//#include "input_iolib_01.h"

#ifndef _input_iolib_01_h
#define _input_iolib_01_h

class_call(parser_read_string(pfc,"columns",&string1,&flag1,errmsg),
           errmsg,errmsg);
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"columns") == 0){
          slen = strlen(pfc->value[index]);
          cmd->columns = (char*) malloc(slen*sizeof(char));
        memcpy(cmd->columns,pfc->value[index],slen+1);
        break;
      }
    }
}

#endif	// ! _input_iolib_01_h
