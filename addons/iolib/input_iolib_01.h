// Use:
//#include "input_iolib_01.h"

#ifndef _input_iolib_01_h
#define _input_iolib_01_h

class_call(parser_read_string(pfc,"columns",&string1,&flag1,errmsg),
           errmsg,errmsg);
gd->columnsFlag=FALSE;
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"columns") == 0){
          cmd->columns = (char*) malloc(MAXLENGTHOFSTRSCMD);
          strcpy(cmd->columns,pfc->value[index]);
          gd->columnsFlag=TRUE;
        break;
      }
    }
}

#endif	// ! _input_iolib_01_h
