// Use:
//#include "input_balls_omp_0357_01.h"

//
// included in file: addons/addons_include/source/addons/class_lib_include_01.h
//

#ifndef _input_balls_omp_0357_01_h
#define _input_balls_omp_0357_01_h
/*
class_call(parser_read_int(pfc,"scanLevel",
                           &param,&flag,errmsg),errmsg,errmsg);
if (flag == TRUE)
  cmd->scanLevel = param;
class_call(parser_read_int(pfc,"scanLevelRoot",
                           &param,&flag,errmsg),errmsg,errmsg);
if (flag == TRUE)
  cmd->scanLevelRoot = param;

class_call(parser_read_string(pfc,"scanLevelMin",&string1,&flag1,errmsg),
           errmsg,errmsg);
if (flag1 == TRUE) {
    for (index=0;index<pfc->size;++index){
      if (strcmp(pfc->name[index],"scanLevelMin") == 0){
          cmd->scanLevelMin = strdup(pfc->value[index]);
        break;
      }
    }
}

class_call(parser_read_int(pfc,"ntosave",
                           &param,&flag,errmsg),errmsg,errmsg);
if (flag == TRUE)
  cmd->ntosave = param;
*/

PARSER_READ(parser_read_int(pfc, "scanLevel", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevel = param;

PARSER_READ(parser_read_int(pfc, "scanLevelRoot", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevelRoot = param;

PARSER_READ(parser_read_int(pfc, "scanLevelMin", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevelMin = param;



#endif	// ! _input_balls_omp_0357_01_h
