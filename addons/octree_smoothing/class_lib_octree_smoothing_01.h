// Use:
//#include "class_lib_octree_smoothing_01.h"

//
// included in file: addons/addons_include/source/addons/class_lib_include_01.h
//

#ifndef _class_lib_octree_smoothing_01_h
#define _class_lib_octree_smoothing_01_h

PARSER_READ(parser_read_int(pfc, "scanLevel", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevel = param;

PARSER_READ(parser_read_int(pfc, "scanLevelRoot", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevelRoot = param;

PARSER_READ(parser_read_int(pfc, "scanLevelMin", &param, &flag, errmsg));
if (flag == TRUE)
  cmd->scanLevelMin = param;

#endif	// ! _class_lib_octree_smoothing_01_h
