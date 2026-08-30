// Use:
//#include "startrun_octree_smoothing_02.h"
//
//
// included in file: addons/addons_include/source/startrun/startrun_include_02.h
//

#ifndef _startrun_octree_smoothing_02_h
#define _startrun_octree_smoothing_02_h

//B
//    cmd->ntosave = GetlParam("ntosave");            // needed?
    cmd->scanLevel = GetiParam("scanLevel");
// Root nodes:
    cmd->scanLevelRoot = GetiParam("scanLevelRoot");
//cmd->scanLevelMin = GetParam("scanLevelMin");
cmd->scanLevelMin = GetiParam("scanLevelMin");
//E

#endif	// ! _startrun_octree_smoothing_02_h
