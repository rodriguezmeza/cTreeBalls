// Use:
//#include "startrun_octree_smoothing_08.h"
//
//
// included in file: addons/addons_include/source/startrun/startrun_include_08.h
//

#ifndef _startrun_octree_smoothing_08_h
#define _startrun_octree_smoothing_08_h

WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"scanLevel",cmd->scanLevel);
// Root nodes:
WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"scanLevelRoot",cmd->scanLevelRoot);
WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI,"scanLevelMin",cmd->scanLevelMin);

#endif	// ! _startrun_octree_smoothing_08_h
