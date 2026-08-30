//=============================================================================
//        1          2          3          4        ^ 5          6          7
// Use:
//#include "octree_smoothing_00.h"

//
// included in file: addons/addons_include/source/tree/treeload_include00.h
//

#ifndef _octree_smoothing_00_h
#define _octree_smoothing_00_h

local void walktree_index_scan_lev(nodeptr, int, int, int);
local INTEGER inodelev;
local INTEGER ibodyleftout;

//B not needed in the public version... they are part of an addon
//local INTEGER Nc1;
//local INTEGER Nc2;
//E


//B Root nodes:
local void walktree_index_scan_lev_root(struct cmdline_data* cmd,
                                        struct  global_data* gd,
                                        nodeptr, int, int);
local INTEGER inodelev_root;
local INTEGER ibodyleftout_root;
//E

local void threadtree_smooth(struct  cmdline_data* cmd,
                             struct  global_data* gd,
                             nodeptr p, nodeptr n, int ifile);

//local INTEGER isel, inosel;
local int smoothBodies(struct  cmdline_data*,
                       struct  global_data* gd, bodyptr, INTEGER);

local int save_nodes(struct  cmdline_data*,
                     struct  global_data* gd, int ifile);
local int save_nodes_root(struct  cmdline_data*,
                          struct  global_data* gd, int ifile);

local char cellsfilePath[MAXLENGTHOFFILES];
local FILE *outcells;

#endif	// ! _octree_smoothing_00_h
