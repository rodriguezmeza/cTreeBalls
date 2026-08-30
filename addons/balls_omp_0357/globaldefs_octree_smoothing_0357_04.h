// Use:
//#include "globaldefs_octree_smoothing_04.h"

#ifndef _globaldefs_octree_smoothing_04_h
#define _globaldefs_octree_smoothing_04_h

//global bodyptr bodytab;
//B
#ifndef BODYTABBF_ON
#define BODYTABBF_ON 1
#endif

global bodyptr bodytabbf;
//global bodyptr bodytabbf;
//E
global bodyptr bodytabsm;
global bodyptr bodytabSel;

global nodeptr *nodetab;

// BALLS
global nodeptr *nodetabscanlev;
// Root nodes:
global nodeptr *nodetabscanlev_root;
//

//B Tree
//global cellptr root;
// BALLS
//global cellptr rootnode;                            // To make treenodes
global bodyptr nodetable;                           // To smooth minimum size
                                                    //  cells
global bodyptr nodetable_root;
//E


#endif	// ! _globaldefs_octree_smoothing_04_h
