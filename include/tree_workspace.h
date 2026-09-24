#ifndef CBALLS_TREE_WORKSPACE_H
#define CBALLS_TREE_WORKSPACE_H
/* Private native-octree construction counters and diagnostic scratch.
 * Sizes preserve the original MAXLEVEL/NbMax numerical contracts. */
#define CBALLS_TREE_MAXLEVEL 32
#define CBALLS_TREE_RADIUS_BINS 33
typedef struct {
    int cellhist[CBALLS_TREE_MAXLEVEL];
    int subnhist[CBALLS_TREE_MAXLEVEL];
    INTEGER NTOT[1];
    INTEGER ip;
    INTEGER cellhistNb[CBALLS_TREE_RADIUS_BINS];
    int cellRadius[CBALLS_TREE_RADIUS_BINS];
    real deltaRadius;
    INTEGER inode;
    char treeinfofilePath[MAXLENGTHOFFILES];
    FILE * outtreeinfo;
    INTEGER inodelevB4;
    INTEGER ibodyleftoutB4;
    INTEGER ncell;
    INTEGER isel;
    INTEGER inosel;
} cballs_tree_workspace;
#endif
