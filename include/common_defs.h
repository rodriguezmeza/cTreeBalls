/*==============================================================================
 HEADER: cmdline_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "common_defs.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4          5          6          7

#ifndef _common_defs_h
#define _common_defs_h

//B Usefule MACROS and other constants:
#define CPUTIME         (second())                  // Gives cputime in seconds

#define PRNUNITOFTIMEUSED   "sec."
#define MAXLENGTHOFSTRSCMD     1024
#define EXTFILES            ".txt"
#define INMB                9.536743116E-7          // 1/(1024*1024)
//E

#define MAXITEMS                100
#define MAXLENGTHOFFILES        1024
#define MAXNSLASHS              5

#define VERBOSENOINFO           0
#define VERBOSEMININFO          1
#define VERBOSENORMALINFO       2
#define VERBOSEDEBUGINFO        3

#endif // ! _common_defs_h

