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
#define MAXLENGTHOFINDIVIDUALFILES        128
#define BUFFERSIZE              2256
#define MAXLENGTHOFFMTFILES     32
#define MAXLENGTHOFREAL         32
#define MAXNSLASHS              20

#define VERBOSENOINFO           0
#define VERBOSEMININFO          1
#define VERBOSENORMALINFO       2
#define VERBOSEDEBUGINFO        3

//B Debug tracking
#define TRACKONLYSEARCH
//#undef TRACKONLYSEARCH

#ifdef DEBUGTRACKING

#ifdef TRACKONLYSEARCH
#define debug_tracking_search(track_step)                              \
  verb_print_debug(1, "Track step (%s)\n", track_step);
#define debug_tracking_s_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %s\n", track_step, extra);
#define debug_tracking_r_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %g\n", track_step, extra);
#define debug_tracking_i_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %d\n", track_step, extra);

#define debug_tracking
#define debug_tracking_s
#define debug_tracking_r
#define debug_tracking_i

#else // ! TRACKONLYSEARCH

#define debug_tracking(track_step)                              \
  verb_print_debug(1, "Track step (%s)\n", track_step);
#define debug_tracking_s(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %s\n", track_step, extra);
#define debug_tracking_r(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %g\n", track_step, extra);
#define debug_tracking_i(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %d\n", track_step, extra);

#define debug_tracking_search(track_step)                              \
  verb_print_debug(1, "Track step (%s)\n", track_step);
#define debug_tracking_s_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %s\n", track_step, extra);
#define debug_tracking_r_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %g\n", track_step, extra);
#define debug_tracking_i_search(track_step, extra)                              \
  verb_print_debug(1, "Track step (%s): %d\n", track_step, extra);

#endif // ! TRACKONLYSEARCH

#else // ! DEBUGTRACKING :: dummies...

#define debug_tracking
#define debug_tracking_s
#define debug_tracking_r
#define debug_tracking_i

#define debug_tracking_search
#define debug_tracking_s_search
#define debug_tracking_r_search
#define debug_tracking_i_search

#endif // ! DEBUGTRACKING
//E

#endif // ! _common_defs_h

