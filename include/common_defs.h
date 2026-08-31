/*==============================================================================
 HEADER: cmdline_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "common_defs.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _common_defs_h
#define _common_defs_h

#include <errno.h>
#include <string.h>

/* octree-balls4-omp always needs the B4 pivot partition, but enabling the
 * addon must not change pivot selection in unrelated search methods. */
#if defined(BALLS4SCANLEV) || defined(OCTREEBALLS4OMP) || defined(OCTREEBALLS4MPI)
#define CBALLS_NEEDS_BALLS4_SCAN 1
#endif

//B Usefule MACROS and other constants:
#define CPUTIME         (second())                  // Gives cputime in seconds

#define PRNUNITOFTIMEUSED   "sec."
#define MAXLENGTHOFSTRSCMD     1024
#define EXTFILES            ".txt"
#define INMB                9.536743116E-7          // 1/(1024*1024)
#define SMALLESTDOUBLE 1.0E-300
#define SMALLESTRMPC 1.0E-00
#define BIGGESTDOUBLE 1.0E+300
// 180*60/Pi
#define RADTOARCMIN   3437.74677
//E
//B octree treeload definitions
//      choose well...
//      If not tdepth and number of cells will grow badly
#define EPSILONLOADBODY 1.0E-4
#define EPSILONFLOATLOADBODY 1.0E-5
//B this is the one working
#define EPSILONHACKCELL 1.0E-16
//E
//B testing this one...
//#define EPSILONHACKCELL 1.0E-12
//E
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


//B added by cBalls
//  be aware that this macro depends on routineName...
#define SET_TAG_NAME(paramtext)                                         \
  do {                                                                  \
    size_t tag_len = strlen(paramtext);                                 \
    if (tag_len >= sizeof(tag[nt])) {                                   \
      snprintf(errmsg, _ERRORMSGSIZE_,                                  \
               "%s: parameter tag '%s' too long\n",                     \
               routineName, paramtext);                                 \
      return FAILURE;                                                   \
    }                                                                   \
    memcpy(tag[nt], paramtext, tag_len + 1);                            \
  } while (0)
//E

//B I/O Macros
#define IPName(param,paramtext)                                 \
  {SET_TAG_NAME(paramtext);                                     \
  addr[nt]=&(param);                                            \
  id[nt++]=INT;}

#define LPName(param,paramtext)                                 \
  {SET_TAG_NAME(paramtext);                                     \
  addr[nt]=&(param);                                            \
  id[nt++]=LONG;}

#define RPName(param,paramtext)                                 \
  {SET_TAG_NAME(paramtext);                                     \
  addr[nt]=&param;                                              \
  id[nt++]=DOUBLE;}

#define BPName(param,paramtext)                                 \
  {SET_TAG_NAME(paramtext);                                     \
  addr[nt]=&param;                                              \
  id[nt++]=BOOLEAN;}

#define SPName(param,paramtext,n)                               \
  {SET_TAG_NAME(paramtext);                                     \
  param=(string) malloc(n);                                     \
  addr[nt]=param;                                               \
  str_size[nt]=(n);                                             \
  id[nt++]=STRING;}
//E I/O Macros

//B Debug tracking
// use: debug_tracking("xx");
#define debug_tracking(track_step)                                      \
  verb_print_debug(1, "Track step (%s)\n", track_step);
// use: debug_tracking_s("xx", string);
#define debug_tracking_s(track_step, extra)                             \
  verb_print_debug(1, "Track step (%s): %s\n", track_step, extra);
// use: debug_tracking_r("xx", real);
#define debug_tracking_r(track_step, extra)                             \
  verb_print_debug(1, "Track step (%s): %g\n", track_step, extra);
// use: debug_tracking_i("xx", integer);
#define debug_tracking_i(track_step, extra)                             \
  verb_print_debug(1, "Track step (%s): %d\n", track_step, extra);
//E

//B generic definitions outside CLASSLIB on
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef SUCCESS
#define SUCCESS 0
#endif

#ifndef FAILURE
#define FAILURE 1
#endif

#ifndef _ERRORMSGSIZE_
#define _ERRORMSGSIZE_ 2048
#endif

#ifndef _BASEPATHSIZE_
#define _BASEPATHSIZE_ 1000
#endif

#ifndef WLCF_ERRORMSG_DEFINED
#define WLCF_ERRORMSG_DEFINED
typedef char ErrorMsg[_ERRORMSGSIZE_];
#endif
//E

#include "cballs_class_call.h"

//B COSMO_FAIL must be equal to cBALLS_FAIL below
#ifdef CLASSLIB
#define COSMO_FAIL(cmd, ...)                                      \
    do {                                                          \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__); \
        return FAILURE;                                           \
    } while (0)

#define COSMO_FAIL_GOTO(cmd, label, ...)                              \
    do {                                                              \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__);  \
        goto label;                                                   \
    } while (0)
#else
#define COSMO_FAIL(cmd, ...) error(__VA_ARGS__)
#define COSMO_FAIL_GOTO(cmd, label, ...) error(__VA_ARGS__)
#endif
//E

#ifdef CLASSLIB
#ifndef cBALLS_FAIL
#define cBALLS_FAIL(cmd, ...)                                      \
    do {                                                          \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__); \
        return FAILURE;                                           \
    } while (0)
#endif
#else
#ifndef cBALLS_FAIL
#define cBALLS_FAIL(cmd, ...) error(__VA_ARGS__)
#endif
#endif

//B not working... check
#ifndef cBALLS_FAIL_GO
#define cBALLS_FAIL_GO(cmd, ...)                                        \
    do {                                                                \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__);    \
    } while (0)
#endif
//E

//B added by cBalls
#ifndef OPEN_OUTPUT_OR_FAIL
#define OPEN_OUTPUT_OR_FAIL(outstr, filename, mode)                    \
    do {                                                               \
        if (stropen_checked((filename), (mode), &(outstr),             \
                            cmd->error_message, _ERRORMSGSIZE_)        \
            == FAILURE)                                                \
            return FAILURE;                                            \
    } while (0)
#endif

#ifndef PRINT_OR_FAIL
#define PRINT_OR_FAIL(expr)        \
    do {                           \
        if ((expr) == FAILURE)     \
            return FAILURE;        \
    } while (0)
#endif
//E

#ifndef WRITE_OUTPUT_OR_FAIL
#define WRITE_OUTPUT_OR_FAIL(outstr, filename, ...)                         \
    do {                                                                    \
        if (fprintf((outstr), __VA_ARGS__) < 0) {                           \
            int _cballs_saved_errno = errno;                                \
            snprintf((cmd)->error_message, _ERRORMSGSIZE_,                  \
                     "%.256s: error writing '%.1024s': %.512s",              \
                     routineName, (filename),                               \
                     strerror(_cballs_saved_errno ? _cballs_saved_errno : EIO)); \
            if ((outstr) != NULL) {                                         \
                fclose(outstr);                                             \
                (outstr) = NULL;                                            \
            }                                                               \
            return FAILURE;                                                 \
        }                                                                   \
    } while (0)
#endif

#ifndef CLOSE_OUTPUT_OR_FAIL
#define CLOSE_OUTPUT_OR_FAIL(outstr, filename)                              \
    do {                                                                    \
        if (cballs_stream_close_checked(cmd, routineName, &(outstr),         \
                                        (filename)) == FAILURE)             \
            return FAILURE;                                                 \
    } while (0)
#endif

#endif // ! _common_defs_h

