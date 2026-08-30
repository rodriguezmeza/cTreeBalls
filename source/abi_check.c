/*==============================================================================
 MODULE: abi_check.c				[wlcov]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	01.05.2026
 Purpose:
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// Application Binary Interface (ABI)
#ifdef CLASSLIB
#include "cballs.h"

size_t sizeof_cmdline_data(void)
{
    return sizeof(struct cmdline_data);
}

size_t sizeof_global_data(void)
{
    return sizeof(struct global_data);
}

#endif

