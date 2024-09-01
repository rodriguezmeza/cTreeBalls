
AddOns include folder

Files of the cTreeBalls tree are:

1. include folder
	- cmdline_defs.h
	- datastruc_defs.h
	- globaldefs.h
	- protdefs.h
2. main
	- main.c
3. source
	- search_omp.c
	- startrun.c
	- testdata.c
	- tpcf_io.c
	- tpcf.c
	- treeload.c
	- treeutils.c
4. addons folder

Any possible AddOn item will be inserted into corresponding files above that are located in the addons/addons_include folder.

Each of the files above has lines like:

#ifdef ADDONS
#include "xxx_include_yyy.h"
#endif

Look for them and think where your AddOn items will be inserted on. Then go to the corresponding  file ("xxx_include.h") and edit it.


Other files and folders are auxiliaries:

1. docs
2. general_lib
3. getparam
4. gsl
5. Makefiles_xxx
6. tests




