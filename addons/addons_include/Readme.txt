
AddOns include folder

Files of the cTreeBalls are:

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

Any possible AddOn item will be inserted in the corresponding to the files above that are located in the addams/addons_include folder.

Each of the files above has lines like:

#ifdef ADDONS
#include "xxx_include.h"
#endif

Look for them and think where your AddOn items will be inserted on. Then go to the corresponding  file ("xxx_include.h") and edit it.


Other files and folders are auxiliaries:

1. class_lib
2. doc
3. general_lib
4. getparam
5. gsl
6. Makefiles
7. patches
8. tests




