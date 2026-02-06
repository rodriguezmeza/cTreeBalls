
If addons for kdtree_cute_box is added:

Edit addons/Makefile_addons_settings and activate cute_box by doing:

KDTREECUTEBOXON = 1

at the end of the file.

1. Include searching method: kdtree-cute-box

If compiled successfully cTreeBalls:

$ cd ./tests/test_cute_box/
$ ../../cballs search=kdtree-cute-box options=params_cola_064.txt

Or choose any other parameter file there are there, params.txt.

Or you can go to the addons/kdtree_cute_box and compile CUTE_BOX:

$ make clean; make




