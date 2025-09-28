
cBalls: 1, 2, 3... for people in a hurry... boys and girls only want plug and play... SNI's under the siege of an evaluation year...

It is assuming that you are in the cTreeBalls directory (where there is a git clon of the code).

First:

$ make clean; make all

Then (after many seconds...)

$ mkdir deleteme; cd deleteme

Run

$ ../cballs

Then you will have an Output folder just created by cBalls filled with all the files resulted of the run. There is a file named "parameters_null-cballs-usedvalues"

$ more Output/parameters_null-cballs-usedvalues

and you will see all the parameters cBalls needs and their default values. A brief description can be found by

$ ../cballs --help

and an extended description is gotten by

$ man ../doc/man/cballs.1



1. Seeing all the parameters:

2. Reading files

3. To test cBalls:

$ ./run_all_tests

or test individual cases:

$ time ../cballs ./In/parameters_to-test_nside256_balls-omp

to plot results (vs reference results):

$ python script/mZetaM_test_kkk_all.py






