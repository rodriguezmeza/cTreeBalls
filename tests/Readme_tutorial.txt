
cBalls: 1, 2, 3... for people in a hurry... boys and girls just want it plug and play... peladito y en la boca... SNI's under the siege of an evaluation year...

(
Some where in your computer $HOME do
$ git clone https://github.com/rodriguezmeza/cTreeBalls.git
)

1. It is assuming that you are in the cTreeBalls directory (where there is a git cloned version of the code).

First:

$ make clean; make all

Then (after too many seconds... a price must be paid to not look for the root administration and get him done some tasks for you...)

$ mkdir deleteme; cd deleteme

Run

$ ../cballs

A distribution of uniformly points in the unit sphere is generated with a gaussian values of convergence field. The (2,3)-point correlation functions are computed and saved as histograms in a folder name Output just created by cBalls. Also you will find a 
file named "parameters_null-cballs-usedvalues". Do

$ more Output/parameters_null-cballs-usedvalues

and you will see all the parameters cBalls needs and their default values. A brief description can be found by

$ ../cballs --help

and an extended description is gotten by

$ man ../doc/man/cballs.1

Read as thoroughly as needed in particular follow the examples at the end for reading files.

2. To test cBalls:

$ cd ../tests
$ ./run_all_tests

or test individual cases:

$ time ../cballs ./In/parameters_to-test_nside256_balls-omp

to plot results (vs reference results):

$ python script/mZetaM_test_kkk_all.py






