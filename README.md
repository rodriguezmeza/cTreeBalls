cTreeBalls: 3 Point Correlation Function computation  {#mainpage}
===========================================================

Author: Mario A. Rodriguez-Meza

with several major inputs from other people, especially xxx, etc.

For download and information, see https://github.com/rodriguezmeza


Compiling cballs and getting started
-----------------------------------

Download the code by cloning it from
https://github.com/rodriguezmeza/cTreeBalls. Go to the cTreeBalls directory
(cd cTreeBalls/) and compile (make clean;
make). If the first compilation attempt fails, you may need to
open the Makefile_machine file and adapt the name of 
the compiler (default: gcc), 
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). The code has been tested with gcc version 10 and would be working with
version 11, 12.

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To check that the code runs, type:

    cd tests
    ../cballs parameters_test_T512_bin

The parameters_test_T512_bin file is a reference input file, containing (and
explaining in the near future) the use of all possible input parameters.

The automatically-generated documentation will be located in

    doc/manual/html/index.html
    doc/manual/cballs_manual.pdf

On top of that, if you wish to modify the code, you will find
comments directly in the files and the modules you may add must go in the folder "addons".

For the moment you may consult man page:

    man doc/cballs.m

or open with a browser the html file: doc/tpcf.html


Python
------

Under development python cballs for importing under python codes or
ipython notebooks. 

You will need to compile not only the code, but also its python wrapper. 
This will be done by typing just
'make all' instead of 'make' (or for speeding up: 'make -j'). 

Plotting utility
----------------

Under development python plotting codes.

Using the code
--------------

You can use cTreeBalls freely, provided that in your publications, you cite
at least the paper `cTreeBalls: 3 Point Correlation Function schemes <http://arxiv.org/abs/xxxx.xxxx>`.

Support
-------

To get support please e-mail the authors.
