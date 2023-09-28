TPCF: 3 Point Correlation Function computation  {#mainpage}
===========================================================

Authors: Mario A. Rodriguez-Meza, et al.

with several major inputs from other people, especially xxx, etc.

For download and information, see https://github.com/cosmoinin


Compiling TPCF and getting started
-----------------------------------

Download the code by cloning it from
https://github.com/cosmoinin/tpcf. Go to the tpcf directory
(cd tpcf/ or src/) and compile (make clean;
make tpcf). You can usually speed up compilation with the option -j:
make -j tpcf. If the first compilation attempt fails, you may need to
open the Makefile and Makefile_machine files and adapt the name of 
the compiler (default: gcc),
of the optimization flag (default: -O4 -ffast-math) and of the OpenMP
flag (default: -fopenmp; this flag is facultative, you are free to
compile without OpenMP if you don't want parallel execution; note that
you need the version 4.2 or higher of gcc to be able to compile with
-fopenmp). Many more details on the TPCF compilation are given on the
wiki page

https://github.com/cosmoinin/tpcf/wiki/Installation

(in particular, for compiling on Mac >= 10.9 despite of the clang
incompatibility with OpenMP).

To check that the code runs, type:

    cd tests
    ../tpcf parameters.in

The parameters.in file is the reference input file, containing and
explaining the use of all possible input parameters.

The automatically-generated documentation will be located in

    doc/manual/html/index.html
    doc/manual/TPCF_manual.pdf

On top of that, if you wish to modify the code, you will find lots of
comments directly in the files.

For the moment you may consult man page:

    man doc/tpcf.m

or open with a browser the html file: doc/tpcf.html


Python
------

Under development python TPCF class for importing under python codes or
ipython notebooks. 

You will need to compile not only the code, but also its python wrapper. 
This can be done by typing just
'make' instead of 'make tpcf' (or for speeding up: 'make -j'). 

Plotting utility
----------------

Under development python plotting codes.

Using the code
--------------

You can use TPCF freely, provided that in your publications, you cite
at least the paper `TPCF: 3 Point Correlation Function schemes <http://arxiv.org/abs/xxxx.xxxx>`.

Support
-------

To get support please e-mail the authors.
