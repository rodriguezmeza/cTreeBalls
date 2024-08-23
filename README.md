# cTreeBalls: Correlation functions computation with Tree/Balls methods

Author: Mario A. Rodriguez-Meza

## Table of Contents

-   [Introduction](#introduction)
-   [Compiling and getting started](#compiling-and-getting-started)
-   [Configuration](#configuration)
-   [Parameters](#parameters)
-   [Python](#python)
-   [License](#license)
-   [Acknowledgements](#acknowledgements)

## Introduction

**c**orrelations function computation with **T**ree/**B**alls methods (short-name **cBalls**) is a C code for computing correlation functions using tree and balls methods. So far can compute 2-point correlation function (2pcf) and 3-point correlation function (3pcf) for counts and scalar fields like convergence, the relevant scalar field in weak lensing.

Complete documentation can be found here: [documentation](docs/_build/html/index.html).

## Compiling and getting started

Download the code by cloning it from https://github.com/rodriguezmeza/cTreeBalls.

Dependencies: cBalls needs gsl version 2.7.1 installed in your system. Go to web page https://www.gnu.org/software/gsl/ for details or ask to your system administrator. Make necessary changes in `Makefile_machine` file and look up for `GSL`. But in case you have problems installing GSL just set it off in `Makefile_settings`, the basic distribution will run okey.

Go to the cTreeBalls directory (`cd cTreeBalls/`) and compile (`make clean; make`). If the first compilation attempt fails, you may need to open the Makefile_machine file and adapt the name of the compiler (default: gcc), of the optimization flag (default: `-O4 -ffast-math`) and of the OpenMP flag (default: `-fopenmp`; this flag is facultative, you are free to compile without OpenMP if you don't want parallel execution; note that you need the version 4.2 or higher of gcc to be able to compile with `-fopenmp`). The code has been tested with gcc version 10 and would be working with version 11, 12. (In particular, for compiling on Mac >= 10.9 despite of the clang incompatibility with OpenMP).

To check that the code runs, if you are in `cTreeBalls` directory, type:

    make clean; make all
    cd tests
    ../cballs

It will run using all default values and a directory named `Output` will be created under `tests`. **cBalls** will save all histograms files and a log file in `Output/tmp`. A file with the parameter values use in the run named `parameters_null-usedvalues` will also be saved. You may use it as a template to create your own parameter files.

If you execute:

    ../cballs options=post-processing script="python scripts/plot2pcf.py"

will do the same but now will plot the 2pcf and save it as a pdf file. Now, let us use a parameter file, execute:

    ../cballs parameters_explained

Directory `Output` is already created then **cBalls** will overwritte all histograms files, parameter file as was run, and the log file. The parameters_explained file is a reference input file, containing (and explaining) the use of all possible input parameters. In this case as can be seen in `parameters_explained` file a catalog of points is read from the file `kappa_nres12_zs9NS256r000.bin`.

To see a plot of the 2pcf, edit parameters_explained and set options to "post-processing", script="python scripts/plot2pcf.py" and execute again:
 
    ../cballs parameters_explained

At the end of the run you will, as before, have a pdf file of the plot.

You may also consult the code´s man page for more detailed information on how to run **cBalls**:

    man ../docs/man/cballs.m


## Configuration

cBalls can be configured by switching on/off several options. Configuration file is `Makefile_setting`.

| Option         | Description                                                                                                                                                   |
|:--------------:|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `DEFDIMENSION`       | `= 3` select dimension of the run: 2 or 3                                                                                                                     |
| `USEGSL`       | `= 0` switch on/off computation using GSL routines. Optional.                                                                                                                      |
| `GSLINTERNAL`     | `= 1` for enabling GSL internal sources<br />(if `= 0`, specify the corresponding compiler flags in `Makefile_machine` file)                                                             |
| `OPENMPMACHINE`     | `= 1` for enabling OpenMP parallelism<br />(Specify the corresponding compiler flag in `Makefile_machine` file)                                                             |
| `SINGLEPON`    | `= 0` for disabling single precision                                                                                                                             |
| `LONGINTON`    | `= 1` for enabling long integers                                                                                                                             |
| `ADDONSON`  | `= 1` for adding more funcionality to the code, like other searching methods, other catalog formats                                                                                                |

**Note**:
After changing `Makefile_settings` in order to have the new settings active in **cBalls** you have to re-compile the code: `make clean; make`. 

## Parameters

The list of available command line parameters can be consulted using the `-h` or `--help` flags:

    ../cballs --help

See also the man page as explained above.


## Python

To install cBalls python module (cballys) go to directory `python` and execute:

    python setup.py build
    python setup.py install --user

Note: make sure you create before cBalls library:

    cd ../
    make clean; make all

To test it go to directory `tests` and run:

    python test_cython_balls.py


## Plotting utilities

Several Jupyter notebooks, written by Abraham Arvizu and Eladio Moreno, are available to process cBalls results. They are in the github repository: 

https://github.com/joar-cafe/CBalls_plots/tree/main/benchmarks

Other python scripts are in directory 'tests/scripts'


## License

**cBalls** is written by Mario A. Rodriguez-Meza, is open source and distributed under the [MIT license](LICENSE.txt). If you use this program in research work that results in publications, please cite the following paper:

Mario A. Rodriguez-Meza et al. 202X, [arXiv:xxxx.xxxxx](https://ui.adsabs.harvard.edu/abs/202XarXivxxxxxxxxxX/abstract)

## Acknowledgements

cBalls use/is based on the following codes or projects:
-   [Zeno](https://home.ifa.hawaii.edu/users/barnes/zeno/index.html)
-   [Gadget-2](https://wwwmpa.mpa-garching.mpg.de/gadget/)
-   [CUTE](https://github.com/damonge/CUTE)
-   [Numerical recipies](https://numerical.recipes/)
-   [GSL](https://www.gnu.org/software/gsl/)
-   [CLASS](https://github.com/lesgourg/class_public)

