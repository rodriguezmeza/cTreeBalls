cTreeBalls: Correlation functions computation with Tree/Balls methods

Author: Mario A. Rodriguez-Meza

For download and information, see https://github.com/rodriguezmeza/cTreeBalls

Introduction
------------

Correlation function computation with Tree/Balls methods (short-name cBalls) is a C code for computing correlation functions using tree and balls methods. So far can compute 2 points correlation function (2pcf) and 3 points correlation function for scalar fields like weak lensing convergence.


Compiling and getting started
-----------------------------

Download the code by cloning it from https://github.com/rodriguezmeza/cTreeBalls 

Dependencies: cBalls optionally needs gsl version 2.7.1 and cfitsio version 4.4.1 installed in your system. Go to its web page `GSL <https://www.gnu.org/software/gsl/>`_ for details or ask to your system administrator. Make necessary changes in ``Makefile_machine`` file and look up for ``GSL``. Can be switched off: USEGSL = 0 in ``Makefile_settings``. I/O cfitsio library is set it OFF. Can be set it ON in ``addons/Makefile_settings``.

Go to the ``cTreeBalls`` directory (``cd cTreeBalls``) and compile (``make clean; make``). If the first compilation attempt fails, you may need to open the Makefile_machine file and adapt the name of the compiler (default: gcc), of the optimization flag (default: ``-O4 -ffast-math``) and of the OpenMP flag (default: -fopenmp; this flag is facultative, you are free to compile without OpenMP if you don't want parallel execution; note that you need the version 4.2 or higher of gcc to be able to compile with ``-fopenmp``). The code has been tested with gcc version 10 and would be working with version 11, 12. (In particular, for compiling on Mac >= 10.9 despite of the clang incompatibility with OpenMP).

To check that the code runs, type::

    cd tests
    ../cballs parameters_explained

The parameters_explained file is a reference input file, containing and explaining the use of all possible input parameters.

On top of that, if you wish to modify the code, you will find comments directly in the files in the folder ``addons``, and the modules you may add must go in this folder.

For the moment you may consult man page::

    man ../docs/cballs.m

or open with a browser the html file: docs/man/cballs.html

See more details about parameters needed by cBalls below (:ref:`PARAMETERS`).

Configuration
-------------

cBalls can be configured by switching on/off several options. Configuration file is ``Makefile_setting``. Search method ``balls-omp`` in folder addons/balls_omp: ``Makefile_settings_balls_omp`` has other options that can be switched on/off.


.. _PARAMETERS:

Parameters
----------

The list of available command line parameters can be consulted using the ``--help`` flag::

    ../cballs --help

or its short version::

    ../cballs -h


Python
------

To install cBalls python module (cballys) go to directory 'python' and execute::

    python setup.py build
    python setup.py install --user

Note: make sure you create cBalls library: in directory cTreeBalls execute::

    make clean; make all

To test it go to directory tests and 'run'::

    python test_cython_balls.py

Note: this interface in Cython was tested in a python environment with python3.7.


Plotting utilities
------------------

Several Jupyter notebooks, written by Abraham Arvizu and Eladio Moreno, are available to process cBalls results. They are in the github repository: 

https://github.com/joar-cafe/CBalls_plots/tree/main/benchmarks

or you can find more scripts in ``tests`` folder.

License
-------

cBalls is written by Mario A. Rodriguez-Meza, and is distributed under the `MIT license <https://github.com/rodriguezmeza/cTreeBalls/blob/main/LICENSE>`_. If you use this program in research work that results in publications, please cite the following paper:

Abraham Arvizu et al. 202X, [arXiv:xxxx.xxxxx](https://ui.adsabs.harvard.edu/abs/202XarXivxxxxxxxxxX/abstract)


Acknowledgements
----------------

cBalls use/is based on the following codes or projects:

* `Zeno <https://home.ifa.hawaii.edu/users/barnes/zeno/index.html>`_
* `Gadget-2 <https://wwwmpa.mpa-garching.mpg.de/gadget/>`_
* `CUTE <https://github.com/damonge/CUTE>`_
* `Numerical recipies <https://numerical.recipes/>`_
* `GSL <https://www.gnu.org/software/gsl/>`_
* `CLASS <https://github.com/lesgourg/class_public>`_
* `CFITSIO <https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_

