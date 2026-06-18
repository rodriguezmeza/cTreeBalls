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

By default **cBalls** reads/writes catalog of points to analyzed files in 4-column format with x, y, z columns first and then value of the convergence field. It has a two line header::

    # nbody NDIM Lx Ly Lz
    # nbody-value NDIM-value Lx Ly Lz - values

Try running::

    ../cballs nbody=6480 o=points_on_sphere testmodel=unit-sphere-random options=stop

In the ``Output`` directory you will have a file: ``points_on_sphere.txt``. View its contents to see the two lines header and the 4 columns structure of data.

Note: in the above example ``points_on_sphere`` was not given an extension. By default cBalls gives to the output files the extension ``.txt``.

On top of that, if you wish to modify the code, you will find comments directly in the files in the folder ``addons``, and the modules you may add must go in this folder. See one of the most simple in ``addons/direct_method``.

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

cBalls python module (cballys) will be installed if you use::

    make clean; make all

Note: if problems, you can go to python directory and execute::

    python setup.py build
    python setup.py install --user

To test it go to directory tests and 'run'::

    python test_cython_balls.py

Note: this interface in Cython was tested in a python environment with ``python3.7``.


Plotting utilities
------------------

Several Jupyter notebooks, written by Abraham Arvizu and Eladio Moreno, are available to process cBalls results. They are in the github repository: 

https://github.com/joar-cafe/CBalls_plots/tree/main/benchmarks

or you can find more scripts in ``tests`` folder.

License
-------

cBalls is written by Mario A. Rodriguez-Meza, and is distributed under the `MIT license <https://github.com/rodriguezmeza/cTreeBalls/blob/main/LICENSE>`_. If you use this program in research work that results in publications, please cite the following paper:

`Abraham Arvizu et al., arXiv:2048.16847 <https://arxiv.org/abs/2408.16847>`_


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
* `HEALPix <https://healpix.sourceforge.io>`_

Also author acknowledges for helpful discussion and testing to the following people:

* Abraham Arvizu
* Alejandro Aviles
* Juan Carlos Hidalgo
* Eladio Moreno
* Gustavo Niz
* Axel Romero Tisnado
* Sofia Samario
