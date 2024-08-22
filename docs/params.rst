
Parameters
========================

This section describes the various parameters cBalls needs for controlling
what the searching process do:

Parameters related to the searching method
-------------------------------------------

:searchMethod: (str or list) [alias: search]
    The searching method to use. Default is ``tree-omp-sincos``. Fastest method so far is ``balls-omp``.

    Use it as::

        searchMethod = balls-omp

    In command line version do not use spaces before and after ``=`` or it won't be parsed correctly. In a parameter file you have more liberty.

:mChebyshev: (positive int) [alias: mcheb]
    The number of multipoles to compute for a 3pcf computation. Number of them are: ``mChebyshev + 1``, because it includes the monopole.

    Use it as::

        mChebyshev = 7
    
    Default value is ``7``. As multipoles comes from a harmonic expansion, we may be interested in computing the whole 3pcf. This process involves a FFT, therefore give ``mChebyshev + 1`` as a power of 2.

Parameters to control the I/O file(s)
-------------------------------------

:infile: (str, default="") [alias: in] File names with points to analyse.

:infileformat: (str, default=columns-ascii) [alias: infmt] Data input files format (columns-ascii, binary or takahasi)

    The input columns for ``columns-ascii`` format:

    - ``x`` = x position of a point in the catalog.
    - ``y`` = y position of a point in the catalog.
    - ``z`` = z position of a point in the catalog.
    - ``kappa`` = The kappa value of the point.

Parameters to control histograms and their output files
-------------------------------------------------------

:useLogHist: (bool, default=true) Which type of binning should be used.
:logHistBinsPD: (float) The minimum separation to include in the histograms.
:rangeN: (float) The maximum separation to include in the histograms.
:rminHist: (float) The minimum separation to include in the histograms.
:sizeHistN: (int) The number of output bins to use.



Set of parameters needed to construct a test model
--------------------------------------------------

:seed: (int, default=123) Random number seed to test run or useful to change a random region in Takahasi simulations.

:testmodel: (str, default=simple-cubic-random) [alias: tstmodel] Test model name to analyse.

:nbody: (int, default=16348) Number of points to test.

:lengthBox: (float, default=10000) [alias: lbox] Length of the box to test.


Miscellaneous parameters
------------------------

:script: (str, default="") Scripts in shell or python that can be run in pre-processing or post-processing.

:stepState: (int, default=10000) number of steps to save a state-run info (pivot number completed in the log file).

:verbose: (int, default=1) [alias: verb] How verbose the code should be during processing.

    - 0 = no output unless there is an error
    - 1 = output warnings
    - 2 = output progress information
    - 3 = output extra debugging lines

:verbose_log: (int, default=1) [alias: verblog] To print messages to a log file ``cballs.log`` in directory ``tmp`` under output directory given by the parameter: ``rootDir`.

    Amount of message information is controlled by the int given.

:numberThreads: (int, default=4) [alias: nthreads] How many OpenMP threads should be used.

    It is needed to switch on OpenMP: ``OPENMPMACHINE = 1`` in ``Makefile_settings`` and recompile cBalls again.

:options: (str, default="") [alias: opt] You may give here various code behavior options.

    Use it as::

        options = str1,str2,str3,...

    where str# is one of the:

    - stop = stop execution before searching process
    - compute-HistN = compute NN encounters and save histogram in a file
    - and-CF = if you use ``compute-HistN`` then you may compute and save the correlation funcion of NN encounters (the equivalent to the radial distribution function in liquids).
    - no-one-ball = during the searching process does not use balls criterion to speed up the code


.. note::

    - It is not necessary to specify all the parameters. You need to give only the ones apropriate to the run. The rest of parameters will use their default values if they are OK with you.

    - When you specify the root output directory using: ``rootDir``, and this is a single directory that will be located in the pwd dir, then do not use ``./`` at the begining of the name or ``/`` at its end.


