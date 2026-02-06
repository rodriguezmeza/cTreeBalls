't" t
.TH cballs 1 "April 2023" UNIX "LSST/CosmoININ PROJECT"
.na
.nh   

.SH NAME
cballs - (2,3)-Point Correlation Function computation code
.SH SYNOPSIS
\fBcballs\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ]
.sp

.SH DESCRIPTION
\fBcballs\fR - computes the 2 or/and 3 point correlation functions using tree/balls methods, complexity O(N Log N).

.SH OPTIONS
All the options have the structure
.sp
$ \fIoption_name\fR=<option_value>
.sp
except option \fB-help\fR. No space between \fIoption_name\fR, '=' and <option_value>
.sp

Options and their possible values are:

.IP "\fB-help\fR" 12
By writting

.sp
$ cballs --help
.sp

you will get the list of all parameters and their default values. An option may have an alias which is a short name of the option. If an option has an alias in the list above it comes after its description surrounded by brackets tagged with 'a:'. For example,

.sp
option_name=<option_value>	... Description ... [a: opt]
.sp
here 'opt' is the short name of the option. In the description of the options below, when an option has an alias it will be noted in the same way but before its description.

.IP "\fBparamfile\fR" 12
is the name file with the values of the input parameters. Overwrite parameters
values below. You may input this filename by only writing:
.sp
$ cballs paramfile=parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start with an "#" or "%". Follow each name option with a blank-space, "=", blank-space and the option value (note the difference when using command line version instead of parameter file version). The order of the option lines does not matter.  Also you may create an example input file by executing
.sp
$ cballs
.sp
This will run the \fBcballs\fR code with default values and when it finish you will have in your running output directory (by default is 'Output', see below in 'rootDir' line) the file "parameters_null-cballs-usedvalues". Now you may edit this file to adapt to your own run parameters. This file can be overwritten so it may be helpful to change this file name to whatever apropriate, and move it elsewhere.

.sp
Note that parameter alias are not allowed in a parameter file. Also behaviour of Mac and Linux plataforms are different when they reach the end of line in a text file. Then would be necessary to add a blank line at the end of the parameter file.

.IP "\fBsearchMethod\fR" 12
[a: search] is the searching method to use. Default is 'balls-omp'. Other search options will be added in the addons folder, like 'octree-ggg-omp' searching method.

.IP "\fBmChebyshev\fR" 12
[a: mcheb] is the total number of chebyshev multipoles to use. (mChebyshev + 1) must be a power of 2. Default value is 7, which means that will be computed 8 multipoles.

.IP "\fBnsmooth\fR" 12
[a: nsm] number of bodies to smooth out (or in a bucket if you are using 'searchMethod=kdtree-omp'). Default is 1, but if you are using 'kdtree-omp' it is convenient to use 'nsmooth=8'.

.IP "\fBrsmooth\fR" 12
[a: rsm] radius of the pivot smoothing neighbourhood. If empty a default is set.

.IP "\fBtheta\fR" 12
control tree search parameter, can be used to increase/decrease computation speed. It has the effect of increasing ('theta<1') or decreasing ('theta>1') the radius of cells. Default value is 1. Experiment with 'theta=1.1' to see the speeding-up of the code. If 'theta=0' is used then computation will be done as it were the direct-2-loops method with complexity O(N^2), useful for testing but not for productions computations.

.IP "\fBusePeriodic\fR" 12
[a: periodic] if false, do not use periodic boundary condition. Useful to compute 2-point correlations functions of periodic boxes steaming from N-body simulations or to analyze halo catalogs from Rockstar for example. \fBcballs\fR is able to read Gadget format or multi-columns format files like the ones given by Rockstar. See below.

.IP "\fBinfile\fR" 12
[a: in]  gives the name of file to input the data to analyse.

.IP "\fBinfileformat\fR" 12
[a: infmt]  gives format of file to input the data to analyse. Common options are: 'columns-ascii', 'binary', 'takahasi'. First one is a file with a header and all the data point info is in column form: positions and convergence field. The binary format is a \fBcballs\fR binary format. Whereas takahasi format is for reading Takahasi simulations data. Other file formats are in the addons. Gadget format is one of them. Other is 'multi-columns-ascii' which can be used in combinations with 'columns' option in order to read files that are Rockstar output. Also, \fBcballs\fR can read fits and healpix files. Note: you need to set properly the format of the point catalog that you are pretending to analyze, because if you mistakenly set wrong format or positions or fields columns, you will obtain wrong results or get a 'segmentation fault' error and the code will stop. Please follow the examples given in the tests folder.

.IP "\fBiCatalog\fR" 12
[a: icats] gives the indexes of point catalogs to analyse. Useful to read several point catalogs to do cross-correlations or to mix several files in one file or to analyse several input file as they were only one.

.sp
Note here that the number of items in \fBinfile\fR, \fBinfileformat\fR and \fBiCatalog\fR have to be the same.

.IP "\fBrootDir\fR" 12
[a: root] gives output dir, where output files will be written. If this folder does not exist it will be created. If it does exist it will be overwritten.

.IP "\fBoutfile\fR" 12
[a: o] will give the name of file to save the data to be analysed. The output is written in column form and a header by default, there is also a binary format. And by default this options is empty then no output file is saved. You can use in combinations with 'options=stop' to just convert input file format that can be use by other applications. It can be use in combination with \fBinfile\fR, \fBinfileformat\fR and \fBiCatalog\fR to read several files and save them as they were only one.

.sp
Note: to this file name, an extension ".txt" is added.

.IP "\fBoufileformat\fR" 12
[a: ofmt]  gives format of file to output the data to be analysed. Common options are: columns-ascii and binary. They follow the \fBinfileformat\fR above.

.sp
Note here that it is recommended to use binary format given that save hard disk space and it is fastest to read.

.IP "\fBthetaL\fR" 12
angle theta left side of the region in radians. To see how many arcmin it is, the factor of conversion is: RADTOARCMIN=3437.74677.

.IP "\fBthetaR\fR" 12
angle theta right side of the region in radians.

.IP "\fBphiL\fR" 12
angle phi left side of the region in radians.

.IP "\fBphiR\fR" 12
angle phi right side of the region in radians.

.sp
Note: phi and theta angles correspond to RA and DEC in the galactic coordinate system. And to convert them to a vector position just use the standard spherical transformation of coordinates.

.sp
Also note that these four parameters can be used as a filter or mask to access only a patch of the full sky region.

.IP "\fBuseLogHist\fR" 12
[a: loghist] if true, it will use log scale in the r-bins of histograms.

.IP "\fBlogHistBinsPD\fR" 12
[a: binspd] Bins per decades in r-bins when 'useLogHist=true'.

.IP "\fBsizeHistN\fR" 12
gives the number of bins to create.

.IP "\fBrangeN\fR" 12
is the maximum separtion of bodies to search for sampling. Default value is 0.0633205 that correspond to 217.68 arcmin.

.IP "\fBrminHist\fR" 12
is the minimum separtion of bodies to search for sampling. Default value is 0.00213811 that correspond to 7.35 arcmin.

.IP "\fBsizeHistPhi\fR" 12
array size for angular histogram.

.IP "\fBhistNNFileName\fR" 12
[a: histNNfname] file name (without extension) to save histograms of NN.

.IP "\fBhistXi2pcfFileName\fR" 12
[a: histXi2pcffname] file name (without extension) to save histograms of Xi2pcf.

.IP "\fBhistZetaFileName\fR" 12
[a: histZfname] prefix of file name to save histograms of matrix Zeta, for the multipoles.

.IP "\fBsuffixOutFiles\fR" 12
[a: suffix] suffix to add to output filenames.

.sp
Note: to all above file names, an extension ".txt" is added.


.IP "\fBseed\fR" 12
gives the seed to init the random number generator. Useful random number seed to test run or useful to change a random region in Takahashi simulations.

.IP "\fBtestmodel\fR" 12
[a: tstmodel] is the name of the model of data to create. By default this option is 'unit-sphere-random' which creates a uniformly random distribution of point over a unit spherical surface.

.IP "\fBnbody\fR" 12
is the total number of bodies to make as a test case. By default a \fBnbody\fR random distribution of points are created in a unit sphere. The random generator is initialize with a seed given by \fBseed\fR option.

.IP "\fBlengthBox\fR" 12
[a: lbox] is the length of the box side. The box is where the bodies reside. When a test simple-cubic or random-cubic models are used then the lenghts of the box are fixed with this parameter. When used with unit-sphere-random then box has lenght of 2.

.IP "\fBpreScript\fR" 12
scripts in shell or python that can be run in pre-processing. Use in combinations with 'options=pre-processing'. Last setting is necessary to activate the script.

.IP "\fBposScript\fR" 12
scripts in shell or python that can be run in post-processing. Use in combinations with 'options=post-processing'. Last setting is necessary to activate the script.

.IP "\fBstepState\fR" 12
gives the frequency to save the the number of pivot computed. Useful to see the state of the run.

.IP "\fBverbose\fR" 12
[a: verb] to print messages to standard output (stdout). There are four levels of information. When is "0", no information is written.

.IP "\fBverbose_log\fR" 12
[a: verblog] to save messages to a log file (in directory "tmp"). When is "0" no information is written.

.IP "\fBnumberThreads\fR" 12
[a: nthreads] To set the number of threads to use (OpenMP).

.IP "\fBcolumns\fR" 12
[a: cols] columns to use as vector position and fields when using multi-columns-ascii or fits formats of point catalog.

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options, like 'smooth-pivot', 'stop', 'pre-processing', 'post-processing', etcetera.

.SH EXAMPLES

Note: It is assumed that working directory is cTreeBalls, as was created when was git cloned. First create the executable, 'cballs':

.sp
$ make clean; make all

.sp
cd tests

.sp
$ time ../cballs nbody=100000 rangeN=0.06 sizeHistN=20 verb=2

.sp
You will have computations written in files in a directory named with option \fBrootDir\fR. Not given here then it is used default 'Output'. If this directory does not exist, it will be created.
Also will be created a directory named "tmp" in this output folder, where a log file will be saved.

.sp
Reading a Takahashi et al. (https://arxiv.org/abs/1706.01472) simulation data file downsized to 256:

.sp
$ time ../cballs infile=./catalogs/Takahasi/allskymap_nres08r081_zs9_mag.fits infmt=fits-healpix options=stop

.sp
To convert a catalog from healpix to ascii do:

.sp
$ ../cballs in=./catalogs/Takahasi/allskymap_nres08r081_zs9_mag.fits infmt=fits-healpix o=takas options=stop

.sp
Note here that \fBin\fR is the alias of \fBinfile\fR and \fBinfmt\fR is the alias of \fBinfileformat\fR.

.sp
In Output directory there will be a file named 'takas.txt', in this case, it is an ASCII file with a two lines header and then the several lines of columns with position '(x,y,z)' and convergence values. You can see its contents:

.sp
$ more Output/takas.txt


.sp
To compute correlations, \fBcballs\fR needs a catalog of bodies which contains at least their positions and values of the field whose correlations we are looking for. With all these values \fBcballs\fR builds up a structure of body info which consist in the general case of: Pos, Kappa, Weight, Mask, Id, Update, among other body information is needed to do the search and compute correlations we are interested in. The searching phase is done quite efficently by building a tree structure of the catalog data.

.sp
We finish with the following examples:

.sp
$ time ../cballs search=octree-ggg-omp infile=./catalogs/Takahasi/allskymap_nres08r081_zs9_mag.fits infmt=fits-healpix

.sp
Plot the results with:

.sp
$ python Xi3pcf_plot_flatten.py

.sp
Now do

.sp
$ time ../cballs search=octree-ggg-omp infile=./allskymap_nres10r081_zs9_mag.fits infmt=fits-healpix options=pre-processing,post-processing,no-normalize-HistZeta,edge-corrections,smooth-pivot preScript="ud_grade < in_ud_grade" posScript="python Xi3pcf_plot_flatten_EE.py; rm -f allskymap_nres10r081_zs9_mag.fits"

.sp
and the results were already plotted with the post-processing action in 'options'. Also the file 'allskymap_nres10r081_zs9_mag.fits' was [removed!

.sp
What does 'no-normalize-HistZeta,edge-corrections,smooth-pivot' mean? 'smooth-pivot' speeds up computation. 'edge-corrections' makes 3-point correlation function corrections if we are using a patch of the sky, in this case catalog given is full-sky, not needed but 'edge-correction' works just fine. 'no-normalize-HistZeta' does not normalize 3-point correlation function histograms as it is needed by 'edge-corrections' option. Note also that these options only work when using 'octree-ggg-omp' searching method.

.sp
Note: python version 3 was used with numpy and matplotlib modules installed. Also, for a long command line as above we may use '\' to split it shorter pieces.

.SH SEE ALSO
fkpt(1) mgpt(1) cballs(1)

.SH COPYRIGHT
Copyright (C) 2023--2026
.br
Mario A. Rodriguez-Meza
.br
