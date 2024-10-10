't" t
.TH cballs 1 "April 2023" UNIX "LSST/CosmoININ PROJECT"
.na
.nh   

.SH NAME
cballs - 3 Point Correlation Function computation code
.SH SYNOPSIS
\fBcballs\fR [ \fIparameter_file_name\fR ] [ \fIoptions\fR ]
.sp

.SH DESCRIPTION
\fBcballs\fR - computes the 3 point correlation function using tree/balls methods, complexity O(N Log N).

.SH OPTIONS
All the options have the structure
.sp
\fIoption_name\fR=<option_value>
.sp
except option \fB-help\fR.
.sp

Options and their possible values are:

.IP "\fB-help\fR" 12
By writting

.sp
cballs --help
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
cballs paramfile=parameters_input_file_name
.sp
Parameter input file may be created by hand with the editor of your choice. Comment lines start with an "%". Follow each name option with a blank space and the option value. The order of the option lines does not matter.  Also you may create an example input file by executing
.sp
cballs
.sp
This will run the \fBcballs\fR code with default values and when it finish you will have in your running output directory the file "parameters_null-cballs-usedvalues". Now you may edit this file to adapt to your own run parameters. This file can be overwritten so it may be helpful to change this file name to whatever apropriate.

.IP "\fBsearchMethod\fR" 12
[a: search] is the searching method to use. Default is 'tree-omp-sincos'. Other search options will be added in the addons folder, like 'balls-omp' searching method.

.IP "\fBmChebyshev\fR" 12
[a: mcheb] is the total number of chebyshev multipoles to use. (mChebyshev + 1) must be a power of 2. Default value is 7, which means that will be computed 8 multipoles.

.IP "\fBnsmooth\fR" 12
[a: nsm] Number of bodies to smooth out (or in a bucket if you are using 'searchMethod=kdtree-omp'). Default is 1, but if you are using kdtree-omp it is convenient to use 'nsmooth=8'.

.IP "\fBrsmooth\fR" 12
[a: rsm] Radius of the pivot smoothing neighbourhood.

.IP "\fBtheta\fR" 12
Control tree search parameter, can be used to increase speed. It has the effect of increasing ('theta<1') or decreasing ('theta>1') the radius of cells. Default value is 1. Experiment with 'theta=1.1' to see the speeding-up of the code.

.IP "\fBcomputeTPCF\fR" 12
[a: tpcf] If true, compute 3pcf.

.IP "\fBusePeriodic\fR" 12
[a: tpcf] If false, do not use periodic boundary condition. Useful to compute 2-point correlations functions of periodic boxes steaming from N-body simulations or to analyze halo catalogs from Rockstar for example. In such a case use in combinations with 'computeTPCF=false'. \fBcballs\fR is able to read Gadget format or multi-columns format files like the ones given by Rockstar. See below.

.IP "\fBinfile\fR" 12
[a: in]  gives the name of file to input the data to analyse.

.IP "\fBinfileformat\fR" 12
[a: infmt]  gives format of file to input the data to analyse. Common options are: 'columns-ascii', 'binary', 'takahasi'. First one is a file with a header and all the data point info is in column form: positions and convergence field. The binary format is a \fBcballs\fR binary format. Whereas takahasi format is for reading Takahasi simulations data. Other file formats are in the addons. Gadget format is ones of them. Other is 'multi-columns-ascii' which can be used in combinations with 'columns' option in order to read files that are Rockstar output. Also, \fBcballs\fR can read fits files. You have to activate this options and recompile the code. Note that in this case you need to edit Makefile_machine to set cfitsio path. Note: you need to set properly the format of the point catalog that you are pretending to analyze, because if you mistakenly set wrong format or positions or fields columns, you will obtain wrong results or get a 'segmentation fault' error and the code will stop. Please follow the examples given in the tests folder.

.IP "\fBiCatalog\fR" 12
[a: icats] gives the indexes of point catalogs to analyse. Useful to read several point catalogs to do cross-correlations.

.IP "\fBrootDir\fR" 12
[a: root] gives output dir, where output files will be written.

.IP "\fBoutfile\fR" 12
[a: o] will give the name of file to save the analysed data. The output is written in column form by default, there is also a binary format. And by default this options is empty then no output file is saved. You can use in combinations with 'options=stop' to just convert input file format that can be use by other applications.

.IP "\fBoufileformat\fR" 12
[a: ofmt]  gives format of file to output the data to be analysed. Common options are: columns-ascii, binary. They follow the \fBinfileformat\fR above.

.IP "\fBuseLogHist\fR" 12
[a: loghist] If true, it will use log scale in the r-bins of histograms.

.IP "\fBlogHistBinsPD\fR" 12
[a: binspd] Bins per decades in r-bins when 'useLogHist=true'.

.IP "\fBsizeHistN\fR" 12
gives the number of bins to create.

.IP "\fBrangeN\fR" 12
is the maximum separtion of bodies to search for sampling.

.IP "\fBrminHist\fR" 12
is the minimum separtion of bodies to search for sampling.

.IP "\fBseed\fR" 12
gives the seed to init the random generator.

.IP "\fBtestmodel\fR" 12
[a: tstmodel] is the name of the model of data to create. By default creates a cubic box of random points. Other option is 'unit-sphere-random' which creats a uniformly random distribution of point over a unit spherical surface.

.IP "\fBnbody\fR" 12
is the total number of bodies to make as a test case. By default a box of size given by \fBlengthBox\fR is created and there \fBnbody\fR random points are generated. The random generator is initialize with a seed given by \fBseed\fR option.

.IP "\fBlengthBox\fR" 12
[a: lbox] is the length of the box side. The box is where the bodies reside.

.IP "\fBscript\fR" 12
    Scripts in shell or python that can be run in pre-processing or post-processing. Use in combinations with 'options=pre-processing' or 'options=post-processing'. Last setting is necessary to activate the script.

.IP "\fBstepState\fR" 12
gives the frequency to save the the number of pivot computed. Useful to see the state of the run.

.IP "\fBverbose\fR" 12
[a: verb] to print messages to stdout. There are three levels of information. When is "0", no information is written.

.IP "\fBverbose_log\fR" 12
[a: verblog] to print messages to a log file (in directory "tmp"). When is "0" no information is written.

.IP "\fBnumberThreads\fR" 12
[a: nthreads] you may give number of threads to compute in parallel.

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options.

.SH EXAMPLES
time ../cballs nbody=100000 rangeN=500 sizeHistN=100 lbox=10000 verb=2

.sp
You will have computation written in a file named with option \fBoutfile\fR in a directory named Output. If this directory does not exist, it will be created.
Also will be created a directory named "tmp" in the folder Output, where a log file will be saved.

.sp
Reading a Takahasi simulation data file downsized to 512:

.sp
time cballs infile=./full_sky_whole_XYZK__zs9_r000_nside512.bin infmt=binary options=stop

.sp
To convert a catalog from binary to ascii do:

.sp
cballs in=./full_sky_whole_XYZK__zs9_r000_nside512.bin infmt=binary o=takas options=stop

.SH SEE ALSO
nplot2d(1) cballs(1)

.SH COPYRIGHT
Copyright (C) 2023
.br
Mario A. Rodriguez-Meza
.br
