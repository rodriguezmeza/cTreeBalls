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
\fBcballs\fR - computes the 3 point correlation function using brute force, complexity O(N^2).

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
cballs -help
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
This will run the \fBcballs\fR code with default values and when it finish you will have in your running directory the file "parameters_null-cballs-usedvalues". Now you may edit this file to adapt to your own plotting parameters. This file can be overwritten so it may be helpful to change this file name to whatever apropriate.

.IP "\fBsearchMethod\fR" 12
[a: search] is the searching method to use. Default is "tree". Other option is to use brute search method: "direct". When the code is compiled with option OPENMPMACHINE=1, then the best options is "treeomp". When using "direct" option, automatically is run the OpenMP version of it.

.IP "\fBinfile\fR" 12
[a: in] will give the name of file to input the analysed data.

.IP "\fBtestmodel\fR" 12
[a: tstmodel] is the name of the model of data to create.

.IP "\fBmChebyshev\fR" 12
[a: mcheb] is the total number of chebyshev multipoles to use.

.IP "\fBoutfile\fR" 12
[a: o] will give the name of file to save the analysed data. The output is written in column form by default.

.IP "\fBnbody\fR" 12
is the total number of bodies to make as a test case. By default a box of size given by \fBlengthBox\fR is created and there \fBnbody\fR random points are generated. The random generator is initialize with a seed given by \fBseed\fR option.

.IP "\fBseed\fR" 12
gives the seed to init the random generator.

.IP "\fBlengthBox\fR" 12
[a: lbox] is the length of the box side. The box is where the bodies reside.

.IP "\fBsizeHistN\fR" 12
gives the number of bins to create.

.IP "\fBrangeN\fR" 12
is the maximum separtion of bodies to search for sampling.

.IP "\fBstatefile\fR" 12
[a: state] is the name of file to save the state of the run. The output is written in binary format. The frequency is given by option \fBstepState\fR.

.IP "\fBstepState\fR" 12
gives the frequency to save the whole state of the run that is save in a file given by \fBstatefile\fR in the "root" directory (usually named "Outputs").

.IP "\fBrestorefile\fR" 12
[a: restore] is the name of file to restore the state of the run.

.IP "\fBverbose\fR" 12
[a: verb] to print messages to stdout. There are three levels of information. When is "0", no information is written.

.IP "\fBverbose_log\fR" 12
[a: verblog] to print messages to a log file (in directory "tmp"). When is "0" no information is written.

.IP "\fBoptions\fR" 12
[a: opt] you may give here various code behavior options.

.SH EXAMPLES
time ../cballs nbody=100000 rangeN=500 sizeHistN=100 lbox=10000 verb=2 mToPlot=2

.sp
You will have computation written in a file named with option \fBoutfile\fR in a directory named Output. If this directory does not exist, it will be created.
Also will be created a directory named "tmp" where a log file will be saved.
You may plot the solutions using \fBnplot2d\fR code:

.sp
nplot2d in=Output/output.txt ws=1 symbolsize=0.1 pj=0 symboltype=4 uc=1:2

.sp
\fBnplot2d\fR is very similar to gnuplot. Then, you can use the one you have available.

.sp
Other files can be plotted as:

.sp
nplot2d in=Output/histN.txt,Output_to-compare-with/histN_direct.txt uc=1:2,1:2

.sp
nplot2d in=Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt uc=1:2,1:3,1:4,1:5,1:6

.sp
Using brute force:

.sp
time cballs o=output options=compute-HistN verb=2 sizeHistN=40 sizeHistTheta=40 rangeN=200 nbody=10000 mToPlot=2 lbox=2000 search=direct-simple-3pcf-omp

.sp
Reading a Takahasi simulation data file:

.sp
time cballs infile=Takahasi/allskymap_nres12r000.zs9.mag.dat infmt=takahasi options=stop lbox=0.2 o=takas_zs9

.sp
The with the output file we run cballs:

.sp
time cballs in=Input/takas_zs9.txt options=compute-HistN verb=2 sizeHistN=40 rangeN=0.02 search=tree-omp

.sp
Other OMP routines:

.sp
time cballs o=output options=compute-HistN verb=2 nbody=64000 search=tree-omp

.sp
time cballs o=output options=compute-HistN verb=2 nbody=64000 search=tree-omp-sincos

.sp
time cballs o=output options=compute-HistN verb=2 nbody=64000 search=direct-omp

.sp
time cballs o=output options=compute-HistN verb=2 nbody=64000 search=direct-omp-sincos

.sp
Using python scripts you can generate several plots:

.sp
cd python

.sp
python plotting_points.py

.sp
python plotting_histN_histXi2pcf.py

.sp
python3 plotting_zetaM.py

.sp
To smooth a data points (reducing its number) do:

.sp
cballs o=output ofmt=columns-ascii options=smooth verb=2 nbody=64000 nsmooth=4

.sp
To convert a catalog from 2D to 3D assuming 2D column format where theta is the first column and phi the second column do:

.sp
cballs in=Input/Taka_nres12r043.zs9_all.txt infmt=columns-ascii-2d-to-3d o=takas_AA_3D options=convert

.SH SEE ALSO
nplot2d(1) cballs(1)

.SH COPYRIGHT
Copyright (C) 2023
.br
Mario A. Rodriguez-Meza et al.
.br
