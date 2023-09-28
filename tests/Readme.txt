

Some command line examples:


./cleansh 
cd /Users/mar/Research/Codigos/NagBody_pkg/NagBody_sources/perturbations/tpcf/
mkdir run
cp -pR tests/* run

Run considered as a reference:

tpcf nbody=100000 rangeN=500 sizeHistN=100 search=direct verb=2 stepState=10000 lbox=10000 > output.log &

time tpcf o=output options=compute-HistN verb=2 nbody=64000 mToPlot=2 search=direct-simple

time tpcf o=output options=compute-HistN verb=2 nbody=64000 mToPlot=2 search=direct-simple

nplot2d in=Output/histN.txt 
nplot2d in=Output/output.txt ws=1 symbolsize=0.1 pj=0 symboltype=4 uc=1:2

nplot2d in=Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt,Output/mhistZeta.txt uc=1:2,1:3,1:4,1:5,1:6

nplot2d in=Output/cputime.txt,Output/cputime.txt,Output/cputime.txt plottype=3 ws=1,0,0 plotjoined=0,1,1 symbolsize=10 uc=1:2,1:4,1:5 xr=3.5:5.2 yr=-2:2


TO COMPARE WITH SIMPLE:

tpcf o=output options=compute-HistN verb=2 nbody=64000 mToPlot=2 search=direct

nplot2d in=Output/histN.txt,Output_simple/histN.txt uc=1:2,1:2

nplot2d in=Output_simple/mhistZeta.txt,Output/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt uc=1:2,1:2,1:3,1:4,1:5,1:6 ws=0,1,0,0,0,0

nplot2d in=Output_simple/mhistZeta.txt,Output/mhistZeta_2.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt uc=1:2,1:2,1:3,1:4,1:5,1:6 ws=0,1,0,0,0,0

nplot2d in=Output_simple/mhistZeta.txt,Output/mhistZeta_1.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt,Output_simple/mhistZeta.txt uc=1:2,1:2,1:3,1:4,1:5,1:6 ws=0,1,0,0,0,0

