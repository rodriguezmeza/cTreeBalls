#
# Check settings in addons/Makefile_addons_settings 'CLASSLIBON' must be:
#   CLASSLIBON = 1
#   also set TPCF = 1 in above file.
#
# If not, set them and execute in cBalls main directory: 'make clean; make all'
#
# Run as:
# python test_cython_octree-ggg-omp.py
#
import os, sys
# Determine the absolute path of the target (cballys) directory
#   these two lines won´t be necessary if cballys is in searching path
target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
sys.path.append(target_directory)
from cballys import cballs
#from cballys import *
Balls = cballs()
Balls.set({'searchMethod':'octree-ggg-omp'})
Balls.set({'infile':'./catalogs/Abraham/kappa_nres12_zs9NS256r000.bin'})
Balls.set({'infileformat':'binary'})
Balls.set({'rangeN':0.0633205,'rminHist':0.00213811,'sizeHistN':20,'numberThreads':16})
Balls.set({'verbose':2,'verbose_log':2})
Balls.set({'rootDir':'Output'})
Balls.set({'options':'compute-HistN,and-CF,out-m-HistZeta'})
#
cputime=Balls.Run()
#
print('cBalls version = ', Balls.getVersion())
print('theta = ', Balls.getTheta())
print('sizeHistN = ', Balls.getsizeHistN())
print('rBins = ', Balls.getrBins())
print('histNN = ', Balls.getHistNN())
print('histCF = ', Balls.getHistCF())
print('histXi2pcf = ', Balls.getHistXi2pcf())
# getHistZetaM_sincos(m, type): m multipole,
#   type: 1 - cos; 2 - sin; 3 - sincos; 4 - cossin
# monopolo cos:
# this line will give segmentation fault if TPCF = 0 in
#   addons/Makefile_addons_settings
print('ZetaM(1,1) = ', Balls.getHistZetaMsincos(1, 1))

#B
print('Searching (Balls.run) cputime=',cputime,' sec.')
