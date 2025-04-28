#
# Check settings in addons/Makefile_addons_settings 'CLASSLIBON' must be:
#   CLASSLIBON = 1
# If not, set it and execute in cBalls main directory: 'make clean; make all'
#   then go to 'python' directory and execute:
#   python setup.py build
#   python setup.py install --user
#
# Run as:
# python test_cython_balls.py
#
from cballys import *
Balls = cballs()
Balls.set({'searchMethod':'balls-omp','scanLevelRoot':0,'scanLevel':11})
Balls.set({'infile':'./catalogs/Abraham/kappa_nres12_zs9NS256r000.bin'})
Balls.set({'infileformat':'binary'})
Balls.set({'rangeN':0.0633205,'rminHist':0.00213811,'sizeHistN':20,'numberThreads':16})
Balls.set({'verbose':2,'verbose_log':2})
Balls.set({'rootDir':'Output'})
Balls.set({'options':'compute-HistN,and-CF,out-m-HistZeta'})
Balls.Run()
#
print('sizeHistN = ', Balls.getsizeHistN())
print('rBins = ', Balls.getrBins())
print('ZetaM(1,1) = ', Balls.getHistZetaM_sincos(1, 1))

