#
# Run as:
# python test_cython_balls.py
#
from cballys import *
Balls = cballs()
Balls.set({'searchMethod':'balls-omp','scanLevelRoot':0,'scanLevel':11})
Balls.set({'infile':'./Takahasi/full_sky_whole_XYZK__zs9_r000_nside512.bin'})
Balls.set({'infileformat':'binary'})
Balls.set({'rangeN':0.0633205,'rminHist':0.00213811,'sizeHistN':20,'numberThreads':16})
Balls.set({'verbose':2,'verbose_log':2})
Balls.set({'rootDir':'Output'})
Balls.set({'options':'compute-HistN,and-CF,out-m-HistZeta'})
Balls.Run()

