#
# Check settings in addons/Makefile_addons_settings:
#   CLASSLIBON = 1
#   TPCFON = 1 only if you require non-zero 3PCF/ZetaM output.
# With TPCFON = 0 or computeTPCF false, PXD ZetaM getters must not segfault;
# they return a zero matrix while MainLoop arrays are live.
#
# If not, set them and execute in cBalls main directory: 'make clean; make all'
#
# Run as:
# python test_cython_octree-ggg-omp.py
#

import os, sys
import numpy as np

# Determine the absolute path of the target (cyballs) directory
#   these two lines won´t be necessary if cballys is in searching path
#target_directory = os.path.abspath('/opt/homebrew/anaconda3/lib/python3.13/site-packages/')
# Append the directory to sys.path
#sys.path.append(target_directory)
from cyballs import cballs
#from cballys import *
Balls = cballs()
Balls.set({'searchMethod':'octree-ggg-omp'})
Balls.set({'infile':'./catalogs/Abraham/kappa_nres12_zs9NS256r000.bin'})
Balls.set({'infileformat':'binary'})
#
# this configuration gives: 'ZetaM contains non-finite values'. Comment it
#Balls.set({'rangeN':0.0633205,'rminHist':0.00213811,'sizeHistN':20,'numberThreads':16})
#
Balls.set({'rangeN':0.0633205,'rminHist':0.00313811,'sizeHistN':20,'numberThreads':16})
#
Balls.set({'verbose':2,'verbose_log':2})
Balls.set({'rootDir':'Output'})
Balls.set({'options':'compute-HistN,and-CF,out-m-HistZeta'})
#

print('cBalls version = ', Balls.getVersion())

cputime = Balls.Run(level=["MainLoop"])

try:
    print('theta = ', Balls.getTheta())
    print('sizeHistN = ', Balls.getsizeHistN())
    print('rBins = ', Balls.getrBins())
    print('histNN = ', Balls.getHistNN())
    print('histCF = ', Balls.getHistCF())
    print('histXi2pcf = ', Balls.getHistXi2pcf())

    # getHistZetaM_sincos(m, type): m multipole,
    #   type: 1 - cos; 2 - sin; 3 - sincos; 4 - cossin
    zeta = Balls.getHistZetaMsincos(1, 1)
    expected_shape = (Balls.getsizeHistN(), Balls.getsizeHistN())

    if zeta.shape != expected_shape:
        raise SystemExit(f"unexpected ZetaM shape: {zeta.shape}, expected {expected_shape}")

    if not np.all(np.isfinite(zeta)):
        raise SystemExit("ZetaM contains non-finite values")

    # monopole:
    print('ZetaM(1,1) = ', zeta)
finally:
    Balls.struct_cleanup()

print('Searching (Balls.run) cputime=', cputime, ' sec.')
