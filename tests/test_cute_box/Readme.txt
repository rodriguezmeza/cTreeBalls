
Note:
If you are using useLogHist = false, then use BALLSON = 0 in addons/Makefile_addons
and in terminal: make clean; make all


In data.txt there is a configuration space sample. Can be visualized with:

$ python data_3d_view.py

If you have cute_box installed, correlation can be done:

$ cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt

To compare with the results obtained using cBalls do:

$ ../../cballs parameters_to-test_data.txt

Compare results with:

$ python plotCF.py


There is also a lpicola n-body simulations results in folder example_cola. Process with cBalls using

$ cballs parameters_to-test_data_cola_064.txt 

and also can be analyzed using CUTE_BOX. For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt


Transform gadget multi-columns-ascii halos to columns-ascii

$ time cballs search=tree-omp-sincos rootDir=Output rangeN=150 rminHist=0.1 sizeHistN=30 nthreads=16 stepState=1000 computeTPCF=false useLogHist=false usePeriodic=true verb=2 in=halos_0.0.ascii infmt=multi-columns-ascii columns=9,10,11 options=stop,only-pos o=halos ofmt=multi-columns-ascii

