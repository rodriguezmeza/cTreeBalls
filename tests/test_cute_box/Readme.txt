
In test_cute_box/data.txt there is a configuration space sample. Can be visualized with:

$ python test_cute_box/data_3d_view.py

If you have cute_box installed, correlation can be done:

$ (PATH_TO_CUTE_BOX)/cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ (PATH_TO_CUTE_BOX)/cute_box params_cola_064.txt

To compare with the results obtained using cBalls do:

$ ../cballs ./In/parameters_to-test_data.txt

Compare results with:

$ python test_cute_box/plotCF.py


There is also a lpicola n-body simulations results in folder example_cola. Process with cBalls using

$ ../cballs/cballs ./In/parameters_to-test_data_cola_064.txt 

and also can be analyzed using CUTE_BOX. For a gadget input (check the params_cola_064.txt for details):

$ (PATH_TO_CUTE_BOX)/cute_box params_cola_064.txt

Transform Rockstar multi-columns-ascii halos to columns-ascii

$ time ../cballs search=octree-sincos-omp rootDir=Output rangeN=150 rminHist=0.1 sizeHistN=30 nthreads=16 stepState=1000 useLogHist=false usePeriodic=true verb=2 in=halos_0.0.ascii infmt=multi-columns-ascii columns=9,10,11 options=stop,only-pos o=halos ofmt=multi-columns-ascii

Note: refer to Rockstar documentation to see how to extract halo catalogs from a Gadget or L-picola simulation.