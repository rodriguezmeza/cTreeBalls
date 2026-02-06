
cBalls has two searching tool to compute correlations for periodic box catalogs:

1. neighbor-boxes-omp
2. kdtree-box-omp

Fastest one is the first and is called using "searchMethod=neighbor-boxes-omp" option. Look for it in the python scripts or in parameter files below.

Computing correlations of N-body simulations like L-picola and halos catalog given by ROCKSTAR

Let us assume we have a halo catalog (extracted with ROCKSTAR from an N-body simulation):

$ time python python/bao_from_halos.py --in catalogs_working/l-picola/r000/rockstar_galaxies/z0p500/halos_0.0.ascii --tag "z0p5" --box 1024 --threads 16 --rmin 0.1 --rmax 150 --nbins 40 --x-col 9 --y-col 10 --z-col 11 --masamin 11 --masamax 13 --out ./Output

At the end you may delete Output/xyz_z0p5.txt which contains all selected halo positions.

CUTE_BOX produce correlations for some data in output_cute_box folder: corr128.dat was the output of data.txt and corr128_cola_064 was the output of L-Picola run in test_cute_box/example_cola. Files in the former are in gadget2 format. For details see CUTE_BOX parameters input files: params.txt and params_cola_064.txt.

In test_cute_box/data.txt is a configuration space sample. Can be visualized with:

$ python test_cute_box/data_3d_view.py

To compare with the results obtained using cBalls do:

$ ../cballs In/parameters_to-test_data.txt

Compare results with:

$ python test_cute_box/plotCF.py


There is also a lpicola n-body simulations results in folder test_cute_box/example_cola. Process them with cBalls using

$ ../cballs In/parameters_to-test_data_cola_064.txt 


Note: others L-Picola halo catalogs (extracted with ROCKSTAR) can be processed in catalogs/l-picola
