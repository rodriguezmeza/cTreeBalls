
cBalls: 1, 2, 3... for people in a hurry... boys and girls just want it plug and play... peladito y en la boca... SNI's under the siege of an evaluation year...

(
Some where in your computer $HOME do
$ git clone https://github.com/rodriguezmeza/cTreeBalls.git
)

1. FIRST STEPS

It is assuming that you are in the cTreeBalls directory (where there is a git cloned version of the code).

First:

$ make clean; make all

Then (after too many seconds... a price must be paid to not look for the root administrator and get him done some tasks for you...)

$ mkdir deleteme; cd deleteme

Run

$ ../cballs

A distribution of uniformly points in the unit sphere is generated with a gaussian values of convergence field. The (2,3)-point correlation functions are computed and saved as histograms in a folder name Output just created by cBalls. Also you will find a 
file named "parameters_null-cballs-usedvalues". Do

$ more Output/parameters_null-cballs-usedvalues

and you will see all the parameters cBalls needs and their default values. A brief description can be found by

$ ../cballs --help

If above is not enough do

$ more ../tests/In/parameters_explained

But if not... and an extended description is gotten by

$ man ../doc/man/cballs.1

Read as thoroughly as needed in particular follow the examples at the end for reading files.


2. TO TEST cBalls

$ cd ../tests
$ ./run_all_tests

or test individual cases:

$ time ../cballs ./In/parameters_to-test_nside256_balls-omp

to plot results (vs reference results):

$ python script/mZetaM_test_kkk_all.py


3. SOME MORE ELABORATE EXAMPLES

$ time ../cballs search=octree-ggg-omp infile=./catalogs/Takahasi/allskymap_nres08r081_zs9_mag.fits infmt=fits-healpix

Plot the results with:

$ python Xi3pcf_plot_flatten.py

Now do (to use ud_grade, it is need to have installed Healpix, go to https://healpix.sourceforge.io/ to download and install it)

$ cp ../tests/in_ud_grade .
$ cp ../tests/Xi3pcf_plot_flatten_EE.py .
$ time ../cballs search=octree-ggg-omp infile=./allskymap_nres10r081_zs9_mag.fits infmt=fits-healpix \
options=pre-processing,post-processing,no-normalize-HistZeta,edge-corrections,smooth-pivot \
preScript="ud_grade < in_ud_grade" posScript="python Xi3pcf_plot_flatten_EE.py; rm -f allskymap_nres10r081_zs9_mag.fits"

and the results were already plotted with the post-processing action in 'options'. Also the file 'allskymap_nres10r081_zs9_mag.fits' was [removed!

The "\" is to connect such a long command line.


4. MIXING SEVERAL INPUT FILES

Get some realizations of Takahashi et al. (http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/)

$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs9.mag.dat
$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs8.mag.dat
$ wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub2/nres12/allskymap_nres12r081.zs10.mag.dat

Convert them to Healpix

$ cp ../tests/catalogs/Takahasi/takahasi_to_fits.py .

Edit this python script in order to have defined one of the above realizations (lines 5 and 6).

$ python takahasi_to_fits.py

Repeat two more times to have all three realization files.

Downsize these Healpix files: edit file in_ud_grade to have the correct lines:

a. 
allskymap_nres12r081.zs9.mag.fits
1024
allskymap_nres10r081_zs9_mag.fits

$ ud_grade < in_ud_grade

b. 
allskymap_nres12r081.zs8.mag.fits
512
allskymap_nres09r081_zs8_mag.fits

$ ud_grade < in_ud_grade

c. 
allskymap_nres12r081.zs10.mag.fits
512
allskymap_nres09r081_zs10_mag.fits

$ ud_grade < in_ud_grade


Then edit a parameter (parameters.txt) file to have the following lines:

searchMethod = octree-ggg-omp
options = all-in-one,GGGCorrelation,no-normalize-HistZeta,edge-corrections,smooth-pivot
infile = allskymap_nres09r081_zs8_mag.fits,allskymap_nres10r081_zs9_mag.fits,allskymap_nres09r081_zs10_mag.fits
infileformat = fits-healpix,fits-healpix,fits-healpix
iCatalogs = 1,2,3

Execute

$ time ../cballs parameters.txt

and plot results:

$ python Xi3pcf_plot_flatten.py

All these command lines process three Takahashi's realizations give them to cBalls, it reads them and creates a single bodies catalog with all of them and then do the neighours searching process and compute (2,3)-point correlation functions and save them as histograms in the default output directory. Results are plotted with the python scripts.


