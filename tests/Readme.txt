
To test the code:

$ time ../cballs options=post-processing script="python scripts/plot2pcf.py"

It will run the code, make a simple cube of random points, and as a post-process plot the 2pcf and save it as a pdf file. All output is put in a folder with name "Output", except output produce by the python script. If "Output" doesn't exist it is created. Next time cBalls is run it is overwritten. The name of the output folder can be change. Look at the help: "../cballs -h" or see parameter file "parameters_explained".

Or you can use a parameter file:

time ../cballs parameters_to-test_nside512_balls-omp

Computation time is about a 45 seconds.

It uses a xyz catalog of points on the celestial sphere with kappa (WL convergence) values. Total number of points is 3145728.

This catalog is a downsizing version of realization Nside1024, zs9, r000 in:

http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/
