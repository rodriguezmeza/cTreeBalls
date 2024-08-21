
In data.txt there is a configuration space sample. Can be visualized with

$ nplot2d in=data.txt plotjoined=0  ws=1

Correlation can be done:

$ cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt
$ nplot2d in=corr128_cola_064.dat

Note: nplot2d is parte of NagBody project:

https://github.com/rodriguezmeza/NagBody_pkg

In this link COLAs n-body codes can be found also.

To compare with the results obtained using cBalls do:

In addons/Makefile_addons set BALLSON = 0 then 

$ make clean; make all

Now execute cBalls:

../cballs parameters_to-test_data-txt_tree-omp-sincos

Compare results with:

$ nplot2d in=test_cute_box/output_cute_box/corr128.txt,Output/histCF.txt uc=1:2,1:2  ws=1,1 plotjoined=1,0 plottype=1

Or use the python script in test_cute_box:

$ python test_cute_box/plotCF.py
