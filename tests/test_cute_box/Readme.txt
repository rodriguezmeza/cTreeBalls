
Note:
If you are using useLogHist = false, then use BALLSON = 0 in addons/Makefile_addons
and in terminal: make clean; make all


In data.txt there is a configuration space sample. Can be visualized with:

$ python data_3d_view.py

If you have cute_box installed, correlation can be done:

$ CUTE_BOX params.txt

To compare with the results obtained using cBalls do:

$ ../../cballs parameters_to-test_data.txt

Compare results with (from tests directory):

$ python plotCF.py


There is also a lpicola n-body simulations results in folder example_cola. Process with cBalls using

$ cballs parameters_to-test_data_cola_064.txt 

and also can be analyzed using CUTE_BOX. For a gadget input (check the params_cola_064.txt for details):

$ CUTE_BOX params_cola_064.txt

