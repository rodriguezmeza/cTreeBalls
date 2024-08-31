
Note:
If you are using useLogHist = false, then use BALLSON = 0 in addons/Makefile_addons
and in terminal: make clean; make all


In data.txt there is a configuration space sample. Can be visualized with:
(Run from tests directory)

$ python test_cute_box/data_3d_view.py

If you have cute_box installed, correlation can be done:

$ cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt

To compare with the results obtained using cBalls do:

../../cballs parameters_to-test_data.txt

Compare results with (from tests directory):

$ python test_cute_box/plotCF.py
