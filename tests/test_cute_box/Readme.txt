
In data.txt there is a configuration space sample. Can be visualized with

$ nplot2d in=data.txt plotjoined=0  ws=1

Correlation can be done:

$ cute_box params.txt

For a gadget input (check the params_cola_064.txt for details):

$ cute_box params_cola_064.txt
$ nplot2d in=corr128_cola_064.dat

Note: nplot2d is parte of NagBody project.
