To use this python wrapper for cballs, you should first compile cballs with 'make all', not just 'make cballs': this is important in order to create libraries. Be sure that you did not remove the -fPIC compiling flag in the Makefile, important for compatibility between OpenMP and the python wrapper.

Then, execute:

$ python setup.py build
$ python setup.py install --user

You can check that these steps work by typing

$ python
>>> from cballys import cballs

If python does not complain, the cballs module has been correctly installed in your python distribution. You can now import it and use its functions from your python codes.  

Notes:
1. GSL is not working so far. Go Makefile_settings (USEGSL) and turn it off.
2. Set NDIM (in ccballys.pxd) as pretend to use cBalls. Set it according to Makefile_settings.