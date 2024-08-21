from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os
import subprocess as sbp
import os.path as osp

# Recover the gcc compiler
GCCPATH_STRING = sbp.Popen(
    ['gcc', '-print-libgcc-file-name'],
    stdout=sbp.PIPE).communicate()[0]
GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()

#liblist = ["cballs","gsl","gslcblas"]
liblist = ["cballs"]
MVEC_STRING = sbp.Popen(
    ['gcc', '-lmvec'],
    stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec","m"]

# define absolute paths
root_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
include_folder = os.path.join(root_folder, "include")
general_libs_folder = os.path.join(root_folder, "general_libs")
getparam_folder = os.path.join(root_folder, "getparam")
gsl_folder = os.path.join(root_folder, "gsl")
class_lib_folder = os.path.join(os.path.join(root_folder, "addons"),"class_lib")
cballys_folder = os.path.join(root_folder, "python")

# Recover the cBalls version
with open(os.path.join(class_lib_folder, 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

# Define cython extension and fix Python version
cballys_ext = Extension("cballys", [os.path.join(cballys_folder, "cballys.pyx")],
                           include_dirs=[nm.get_include(), include_folder, class_lib_folder, general_libs_folder, getparam_folder],
                           libraries=liblist,
                           library_dirs=[root_folder, GCCPATH],
                           extra_link_args=['-lgomp']
#                           extra_link_args=['-lgomp -lgsl -lgslcblas']
                       )
import sys
cballys_ext.cython_directives = {'language_level': "3" if sys.version_info.major>=3 else "2"}

setup(
    name='cballys',
    version=VERSION,
    description='Python interface to the nPCF code cballs',
#    url='http://github.com/rodriguezmeza/cTreeBalls.git',
    cmdclass={'build_ext': build_ext},
    ext_modules=[cballys_ext],
)
