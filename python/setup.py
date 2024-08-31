#
# Test it with python 3.7
#   there are some issues with python 3.9...
#
#https://blog.ganssle.io/articles/2021/10/setup-py-deprecated.html
#https://build.pypa.io/en/stable/
#
# Should be installed with 'pip install .' ... check it!
#
# pip install setuptools --upgrade
# pip install numpy --upgrade
#
# Could be useful to add
# pipVersion = pkg_resources.require("pip")[0].version
#setuptoolsVersion = pkg_resources.require("setuptools")[0].version
#
#print("\n PIP Version", pipVersion, "\n")
#print("\n Setuptools Version", setuptoolsVersion, "\n")
#

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
pxd_folder = os.path.join(os.path.join(root_folder, "addons"),"pxd")
cballys_folder = os.path.join(root_folder, "python")

# Recover the cBalls version
with open(os.path.join(class_lib_folder, 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

# Look for better way to find cfitsio path lib...
#   looking for cfitsio libs are OFF by default...
#   find its path in your system, edit line below and uncomment it
# Define cython extension and fix Python version
cballys_ext = Extension("cballys", [os.path.join(cballys_folder, "cballys.pyx")],
                           include_dirs=[nm.get_include(), include_folder, class_lib_folder, general_libs_folder, getparam_folder, pxd_folder],
                           libraries=liblist,
                           library_dirs=[root_folder, GCCPATH,
#                            '/Users/mar/NagBody_pkg/local/cfitsio/lib/'
                           ],
                           extra_link_args=['-lgomp']
#                           extra_link_args=['-lgomp', '-lcfitsio']
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
