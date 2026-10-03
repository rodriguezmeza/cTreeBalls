gsl INSTALLATION

1. The way to install the package is to follow the steps:

If necessary, define the following environment variables:

For a standard Linux and Mac OS X:
export CC=gcc
export CXX=g++
export F77=gfortran
export FC=gfortran
export F90=gfortran

Intel compilers:
export CC=icc
export CXX=icc
export F77=ifort
export FC=ifort
export F90=ifort

2. Configure, make and install:

./configure --prefix=$HOME/NagBody_pkg/local/gsl 2>&1 | tee configure_gsl.log
make 2>&1 | tee make_gsl.log
make check 2>&1 | tee check_gsl.log
make install 2>&1 | tee install_gsl.log

3. Clean directory:

make distclean

4. We assume you are using bash shell. Modify env_config/nagbodyrc.sh to set the appropriate environment variable: make it to contain the following lines:

export PATH=${HOME}/NagBody_pkg/local/gsl/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/gsl/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/gsl/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/gsl/lib:${DYLD_LIBRARY_PATH}

Then 

source env_config/nagbodyrc.sh.

That's it! 

