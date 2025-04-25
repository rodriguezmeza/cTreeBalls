# this code cannot be run directly
# do 'source env_config/nagbodyrc.sh' from your sh shell or put it in your profile
#
#

alias nmake='make -f NagBody'

#
#To solve some configuration problems related to zsh visit:
#https://scriptingosx.com/2019/06/moving-to-zsh-part-2-configuration-files/
#
#
# NagBody Defintion of variables ...

# To activate nnbodykit environment (added to load notebooks in jupyter, does not work)
# (For Anaconda2)
#source /anaconda2/bin/activate /anaconda2/envs/nbodykit-env/

# NagBody:
export PATH=${HOME}/NagBody_pkg/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/man:${MANPATH}

# CFITSIO (NEEDED BY CAMB and Healpix):
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/cfitsio/lib:${LD_LIBRARY_PATH}

# FFTW-2 (Needed by gadget type of codes):
# INSTALADO CON PORT ... fftw, fftw-single, falto que se incluya con openmpi ...
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw2/lib:${DYLD_LIBRARY_PATH}

# FFTW3 (Needed by cola type of codes):
# Instalar con port ... fftw-3 ... tambien la precision simple fftw-3-single
export PATH=${HOME}/NagBody_pkg/local/fftw3/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/fftw3/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/fftw3/lib/pkgconfig:${PKG_CONFIG_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${DYLD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/fftw3/lib:${LD_LIBRARY_PATH}

# GSL (Needed by gadget and cola type of codes):
# INSTALADO CON PORT ...
export PATH=${HOME}/NagBody_pkg/local/gsl/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/gsl/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/gsl/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/gsl/lib:${LD_LIBRARY_PATH}
#

# Healpix (Manual installation in $(HOME)/NagBody_pkg/local/Healpix_3.30):
# If this env configuration was not done in .bash_profile use these line:
#
# modifications by HEALPixAutoConf 3.30
# (Note: see this particular user path, change it properly as it si done in the following line)
#[ -r /Users/mar/.healpix/3_30_Darwin/config ] && . /Users/mar/.healpix/3_30_Darwin/config
#[ -r ${HOME}/.healpix/3_30_Darwin/config ] && . ${HOME}/.healpix/3_30_Darwin/config
#
#
# Modification suggested by Healpix in the configuration process:
#  Where the file /global/homes/m/mrodrig/.healpix/3_81_Linux/config contains:
# configuration for Healpix 3.81
#HEALPIX=/global/homes/m/mrodrig/NagBody_pkg/local/Healpix_3.81 ; export HEALPIX
#HPX_CONF_DIR=/global/homes/m/mrodrig/.healpix/3_81_Linux
#if [ -r ${HPX_CONF_DIR}/idl.sh ] ; then . ${HPX_CONF_DIR}/idl.sh ; fi
#if [ -r ${HPX_CONF_DIR}/gdl.sh ] ; then . ${HPX_CONF_DIR}/gdl.sh ; fi
#if [ -r ${HPX_CONF_DIR}/fl.sh ]  ; then . ${HPX_CONF_DIR}/fl.sh  ; fi
#if [ -r ${HPX_CONF_DIR}/f90.sh ] ; then . ${HPX_CONF_DIR}/f90.sh ; fi
#if [ -r ${HPX_CONF_DIR}/cpp.sh ] ; then . ${HPX_CONF_DIR}/cpp.sh ; fi
#if [ -r ${HPX_CONF_DIR}/c.sh ] ;   then . ${HPX_CONF_DIR}/c.sh ;   fi
#
#
# modifications by HEALPixAutoConf 3.81
[ -r /Users/mar/.healpix/3_81_Darwin/config ] && . /Users/mar/.healpix/3_81_Darwin/config
# Be aware that $HOME is different in other platforms. Change accordingly.


# LAPACK (can be installed with ports en Mac):
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/lapack/lib:${LD_LIBRARY_PATH}

###################################################################
# PLPlot (ver. 5.15.0):
export PATH=${HOME}/NagBody_pkg/local/plplot/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/plplot/share/man:${MANPATH}
export PKG_CONFIG_PATH=${HOME}/NagBody_pkg/local/plplot/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib/plplot5.15.0/driversd:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/plplot/lib/plplot5.15.0/driversd:${DYLD_LIBRARY_PATH}

###################################################################
# HDF5:
# INSTALADO CON PORT ... hdf5-18 ...
export PATH=${HOME}/NagBody_pkg/local/hdf5/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/hdf5/share/man:${MANPATH}
export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/hdf5/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/hdf5/lib:${DYLD_LIBRARY_PATH}

# PHDF5 (Parallel):
#export PATH=${HOME}/NagBody_pkg/local/phdf5/bin:${PATH}
#export MANPATH=${HOME}/NagBody_pkg/local/phdf5/share/man:${MANPATH}
#export LD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/phdf5/lib:${LD_LIBRARY_PATH}
#export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/phdf5/lib:${DYLD_LIBRARY_PATH}

# GADGETVIEWER:
export PATH=${HOME}/NagBody_pkg/local/gadgetviewer/bin:${PATH}

###################################################################
# Planck 2018:
# To activate Planck 2018 for CLASS:
#source $HOME/NagBody_pkg/NagBody_sources/class/planck2018/src/code/plc_3.0/plc-3.01/bin/clik_profile.sh
source $HOME/NagBody_pkg/bin/clik_profile.sh



###################################################################
# Python environment (ver 3.9 and 3.10): Useful for TreeCorr
#export PYTHONPATH=${HOME}/.local/lib/python3.9/site-packages:${PYTHONPATH}
#export PYTHONPATH=${HOME}/NagBody_pkg/local/TreeCorr/lib/python3.10/site-packages:${PYTHONPATH}
#export PATH=${HOME}/NagBody_pkg/local/TreeCorr/bin:${PATH}

###################################################################
# Open MPI (ver 4.1):
export PATH=${HOME}/NagBody_pkg/local/openmpi/bin:${PATH}
export MANPATH=${HOME}/NagBody_pkg/local/openmpi/share/man:${MANPATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/openmpi/lib:${DYLD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${HOME}/NagBody_pkg/local/openmpi/lib/openmpi:${DYLD_LIBRARY_PATH}

