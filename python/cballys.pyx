"""
.. module:: cballys
    :synopsis: Python wrapper around cTreeBalls
.. moduleauthor:: Mario A. Rodriguez-Meza <marioalberto.rodriguezmeza@gmail.com>

.. based on Julien Lesgourges' CLASS

This module defines a class called cballs.

# MAR 27.04.2023: TODO:

"""

import numpy as np
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython

import time

import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

ctypedef np.float64_t DTYPE_t
ctypedef np.int32_t DTYPE_i

from ccballys cimport *

class cBallsError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in cballs: ' + self.message

class cBallsSevereError(cBallsError):
    """
    Raised when cballs failed to understand one or more input parameters.

    This case would not raise any problem in cballs default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong computation would be selected.
    """
    pass

class cBallsSevereErrorDummy():
    """
    Raised when cballs failed to understand one or more input parameters.

    This case would not raise any problem in cballs default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong computation would be selected.
    """
    pass

class cBallsComputationError(cBallsError):
    """
    Raised when cballs could not compute the computation at this point.

    """
    pass

cdef class cballs:
    """
    cballs wrapping, creates the glue between C and python

    The actual cballs wrapping, the only class we will call from Python
    (indeed the only one we will import, with the command:
    from cballys import cballs

    """
    cdef cmdline_data cmd
    cdef global_data gd
    cdef file_content fc

    cdef int nthreads
    cdef double cputime

    cdef int computed # Flag to see if cballys has already computed with the given pars
    cdef int allocated # Flag to see if cballys structs are allocated already
    cdef object _pars # Dictionary of the parameters
    cdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    def set_default(self):
        _pars = {
            "searchMethod":"octree-ggg-omp",
            "rminHist":0.00213811,
            "rangeN":0.0633205,
            "sizeHistN":20,
            "rootDir":"Output",
            "infileformat":"fits-healpix",
            "nsmooth":8,
            "rsmooth":"\0",
            "theta":1.05,
            "iCatalogs":"1",
            "columns":"1,2,3,4",
            "stepState":1000000,
            "numberThreads":16,
            "verbose":0,
            "verbose_log":0,
            "options":"compute-HistN",
            }
        self.set(**_pars)

    def __cinit__(self, default=True):
        cdef char* dumc
        self.allocated = False
        self.computed = False
        self._pars = {}
        self.fc.size=0
        self.fc.filename = <char*>malloc(sizeof(char)*30)
        assert(self.fc.filename!=NULL)
        dumc = "NOFILE"
        sprintf(self.fc.filename,"%s",dumc)
        self.ncp = set()
        if default: self.set_default()

    def __dealloc__(self):
        if self.allocated:
          self.struct_cleanup()
        self.clean()
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
            free(self.fc.filename)

    def set(self,*pars,**kars):
        oldpars = self._pars.copy()
        if len(pars)==1:
            self._pars.update(dict(pars[0]))
        elif len(pars)!=0:
            raise cBallsSevereError("bad call")
        self._pars.update(kars)
        if viewdictitems(self._pars) <= viewdictitems(oldpars):
          return # Don't change the computed states, if the new dict was already contained in the previous dict
        self.computed=False
        return True

    def clean(self):
        self._pars = {}
        self.computed = False

    def _fillparfile(self):
        cdef char* dumc

        if self.fc.size!=0:
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
        self.fc.size = len(self._pars)
        self.fc.name = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.name!=NULL)

        self.fc.value = <FileArg*> malloc(sizeof(FileArg)*len(self._pars))
        assert(self.fc.value!=NULL)

        self.fc.read = <short*> malloc(sizeof(short)*len(self._pars))
        assert(self.fc.read!=NULL)

        i = 0
        for kk in self._pars:

            dumcp = kk.strip().encode()
            dumc = dumcp
            sprintf(self.fc.name[i],"%s",dumc)
            dumcp = str(self._pars[kk]).strip().encode()
            dumc = dumcp
            sprintf(self.fc.value[i],"%s",dumc)
            self.fc.read[i] = FALSE
            i+=1

    def struct_cleanup(self):
        if(self.allocated != True):
          return
#        print('Flags (struct_cleanpup): ',self.gd.tree_allocated,self.gd.gd_allocated_2,
#                        self.gd.bodytable_allocated,self.gd.histograms_allocated,
#                        self.gd.gd_allocated,self.gd.cmd_allocated)
#B        if "MainLoop" in self.ncp:
#B this memory freeing may cause segmentation fault (core dumped) in linux
#   when using cballys in a python loop. Must be analysed and corrected
#        if self.gd.tree_allocated:
#            EndRun_FreeMemory_tree(&(self.cmd), &(self.gd))
#E
        if self.gd.gd_allocated_2:
            EndRun_FreeMemory_gd_2(&(self.cmd), &(self.gd))
        if self.gd.bodytable_allocated:
            EndRun_FreeMemory_bodytable(&(self.cmd), &(self.gd))
        if self.gd.histograms_allocated:
            EndRun_FreeMemory_histograms(&(self.cmd), &(self.gd))
        if self.gd.gd_allocated:
            EndRun_FreeMemory_gd(&(self.cmd), &(self.gd))
        if self.gd.cmd_allocated:
            EndRun_FreeMemory_cmd(&(self.cmd), &(self.gd))
#E
        self.ncp = set()

        self.allocated = False
        self.computed = False

    def clean_all(self):
        self.struct_cleanup()
        self.clean()

    def _check_task_dependency(self, level):
        """
        Fill the level list with all the needed modules

        .. warning::

            the ordering of modules is obviously dependent on CLASS module order
            in the main.c file. This has to be updated in case of a change to
            this file.

        Parameters
        ----------

        level : list
            list of strings, containing initially only the last module required.
            For instance, to recover all the modules, the input should be
            ['lensing']

        """
        if "EndRun" in level:
            if "MainLoop" not in level:
                level.append("MainLoop")
#        if "EndRun" in level:
#            if "EndRun_FreeMemory" not in level:
#                level.append("EndRun_FreeMemory")
#        if "EndRun_FreeMemory" in level:
#            if "MainLoop" not in level:
#                level.append("MainLoop")
        if "MainLoop" in level:
            if "StartRun_Common" not in level:
                level.append("SetNumberThreads")
        if "SetNumberThreads" in level:
            if "PrintParameterFile" not in level:
                level.append("PrintParameterFile")
        if "PrintParameterFile" in level:
            if "StartRun_Common" not in level:
                level.append("StartRun_Common")
        if len(level)!=0 :
            if "input" not in level:
                level.append("input")
        return level

    def _pars_check(self, key, value, contains=False, add=""):
        val = ""
        if key in self._pars:
            val = self._pars[key]
            if contains:
                if value in val:
                    return True
            else:
                if value==val:
                    return True
        if add:
            sep = " "
            if isinstance(add,str):
                sep = add

            if contains and val:
                    self.set({key:val+sep+value})
            else:
                self.set({key:value})
            return True
        return False


#    def Run(self, level=["EndRun_FreeMemory"]):
    def Run(self, level=["MainLoop"]):
        """
        Run(level=["MainLoop"])

        Main function, execute all routines in cballs.

        Parameters
        ----------
        level : list
                list of the last module desired.

        .. warning::

            level default value

        """
        cdef ErrorMsg errmsg

#
# tracking...
#        print('INIT: computed, allocated: ', self.computed, self.allocated)
#
        level = self._check_task_dependency(level)

        if self.computed and self.ncp.issuperset(level):
            return

        if self.allocated:
#            self.gd.tree_allocated=self.getTreeAllocated()
#            self.gd.gd_allocated_2=self.getAllocated2()
#            self.gd.bodytable_allocated=self.getBodytableAllocated()
#            self.gd.histograms_allocated=self.getHistogramsAllocated()
#            self.gd.gd_allocated=self.getGDAllocated()
#            self.gd.cmd_allocated=self.getCMDAllocated()
#            print('Flags (allocated): ',self.gd.tree_allocated,self.gd.gd_allocated_2,
#                        self.gd.bodytable_allocated,self.gd.histograms_allocated,
#                        self.gd.gd_allocated,self.gd.cmd_allocated)
            self.struct_cleanup()

        self.computed = False

        self._fillparfile()

        self.ncp=set()
        self.allocated = True

        if "input" in level:
#
# tracking...
#            print('Track step: input... (0000)')
#
            if input_read_from_file(&self.cmd, &self.gd, &self.fc, errmsg) == FAILURE:
                raise cBallsSevereError(errmsg)
            self.ncp.add("input")
            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == FALSE:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i].decode())
            if problem_flag:
                raise cBallsSevereError(
                    "cballs did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        if "StartRun_Common" in level:
#
# tracking...
#            print('Track step: StartRun_Common... (000)')
#
            if StartRun_Common(&(self.cmd), &(self.gd)) == FAILURE:
                self.struct_cleanup()
                raise cBallsComputationError(self.op.error_message)
#
# tracking...
#            print('Track step: StartRun_Common... (018)')
#
            self.ncp.add("StartRun_Common")

        if "PrintParameterFile" in level:
#
# tracking...
#            print('Track step: PrintParameterFile...')
#
            if PrintParameterFile(&(self.cmd), &(self.gd), "cballys_param.txt") == FAILURE:
                self.struct_cleanup()
                raise cBallsComputationError(self.op.error_message)
            self.ncp.add("PrintParameterFile")

        if "SetNumberThreads" in level:
            if SetNumberThreads(&(self.cmd)) == FAILURE:
                self.struct_cleanup()
                raise cBallsComputationError(self.op.error_message)
            self.ncp.add("SetNumberThreads")
            self.nthreads=self.getNThreads()

# Consider process_time() as another good option...
        if "MainLoop" in level:
#            print('rootDir: ',self.getRootDir())

#            print('Flags (cleanpup): ',self.gd.tree_allocated,self.gd.gd_allocated_2,
#                        self.gd.bodytable_allocated,self.gd.histograms_allocated,
#                        self.gd.gd_allocated,self.gd.cmd_allocated)
            start_wall_time_p = time.process_time()
            if MainLoop(&(self.cmd), &(self.gd)) == FAILURE:
                self.struct_cleanup()
                raise cBallsComputationError(self.op.error_message)
            self.ncp.add("MainLoop")
            end_wall_time_p = time.process_time()
#            end_wall_time = time.perf_counter()
#            self.cputime = end_wall_time - start_wall_time
#            print('tree_allocated',self.gd.tree_allocated)
#            print('gd_allocated_2',self.gd.gd_allocated_2)
#            print('bodytable_allocated',self.gd.bodytable_allocated)
#            print('histograms_allocated',self.gd.histograms_allocated)
#            print('gd_allocated',self.gd.gd_allocated)
#            print('cmd_allocated',self.gd.cmd_allocated)
            self.gd.tree_allocated=self.getTreeAllocated()
            self.gd.gd_allocated_2=self.getAllocated2()
            self.gd.bodytable_allocated=self.getBodytableAllocated()
            self.gd.histograms_allocated=self.getHistogramsAllocated()
            self.gd.gd_allocated=self.getGDAllocated()
            self.gd.cmd_allocated=self.getCMDAllocated()
#            print('Flags (MainLoop): ',self.gd.tree_allocated,self.gd.gd_allocated_2,
#                        self.gd.bodytable_allocated,self.gd.histograms_allocated,
#                        self.gd.gd_allocated,self.gd.cmd_allocated)
            self.cputime = (end_wall_time_p - start_wall_time_p)/self.nthreads
#
# tracking...
#        print('Track step: After MainLoop... (022a)')
#

#        if "EndRun_FreeMemory" in level:
#            if EndRun_FreeMemory(&(self.cmd), &(self.gd)) == FAILURE:
#                self.struct_cleanup()
#                raise cBallsComputationError(self.op.error_message)
#            self.ncp.add("EndRun_FreeMemory")


        self.computed = True

#        print('tree_allocated',self.gd.tree_allocated)
#        print('gd_allocated_2',self.gd.gd_allocated_2)
#        print('bodytable_allocated',self.gd.bodytable_allocated)
#        print('histograms_allocated',self.gd.histograms_allocated)
#        print('gd_allocated',self.gd.gd_allocated)
#        print('cmd_allocated',self.gd.cmd_allocated)

#
# tracking...
#        print('END: computed, allocated: ', self.computed, self.allocated)
#
        return self.cputime

#
#B Interfaces to PXD functions
#
#B flags
#
    def getTreeAllocated(self):
        cdef short value
        cdef short out_value
        if get_tree_allocated(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getAllocated2(self):
        cdef short value
        cdef short out_value
        if get_allocated_2(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getBodytableAllocated(self):
        cdef short value
        cdef short out_value
        if get_bodytable_allocated(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getHistogramsAllocated(self):
        cdef short value
        cdef short out_value
        if get_histograms_allocated(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getGDAllocated(self):
        cdef short value
        cdef short out_value
        if get_gd_allocated(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getCMDAllocated(self):
        cdef short value
        cdef short out_value
        if get_cmd_allocated(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

#E

#B parameters
#
    def getNThreads(self):
        cdef int value
        cdef int out_value
        if get_nthreads(&self.cmd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_value = value
        return out_value

    def getnMonopoles(self):
        cdef int value
        if get_nmonopoles(&self.cmd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        return value

    def getTheta(self):
        cdef double value
        if get_theta(&self.cmd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        return value

    def getrsmooth(self):
        cdef double value
        cdef double out_value
        if get_rsmooth(&self.gd,&value)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_theta = value
        return out_value

    def getCPUTime(self):
        cdef double cputime
        cdef double out_cputime
        if get_cputime(&self.gd,&cputime)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_cputime = cputime
        return out_cputime

    def getsizeHistN(self):
        cdef int sizeHistN
        cdef int out_sizeHistN
        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()
        out_sizeHistN = sizeHistN
        return out_sizeHistN

    def getVersion(self):
        cdef char *param
        cdef char *out_param
        cdef int buffer_size = 20 # Allocate enough space
        # Allocate memory for the destination string
        param = <char*> malloc(buffer_size * sizeof(char))
        out_param = <char*> malloc(buffer_size * sizeof(char))
        if out_param is NULL:
            raise MemoryError("Failed to allocate memory")

        if get_version(&self.cmd,param)== FAILURE:
            raise cBallsSevereErrorDummy()

#        print(f"param: {<bytes>param}.decode('utf-8')")

        # Copy the string using strcpy (be careful with buffer overflows)
        strcpy(out_param, param)
        return out_param

# in common_defs.h:
#define MAXLENGTHOFFILES        1024
#
    def getRootDir(self):
        cdef char *value
#        cdef char *out_value
        cdef int buffer_size = 1024 # Allocate enough space
        # Allocate memory for the destination string
        value = <char*> malloc(buffer_size * sizeof(char))
#        out_value = <char*> malloc(buffer_size * sizeof(char))
#        if out_value is NULL:
#            raise MemoryError("Failed to allocate memory")

        if get_rootDir(&self.cmd,value)== FAILURE:
            raise cBallsSevereErrorDummy()

#        print(f"param: {<bytes>param}.decode('utf-8')")

        # Copy the string using strcpy (be careful with buffer overflows)
#        strcpy(out_value, value)
        return value

#E parameters

#B histograms
    def getrBins(self):
        cdef int sizeHistN
        cdef int index_r

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()

        cdef np.ndarray[DTYPE_t, ndim=1] out_rBins = np.zeros(sizeHistN,'float64')

        if get_rBins(&self.cmd, &self.gd)==FAILURE:
            raise cBallsSevereErrorDummy()

        for index_r in range(sizeHistN):
            out_rBins[index_r] = self.gd.rBins[index_r+1]

        return out_rBins

    def getHistNN(self):
        cdef int sizeHistN
        cdef int index_r

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()

        cdef np.ndarray[DTYPE_t, ndim=1] out_HistNN = np.zeros(sizeHistN,'float64')

        if get_HistNN(&self.cmd, &self.gd)==FAILURE:
            raise cBallsSevereErrorDummy()

        for index_r in range(sizeHistN):
            out_HistNN[index_r] = self.gd.vecPXD[index_r+1]

        return out_HistNN

    def getHistCF(self):
        cdef int sizeHistN
        cdef int index_r

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistCF(&self.cmd, &self.gd)==FAILURE:
            raise cBallsSevereErrorDummy()

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist

    def getHistXi2pcf(self):
        cdef int sizeHistN
        cdef int index_r

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistXi2pcf(&self.cmd, &self.gd)==FAILURE:
            raise cBallsSevereErrorDummy()

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist

    def getHistZetaMsincos(self, int m, int type):
        cdef int sizeHistN
        cdef int index_r
        cdef int sizesqr
        cdef ErrorMsg errmsg

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise cBallsSevereErrorDummy()
        
        sizesqr = sizeHistN*sizeHistN

        rows = sizeHistN
        cols = sizeHistN
        cdef np.ndarray[np.float64_t, ndim=2] matrix = np.zeros((rows, cols), dtype=np.float64)
        
        if get_HistZetaMsincos(&self.cmd, &self.gd, m, type, errmsg)==FAILURE:
            raise cBallsSevereError(errmsg)

        # You can then populate the matrix
        for i in range(1,rows):
            for j in range(1,cols):
                matrix[i, j] = self.gd.matPXD[i][j]

        return matrix
#
#E histograms

#
#E Interfaces to PXD functions
#
