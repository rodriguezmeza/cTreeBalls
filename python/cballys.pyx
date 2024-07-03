"""
.. module:: cballys
    :synopsis: Python wrapper around cTreeBalls
.. moduleauthor:: Mario A. Rodriguez-Meza <marioalberto.rodriguezmeza@gmail.com>

This module defines a class called cballs.

# MAR 27.04.2023: TODO:

"""
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython

import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

from ccballys cimport *

class CosmoError(Exception):
    def __init__(self, message=""):
        self.message = message.decode() if isinstance(message,bytes) else message

    def __str__(self):
        return '\n\nError in cballs: ' + self.message


class CosmoSevereError(CosmoError):
    """
    Raised when cballs failed to understand one or more input parameters.

    This case would not raise any problem in cballs default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood, otherwise the wrong computation would be selected.
    """
    pass

class CosmoComputationError(CosmoError):
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

    cdef int computed # Flag to see if cballys has already computed with the given pars
    cdef int allocated # Flag to see if cballys structs are allocated already
    cdef object _pars # Dictionary of the parameters
    cdef object ncp   # Keeps track of the structures initialized, in view of cleaning.

    property pars:
        def __get__(self):
            return self._pars
    property state:
        def __get__(self):
            return True
    property theta:
        def __get__(self):
            return self.cmd.theta
#    property histN:
#        def __get__(self):
#            return self.gd.histN

    def set_default(self):
        _pars = {
            "searchMethod":"tree-omp-sincos",}
        self.set(**_pars)

    def __cinit__(self, default=False):
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
            raise CosmoSevereError("bad call")
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
            self.fc.read[i] = _FALSE_
            i+=1

    def struct_cleanup(self):
        if(self.allocated != True):
          return
# Add necesary calls to free memory (bodytab, histograms...)
#        if "thermodynamics" in self.ncp:
#            thermodynamics_free(&self.th)
#        if "background" in self.ncp:
#            background_free(&self.ba)
        self.allocated = False
        self.computed = False

    def _check_task_dependency(self, level):
        """
        Fill the level list with all the needed modules

        .. warning::

            the ordering of modules is obviously dependent on cBalls module order
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
        if "MainLoop" in level:
            if "startrun_Common" not in level:
                level.append("set_number_threads")
        if "set_number_threads" in level:
            if "startrun_Common" not in level:
                level.append("startrun_Common")
        if "startrun_Common" in level:
            if "PrintParameterFile" not in level:
                level.append("PrintParameterFile")
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


    def Run(self, level=["EndRun"]):
        """
        compute()

        Main function, execute all routines in cballs.

        Parameters
        ----------
        level : list
                list of the parameters.

        .. warning::

            level default value

        """
        cdef ErrorMsg errmsg

        level = self._check_task_dependency(level)

        if self.computed and self.ncp.issuperset(level):
            return

        if self.allocated:
            self.struct_cleanup()

        self.computed = False

        self._fillparfile()

        self.ncp=set()
        self.allocated = True

        if "input" in level:
            if input_read_from_file(&self.cmd, &self.fc, errmsg) == _FAILURE_:
                raise CosmoSevereError(errmsg)
            self.ncp.add("input")
            problem_flag = False
            problematic_parameters = []
            for i in range(self.fc.size):
                if self.fc.read[i] == _FALSE_:
                    problem_flag = True
                    problematic_parameters.append(self.fc.name[i].decode())
            if problem_flag:
                raise CosmoSevereError(
                    "cballs did not read input parameter(s): %s\n" % ', '.join(
                    problematic_parameters))

        if "startrun_Common" in level:
            if startrun_Common(&(self.cmd), &(self.gd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.op.error_message)
            self.ncp.add("startrun_Common")

        if "PrintParameterFile" in level:
            if PrintParameterFile(&(self.cmd), "cballys_param.txt") == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.op.error_message)
            self.ncp.add("PrintParameterFile")

        if "set_number_threads" in level:
            if set_number_threads(&(self.cmd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.op.error_message)
            self.ncp.add("set_number_threads")

        if "MainLoop" in level:
            if MainLoop(&(self.cmd), &(self.gd)) == _FAILURE_:
                self.struct_cleanup()
                raise CosmoComputationError(self.op.error_message)
            self.ncp.add("MainLoop")


        self.computed = True

        return


    def theta(self):
        return self.cmd.theta

    def sizeHistN(self):
        return self.cmd.sizeHistN

#    def histN(self):
#        return self.gd.histN
