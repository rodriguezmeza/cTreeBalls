"""
.. module:: cyballs
    :synopsis: Python wrapper around cTreeBalls
.. moduleauthor:: Mario A. Rodriguez-Meza <marioalberto.rodriguezmeza@gmail.com>

.. based on Julien Lesgourges' CLASS

This module defines a class called cballs.

# MAR 15.02.2026:

"""

from math import exp,log
import numpy as np
from os.path import abspath, dirname
cimport numpy as np
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
import cython
cimport cython
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

import time

import sys
def viewdictitems(d):
    if sys.version_info >= (3,0):
        return d.items()
    else:
        return d.viewitems()

ctypedef np.float64_t DTYPE_t
ctypedef np.int32_t DTYPE_i

from libc.stdio cimport snprintf
from libc.stddef cimport size_t

class CallableFloat(float):
  def __call__(self):
    return self


# Import the .pxd containing definitions
from ccyballs cimport *

DEF _MAXTITLESTRINGLENGTH_ = 8000

__version__ = _VERSION_.decode("utf-8")


def search_method_id(method):
    """Return the compiled search-method id, or ``-1`` when unavailable."""
    cdef bytes encoded

    if not isinstance(method, str):
        raise TypeError("method must be a string")
    if "\0" in method:
        raise ValueError("method contains an embedded NUL byte")
    encoded = method.encode("utf-8")
    return cballs_search_method_id(<const char *>encoded)

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
    understood.
    """
    pass

class CosmoSevereErrorDummy(CosmoSevereError):
    """
    Raised when cballs failed to understand one or more input parameters.

    This case would not raise any problem in cballs default behaviour. However,
    for parameter extraction, one has to be sure that all input parameters were
    understood.
    """
    pass

class CosmoComputationError(CosmoError):
    """
    Raised when cballs could not compute at this point.

    This will be caught by the parameter extraction code to give an extremely
    unlikely value to this point
    """
    pass

cdef inline void safe_copy_cstr(char *dest, size_t dest_size, bytes value, str label) except *:
    cdef size_t value_len = len(value)

    if value_len >= dest_size:
        raise CosmoSevereError(
            f"{label} is too long: {value_len} bytes, max is {dest_size - 1}"
        )

    if b"\0" in value:
        raise CosmoSevereError(
            f"{label} contains an embedded NUL byte"
        )

    memcpy(dest, <const char *> value, value_len)
    dest[value_len] = 0

cdef class cballs:
    """
    cballs wrapping, creates the glue between C and python

    The actual cballs wrapping, the only class we will call from Python
    (indeed the only one we will import, with the command:
    from cyballs import cballs

    """
    # List of used structures, defined in the header file. They have to be
    # "cdefined", because they correspond to C structures
    cdef file_content fc
    cdef cmdline_data cmd
    cdef global_data gd
    cdef cballs_runtime_state *runtime_state

    cdef ErrorMsg error_message

    cdef int nthreads
    cdef double cputime

    cdef int computed
    cdef int allocated
    cdef object _pars
    cdef object _memory_catalogs
    cdef object ncp

    cdef char path_to_this[1000]

    _levellist = ["input","StartRun_Common","PrintParameterFile","SetNumberThreads","Initial", "MainLoop", "EndRun"]

    # Special properties
    @property
    def pars(self):
      return self._pars
    @property
    def state(self):
      return True

#B definition for abi useful check

    cdef void _check_abi(self) except *:
        cdef size_t c_cmd_size = sizeof_cmdline_data()
        cdef size_t c_gd_size = sizeof_global_data()

        if sizeof(cmdline_data) != c_cmd_size:
            raise CosmoSevereError(
                "ABI mismatch for cmdline_data: "
                f"Cython sees {sizeof(cmdline_data)} bytes, "
                f"C library sees {c_cmd_size} bytes"
            )

        if sizeof(global_data) != c_gd_size:
            raise CosmoSevereError(
                "ABI mismatch for global_data: "
                f"Cython sees {sizeof(global_data)} bytes, "
                f"C library sees {c_gd_size} bytes"
            )

    cdef void _activate_runtime(self) except *:
        if (self.runtime_state == NULL or
                cballs_runtime_activate(self.runtime_state) == FAILURE):
            raise CosmoSevereError("cballs runtime state is unavailable")
#E

    # Now we can start with the actual code describing the cballs

    # One important consequence of set_default definition:
    #   if the default searchMethod, infileformat, theta, etc.
    #   are wrong for the active build,
    #   fix them in C input_default_params() / input_default_params.h,
    #   not in cyballs.pyx.
    def set_default(self, **overrides):
        """
        Reset to the compiled C defaults.

        Python only stores user overrides in self._pars. The real defaults come
        from input_default_params(), the same routine used by C input parsing.
        """
        self._activate_runtime()

        if self.allocated:
            self.struct_cleanup()

        self.computed = False
        self.ncp = set()
        self._pars = {}
        self._memory_catalogs = []

        memset(&self.cmd, 0, sizeof(cmdline_data))
        memset(&self.gd, 0, sizeof(global_data))
        self.gd.startrun_cputime = False

        if input_default_params(&self.cmd) == FAILURE:
            raise CosmoSevereError(
                (<char *> self.cmd.error_message).decode("utf-8", "replace")
            )

        if overrides:
            self.set(**overrides)

        return True

    def __cinit__(self, default=True):
        self.runtime_state = NULL
        memset(&self.cmd, 0, sizeof(cmdline_data))
        memset(&self.gd, 0, sizeof(global_data))
        memset(&self.fc, 0, sizeof(file_content))

        self.allocated = False
        self.computed = False
        self._pars = {}
        self._memory_catalogs = []
        self.gd.startrun_cputime = False

        self._check_abi()

        self.runtime_state = cballs_runtime_create()
        if self.runtime_state == NULL:
            raise CosmoSevereError("not enough memory allocating C runtime state")

        self.fc.filename = <char*> malloc(sizeof(char) * 30)
        if self.fc.filename == NULL:
            raise CosmoSevereError("not enough memory allocating fc.filename")

        dumc = "NOFILE"
        safe_copy_cstr(self.fc.filename, 30, b"NOFILE", "fc.filename")
        self.ncp = set()
        if default: self.set_default()
        try:
          import importlib.resources
          resource_path = abspath(importlib.resources.files('cyballs'))
        except ImportError as ie:
          resource_path = dirname(abspath(__file__))
        path_to_this_as_bytes = resource_path.encode()
        safe_copy_cstr(self.path_to_this, sizeof(self.path_to_this),
               path_to_this_as_bytes, "path_to_this")

    def __dealloc__(self):
        if self.runtime_state != NULL:
            cballs_runtime_activate(self.runtime_state)
        if self.allocated:
            cballs_end_run_free_memory_guarded(&self.cmd, &self.gd)
            self.allocated = False
            self.computed = False
        # Reset all the fc to zero if its not already done
        if self.fc.size !=0:
            self.fc.size=0
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
        free(self.fc.filename)
        if self.runtime_state != NULL:
            cballs_runtime_destroy(self.runtime_state)
            self.runtime_state = NULL

    # Set up the dictionary
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

    # Create an equivalent of the parameter file. Non specified values will be
    # taken at their default (in cballs)
    def _fillparfile(self):
        cdef int new_size
        cdef int i
        cdef int memory_parameter_count
        cdef FileArg *new_name = NULL
        cdef FileArg *new_value = NULL
        cdef short *new_read = NULL

        if self.fc.size != 0:
            free(self.fc.name)
            free(self.fc.value)
            free(self.fc.read)
            self.fc.name = NULL
            self.fc.value = NULL
            self.fc.read = NULL
            self.fc.size = 0

        if self._memory_catalogs and (
                'infile' in self._pars or 'infileformat' in self._pars):
            raise CosmoSevereError(
                "in-memory catalogs cannot be combined with infile or infileformat"
            )

        memory_parameter_count = 2 if self._memory_catalogs else 0
        new_size = len(self._pars) + memory_parameter_count
        if 'base_path' not in self._pars:
            new_size += 1

        try:
            new_name = <FileArg*> malloc(sizeof(FileArg) * new_size)
            if new_name == NULL:
                raise CosmoSevereError("not enough memory allocating fc.name")

            new_value = <FileArg*> malloc(sizeof(FileArg) * new_size)
            if new_value == NULL:
                raise CosmoSevereError("not enough memory allocating fc.value")

            new_read = <short*> malloc(sizeof(short) * new_size)
            if new_read == NULL:
                raise CosmoSevereError("not enough memory allocating fc.read")

            i = 0
            for kk in self._pars:
                dumcp = kk.strip().encode()
                safe_copy_cstr(<char *> new_name[i], sizeof(FileArg), dumcp, "parameter name")

                dumcp = str(self._pars[kk]).strip().encode()
                safe_copy_cstr(<char *> new_value[i], sizeof(FileArg), dumcp, "parameter value")

                new_read[i] = FALSE
                i += 1

            if self._memory_catalogs:
                memory_names = ",".join(
                    "m%d" % catalog
                    for catalog in range(len(self._memory_catalogs))
                ).encode()
                memory_formats = ",".join(
                    "lya-ascii" if self._memory_catalogs[catalog][6] is not None else "memory"
                    for catalog in range(len(self._memory_catalogs))
                ).encode()

                safe_copy_cstr(<char *> new_name[i], sizeof(FileArg),
                               b"infile", "parameter name")
                safe_copy_cstr(<char *> new_value[i], sizeof(FileArg),
                               memory_names, "in-memory catalog descriptors")
                new_read[i] = FALSE
                i += 1

                safe_copy_cstr(<char *> new_name[i], sizeof(FileArg),
                               b"infileformat", "parameter name")
                safe_copy_cstr(<char *> new_value[i], sizeof(FileArg),
                               memory_formats, "in-memory catalog formats")
                new_read[i] = FALSE
                i += 1

            if 'base_path' not in self._pars:
                safe_copy_cstr(<char *> new_name[i], sizeof(FileArg), b"base_path", "parameter name")
                safe_copy_cstr(
                    <char *> new_value[i],
                    sizeof(FileArg),
                    (<char *> self.path_to_this)[:strlen(self.path_to_this)],
                    "base_path"
                )
                new_read[i] = FALSE

            self.fc.name = new_name
            self.fc.value = new_value
            self.fc.read = new_read
            self.fc.size = new_size

        except Exception:
            if new_name != NULL:
                free(new_name)
            if new_value != NULL:
                free(new_value)
            if new_read != NULL:
                free(new_read)
            raise

    # Called at the end of a run, to free memory
    def struct_cleanup(self):
        if(self.allocated != True):
            return

        self._activate_runtime()

        if cballs_end_run_free_memory_guarded(&(self.cmd), &(self.gd)) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        self.ncp = set()
        self.allocated = False
        self.computed = False
        

    def clean_all(self):
        self.struct_cleanup()
        self.clean()
        self._memory_catalogs = []

    cdef object _prepare_catalog_vector(self, object values,
                                        Py_ssize_t nbody, str name):
        cdef object result

        if values is None:
            return None
        try:
            result = np.ascontiguousarray(values, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise CosmoSevereError(
                f"{name} must be convertible to float64"
            ) from exc
        if result.ndim != 1 or result.shape[0] != nbody:
            raise CosmoSevereError(
                f"{name} must have shape ({nbody},), got {result.shape}"
            )
        if not np.isfinite(result).all():
            raise CosmoSevereError(f"{name} contains non-finite values")
        return result

    def set_catalog(self, positions, kappa=None, weights=None, mask=None,
                    gamma1=None, gamma2=None, int catalog=0):
        """Register one NumPy catalog for the next run.

        Parameters
        ----------
        positions : array_like, shape (N, compiled_ndim)
            Cartesian positions. They are copied into C-owned storage when
            :meth:`Run` starts, so the input array is never modified.
        kappa, weights : array_like, shape (N,), optional
            Scalar field and statistical weights. Omitted arrays default to 1.
        mask : array_like, shape (N,), optional
            Boolean or 0/1 validity mask. Omitted values are all valid.
        gamma1, gamma2 : array_like, shape (N,), optional
            Shear components for a THREEPCFSHEAR build. Supply both or neither.
        catalog : int, optional
            Zero-based catalog slot. Catalogs must be contiguous from zero.
        """
        cdef object positions_array
        cdef object kappa_array
        cdef object weights_array
        cdef object mask_array
        cdef object gamma1_array
        cdef object gamma2_array
        cdef Py_ssize_t nbody
        cdef int ndim = cballs_compiled_ndim()

        if catalog < 0:
            raise CosmoSevereError("catalog must be non-negative")
        if catalog >= cballs_max_memory_catalogs():
            raise CosmoSevereError(
                f"catalog must be below {cballs_max_memory_catalogs()}"
            )
        if catalog > len(self._memory_catalogs):
            raise CosmoSevereError(
                "catalogs must be registered contiguously from index 0"
            )
        if (gamma1 is None) != (gamma2 is None):
            raise CosmoSevereError(
                "gamma1 and gamma2 must either both be supplied or both omitted"
            )

        try:
            positions_array = np.ascontiguousarray(positions, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise CosmoSevereError(
                "positions must be convertible to a float64 array"
            ) from exc
        if positions_array.ndim != 2 or positions_array.shape[1] != ndim:
            raise CosmoSevereError(
                f"positions must have shape (N, {ndim}), got {positions_array.shape}"
            )
        nbody = positions_array.shape[0]
        if nbody < 3:
            raise CosmoSevereError("an in-memory catalog needs at least 3 bodies")
        if not np.isfinite(positions_array).all():
            raise CosmoSevereError("positions contains non-finite values")

        kappa_array = self._prepare_catalog_vector(kappa, nbody, "kappa")
        weights_array = self._prepare_catalog_vector(weights, nbody, "weights")
        gamma1_array = self._prepare_catalog_vector(gamma1, nbody, "gamma1")
        gamma2_array = self._prepare_catalog_vector(gamma2, nbody, "gamma2")

        if mask is None:
            mask_array = None
        else:
            try:
                mask_array = np.ascontiguousarray(mask)
            except (TypeError, ValueError) as exc:
                raise CosmoSevereError(
                    "mask must be convertible to a one-dimensional array"
                ) from exc
            if mask_array.ndim != 1 or mask_array.shape[0] != nbody:
                raise CosmoSevereError(
                    f"mask must have shape ({nbody},), got {mask_array.shape}"
                )
            if not np.all((mask_array == 0) | (mask_array == 1)):
                raise CosmoSevereError("mask values must be boolean or 0/1")
            mask_array = np.ascontiguousarray(mask_array, dtype=np.uint8)

        if self.allocated:
            self.struct_cleanup()
        entry = (positions_array, kappa_array, weights_array, mask_array,
                 gamma1_array, gamma2_array, None)
        if catalog == len(self._memory_catalogs):
            self._memory_catalogs.append(entry)
        else:
            self._memory_catalogs[catalog] = entry
        self.computed = False
        self.ncp = set()
        return True

    def set_forest_catalog(self, positions, delta, weights, forest_ids):
        """Register observer-centered Ly-alpha pixels without reading a file.

        Coordinates are comoving Cartesian distances. IDs must be integer
        arrays, never floating point (DESI IDs can exceed 2**53). Equal IDs
        denote the same quasar. IDs are compacted internally to fit INTEGER;
        original arrays remain unchanged. Requires a 3D Ly-alpha addon build.
        """
        ids = np.asarray(forest_ids)
        pos = np.asarray(positions)
        if ids.dtype.kind not in "iu" or ids.ndim != 1:
            raise CosmoSevereError("forest_ids must be a one-dimensional integer array")
        if pos.ndim != 2 or pos.shape[1] != 3 or ids.size != pos.shape[0]:
            raise CosmoSevereError("positions must be (N, 3) with N forest IDs")
        distance = np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])
        if not np.all(np.isfinite(distance) & (distance > 0)):
            raise CosmoSevereError("forest coordinates need positive finite observer distances")
        w = np.asarray(weights)
        if w.ndim != 1 or w.size != ids.size or not np.all(np.isfinite(w) & (w >= 0)):
            raise CosmoSevereError("forest weights must be finite non-negative values")
        if len(self._memory_catalogs) > 1:
            raise CosmoSevereError("forest input requires exactly one catalog")
        dense_ids = np.ascontiguousarray(np.unique(ids, return_inverse=True)[1],
                                        dtype=np.int64)
        self.set_catalog(positions, kappa=delta, weights=weights)
        self._memory_catalogs[0] = self._memory_catalogs[0][:6] + (dense_ids,)
        return True

    def clear_catalogs(self):
        """Remove registered in-memory catalogs and release any active run."""
        if self.allocated:
            self.struct_cleanup()
        self._memory_catalogs = []
        self.computed = False
        self.ncp = set()

    @property
    def catalog_count(self):
        """Number of in-memory catalogs registered on this object."""
        return len(self._memory_catalogs)

    cdef void _load_memory_catalogs(self) except *:
        cdef int ifile
        cdef int status
        cdef object entry
        cdef np.ndarray positions_array
        cdef np.ndarray kappa_array
        cdef np.ndarray weights_array
        cdef np.ndarray mask_array
        cdef np.ndarray gamma1_array
        cdef np.ndarray gamma2_array
        cdef np.ndarray forest_ids_array
        cdef const double *kappa_ptr
        cdef const double *weights_ptr
        cdef const unsigned char *mask_ptr
        cdef const double *gamma1_ptr
        cdef const double *gamma2_ptr

        for ifile in range(len(self._memory_catalogs)):
            entry = self._memory_catalogs[ifile]
            positions_array = entry[0]
            kappa_ptr = NULL
            weights_ptr = NULL
            mask_ptr = NULL
            gamma1_ptr = NULL
            gamma2_ptr = NULL

            if entry[1] is not None:
                kappa_array = entry[1]
                kappa_ptr = <const double *>np.PyArray_DATA(kappa_array)
            if entry[2] is not None:
                weights_array = entry[2]
                weights_ptr = <const double *>np.PyArray_DATA(weights_array)
            if entry[3] is not None:
                mask_array = entry[3]
                mask_ptr = <const unsigned char *>np.PyArray_DATA(mask_array)
            if entry[4] is not None:
                gamma1_array = entry[4]
                gamma1_ptr = <const double *>np.PyArray_DATA(gamma1_array)
                gamma2_array = entry[5]
                gamma2_ptr = <const double *>np.PyArray_DATA(gamma2_array)

            status = cballs_load_memory_catalog(
                &self.cmd, &self.gd, ifile,
                <const double *>np.PyArray_DATA(positions_array),
                <size_t>positions_array.shape[0],
                <int>positions_array.shape[1],
                kappa_ptr, weights_ptr, mask_ptr, gamma1_ptr, gamma2_ptr)
            if status == FAILURE:
                raise CosmoSevereError(
                    (<char *>self.cmd.error_message).decode("utf-8", "replace")
                )
            if entry[6] is not None:
                forest_ids_array = entry[6]
                status = cballs_set_memory_forest_ids(
                    &self.cmd, &self.gd, ifile,
                    <const int64_t *>np.PyArray_DATA(forest_ids_array),
                    <size_t>positions_array.shape[0])
                if status == FAILURE:
                    raise CosmoSevereError(
                        (<char *>self.cmd.error_message).decode("utf-8", "replace")
                    )

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
        # If it's a string only, treat as a list
        if isinstance(level, str):
          level=[level]
        # For each item in the list
        levelset = set()
        for item in level:
          # If the item is not in the list of allowed levels, make error message
          if item not in self._levellist:
            raise CosmoSevereError("Unknown computation level: '{}'".format(item))
          # Otherwise, add to list of levels up to and including the specified level
          levelset.update(self._levellist[:self._levellist.index(item)+1])
        return levelset

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

    def Run(self, level=["MainLoop"]):
        """
        Run(level=["MainLoop"])

        Main function, execute all the methods for all desired modules.
        This is called in Python, and this ensures that the cballs instance
        of this class contains all the relevant quantities. Then, one can deduce
        integral quantities, etc...

        Parameters
        ----------
        level : list
                list of the last module desired. The internal function
                _check_task_dependency will then add to this list all the
                necessary modules to compute in order to initialize this last
                one. The default last module is "MainLoop" so Python getters can inspect live C arrays.

        .. warning::

            level default value should be left as an array (it was creating
            problem when casting as a set later on, in _check_task_dependency)

        """
        cdef ErrorMsg errmsg

        # Append to the list level all the modules necessary to compute.
        level = self._check_task_dependency(level)

        # Check if this function ran before (self.computed should be true), and
        # if no other modules were requested, i.e. if self.ncp contains (or is
        # equivalent to) level. If it is the case, simply stop the execution of
        # the function.
        if self.computed and self.ncp.issuperset(level):
            return self.cputime

        self._activate_runtime()

        # Check if already allocated to prevent memory leaks
        if self.allocated:
            self.struct_cleanup()

        # Otherwise, proceed with the normal computation.
        self.computed = False

        # Equivalent of writing a parameter file
        self._fillparfile()

        # self.ncp will contain the list of computed modules (under the form of
        # a set, instead of a python list)
        self.ncp=set()

#B correction
        # Up until the empty set, all modules are allocated
        # (And then we successively keep track of the ones we allocate additionally)
        self.allocated = True

        try:
            # --------------------------------------------------------------------
            # Check the presence for all cballs modules in the list 'level'. If a
            # module is found in level, execute its routine.
            # --------------------------------------------------------------------
            # The input module should raise a CosmoSevereError, because
            # non-understood parameters asked to the wrapper is a problematic
            # situation.
            if "input" in level:
                if input_read_from_file_guarded(&self.cmd, &self.gd, &self.fc, errmsg) == FAILURE:
                    raise CosmoSevereError(errmsg)
                self.ncp.add("input")

                problem_flag = False
                problematic_parameters = []
                for i in range(self.fc.size):
                    if self.fc.read[i] == FALSE:
                        problem_flag = True
                        problematic_parameters.append((<char *> self.fc.name[i]).decode("utf-8"))

                if problem_flag:
                    raise CosmoSevereError(
                        "cballs did not read input parameter(s): %s\n" %
                        ', '.join(problematic_parameters)
                    )

                if self._memory_catalogs:
                    self._load_memory_catalogs()

            # The following list of computation is straightforward. If the "_init"
            # methods fail, call `struct_cleanup` and raise a CosmoComputationError
            # with the error message from the faulty module of CLASS.
            if "StartRun_Common" in level:
                if cballs_start_run_common_guarded(&(self.cmd), &(self.gd)) == FAILURE:
                    raise CosmoComputationError((<char *> self.cmd.error_message).decode("utf-8", "replace"))
                self.ncp.add("StartRun_Common")

            # keep the rest of the C stages here too

            if "PrintParameterFile" in level:
                if cballs_print_parameter_file_guarded(&(self.cmd), &(self.gd), "cyballs_param.txt") == FAILURE:
                    raise CosmoComputationError((<char *> self.cmd.error_message).decode("utf-8", "replace"))
                self.ncp.add("PrintParameterFile")

            if "SetNumberThreads" in level:
                if cballs_set_number_threads_guarded(&(self.cmd)) == FAILURE:
                    raise CosmoComputationError((<char *> self.cmd.error_message).decode("utf-8", "replace"))
                self.ncp.add("SetNumberThreads")
                self.nthreads=self.getNThreads()

            if "MainLoop" in level:
                start_wall_time_p = time.process_time()
                if cballs_main_loop_guarded(&(self.cmd), &(self.gd)) == FAILURE:
                    raise CosmoComputationError((<char *> self.cmd.error_message).decode("utf-8", "replace"))
                self.ncp.add("MainLoop")
                end_wall_time_p = time.process_time()
                self.cputime = (end_wall_time_p - start_wall_time_p)/self.nthreads

            if "EndRun" in level:
                if cballs_end_run_guarded(&(self.cmd), &(self.gd)) == FAILURE:
                    raise CosmoComputationError((<char *> self.cmd.error_message).decode("utf-8", "replace"))
                self.ncp.add("EndRun")
                self.allocated = False

        except Exception:
            self.struct_cleanup()
            raise
#E

        self.computed = True

        return self.cputime

        # At this point, the cballs instance contains everything needed. The
        # following functions are only to output the desired numbers

#
#
#B cballs definitions
#

    def abi_sizes(self):
        return {
            "cmdline_data": sizeof_cmdline_data(),
            "global_data": sizeof_global_data(),
        }

    def _set_allocation_failure_after_for_tests(self, long successful_allocations=0):
        """Inject one C allocation failure; intended only for regression tests."""
        cballs_test_fail_allocation_after(successful_allocations)

    def _reset_allocation_failure_for_tests(self):
        """Disable the C allocation-failure regression hook."""
        cballs_test_reset_allocation_failure()

    def _runtime_bodytable_address(self, int ifile=0):
        """Return an internal allocation address for ownership regression tests."""
        return <size_t> cballs_runtime_bodytable_at(self.runtime_state, ifile)

#
#B Interfaces to PXD functions
#
#B flags
#E

#B parameters

# in common_defs.h:
#
    def getNThreads(self):
        cdef int value
        if get_nthreads(&self.cmd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value

#E parameters

#B FROM OLD version

    def getTreeAllocated(self):
        cdef short value
        cdef short out_value
        if get_tree_allocated(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getAllocated2(self):
        cdef short value
        cdef short out_value
        if get_allocated_2(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getBodytableAllocated(self):
        cdef short value
        cdef short out_value
        if get_bodytable_allocated(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getHistogramsAllocated(self):
        cdef short value
        cdef short out_value
        if get_histograms_allocated(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getGDAllocated(self):
        cdef short value
        cdef short out_value
        if get_gd_allocated(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getCMDAllocated(self):
        cdef short value
        cdef short out_value
        if get_cmd_allocated(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

#E

#B parameters
#
    #B from cBalls
    def getScanLevel(self):
        cdef int value
        if get_scanLevel(&self.cmd, &value) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value

    def getScanLevelRoot(self):
        cdef int value
        if get_scanLevelRoot(&self.cmd, &value) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value

    def getScanLevelMin(self):
        cdef int value
        if get_scanLevelMin(&self.cmd, &value) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value
    #

    def getnMultipoles(self):
        cdef int value
        if get_nmultipoles(&self.cmd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value

    def getTheta(self):
        cdef double value
        if get_theta(&self.cmd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value

    def getrsmooth(self):
        cdef double value
        cdef double out_value
        if get_rsmooth(&self.gd,&value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_value = value
        return out_value

    def getCPUTime(self):
        cdef double cputime
        cdef double out_cputime
        if get_cputime(&self.gd,&cputime)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_cputime = cputime
        return out_cputime

    def getsizeHistN(self):
        cdef int sizeHistN
        cdef int out_sizeHistN
        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        out_sizeHistN = sizeHistN
        return out_sizeHistN

#B
    def getVersion(self):
        cdef char param[64]

        if get_version(&self.cmd, param) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        
        return (<char *>param).decode("utf-8", "replace")
#E

    def getRootDir(self):
        cdef char value[1024]

        if get_rootDir(&self.cmd, value) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        return (<char *>value).decode("utf-8", "replace")

#E parameters

#B parameters
    def getNBody(self):
        cdef INTEGER value
        if get_nbody(&self.cmd, &self.gd, &value)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        return value
#E parameters

#B histograms

#B added by cBalls
    cdef void _require_live_histograms(self) except *:
        cdef short value

        if self.allocated != True or "MainLoop" not in self.ncp or "EndRun" in self.ncp:
            raise CosmoSevereError(
                'PXD histogram getters require live arrays; call Run(level=["MainLoop"]) before reading them, and clean after reading.'
            )

        if get_histograms_allocated(&self.gd, &value) == FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        if value != True:
            raise CosmoSevereError(
            'C histogram arrays are not allocated; call Run(level=["MainLoop"]) before histogram getters.'
            )

    cdef void _require_shear_results(self) except *:
        self._require_live_histograms()

        if (self.gd.histShearXiPlusRe == NULL or
                self.gd.histShearXiPlusIm == NULL or
                self.gd.histShearXiMinusRe == NULL or
                self.gd.histShearXiMinusIm == NULL or
                self.gd.histShearXiWeight == NULL or
                self.gd.histShearGammaNumeratorRe == NULL or
                self.gd.histShearGammaNumeratorIm == NULL or
                self.gd.histShearGammaMultipoleRe == NULL or
                self.gd.histShearGammaMultipoleIm == NULL or
                self.gd.histShearDenominatorRe == NULL or
                self.gd.histShearDenominatorIm == NULL or
                self.gd.histShearGammaRe == NULL or
                self.gd.histShearGammaIm == NULL):
            raise CosmoSevereError(
                'Shear results are unavailable; run searchMethod="octree-shear-omp" first.'
            )

        if (self.cmd.sizeHistN <= 0 or self.gd.shearMultipoleMax < 0 or
                self.gd.shearAngularBins <= 0 or
                self.gd.shearAngularBins != self.cmd.sizeHistPhi):
            raise CosmoSevereError('Shear result dimensions are inconsistent.')
#E

    def getrBins(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_rBins = np.zeros(sizeHistN,'float64')

        if get_rBins(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_rBins[index_r] = self.gd.rBins[index_r+1]

        return out_rBins

    def getHistNN(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_HistNN = np.zeros(sizeHistN,'float64')

        if get_HistNN(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_HistNN[index_r] = self.gd.vecPXD[index_r+1]

        return out_HistNN

    def getHistCF(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistCF(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist

    def getHistXi2pcf(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistXi2pcf(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist

#B cross
    def getHistXi2pcf12(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistXi2pcf12(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist

    def getHistXi2pcf13(self):
        cdef int sizeHistN
        cdef int index_r

        self._require_live_histograms()

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        cdef np.ndarray[DTYPE_t, ndim=1] out_hist = np.zeros(sizeHistN,'float64')

        if get_HistXi2pcf13(&self.cmd, &self.gd)==FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        for index_r in range(sizeHistN):
            out_hist[index_r] = self.gd.vecPXD[index_r+1]

        return out_hist
#E

    def getShearXiPlus(self):
        """Return the weighted flat-sky shear :math:`\\xi_+` bins."""
        cdef int bins
        cdef int index
        cdef np.ndarray[np.float64_t, ndim=1] real_part
        cdef np.ndarray[np.float64_t, ndim=1] imag_part

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        real_part = np.empty(bins, dtype=np.float64)
        imag_part = np.empty(bins, dtype=np.float64)
        for index in range(bins):
            real_part[index] = self.gd.histShearXiPlusRe[index]
            imag_part[index] = self.gd.histShearXiPlusIm[index]
        return real_part + 1j*imag_part

    def getShearXiMinus(self):
        """Return the weighted flat-sky shear :math:`\\xi_-` bins."""
        cdef int bins
        cdef int index
        cdef np.ndarray[np.float64_t, ndim=1] real_part
        cdef np.ndarray[np.float64_t, ndim=1] imag_part

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        real_part = np.empty(bins, dtype=np.float64)
        imag_part = np.empty(bins, dtype=np.float64)
        for index in range(bins):
            real_part[index] = self.gd.histShearXiMinusRe[index]
            imag_part[index] = self.gd.histShearXiMinusIm[index]
        return real_part + 1j*imag_part

    def getShearXiWeight(self):
        """Return the pair-weight denominator used for both shear 2PCFs."""
        cdef int bins
        cdef int index
        cdef np.ndarray[np.float64_t, ndim=1] result

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        result = np.empty(bins, dtype=np.float64)
        for index in range(bins):
            result[index] = self.gd.histShearXiWeight[index]
        return result

    cdef object _copy_shear_complex(self, double *source_re,
                                    double *source_im, Py_ssize_t count,
                                    tuple shape):
        cdef Py_ssize_t index
        cdef np.ndarray[np.float64_t, ndim=1] real_part
        cdef np.ndarray[np.float64_t, ndim=1] imag_part

        if source_re == NULL or source_im == NULL:
            raise CosmoSevereError('Requested shear result was released.')
        real_part = np.empty(count, dtype=np.float64)
        imag_part = np.empty(count, dtype=np.float64)
        for index in range(count):
            real_part[index] = source_re[index]
            imag_part[index] = source_im[index]
        return (real_part + 1j*imag_part).reshape(shape)

    def getShearUpsilonMultipoles(self):
        """Backward-compatible alias for ``getShearUpsilonXMultipoles()``."""
        cdef int bins
        cdef int multipoles
        cdef Py_ssize_t count

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        multipoles = 2*self.gd.shearMultipoleMax + 1
        count = <Py_ssize_t>4*multipoles*bins*bins
        return self._copy_shear_complex(
            self.gd.histShearGammaNumeratorRe,
            self.gd.histShearGammaNumeratorIm,
            count, (4, multipoles, bins, bins))

    def getShearUpsilonXMultipoles(self):
        """Return raw Porth x-projection 3PCF numerator multipoles."""
        return self.getShearUpsilonMultipoles()

    def getShearGammaMultipoles(self):
        """Backward-compatible alias for ``getShearGammaXMultipoles()``.

        The first axis contains :math:`\\Gamma_0` through
        :math:`\\Gamma_3`; the second contains orders ``-nmax..nmax``.
        """
        cdef int bins
        cdef int multipoles
        cdef Py_ssize_t count

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        multipoles = 2*self.gd.shearMultipoleMax + 1
        count = <Py_ssize_t>4*multipoles*bins*bins
        return self._copy_shear_complex(
            self.gd.histShearGammaMultipoleRe,
            self.gd.histShearGammaMultipoleIm,
            count, (4, multipoles, bins, bins))

    def getShearGammaXMultipoles(self):
        """Return edge-corrected Porth x-projection 3PCF multipoles.

        The first axis contains :math:`\\Gamma_0^{\\times}` through
        :math:`\\Gamma_3^{\\times}`; the second contains orders
        ``-nmax..nmax``.
        """
        return self.getShearGammaMultipoles()

    def getShearWindowMultipoles(self):
        """Return window multipoles ``N_n`` with shape ``(4*nmax+1, B, B)``."""
        cdef int bins
        cdef int multipoles
        cdef Py_ssize_t count

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        multipoles = 4*self.gd.shearMultipoleMax + 1
        count = <Py_ssize_t>multipoles*bins*bins
        return self._copy_shear_complex(
            self.gd.histShearDenominatorRe,
            self.gd.histShearDenominatorIm,
            count, (multipoles, bins, bins))

    def getShearGamma(self):
        """Backward-compatible alias for ``getShearGammaX()``."""
        cdef int bins
        cdef int phi_bins
        cdef Py_ssize_t count

        self._require_shear_results()
        bins = self.cmd.sizeHistN
        phi_bins = self.gd.shearAngularBins
        count = <Py_ssize_t>4*phi_bins*bins*bins
        return self._copy_shear_complex(
            self.gd.histShearGammaRe, self.gd.histShearGammaIm,
            count, (4, phi_bins, bins, bins))

    def getShearGammaX(self):
        """Return angular Porth x-projection 3PCFs as ``(4, P, B, B)``."""
        return self.getShearGamma()

    def getShearMultipoleOrders(self):
        """Return the order axis corresponding to shear 3PCF multipoles."""
        self._require_shear_results()
        return np.arange(-self.gd.shearMultipoleMax,
                         self.gd.shearMultipoleMax + 1, dtype=np.int32)

    def getShearPhiCenters(self):
        """Return angular-bin centers in radians on ``[-pi, pi)``."""
        cdef int phi_bins

        self._require_shear_results()
        phi_bins = self.gd.shearAngularBins
        return -np.pi + (np.arange(phi_bins, dtype=np.float64) + 0.5) \
            *(2.0*np.pi/phi_bins)


    def getHistZetaMsincos(self, int m, int type):
        self._require_live_histograms()

        cdef int sizeHistN
        cdef int index_r
        cdef int sizesqr
        cdef short computeTPCF
        cdef ErrorMsg errmsg

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        
        sizesqr = sizeHistN*sizeHistN

        rows = sizeHistN
        cols = sizeHistN
        cdef np.ndarray[np.float64_t, ndim=2] matrix = np.zeros((rows, cols), dtype=np.float64)

        if get_computeTPCF(&self.cmd, &self.gd, &computeTPCF)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        if computeTPCF==0:
            return matrix

        if get_HistZetaMsincos(&self.cmd, &self.gd, m, type, errmsg)==FAILURE:
            raise CosmoSevereError(errmsg)

        # You can then populate the matrix
        for i in range(0,rows):
            for j in range(0,cols):
                matrix[i, j] = self.gd.matPXD[i][j]

        return matrix

    def getHistZetaM_EE(self, int m):
        self._require_live_histograms()

        cdef int sizeHistN
        cdef int index_r
        cdef int sizesqr
        cdef short computeTPCF
        cdef ErrorMsg errmsg

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        
        sizesqr = sizeHistN*sizeHistN

        rows = sizeHistN
        cols = sizeHistN
        cdef np.ndarray[np.float64_t, ndim=2] matrix = np.zeros((rows, cols), dtype=np.float64)

        if get_computeTPCF(&self.cmd, &self.gd, &computeTPCF)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        if computeTPCF==0:
            return matrix

        if get_HistZetaM_EE(&self.cmd, &self.gd, m, errmsg)==FAILURE:
            raise CosmoSevereError(errmsg)

        # You can then populate the matrix
        for i in range(0,rows):
            for j in range(0,cols):
                matrix[i, j] = self.gd.matPXD[i][j]

        return matrix

    def getHistZetaM_EE_Im(self, int m):
        self._require_live_histograms()

        cdef int sizeHistN
        cdef short computeTPCF
        cdef ErrorMsg errmsg

        if get_sizeHistN(&self.cmd,&sizeHistN)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))

        rows = sizeHistN
        cols = sizeHistN
        cdef np.ndarray[np.float64_t, ndim=2] matrix = np.zeros(
            (rows, cols), dtype=np.float64)

        if get_computeTPCF(&self.cmd, &self.gd, &computeTPCF)== FAILURE:
            raise CosmoSevereErrorDummy((<char *> self.cmd.error_message).decode("utf-8", "replace"))
        if computeTPCF==0:
            return matrix

        if get_HistZetaM_EE_Im(&self.cmd, &self.gd, m, errmsg)==FAILURE:
            raise CosmoSevereError(errmsg)

        for i in range(0,rows):
            for j in range(0,cols):
                matrix[i, j] = self.gd.matPXD[i][j]
        return matrix

    def getHistZetaM_EE_complex(self, int m):
        return self.getHistZetaM_EE(m) + 1j*self.getHistZetaM_EE_Im(m)
#
#E histograms

#E FROM OLD version

#
#E Interfaces to PXD functions
#

#E cballs definitions
