#!/usr/bin/env python
# coding: utf-8


### Module to parse the cballs code
import numpy as np
import sys, platform, os, subprocess

### Modules to compute...
import scipy

path_to_cballs='../../../cTreeBalls/'


def run_cballs(proc = os.getpid(),
            searchMethod = "tree-omp-sincos",
            scanLevelRoot = 0,
            options  = "compute-HistN,and-CF,out-m-HistZeta",
            infile = "../../Eladio/full_sky_whole_XYZK__zs9_r000_nside512.txt",
            rangeN = 0.0633205,
            rminHist = 0.00213811,
            sizeHistN = 20,
            stepState = 1000,
            verbose = 2,
            verbose_log = 2,
            numberThreads = 16
            ):
    
    '''
    Runs the cballs code.

    Args:
        proc: number of process,
        searchMethod: name of the search method,
        scanLevelRoot: where to scan the tree,
        options: options to change code behavior,
        infile: name of the data catalog file,
        rangeN: radial range of histograms,
        rminHist: radial minimum of histograms,
        sizeHistN: number of bins in the histograms,
        stepState: number of pivots completed,
        verbose: verbosity level at the stdout.
        verbose_log: verbosity level at log file.
        numberThreads: cosmological model,


    Returns:
        rbins array.
    '''


    args_exe = [path_to_cballs+'cballs',
                'searchMethod=%s'%searchMethod,
                'scanLevelRoot=%d'%scanLevelRoot,
                'options=%s'%options,
                'infile=%s'%infile,
                'rangeN=%f'%rangeN,
                'rminHist=%f'%rminHist,
                'sizeHistN=%d'%sizeHistN,
                'stepState=%d'%stepState,
                'verbose=%d'%verbose,
                'verbose_log=%d'%verbose_log,
                'numberThreads=%d'%numberThreads
                ]
    print(args_exe)
    
    res_proc = subprocess.call(args_exe)
    
    ### Part to read data output files
    name_read = 'Output/rbins'+'.txt'
    data = np.loadtxt(name_read, skiprows=0)

    return({'r':data[0]})


