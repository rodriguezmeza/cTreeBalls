---
title: 'cTreeBalls: a fast 3-point correlation function code for clustering measurements'
tags:
  - Python
  - C
  - Cython
  - cosmology
  - astronomy
authors:
  - name: Mario A. Rodriguez-Meza
    orcid: 0000-0003-1160-1488
    affiliation: 1       # (add quotes to multiple affiliations)
    corresponding: true
  - name: Eladio A. Moreno-Alcala
    orcid: 0000-0002-5400-2584
    affiliation: 2       # (add quotes to multiple affiliations)
    corresponding: true
  - name: Alejandro Aviles
    orcid: 0000-0001-5998-3986
    affiliation: 3       # (add quotes to multiple affiliations)
    corresponding: true
  - name: Gustavo Niz
    orcid: 0000-0002-1544-8946
    affiliation: 2       # (add quotes to multiple affiliations)
    corresponding: true
affiliations:
  - name: Departamento de Física, Instituto Nacional de Investigaciones Nucleares,
          Apartado Postal 18-1027,
          Col. Escandón, Ciudad de México, 11801, México
    index: 1
  - name: Departamento de Ciencias e Ingenierías, Universidad de Guanajuato,
          37150, León, Guanajuato, México,
    index: 2
  - name: Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México,
          62201, Cuernavaca, Morelos, México
    index: 3
date: 20 February 2026
bibliography: paper.bib


---

# Summary

`cTreeBalls`[^1] (cBalls for short) is a Python/C package useful to measure 
(2,3)-point clustering statistics. `cBalls` can efficiently calculate 
3 point correlations of more than 200 million data points (a full sky 
simulation with Nside=4096) in  less than 10 minutes on a single high-performance 
computing node, enabling a feasible analysis for the upcoming LSST data. 
It builds upon octree and kd-tree algorithms [-@Barnes:1986;-@Bentley:1975],
and supplies a user-friendly interface with flexible input/output (I/O) of
catalogue data and measurement results, with the built program
configurable through external parameter files and tracked through enhanced
logging and warning/exception handling. For completeness and
complementarity, methods for measuring two-point clustering statistics 
for periodic boxes are also included in the package.


[^1]: [github.com/rodriguezmeza/cTreeBalls](https://github.com/rodriguezmeza/cTreeBalls)


# Statement of need

Correlation functions are frequently used to extract relevant information about the large scale
structure of the universe, which in turn can be used to discriminate between cosmological models.
A common practice is to employ two point statistics, such as the two point correlation function
(2PCF) or its Fourier space counterpart, the power spectrum. A random gaussian field is fully
characterized by the two point statistics, however, the non-linear cosmic evolution of matter
produces a measurable deviation from this random gaussian field assumption. Equivalently,
primordial non-gaussianities may also imprint a non-gaussian signal. It is therefore crucial to
use alternative methods to extract the non-gaussian information that is not captured by the
two point statistics. Higher order point correlations are frequently used as statistical tools to
extract such information. The three point correlation function (3PCF), or the bispectrum in
Fourier space, is the next order term in the hierarchy of n-point correlation functions based
on signal to noise ratio on sufficiently large scales. Just as the 2PCF describes characteristic
length-scales of data, the 3PCF encodes all the spatial triangles in that distribution. As a
result, an isotropic 3PCF depends on three scales that determine the side lengths of these
triangles.


Specifically, `cBalls` can calculate:

  * three-point clustering statistics, namely multipoles of the three-point 
    correlation function
    (3PCF) in configuration space for projected fields;

  * two-point clustering statistics, namely two-point correlation function 
    (2PCF), for projected field and for a cubic-box simulation.



# Implementation

A binned estimator to study the 3PCF signal of the weak lensing convergence field was
introduced in Arvizu et al., 2025. Extract this kind of signal is prohibitively CPU time
consuming when we are dealing with a vast amount of sources. Brute force algorithms scales
with the number of galaxies as O(N 3). We can speed up if we introduce a tree structure
(kd-tree or octree) for the spatial part of the catalog map. We are using mainly octree
structures (Barnes & Hut, 1986), but kd-tree is also used (Bentley, 1975). Now, the real
improvement comes when we use a harmonic decomposition of the weak lensing convergence
field and impose homogeneity and isotropy, and using the fact the field is a spin-0 field (a real
number). In such a case complexity of the code goes as O(N log N ) as it does a 2PCF code.
Multipoles are the estimators we are looking for that characterize the three-point statistics.



# Features

`cBalls` has many features that makes it appropriate to be integrated in massive workflow like
TXPipe (Prat et al., 2023) used in Vera C. Rubin Observatory with its Legacy Survey of Space
and Time (LSST) collaboration[^3]:

* The code has a Python interface (cyballs) that was implemented using Cython to bind
the C searching engine. Python scripts are given in such a way users can run the code
in a linux terminal or can be run using a Jupyter notebook. The C code can also be
compiled and executed independently of Python.

* Measurement pipelines can be configured through external parameter files
  (in the YAML format for the Python program), cleanly separating user
  inputs from the program itself. Alternatively, measurement parameters can
  be set for individual Python methods without the use of a parameter file.

* Parameters files can be given as standard ASCII files or in command line. If using a
Jupyter Notebook or a Python script, parameters are given very easily as a very simple
dictionary without the need of a parameter file. Before computing correlations parameter
values are checked to secure the computation.

[^3]: http://www.lsst.org


# Performance

Brute force algorithms to compute 3-point correlation function scale as $O(N^3)$. 
Tree methods can improve algorithms performance. However, as explained above, 
CPU time can be dramatically reduced if we use a harmonic basis decomposition, 
a sharing neighbor list and a tighter oc-ball-tree structure. 
Resulting that computing this higher order statistics will scale with the 
number of objects similar to computing 2-point correlation function. 
Full sky analysis for an Nside = 8192, around 800 millions of pixels and 
for maximal searching radius distance of 200 arcmin takes around of 40 CPU minutes 
(wall-clock) on 128 threads of a single Perlmutter-NERSC node.[^4] 
For this vast amount of pixels the CPU time consumption comes
mainly from scanning the three and updating the histograms by summing up the neighbor lists.
For each multipole component, we need to compute its corresponding Chebyshev polynomials
which translate toa naïve $O(mmaxN log N )$ scaling, where mmax is the maximal multipole
number. In summary, we found ([@Arvizu:2025]), for a set of `HEALPix` Nsides parameters
256, 512, 1024, 2048 and 4096, that the power scaling law goes as N 1.1. The scaling with
maximum number of multipoles shows a linear scaling $t(m) = 1.85m + 63.7$ seconds for Nside
$= 1024$. We notice that the slope of $t(m)$ is close to 2, which is the number of operations
needed in the recursive relations of the Chebyshev polynomials. However, the dependence on
$m$ is moderate, being $\sim 20$ % slower to compute up to $m = 8$ than up to $m = 1$.




[^4]: https://docs.nersc.gov/systems/perlmutter/architecture



# Acknowledgements

Author acknowledges support by CONAHCyT CBF2023-2024-162  and  CBF2023-2024-589. 
Also the author acknowledges for the computational resources of the 
DiRAC@Durham facility managed by the Institute for Computational Cosmology on behalf 
of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The equipment was funded by BEIS 
capital funding via STFC capital grants ST/K00042X/1, ST/P002293/1, ST/R002371/1 and ST/S002502/1, 
Durham University and STFC operations grant ST/R000832/1. DiRAC is part of the 
National e-Infrastructure in the U.K.


# References
