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
    affiliation: 1
  - name: Eladio Moreno
    orcid: 0000-0002-5400-2584
    affiliation: 2
    corresponding: true
  - name: Alejandro Aviles
    orcid: 0000-0001-5998-3986
    affiliation: 3
  - name: Gustavo Niz
    orcid: 0000-0002-1544-8946
    affiliation: 2
affiliations:
  - name: Departamento de Física, Instituto Nacional de Investigaciones Nucleares, Apartado Postal 18-1027, Col. Escandón, Ciudad de México,11801, México
    index: 1
  - name: Departamento de Ciencias e Ingenierías, Universidad de Guanajuato, 37150, León, Guanajuato, México
    index: 2
  - name: Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México, 62210, Cuernavaca, Morelos, México
    index: 3
date: 9 April 2026
bibliography: paper.bib
---

# Summary

`cTreeBalls`[^1] (`cBalls` for short) is a Python/C package useful to measure (2,3)-point clustering statistics. `cBalls` can efficiently calculate 3 point correlations of more than  200 million of healpix pixels (a full sky simulation with Nside=4096) in  less than 10 minutes on a single high-performance computing node, enabling a feasible analysis for the upcoming LSST data. It builds upon octree [@Barnes:1986] and kd-tree algorithms [@Bentley:1975] and supplies a user-friendly interface with flexible input/output (I/O) of catalogue data and measurement results, with the built program configurable through external parameter files and tracked through enhanced logging and warning/exception handling. For completeness and complementarity, methods for measuring two-point clustering statistics  for periodic boxes are also included in the package. `cTreeBalls` was developed for its use in the Dark Energy Science Collaboration (DESC) of the Rubin Observatory Legacy Survey of Space and Time
(LSST).[^2]

[^1]: http://github.com/rodriguezmeza/cTreeBalls.git/

[^2]: http://www.lsst.org/

# Statement of need



Correlation functions are frequently used to extract relevant information about the large scale structure of the universe, which in turn can be used to discriminate between cosmological models. A common practice is to employ two point statistics, such as the two point correlation function (2PCF) or its Fourier space counterpart, the power spectrum. A random gaussian field is fully characterized by the two point statistics, however, the non-linear cosmic evolution of matter produces a measurable deviation from this random gaussian field assumption. Equivalently, primordial non-gaussianities may also imprint a non-gaussian signal. Therefore, it is crucial to use alternative methods to extract the non-gaussian information that is not captured by the two point statistics. Higher order point correlations are frequently used as statistical tools to extract such information. The three point correlation function (3PCF), or the bispectrum in Fourier space, is the next order term in the hierarchy of n-point correlation functions based on signal to noise ratio on sufficiently large scales. Just as the 2PCF describes characteristic length-scales of data, the 3PCF encodes all the spatial triangles in that distribution. As a result, an isotropic 3PCF depends on three scales that determine the side lengths of these triangles.

In this work, we provide an infrastructure to calculate 2PCF and 3PCF for scalar quantities on the sphere, such as the weak lensing convergence field or the CMB lensing field. A naïve estimator of the 3PCF requires $O(N^3)$ evaluations, where $N$ is the number of discrete data points or pixels on a continuum field. For present galaxy surveys, such as the one that will be released by the LSST [@Vera:2019], this scaling results in a prohibited amount of time. Using tree structures in the search for neighboring data points can reduce the scaling to $O(N^2\log N)$, which is still not viable. However, one can map the information of the 3PCF onto an orthogonal basis whose coefficients scale as pair counts, $O(N\log N)$, hoping that only a small and finite number of eigenvectors on that basis capture most of the information encoded in the 3PCF. This is the case of a harmonic decomposition of one of the angles in the triangle configurations of the 3PCF for the weak lensing convergence [@Arvizu:2025]. This expansion was first proposed by Szapudi [@Szapudi:2004] and has since been applied to different algorithms and codes (e.g. [@Zheng:2004], [@Slepian:2018], [@Philcox:2022]).


Following above ideas
we have developed a code that can efficiently calculate 3 point correlations of more than 200 million pixels (a full sky simulation with Nside=4096) in less than 10 minutes on a single high-performance computing node. Therefore we will be able to analyse the upcoming LSST data that is planned to be of more than 20 billions of objects.

Specifically, `cBalls` can calculate:

- three-point clustering statistics, namely multipoles of the three-point correlation function (3PCF) in configuration space for projected fields;
- two-point clustering statistics, namely the two-point correlation function (2PCF), for projected fields or cubic-box simulation mocks.


There are existing software packages that use different decomposition of the three-point clustering statistics or that
use different tree schemes, for instance:

- `TreeCorr` [@Jarvis:2004] computes (2,3)-point correlation function for counts, convergence and shear weak lensing fields. It uses a ball-tree scheme.

 
 
- `Triumvirate` [@Wang:2023] computes three point statistics in Fourier and configuration spaces using a tri-polar spherical harmonic decomposition [@Sugiyama:2019].


In particular, the last one  uses a different decomposition of three-point clustering statistics and may have different constraining power on different cosmological parameters. Therefore, `cBalls` fulfills complementary needs in current galaxy clustering analyses.

# Implementation

A binned estimator to study the 3PCF signal of the weak lensing convergence field was introduced in [@Arvizu:2025]. Extract this kind of signal is prohibitively CPU time consuming when we are dealing with a vast amount of sources. Brute force algorithms scales with the number of galaxies as $O(N^3)$. We can speed up if we introduce a tree structure (kd-tree or octree) for the spatial part of the catalog map. We are using mainly octree structures [@Barnes:1986], but kd-tree is also used [@Bentley:1975]. Now, the real improvement comes when we use a harmonic decomposition of the  weak lensing convergence field and impose homogeneity and isotropy, and using the fact the field is a spin-0 field (a real number). In such a case complexity of the code goes as $O(N \log N)$ as it does a 2PCF code. Multipoles are the estimators we are looking for that characterize the three-point statistics.

Then, computation loop goes as: pick up a galaxy from the survey catalog, this is called the pivot, a vertex of a triangle, find all its neighbours by scanning the tree in the spatial range of search and fill the multipoles bins. Here, we need to evaluate $\sin$ and $\cos$ functions, we instead turn this into evaluation of Chebyshev's polynomials which is much faster than evaluating trigonometric functions.
 

Another speed improvement comes from the following idea:
spatial closer pivots share neighbor list. We just need to define and quantify "closeness". At least there are four length scales. First, the length scale of the root in the octree hierarchy, the cube that contain all the galaxies to be analyzed. Second, the scale length of the smallest cell in this hierarchy. Third, the radial length of each bin in the histograms that contain all correlation information. And fourth, the mean separation between galaxies. `cBalls` can compute the pivot radius of the neighborhood automatically using the latter scale length. But user can give one if needed.

Furthermore, octree cells can be tightened (by pruning) and create a kind of *oc-ball-tree*. Therefore, when scanning the tree, we only visit regions which are filled up with galaxies, instead of visiting cubic cells which may have not negligible void space. 

# Features

`cBalls` has many features that makes it appropriate to be integrated
in massive workflow like `TXPipe` [@Prat:2023] used in Vera C. Rubin Observatory with its Legacy Survey of Space and Time (LSST) collaboration[^3]:

- The code has a Python interface (`cyballs`) that was implemented using Cython to bind the C searching engine. Python scripts are given in such a way users can run the code in a linux terminal or can be run using a Jupyter notebook. The C code can also be compiled and executed independently of Python.
- Parameters files can be given as standard `ASCII` files or in command line. If using a Jupyter Notebook or a Python script, parameters are given very easily as a very simple dictionary without the need of a parameter file. Before computing correlations parameter values are checked to secure the computation.
- Catalog map files can be given in several convenient formats, such as `CFITSIO`, `HEALPix`, `GADGET-2` [@Springel:2005] and  `ASCII` columns format, like `ROCKSTAR` [@Behroozi:2013] halos catalogs. Also, Takahashi *et al*. [@Takahashi:2017]  weak lensing realizations can be read by the code without being transformed to `HEALPix` format.  Several map files can be read in order to compute cross-correlations function or to apply a mask file that is useful for current surveys that are unable to scan the full sky.
- Searching engines are parallelised with OpenMP threading, i.e., pivot pixels or galaxies are distributed amongst multiple CPU threads.
- Clustering statistics can be done on a cubic box with periodic boundary condition (only 2pcf) or on a unit sphere under the Limber approximation for projected scalar fields.
- Second and third order statistics can be done in a single run or they can be done independently.  
- Edge corrections are implemented [@Slepian:2015] in the case of computing 3-pcf multipoles over volumes that are not the full sky. 
- A four levels logger is provided for runtime tracking. Running information can be printed on the screen or saved in a log file. CPU time and memory usage can be reported at each main stages of the computation.

[^3]: http://www.lsst.org/

# Performance

Brute force algorithms to compute 3-point correlation function scale as $O(N^3)$. Tree methods can improve algorithms performance. However, as explained above, CPU time can be dramatically reduced if we use a harmonic basis decomposition, a sharing neighbor list and a tighter oc-ball-tree structure.  Resulting that computing this higher order statistics will scale with the number of objects similar to computing 2-point correlation function. Full sky analysis for an Nside $=8192$, around 800 millions of pixels and for maximal searching radius distance of 200 arcmin takes around of 40 CPU minutes (wall-clock) on 128 threads of a single Perlmutter-NERSC node[^4].
For this vast amount of pixels the CPU time consumption comes mainly from scanning the three and updating the histograms by summing up the neighbor lists. For each multipole component, we need to compute its corresponding Chebyshev polynomials which translate toa naïve $O(m_{max}N\log N)$ scaling, where $m_{max}$ is the maximal multipole number. In summary, we found [@Arvizu:2025], for a set of `HEALPix` Nsides parameters 256, 512, 1024, 2048 and 4096, that the power scaling law goes as $N^{1.1}$. The scaling with maximum number of multipoles shows a linear scaling $t(m)=1.85 +63.7$ seconds for Nside $=1024$. We notice that the slope of $t(m)$ is close to 2, which is the number of operations needed in the recursive relations of the Chebyshev polynomials. However, the dependence on $m$ is moderate, being $\sim 20$ \% slower to compute up to $m=8$ than up to $m=1$.

[^4]: https:/docs.nersc.gov/systems/perlmutter/architecture


# Acknowledgements

The authors acknowledge support by SECIHTI (previously CONAHCyT) grants CBF2023-2024-162,  CBF2023-2024-589 and BF-2025-I-2795. AA also acknowledges DGAPA-PAPIIT IA101825. EM and GN acknowledge the support of the DAIP-UG grant CIIC-254/2026 and the computational resources of the DCI-UG DataLab. MARM also acknowledges the computational resources of the DiRAC@Durham facility managed by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). 

For the purpose of open access, the authors have applied a Creative Commons Attribution (CC BY) license to any Author Accepted Manuscript version arising from this submission. 

# References
