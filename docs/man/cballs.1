.TH CBALLS 1 "September 2026" "cTreeBalls 1.1.0" "User Commands"
.SH NAME
cballs \- compute two- and three-point correlation functions
.SH SYNOPSIS
.B cballs
.RI [ parameter-file ]
.RI [ name=value ... ]
.SH DESCRIPTION
.B cballs
computes correlation functions for scalar fields, point counts, full-sky
spin-2 shear, Lyman-alpha forests, and physical three-dimensional catalogs.
Available search methods are selected at compile time.
.PP
Command-line assignments must not contain spaces around the equals sign.
Parameter files use one name/value assignment per line and may contain comments.
.SH DISCOVERY
.TP
.B options=make-info
Print the Makefile settings, resolved toolchain, precision, library, OpenMP,
and MPI profile embedded in the executable.
.TP
.B options=print-options
Print registered runtime option names, scope, and a short description.
.TP
.B options=print-search-methods
Print only search methods registered in this executable, with geometry,
statistics, build switch, and usage notes.
.SH ACTIVE SEARCH METHODS
The maintained profile enables these families:
.TP
.B octree-sincos-omp
Core OpenMP octree method for scalar 2PCF and sine/cosine 3PCF multipoles.
.TP
.B kdtree-2balls-omp, kdtree-2balls-mpi
Median KD-tree dual-node 2PCF and LogMultipole 3PCF.
.TP
.B balltree-2balls-omp, balltree-2balls-mpi
PCA ball-tree dual-node 2PCF and LogMultipole 3PCF.
.TP
.B octree-2balls-omp, octree-2balls-mpi
Native-octree dual-node 2PCF and LogMultipole 3PCF.
.TP
.B octree-shear-sphere-2balls-omp
Full-sky spin-2 correlation functions over the native octree.
.TP
.B kdtree-shear-sphere-2balls-omp
Full-sky spin-2 correlation functions over a median KD tree.
.TP
.B balltree-shear-sphere-2balls-omp
Full-sky spin-2 correlation functions over a PCA ball tree.
.TP
.B lya-*-omp, lya-*-mpi
Anisotropic, radial, and radial interval-tree forest estimators. Run
options=print-search-methods for the complete names.
.TP
.B octree-3pcf-3d-omp, octree-3pcf-3d-mpi
Physical-3D Legendre multipoles and data/random survey-window estimation.
.TP
.B kdtree-box-omp, neighbor-boxes-omp
Periodic Cartesian 2PCF methods.
.SH COMMON PARAMETERS
.TP
.BI searchMethod= name
Select a method printed by options=print-search-methods. Alias: search.
.TP
.BI infile= path
Input catalog. Alias: in.
.TP
.BI infileformat= format
Catalog format. Alias: infmt. Common values include columns-ascii, binary,
fits-healpix, gadget, and lya-ascii when their input add-ons are enabled.
.TP
.BI iCatalogs= list
One-based catalog selectors used by auto- and cross-correlation methods.
.TP
.BI rootDir= path
Output directory. Alias: root.
.TP
.BI rminHist= value
Minimum measured separation.
.TP
.BI rangeN= value
Maximum measured separation.
.TP
.BI sizeHistN= bins
Number of radial bins.
.TP
.BI mChebyshev= order
Largest positive 3PCF multipole order.
.TP
.BI useLogHist= bool
Choose logarithmic or linear radial bins.
.TP
.BI numberThreads= count
OpenMP threads per process or per MPI rank.
.TP
.BI options= list
Comma-separated behavior controls.
.SH SEARCH CONTROLS
.TP
.B only-2pcf, only-3pcf
Run only the requested compiled correlation order.
.TP
.B no-two-balls
Disable dual-node cell aggregation and use the exact body-pair limit.
.TP
.B dual-node-bin-slop
Enable bin-position-aware Log/Linear node acceptance.
.TP
.B no-smooth-pivot
Disable the build-default pivot smoothing on methods that support it.
.TP
.B read-mask
Apply a supported companion or in-memory binary mask.
.TP
.B edge-corrections,no-normalize-HistZeta
Compute complex scalar or shear 3PCF window correction. This requires 3PCF.
.TP
.B weights-norm
Use catalog weights in signal and normalization moments.
.SH PYTHON
The compiled extension is imported as
.BR cyballs .
The convergence, shear, and forest drivers are under
.IR tests/python .
They retain one NumPy catalog across selected engines and write timing and
comparison summaries.
.SH MPI
Use one MPI implementation for MPICC, the runtime launcher, the extension, and
mpi4py. numberThreads is per rank. All ranks enter the same run and cleanup
sequence; rank zero writes results.
.SH EXAMPLES
.PP
Inspect the executable:
.PP
.nf
  cballs options=make-info
  cballs options=print-options
  cballs options=print-search-methods
.fi
.PP
Run a scalar angular 2PCF:
.PP
.nf
  cballs search=octree-2balls-omp infile=catalog.fits \
    infmt=fits-healpix options=only-2pcf rootDir=Output
.fi
.PP
Run masked edge-corrected 3PCF:
.PP
.nf
  cballs search=kdtree-2balls-omp infile=map.fits,mask.fits \
    infmt=fits-healpix,fits-healpix \
    options=read-mask,only-3pcf,edge-corrections,no-normalize-HistZeta \
    rootDir=Output_edge
.fi
.PP
Run an MPI forest estimator:
.PP
.nf
  mpiexec -n 2 cballs addons/lya_forest_mpi/parameters.ini
.fi
.SH FILES
.TP
.I addons/Makefile_addons_settings
Compile-time addon switches.
.TP
.I tests/python
Catalog drivers, plots, and README examples.
.TP
.I docs
Sphinx documentation.
.SH SEE ALSO
The project README and Sphinx guide document input schemas, normalization,
scientific conventions, and validation tests.
.SH COPYRIGHT
Copyright 2023-2026 Mario A. Rodriguez-Meza. Distributed under the MIT license.
