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
.TP
.B options=build-fingerprint
Print the source/profile/toolchain JSON identity. The Python extension exposes
the corresponding identity through cyballs.build_info().
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
.B lya-los-tree-2pcf-omp, lya-los-tree-3pcf-omp, lya-los-tree-2pcf-3pcf-omp
Exact anisotropic 3D forest estimators using octree forest discovery and
per-LOS radial trees. Transverse distance is retained. These methods have no
MPI counterpart and are enabled by LYAFORESTOMPON=1.
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
This alone does not make the native octree 3PCF pivot scan exact.
.TP
.B no-one-ball,no-two-balls,no-smooth-pivot
Exact unsmoothed body-level reference for scalar and spherical shear two-ball
methods. BALLS4SCANLEVON may remain enabled: it controls scheduling.
.TP
.B legacy-one-ball
Select the privately linked compatibility kernel behind an active two-ball
method. Its smoothing behavior differs from native octree dual node traversal.
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
.TP
.B shear-pivot-reuse
Opt-in native octree or ball-tree OpenMP 3PCF aggregate-pivot reuse with
inherited unresolved neighbors and bounded spherical transport. Requires
BALLS4SCANLEVON=1, no-smooth-pivot, positive theta and full pivot coverage.
KD-tree and compatibility mode reject this option. The independent 2PCF path
is unchanged. Validate against exact output before production use.
.TP
.B no-balltree-shear-member-cache
Disable the transient source-frame construction cache (default cap 256 MiB).
This diagnostic preserves direct member moments and conservative bounds.
.TP
.BI stepState= count
With verbosity enabled, report completed body-pivot progress in the scalar
octree, KD-tree and ball-tree OpenMP two-ball 3PCF scans.
.SH ENVIRONMENT
.TP
.B CBALLS_SHEAR_PIVOT_TOL
Finite phase budget in radians from 0 to 3, default 0.1. Not a relative-error
guarantee. With reuse enabled, this budget controls 3PCF acceptance; theta's
magnitude does not. Zero disables reuse, not ordinary neighbor approximation.
.TP
.B CBALLS_SHEAR_PROFILE
Set to 1 for per-thread phase timers and SHEAR_REUSE counters.
.SH PYTHON
The compiled extension is imported as
.BR cyballs .
The convergence, shear, and forest drivers are under
.IR tests/python .
They retain one NumPy catalog across selected engines and write timing and
comparison summaries.
Their timing tables separate setup, complete Python compute calls, and native
MainLoop time. mainloop_wall_s and mainloop_cpu_s exclude Python provenance
capture; compute columns include it. MPI wall time is the maximum
across participating ranks and CPU time is their sum. The forest driver reads
DESI and eBOSS/PICCA delta FITS, compares compatible estimator families, and
supports model distortion/covariance analysis. Consult its README for input
contracts and separately scoped external reference timings.
.SH MPI
Use one MPI implementation for MPICC, the runtime launcher, the extension, and
mpi4py. numberThreads is per rank. All ranks enter the same run and cleanup
sequence; rank zero writes results.
.SH CAPABILITY CATALOGUE
capabilities/engines.json declares the public active search methods, aliases
and regression ownership. Run scripts/generate_capabilities.py after editing
it; builds reject stale generated registration files. print-search-methods
reports only entries enabled in the current executable. Inactive addons
and addons/python_env are not part of the public distribution.
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
