.TH CBALLS 1 "October 3, 2026" "cTreeBalls 1.1.0" "User Commands"
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
Parameter files use one name/value assignment per line and may contain
comments.
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
The maintained Makefile profile enables exactly these 41 methods.
.TP
.B octree-sincos-omp
2D/3D scalar; octree; OpenMP. standard 2PCF and sine/cosine 3PCF multipoles.
.TP
.B kdtree-2balls-omp
scalar; median k-d tree; dual-node OpenMP scan. 2PCF and LogMultipole angular
3PCF, masks, and complex edge correction.
.TP
.B kdtree-2balls-mpi
scalar; median k-d tree; deterministic MPI+OpenMP dual-node scan. distributed
2PCF and LogMultipole angular 3PCF, masks, and complex edge correction.
.TP
.B balltree-2balls-omp
scalar; FCFC PCA ball tree; dual-node dual/triple-node OpenMP scan. 2PCF and
angular-multipole 3PCF with auto- and cross-catalog support.
.TP
.B balltree-2balls-mpi
scalar; FCFC PCA ball tree; deterministic MPI+OpenMP dual/triple-node scan.
distributed 2PCF and angular-multipole 3PCF with auto- and cross-catalog
support.
.TP
.B octree-2balls-omp
scalar; native octree; dual-node 2PCF and LogMultipole 3PCF. 2PCF and
angular-multipole 3PCF with auto- and cross-catalog support.
.TP
.B octree-2balls-mpi
scalar; native octree; deterministic MPI+OpenMP frontier. distributed 2PCF and
LogMultipole angular-multipole 3PCF.
.TP
.B octree-shear-sphere-2balls-omp
full-sky spin-2 shear; dual-node dual-node octree; OpenMP. xi+/xi- and four
Gamma-x shear 3PCF components.
.TP
.B octree-shear-sphere-2balls-mpi
full-sky spin-2 shear; native octree; adaptive MPI+OpenMP frontier.
distributed xi+/xi- and four Gamma-x shear 3PCF components.
.TP
.B kdtree-shear-sphere-2balls-omp
full-sky spin-2 shear; dual-node median KD tree; OpenMP. xi+/xi- and four
Gamma-x shear 3PCF components.
.TP
.B kdtree-shear-sphere-2balls-mpi
full-sky spin-2 shear; median KD tree; adaptive MPI+OpenMP frontier.
distributed xi+/xi- and four Gamma-x shear 3PCF components.
.TP
.B balltree-shear-sphere-2balls-omp
full-sky spin-2 shear; dual-node FCFC PCA ball tree; OpenMP. xi+/xi- and four
Gamma-x shear 3PCF components.
.TP
.B balltree-shear-sphere-2balls-mpi
full-sky spin-2 shear; FCFC PCA ball tree; adaptive MPI+OpenMP frontier.
distributed xi+/xi- and four Gamma-x shear 3PCF components.
.TP
.B kdtree-box-omp
periodic Cartesian box; k-d tree; OpenMP. box 2PCF.
.TP
.B neighbor-boxes-omp
periodic Cartesian box; linked boxes; OpenMP. periodic pair counts and
unweighted density correlation function.
.TP
.B octree-3pcf-3d-omp
3D scalar; exact octree leaves; OpenMP. spherical-harmonic 2PCF/3PCF
multipoles.
.TP
.B octree-3pcf-3d-mpi
3D scalar; exact octree leaves; MPI+OpenMP pivot blocks. spherical-harmonic
2PCF/3PCF; data/random survey estimator and edge correction.
.TP
.B lya-2pcf-omp
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
OpenMP. weighted anisotropic 2PCF.
.TP
.B lya-3pcf-omp
3D observer-centered forest pixels; exact default, opt-in bounded cell-bin
approximation; OpenMP. weighted five-dimensional 3PCF.
.TP
.B lya-2pcf-3pcf-omp
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
OpenMP. weighted 2PCF and 3PCF; shared discovery or independently optimized
passes.
.TP
.B lya-los-tree-2pcf-omp
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
OpenMP. exact weighted anisotropic 2PCF.
.TP
.B lya-los-tree-3pcf-omp
3D observer-centered forest pixels; exact default, opt-in bounded cell-bin
approximation; OpenMP. exact weighted five-dimensional 3PCF.
.TP
.B lya-los-tree-2pcf-3pcf-omp
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
OpenMP. weighted 2PCF and 3PCF; shared discovery or independently optimized
passes.
.TP
.B lya-1d-2pcf-omp
radial Lyman-alpha pixels; sorted 1D range scan; OpenMP. weighted radial-only
2PCF.
.TP
.B lya-1d-3pcf-omp
radial Lyman-alpha pixels; sorted 1D range scan; OpenMP. weighted radial-only
3PCF.
.TP
.B lya-1d-2pcf-3pcf-omp
radial Lyman-alpha pixels; sorted 1D range scan; OpenMP. weighted radial-only
2PCF and 3PCF.
.TP
.B lya-1d-tree-2pcf-omp
radial Lyman-alpha pixels; exact 1D interval tree; OpenMP. weighted
radial-only 2PCF.
.TP
.B lya-1d-tree-3pcf-omp
radial Lyman-alpha pixels; exact 1D interval tree; OpenMP. weighted
radial-only 3PCF.
.TP
.B lya-1d-tree-same-los-2pcf-omp
per-LOS radial Lyman-alpha pixels; exact 1D interval trees; OpenMP. equal-LOS
mean of weighted within-forest radial 2PCFs.
.TP
.B lya-2pcf-mpi
Lyman-alpha pixels; 3D octree; MPI+OpenMP. weighted anisotropic 2PCF.
.TP
.B lya-3pcf-mpi
Lyman-alpha pixels; 3D octree; MPI+OpenMP. weighted five-dimensional 3PCF.
.TP
.B lya-2pcf-3pcf-mpi
Lyman-alpha pixels; 3D octree; MPI+OpenMP. weighted 2PCF and 3PCF.
.TP
.B lya-1d-2pcf-mpi
Lyman-alpha pixels; radial range scan; MPI+OpenMP. weighted radial-only 2PCF.
.TP
.B lya-1d-3pcf-mpi
Lyman-alpha pixels; radial range scan; MPI+OpenMP. weighted radial-only 3PCF.
.TP
.B lya-1d-2pcf-3pcf-mpi
Lyman-alpha pixels; radial range scan; MPI+OpenMP. weighted radial-only 2PCF
and 3PCF.
.TP
.B lya-1d-tree-2pcf-mpi
Lyman-alpha pixels; radial interval tree; MPI+OpenMP. weighted radial-only
2PCF with same-forest subtraction.
.TP
.B lya-1d-tree-3pcf-mpi
Lyman-alpha pixels; radial interval tree; MPI+OpenMP. weighted radial-only
3PCF with exact three-forest exclusion.
.TP
.B lya-anisotropic-multipole-3pcf-omp
3D observer-centered forest pixels; exact radial and LOS polar bins; LOS-tree
discovery; OpenMP. anisotropic Legendre raw moments; approximate mu
reconstruction; optional exact five-dimensional bins.
.TP
.B lya-los-tree-2pcf-mpi
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
replicated catalog MPI + OpenMP. exact weighted anisotropic 2PCF.
.TP
.B lya-los-tree-3pcf-mpi
3D observer-centered forest pixels; exact default, opt-in bounded cell-bin
approximation; replicated catalog MPI + OpenMP. exact weighted
five-dimensional 3PCF.
.TP
.B lya-los-tree-2pcf-3pcf-mpi
3D observer-centered forest pixels; exact default, opt-in forest cell pairs;
replicated catalog MPI + OpenMP. weighted 2PCF and 3PCF; shared discovery or
independently optimized passes.
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
.B dual-node-bin-theta
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
Opt-in native octree, KD-tree or ball-tree OpenMP 3PCF aggregate-pivot reuse
with
inherited unresolved neighbors and bounded spherical transport. Requires
BALLS4SCANLEVON=1, no-smooth-pivot, positive theta and full pivot coverage.
MPI and compatibility mode reject this option. The independent 2PCF path
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
.B CBALLS_SHEAR_BIN_THETA
Internal radial-bin assignment slop in [0,1] bin widths, default 0.
Used only with shear-pivot-reuse. Radial range cutoffs remain strict.
Resolved radial pairs complete at parent cells; descendants inherit completion
masks. Use nsmooth=1 for a finer binary pivot hierarchy. Qualify raw and
window-corrected multipoles on the actual catalog.
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
Copyright 2023-2026 Mario A. Rodriguez-Meza. Distributed under the MIT
license.

.SH SCALAR HIERARCHICAL REUSE
The scalar octree-2balls, kdtree-2balls and balltree-2balls OpenMP and MPI
engines accept options=scalar-pivot-reuse,no-smooth-pivot for 3PCF.
CBALLS_SCALAR_PIVOT_TOL is a runtime phase budget in [0,3] radians (default
0.1).
CBALLS_SCALAR_BIN_THETA is a runtime internal-bin allowance in [0,1] bin
widths
(default 0). Both must be finite and agree on all MPI ranks. Outer radial
cuts remain strict. Positive theta is required. Zero phase budget, exact
controls, smoothing, periodic runs and only-2pcf retain the original path.
The phase budget is not a relative coefficient-error guarantee.
See docs/SCALAR_HIERARCHICAL_REUSE.md for qualification and benchmark
commands.

.SH LYMAN-ALPHA HIERARCHICAL REUSE
lya2Kernel=1 enables certified forest pairs and radial range moments.
lya3Kernel=5 selects adaptive radial/polar moment reuse for the
five-dimensional
histogram, with exact segment fallback for sparse pivot cells. The original
3D OpenMP/MPI names support these controls; kernels 3/4 also support MPI.
Zero slop preserves the estimator. See docs/LYA_HIERARCHICAL_REUSE.md.

.SH SOURCE LAYOUT
Python benchmark, analysis and regression scripts are under tests/python.
The python directory contains only Cython binding sources. Shared builders and
compatibility kernels required by enabled addons are in support. The branch
excludes inactive standalone addons and addons/python_env.
