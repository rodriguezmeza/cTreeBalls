octree-ggg-omp

Enable OCTREEGGGOMPON=1. Use logarithmic bins (useLogHist=true).
TWOPCFON=1 enables pairs and TPCFON=1 enables scalar angular multipoles.
Use only-2pcf or only-3pcf to select runtime work when supported by the build.
no-one-ball disables neighbor-cell aggregation for an exact reference.

For 3D sky catalogs, keep observer-centered positions and use chord-distance
bins. Tangent-plane bearings define the scalar Fourier modes. Acceptance
bounds projected angular error; the repaired C writer publishes the computed
histograms. See docs/3pcf.rst for the precise estimator and component convention.

For raw comparison use no-normalize-HistZeta,weights-norm. Raw moments exclude
repeated neighbors inside the native kernel. Old affected multipoles must be
recomputed. SMOOTHPIVOTON only compiles capability; smooth-pivot is explicit
runtime opt-in. rsmooth alone does not activate smoothing.

read-mask excludes invalid pivots and neighbors. Complex edge correction uses
edge-corrections,no-normalize-HistZeta, with NMultipolesON=1 and
NONORMHISTON=1. Signal modes run through mChebyshev; corrected windows extend
through 2*mChebyshev. Empty/singular systems return zero.
Without correction, ggg-full-window explicitly retains the extended window.
Use a fresh output directory when changing the window policy.

GGG_OMP_PIVOT_CHUNK_SIZE sets the deterministic pivot-chunk size at build time.
ggg-profile reports worker work/wait/merge wall intervals when verbosity is
enabled. These sums are not process CPU time. Validate accuracy before timing.

OCTREEGGGMPION=1 enables the shared-kernel MPI sibling, octree-ggg-mpi.
See addons/octree_ggg_mpi/README.md for scheduling and rank ownership.
See python/README_kappa_corr_all_engines.md for one-read FITS/NPZ workflows.

Checks, from the source root with matching binaries:
    python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py
    make test-octree-ggg-fast-path
