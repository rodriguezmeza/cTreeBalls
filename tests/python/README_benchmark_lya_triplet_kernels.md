# Exact and anisotropic multipole Ly-alpha 3PCF kernels

`benchmark_lya_triplet_kernels.py` compares the same selected pixels, weights,
forest exclusions, radial bins and polar bins. It retains every numerical
product, command, catalog, binary hash, process timing and per-process peak RSS.
It uses the same FITS reader as `lya_corr_all_engines.py`.

From the source checkout, after building `cballs`:

```bash
python tests/python/benchmark_lya_triplet_kernels.py \
  --fits /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --max-forests 1000 --pixel-stride 10 --threads 16 \
  --r3-max 160 --r3-bins 4 --theta-bins 4 --mu-bins 4 \
  --warmups 1 --repeats 3 --lmax 4 8 16 32 \
  --outdir results/lya-triplet-calibration
```

The output directory must be new. `--catalog saved.npz`, `--ascii pixels.txt`
and `--synthetic` are also supported. Omit values after `--lmax` to test only
exact kernels. A full reference run can still take minutes on a dense catalog.
All runs are sequential fresh processes. Native process wall/CPU timing includes
reading, tree construction, search, and output; Python FITS preprocessing is
outside that timing. Native search CPU is recorded separately. Peak RSS covers
the child process, not the Python FITS reader, and is unavailable without wait4.
Compare identical builds, hosts, thread affinity and workload for timing claims.

## Exact paths

The existing `lya-3pcf-omp`, `lya-2pcf-3pcf-omp`, corresponding LOS-tree methods,
and original 3D MPI methods retain their estimator and output format.
The kernel selection parameter is:

| `lya3Kernel` | Implementation |
|---|---|
| 0 (default) | Exact per-forest segment aggregation, tiled pixel fallback |
| 1 | Retained individual-neighbor-pair reference loop |
| 2 | Tiled direct loop, without segment aggregation |

OpenMP 3PCF work publishes fixed blocks in catalog order: automatic selection
uses eight pivots for runs with <=16384 pivots, otherwise 64. Full-catalog timing
showed that using eight everywhere increased overhead, despite helping the
small dense-pivot pilot. `lya3PivotBlock=1..4096` overrides the selection; the
benchmark exposes `--pivot-block`. Selection never depends on thread count.
Pair-only work retains 64; original MPI rank ownership retains one. Scheduling
and segment sums change floating summation grouping, so agreement is to
roundoff, not bit-for-bit with the previous binary. A fixed kernel's output is
deterministic across OpenMP thread counts. MPI reduction can change rounding.

Each pivot's neighbors are grouped by **forest ID, radial bin and LOS polar
bin**, sorted by measured radius, and split into a binary hierarchy. Summaries
store component bounds of the actual pixel directions and sums of `w` and
`w*delta`. For two nodes from distinct forests, a conservative interval for
`u_q dot u_r` must lie wholly within one mu bin before their products are
aggregated. Otherwise nodes are split, then individual pixels are evaluated.
The direct leaf fallback batches exact geometry and, for up to 64 mu bins,
combines deposits in a small local histogram before updating the worker grid.
Larger mu grids use the general tiled loop.
Noncollinear pixels sharing a forest ID remain supported. Both leg orders and
zero-weight geometric triplet counts are preserved. No third-side cutoff is
introduced. No line-of-sight approximation is made.

Extreme weight/delta magnitudes use leaf arithmetic. Ambiguous boundary bounds
also descend. Bounds assume normal IEEE arithmetic; aggregation is disabled
under fast-math or a non-double arithmetic profile. Default builds use
`-fno-fast-math`. New scratch sizes are checked and included in the resource
preflight; the memory budget remains a planned-allocation limit, not a process
RSS cap. `theta` does not tune any of these forest estimators.

## Explicit anisotropic multipoles

`searchMethod=lya-anisotropic-multipole-3pcf-omp` (ID 209) uses LOS-tree discovery
and `lya3LMax=0..32` (default 8). This is an experimental, opt-in method.
It keeps **both radial bins and both polar angles relative to the pivot LOS**.
It is different from the isotropic `octree-3pcf-3d-omp` estimator.

For a pivot p and radial/polar cell a, define real normalized harmonics H with
`sum_m H_lm(u) H_lm(v) = P_l(u dot v)`. A forest's numerator moments are
`A_lm(f,a) = sum_(q in f,a) w_q delta_q H_lm(u_pq)`; denominator moments replace
`w_q delta_q` by `w_q`. The pivot forest is excluded during discovery.

The raw numerator is

```
T_l(a,b) = w_p delta_p sum_(f != g) sum_m A_lm(f,a) A_lm(g,b).
```

The denominator uses weight-only moments and pivot factor `w_p`. The implementation
uses an occupancy-selected forest-block hierarchy or optimized prefix,
and deposits both leg orders. `--multipole-kernel 0` is automatic (default),
`1` retains the original prefix reference, and `2` forces four-forest groups.
See [the reuse guide](../../docs/LYA_MULTIPOLE_HIERARCHICAL_REUSE.md). This is algebraically the total moment product minus same-forest
products, without subtracting two nearly equal, large auto products. It removes
both coincident neighbors and distinct pixels in the same forest. No pairwise
mu binning is performed in this kernel. Selected multipole moments themselves
are exact up to floating arithmetic; truncation only enters reconstruction.

For a mu bin [a,b], reconstruct numerator and denominator **separately** using

```
c_l = (2*l+1)/2 * integral_a^b P_l(mu) dmu
N_bin ~= sum_(l=0..L) c_l T_l^N
D_bin ~= sum_(l=0..L) c_l T_l^D
zeta_bin ~= N_bin / D_bin.
```

Two distinct files avoid confusing approximate histograms with exact output:

* `histZetaM_lya_multipoles.txt`: `b1 b2 t1 t2 ell numerator denominator`.
  These are signed raw Legendre moments, not normalized correlation ratios.
* `histZetaM_lya5d_multipole.txt`: the usual 13-column five-dimensional layout,
  explicitly marked approximate. All bins are retained. Nonpositive reconstructed
  denominators produce NaN correlation; signed raw values are never clipped.

Finite-order top-hat reconstruction can ring, leak into empty bins, and have
large relative errors near correlation zero crossings. Increasing L does not
guarantee monotonically improving error in every individual bin. `L=32` is
not a universal accuracy guarantee. At an interior mu-bin edge, the infinite
Legendre series itself converges to the midpoint of a top-hat jump; atoms exactly
on that edge need direct treatment to reproduce half-open bin conventions. Use the raw moments directly if they are
the intended observable, or calibrate reconstructed bins on representative data.
Costs scale roughly with neighbor count times (L+1)^2, plus products between
occupied radial/polar cells across forests. High L, many populated cells, or
sparse neighborhoods can make this slower than the exact segment method.

## Acceptance and regression evidence

The benchmark checks all exact products against kernel 1, including raw sums,
with `rtol=1e-9, atol=3e-11` for large accumulated catalogs. Small independent
oracle tests use tighter tolerances. Every repeated output must be identical.

For approximate reconstruction, `summary.json` reports raw-sum relative L2,
maximum and 95th-percentile correlation relative errors, invalid denominators,
empty-bin leakage and excluded low-signal bins. By default the acceptance test
requires <=5% relative error in **every occupied reference bin with
abs(zeta)>1e-6**, with no invalid eligible bins. The floor is explicit and
configurable; excluded bins are counted and must not be represented as passing.
`--require-accepted-multipole` makes lack of an accepted L a failing exit.
The script never silently substitutes approximate output in exact comparisons.

Run native regressions with:

```bash
CBALLS="$PWD/cballs" python -m pytest -q \
  tests/python/test_lya_triplet_acceleration.py \
  tests/python/test_lya_forest_los_tree.py
python tests/python/test_lya_forest_omp.py
python scripts/generate_capabilities.py --check
```

The capability manifest declares the new engine, shared-file ownership and
oracle tests. The active release gate retains independent multipole oracle cases
and compares one versus multiple threads. The standard all-engines script keeps
its exact-estimator comparisons separate; use this calibration script for the
new anisotropic multipole engine. Cython accepts its name and parameters through
`set()` and `set_forest_catalog()`; `Run()` writes the files above and
`getRunMetadata()` records the kernel, L and approximate reconstruction flag.
Rebuild both cballs and cyballs after updating native command fields.

## Explicit cell-geometry calibration

See [the cell calibration guide](README_benchmark_lya_cell_approximation.md) for mu-only slop, persistent forest nodes, pivot-cell aggregation, and retained accuracy/timing/memory evidence. All slops default to zero.

## Hierarchical radial moments

`lya2Kernel=1` and `lya3Kernel=5` select certified pair-range sums and adaptive
radial/polar moment combinations for the original 3D OpenMP/MPI methods.
Combined runs share the forest tree; zero slop preserves the histogram estimator.
Persistent kernels 3/4/5 also support MPI. The all-engines driver accepts
`--lya2-kernel 1 --lya3-kernel 5 --lya3-pivot-cell-max 8`; these controls apply
only to eligible 3D Ly-alpha statistics. Default kernels remain 0.
See the [implementation and benchmark guide](../../docs/LYA_HIERARCHICAL_REUSE.md).


## Exact mu output alongside raw moments

`--mu-mode exact` sets `lya3MuMode=1` on the anisotropic multipole runs. They
reuse discovery/sorting and additionally calculate exact hard angular bins.
The driver compares `_lya5d.txt` against the direct reference, checks raw sums,
and records the selected `mu_product` and `filename` with every comparison.
The original finite-L `_lya5d_multipole.txt` is also evaluated separately under
`reconstruction_comparisons`; its accuracy can still fail.

```sh
python tests/python/benchmark_lya_triplet_kernels.py \
  --synthetic --mu-mode exact --lmax 4 8 16 32 --threads 2 \
  --warmups 1 --repeats 3 --require-accepted-multipole \
  --outdir results/lya-exact-mu-calibration
```

Default `--mu-mode finite-multipole` preserves the previous approximate-only
test. `--require-accepted-multipole` refers to the **selected** mu product of a
multipole run, and exact-mode timing includes both multipoles and exact bins.
It must not be quoted as moment-only performance or a successful finite-L
reconstruction. See [the detailed guide](../../docs/LYA_MULTIPOLE_RECONSTRUCTION.md)
for the finite-moment limitation, output metadata, memory and API contracts.
