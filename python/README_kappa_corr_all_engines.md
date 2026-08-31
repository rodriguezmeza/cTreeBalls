# In-memory kappa correlations

`kappa_corr_all_engines.py` reads and converts a convergence catalog once,
registers the resulting NumPy arrays with one `cyballs.cballs` object, and runs
any number of scalar angular search methods without rereading the catalog.
This retains prepared arrays and copies them into C body storage per run;
it is not a zero-copy API. See [the Sphinx method guide](../docs/search_methods.rst)
and [scalar numerical contract](../docs/3pcf.rst).

## Coordinates, bins, and raw comparisons

Keep unit-sphere positions centered on the observer. The repaired angular
engines use tangent bearings and Euclidean chord bins, not chord interior
angles. Rotate the whole catalog without translating it.

The current driver's legacy `--theta-scale degree` converts the supplied
limits to radians and forwards those numbers as `rminHist/rangeN`;
`--theta-scale radian` forwards them directly. It does **not** apply the
exact conversion `r=2*sin(alpha/2)`. This is a small-angle convention.
For exact wide-angle limits, compute chord limits and pass them with the
direct-forwarding `radian` setting, or use the in-memory cyballs API with
explicit chord-valued parameters. Actual angle is `2*asin(r/2)`.
The general CPU benchmark script, by contrast, converts degrees to chords.

For comparable raw scalar 3PCF across engines, add:

```sh
--more-options no-normalize-HistZeta,weights-norm
```

This selects weighted ordered distinct triplets. Native KD-tree, legacy
balltree, and BALLS4 kernels now remove repeated neighbors; do not subtract
them again. Ordinary non-edge driver runs do not force a shared raw contract.
Reconstruct `coscos+sinsin + 1j*(sincos-cossin)` before comparing higher
modes. Old affected multipoles must be recomputed with rebuilt C/cyballs.
Compiling `SMOOTHPIVOTON=1` alone does not enable smoothing. Raw KD/legacy
balltree runs reject `smooth-pivot`; BALLS4 never supports it.

List the methods supported by both the current `cballs` executable and the
installed `cyballs` extension:

```sh
python3 python/kappa_corr_all_engines.py --list-engines
```

Run selected OpenMP methods on one FITS map:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits \
  --engine octree-ggg-omp \
  --engine kdtree-omp \
  --engine balltree-2balls-omp \
  --threads 16 \
  --outdir Output_kappa
```

Use `--engine all-omp`, `--engine all-mpi`, or `--engine all` to select all
available engines in a group. MPI engines, including one-rank suites, need
`mpi4py` built against
the same MPI implementation used for cTreeBalls:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits \
  --engine all-mpi \
  --mpi-ranks 4 \
  --threads 8 \
  --outdir Output_kappa_mpi
```

The former `octree-kkk-balls4-omp` engine is now `octree-balls4-omp`.
It forces logarithmic bins and rejects `--more-options smooth-pivot`,
independently of the global `SMOOTHPIVOTON` build setting.

The script may also be launched directly by MPI:

```sh
mpiexec -n 4 python3 python/kappa_corr_all_engines.py \
  --fits map.fits --engine octree-ggg-mpi --threads 8
```

An existing MPI launch is reused without nesting another mpiexec. Python owns
MPI for the entire suite, so cleaning up one engine does not finalize MPI before
the next engine. Input failures on rank 0 are broadcast before catalog transfer,
allowing every rank to exit instead of waiting indefinitely.

## Masked maps

`octree-2balls-omp` and `octree-2balls-mpi` support masked 2PCF and 3PCF.
Use a rebuilt `cyballs` extension containing the native octree mask support:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits --mask mask.fits \
  --engine octree-2balls-omp --threads 16 \
  --outdir Output_kappa_mask
```

For its MPI variant:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits --mask mask.fits \
  --engine octree-2balls-mpi --mpi-ranks 4 --threads 8 \
  --outdir Output_kappa_mask_mpi
```

The script automatically adds `read-mask`. Pixels with mask values greater
than `--mask-threshold` (default `0`) are retained; nonfinite mask values are
excluded. The map and mask are each read once on rank 0, then the in-memory
catalog and mask are shared across the requested engines and broadcast to MPI
ranks. The field mean, unless disabled, is computed using retained pixels.

Select mask-capable engines explicitly: `octree-2balls-omp`,
`octree-2balls-mpi`, `octree-ggg-omp`, `octree-ggg-mpi`,
`octree-balls4-omp`, or `octree-balls4-mpi`. Groups such as
`all-omp` automatically filter out methods without mask support, including
when the mask is supplied inside an NPZ catalog. An explicitly requested
incompatible method is rejected.
`--list-engines` displays mask support next to each method.

An NPZ catalog may instead contain a binary `mask` array. For an existing
notebook catalog, use `KappaCatalog(positions=xyz, kappa=kappa, mask=valid)`
with `valid` boolean or 0/1, and select one of the methods above. Masking excludes
invalid pivots and neighbors from the autocorrelation; it does not itself apply
a survey-window edge correction. To apply angular 3PCF window deconvolution,
select `--type edge_effects` as described below.

## Edge-corrected 3PCF

Enable complex scalar angular 3PCF edge correction with `--edge-corrections`.
The compatible engines are:

- `octree-ggg-omp` and `octree-ggg-mpi`;
- `octree-balls4-omp` and `octree-balls4-mpi` (logarithmic bins);
- `octree-2balls-omp` and `octree-2balls-mpi`;
- `balltree-2balls-omp` and `balltree-2balls-mpi`;
- `balltree-2balls-omp_3pcf` and `balltree-2balls-mpi_3pcf`.

Run every available compatible OpenMP engine on an unmasked map:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits --engine all-omp --edge-corrections \
  --threads 16 --outdir Output_kappa_edge
```

For a masked map, the same group selects only its mask-capable methods:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits --mask mask.fits \
  --engine all-omp --threads 16 \
  --edge-corrections --outdir Output_kappa_mask_edge
```

For MPI, use `--engine all-mpi --mpi-ranks 4 --edge-corrections`, or select a
specific method such as `octree-2balls-mpi`. `--engine all --edge-corrections`
includes both OpenMP and MPI engines; add `--mpi-ranks` for multiple ranks.
Only methods present in both the executable and extension are selected.
Explicit unsupported engine names fail instead of silently disabling correction.

`--type edge_effects` remains supported. Passing
`--more-options edge-corrections` also selects corrected output, even without
`--type edge_effects`. These three spellings are equivalent. The driver adds
`edge-corrections,no-normalize-HistZeta` automatically. It saves real,
imaginary, and complex corrected modes as `zeta_edge_N`, `zeta_edge_im_N`,
and `zeta_edge_complex_N` in `histograms.npz`. `N=1` is the monopole.
The existing `histZetaM_EE_N.txt` files contain the real part; the companion
`histZetaM_EE_Im_N.txt` files retain the imaginary part.
Missing, wrong-shaped, or nonfinite corrected arrays fail the engine run;
there is no fallback to raw multipoles. Each `result.json` and the suite
`summary.json` record whether edge correction was requested.

The build must enable `TPCFON=1`. GGG OpenMP/MPI additionally require
`NMultipolesON=1` and `NONORMHISTON=1`; the two-ball engines do not need these
two extra settings. `--list-engines` reports algorithm support, while a missing
build prerequisite is reported by the native engine. Legacy real-only or
file-only edge routines in `octree-ggg`, `octree-kkk-omp`, and triangle
variants are not exposed as complex corrected products by this driver.

For BALLS4, enable `OCTREEBALLS4OMPON=1` and/or `OCTREEBALLS4MPION=1`.
Its new edge path does not require `NMultipolesON` or `NONORMHISTON`.
Both engines accept masks and reject `smooth-pivot`. Example:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits --mask mask.fits --engine octree-balls4-mpi \
  --mpi-ranks 2 --threads 8 --edge-corrections \
  --more-options weights-norm --outdir Output_balls4_mask_edge
```

Use `--engine octree-balls4-omp` without `--mpi-ranks` for OpenMP alone.
The driver reads the catalog once and shares it with cyballs; MPI broadcasts
the in-memory catalog from rank 0. Edge runs use body pivots within the B4
partition and distinct-neighbor angular moments. Either `no-one-ball` or
`no-two-balls` disables neighbor-cell aggregation in this path.
Uncorrected `no-normalize-HistZeta` runs now use the same body-pivot
distinct-neighbor kernel. Only the legacy normalized non-edge path retains
the original cell-pivot estimator.

The signal uses modes through `--multipoles`; the selection-window moments
extend through twice that order and use the same accepted neighbors and
weights. Empty or singular window systems produce zero, with counts reported
at native verbosity 2. Use `--more-options no-two-balls` for an exact small-data
two-ball check (`no-one-ball` for GGG). Raw numerator and corrected
correlation are different products. Add `--more-options weights-norm`
to use catalog weights in the signal and window moments.
Narrow windows can amplify noise and roundoff even when nonsingular; check
stability with binning, multipole order, and traversal tolerance.
This corrects angular mode mixing for sampled scalar fields, not galaxy
number-count estimators requiring independent random catalogs. The 2PCF
weighted pair average is unchanged, and `only-2pcf` cannot be combined with
edge correction.

## Notebook usage

For a notebook or another Python program, import the catalog and runner types:

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path("python").resolve()))
from kappa_corr_all_engines import KappaCatalog, RunConfig, run_engine_suite

catalog = KappaCatalog(positions=xyz, kappa=kappa, weights=weights)
config = RunConfig(
    engines=("octree-ggg-omp", "kdtree-omp"),
    output_dir=Path("Output_kappa"),
    threads=16,
)
results = run_engine_suite(catalog, config)
```

For corrected notebook results, select compatible engines and set
`edge_corrections=True` (or `result_type="edge_effects"`):

```python
config = RunConfig(
    engines=("octree-ggg-omp", "octree-2balls-omp"),
    output_dir=Path("Output_kappa_edge"),
    edge_corrections=True,
    threads=16,
)
results = run_engine_suite(catalog, config)
zeta = results["octree-2balls-omp"]["zeta_edge_complex_1"]
```

`positions` must have shape `(N, 3)` and contain Cartesian unit vectors.
`kappa`, optional `weights`, and optional `mask` must have shape `(N,)`.

This driver intentionally excludes shear, Ly-alpha, cross-catalog, periodic-box,
and physical-3D methods. Those engines require different catalog fields or
estimator semantics and cannot correctly consume a scalar angular kappa map.
Use `python/lya_corr_all_engines.py` for DESI Ly-alpha catalogs, including all
radial/3D OpenMP and MPI variants through `set_forest_catalog()`.
The physical-3D aliases `octree-ggg-3d-omp` and `octree-ggg-3d-mpi` are
recognized as aliases, not incorrectly reported as C/Cython build mismatches.

Each engine writes `histograms.npz`, `result.json`, and compatible text files
under `<outdir>/<engine>/python/`. The top-level `summary.json` records timings,
failures, MPI rank count, and the one-read/one-registration contract.
