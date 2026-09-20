# Ly-alpha forest all-engines driver

Run this driver from `tests/python`, not the Cython-only `python/` directory.
Native `timing_report.txt` and `summary.json["timings"]` separate parameter/thread
setup and MainLoop wall/CPU time. MPI wall time is the maximum rank time and CPU
is summed across participating ranks. Per-rank measurements, thread counts and
native parameters are retained. Catalog registration and Python analysis are
excluded; requested native histogram output is included. The historical
`wall_seconds` field also includes cleanup and synchronization. Upstream
reference timings retain their separately documented scope and are not mixed
into the native timing table.

`lya_corr_all_engines.py` reads a forest catalog once and runs the active
cTreeBalls forest and physical-3D multipole methods. `LYAFORESTOMPON=1`
provides twelve OpenMP names; `LYAFORESTMPION=1` provides eight MPI counterparts.
`OCTREE3PCF3DOMPON=1` and `OCTREE3PCF3DMPION=1` add two physical-3D methods.

List the methods compiled into the current extension:

```sh
python3 tests/python/lya_corr_all_engines.py --list-engines
```

## Input

NPZ catalogs contain `positions`, `delta`, `weights`, and integer
`forest_ids`. Six-column ASCII input is `x y z delta weight forest_id`.
`--fits` auto-detects DESI image-layout or eBOSS/PICCA forest-HDU delta FITS
files, including `.fits.gz` and quoted globs. Use `--fits-layout desi` or
`--fits-layout eboss` to require a format; do not mix layouts in one run.
The driver converts absorption redshift to comoving Mpc/h using the explicitly
reported fiducial cosmology. See the eBOSS and preprocessing details below.

## Examples

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz --engine all-omp --statistics 2pcf \
  --threads 8 --output Output_lya_2pcf
```

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz \
  --engine lya-1d-tree-3pcf-omp octree-3pcf-3d-omp \
  --statistics 3pcf --threads 8 --output Output_lya_3pcf
```

The radial and physical-3D estimators have different geometry and are timed
and reported together, not asserted to be numerically identical.

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz \
  --engine lya-2pcf-mpi lya-1d-tree-2pcf-mpi \
  --statistics 2pcf --mpi-ranks 2 --threads 4 \
  --output Output_lya_mpi
```

Use repeated `--mpi-extra-arg` arguments, for example
`--mpi-extra-arg='--bind-to none'`. The launcher must match the MPI libraries
used by cyballs and mpi4py. OpenMP-only methods run on rank zero, MPI methods
on every rank, and external references/analysis once on rank zero.

## Scientific contracts

The three `lya-los-tree-*` names use a 3D octree to discover neighboring
forests and then query each forest's 1D radial tree. They preserve the
**anisotropic 3D** estimator, unlike `lya-1d-tree-*`. `all-3d` and `all-tree`
include these new methods; `all-1d` does not. They are OpenMP-only at present.
Rebuild the Cython extension before using newly registered names.

For a like-for-like 3D comparison on the same retained catalog:

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz --statistics both --threads 8 \
  --engine lya-2pcf-3pcf-omp lya-los-tree-2pcf-3pcf-omp \
  --output Output_lya_los_comparison
```

The anisotropic methods bin parallel and transverse separation and exclude
same-forest pairs. Radial methods use signed or absolute line-of-sight lags as
documented by each engine. `lya-1d-tree-same-los-2pcf-omp` is deliberately
different: it accepts pairs within one forest, normalizes each occupied
forest/bin, and then averages forests equally.

Three-point forest methods require three distinct forest IDs. The physical-3D
multipole methods use `exclude-all-same-los`; they are not interchangeable with
the anisotropic five-dimensional estimator.

`summary.json` records catalog provenance, timings, products, comparisons, and
plot paths. Empty denominator bins publish finite zero in native output and
appear as missing data in plots.

## eBOSS and preprocessing

eBOSS tables need `LOGLAM` (log10 wavelength in Angstrom), `DELTA` (or
`DELTA_BLIND`), and `WEIGHT`. Each forest header needs `RA`, `DEC`, and an
integer `THING_ID`, `LOS_ID`, or `TARGETID`. IDs retain signed 64-bit precision;
duplicate IDs across HDUs/files are rejected. Header angles default to radians;
`--eboss-angle-unit deg` is an explicit alternative. DESI angles use their FITS
column units. **Inputs must be extracted delta products, not raw flux spectra.**
Continuum/delta extraction must already have been performed.

Absorption redshift is `lambda/1215.67 - 1`. The reported flat LambdaCDM model
(`--omega-m 0.315 --h 0.674`, radiation disabled) converts it to comoving Mpc/h.
All engines, including the external reference, receive the same retained
positions, deltas, weights and forest IDs. `--max-forests`, `--z-min`,
`--z-max`, and `--pixel-stride` define subsets; stride is not rebinning.

By default, deltas and positive finite weights are unchanged. Optional
`--project-delta` removes each retained forest's weighted mean and weighted
log-wavelength slope, following equation (2.7) of arXiv:2507.00129. Optional
`--redshift-weight-exponent E --weight-z-ref Z` multiplies input weights by
`((1+z)/(1+Z))**E` before projection. For the preprocessing convention in the
supplied `delta_reader_eboss.py`, use:

```sh
--omega-m 0.3153 --project-delta --redshift-weight-exponent 1.9 --weight-z-ref 2.25
```

Do not repeat weighting/projection already present in your products. Projection
here applies after selection/stride, not to discarded pixels. Matching the
original full-forest preprocessing requires retaining those pixels. Distance
integration uses Astropy, not the upstream interpolated distance table.
Selection, cosmology, blinding and preprocessing are recorded in `summary.json`
and any `--save-catalog` NPZ.

## Compare against lya2pcf

`--lya2pcf-source PATH` enables the actual upstream Numba CPU pair kernel from
that checkout. It does not modify the reference source, import its command-line
readers, require CUDA, or load pickled forest objects. Each reference binning
runs in a fresh process because Numba captures module globals on compilation.
NumPy, SciPy, Numba and the upstream dependencies must be installed in the driver
Python environment; covariance also needs Healpy and the upstream
`post_processing.py` dependencies, including fitsio.

A small eBOSS validation of every available native method, including MPI:

```sh
python3 tests/python/lya_corr_all_engines.py \
  --fits '/path/to/eboss/deltas/*.fits.gz' --fits-layout eboss \
  --max-forests 64 --pixel-stride 8 \
  --engine all --statistics both --mpi-ranks 2 --threads 4 \
  --lya2pcf-source /Users/mar/Documents/Codex/lyman_alpha/lya2pcf_2026-03-14 \
  --rp-max 200 --rt-max 200 --rp-bins 50 --rt-bins 50 \
  --r3-max 20 --r3-bins 4 --theta-bins 4 --mu-bins 4 \
  --reference-covariance --reference-nside 32 \
  --rtol 1e-8 --atol 1e-12 --relative-floor 1e-12 \
  --fail-on-mismatch --output Output_lya_eboss_validation
```

For a non-MPI run, replace `--engine all --mpi-ranks 2` with
`--engine all-omp`. `--statistics 2pcf` selects standalone 2PCF names;
`both` also includes combined methods, comparing both outputs. Only compiled
methods are selected by `all`; inspect `--list-engines` first. Use a fresh
output directory each time.

| Reference name | Valid comparison |
| --- | --- |
| `lya2pcf-cpu` | Anisotropic 3D 2PCF, including `lya-los-tree-*` |
| `lya2pcf-kernel-radial` | Cross-forest radial-only 2PCF: collinear sightlines and one transverse bin |
| `lya2pcf-kernel-same-los` | Within-forest 2PCF: disjoint pixel slices, no self pairs, equal-forest normalization |

The last two are explicitly adapted kernel checks, not native 1D features of
lya2pcf. There is no external 3PCF implementation in this reference checkout;
3PCFs are compared among native engines, separately by estimator family.
The anisotropic reference requires each forest ID to lie on one observer
sightline. Neighbor enumeration uses an angular tree and orders equal-RA
forests by integer ID, avoiding the upstream strict-RA tie omission. The pair
kernel, including its small-angle approximation, is unchanged. Pairs extremely
close to bin boundaries can consequently differ from native geometry; such
differences are reported rather than hidden.

Start small. Radial-only cross-forest references have no transverse cut and
can require all forest pairs. The same-LOS adapter is a correctness reference,
not an optimized benchmark. Full 3PCF catalogs can also be expensive. For a
production 3D 2PCF comparison select, for example,
`--engine lya-2pcf-omp lya-los-tree-2pcf-omp --statistics 2pcf`.
`--reference-timeout SECONDS` bounds each worker; default is no timeout.
Worker progress and tracebacks are retained in `lya2pcf_reference/*.log`.

## Comparisons and plots

All pairs within each compatible product family get CSVs with bin indices,
reference, candidate, signed difference and relative difference. `summary.json`
reports errors in correlation, numerator and denominator, occupied-bin
disagreements, and tolerance checks. The criterion is
`abs(candidate-reference) <= atol + rtol*abs(reference)`; occupancy must match.
`--fail-on-mismatch` exits nonzero after writing the comparison report.
Relative errors below `--relative-floor` are NaN in CSVs and omitted from
plots, not presented as zero error. No projection/interpolation is used for
bin-by-bin equality checks.

Plots include estimator maps, flattened 3PCFs, per-bin absolute and percentage
differences, `r^2 xi(r_parallel,r_transverse)` maps, and angular wedges analogous
to `two_point_analysis.ipynb`. `--no-plots` keeps numerical outputs only;
`--no-analysis` omits the additional 3D 2PCF archives and wedges.
The three-panel comparison figures include the full, unprojected histogram
bins. Separate multipole radial maps show only `ell=0`, not a sum of unrelated
multipoles; the numerical comparisons retain every multipole.

Wedges default to mu boundaries `0,0.5,0.8,0.95,1`, 50 radial bins, and maximum
`min(rp_max,rt_max)`. Configure `--wedge-mu-edges`, `--wedge-bins`, and
`--wedge-r-max`. Bin-area overlap uses `--wedge-subsamples 10` subpixels per
dimension; use `100` for the notebook's sampling resolution. With covariance,
wedges use inverse diagonal variance and propagate `W @ C @ W.T`. Without
covariance they use pair-weight sums and do not invent statistical error bars.
Empty wedges are NaN.
Bins with zero estimated variance are excluded from inverse-variance wedges,
not given infinite weight. Very small or degenerate samples may consequently
show no usable wedge bins even though some measured correlations are populated.

`--reference-covariance` computes the upstream unsmoothed weighted-subsampling
covariance from HEALPix partial histograms. Pairs belong to the sky region of
their first forest in RA/ID order; `--reference-nside` defaults to 32. At least
two occupied regions are necessary. This is a shared reference covariance,
not independently estimated by each native engine, and not automatically a
well-conditioned precision matrix. No covariance smoothing or 3PCF uncertainty
model is silently applied. Alternatively supply `--covariance covariance.npy`
(or NPZ/FITS). `--max-hist-mib` bounds dense analysis/reference workspaces;
tree/catalog storage is additional.

Each `analysis_2pcf/ENGINE.npz` contains correlation, raw sums, bin metadata,
wedge centers, wedge correlations/errors, plus covariance when available.
These archives can be loaded directly in notebooks using
`np.load(..., allow_pickle=False)`.
Keep the sibling helpers `lya_fits.py`, `lya_reference.py`, and `lya_analysis.py`
with the driver when copying this workflow to another installation.

## Distortion matrix

The matrix in equations (2.5)-(2.7) is a forward model:
`xi_observed_model = D @ xi_unprojected_model`. It is not an inverse correction
to data and must not be applied to already distorted measurements. Add these
options to a 3D 2PCF run:

```sh
--distortion-matrix /path/to/distortion.npy \
--model-correlation /path/to/unprojected_model.npy \
--covariance /path/to/covariance.npy
```

Supported inputs are NPY; NPZ keys `distortion`, `model` (or `correlation`),
`covariance`; or FITS table columns `DM`, `DA`, `CO`. Currently matrices must
be square on the same non-negative grid as the measurement: flatten with
parallel-bin major, transverse-bin minor. Models must be finite, including
unoccupied input bins. Dimensions and supplied NPZ/FITS bin metadata are
checked. NPY has no metadata, so matching bin order/ranges, cosmology,
selection, weights and continuum convention is the user's responsibility.
All-NaN distortion rows are allowed only for unoccupied output bins; partially
invalid rows are rejected.

Analysis saves the forward-distorted model and residuals, overlays the model
on wedges, and plots separate absolute/relative model residuals. It does not
change any native or reference measurement. Computing a new distortion matrix
remains an upstream CUDA workflow (`distortion.py`); this driver reads that
output and does not claim a new CPU matrix estimator.

## Timing and tests

Native times include C startup, tree/search, output and cleanup, excluding
FITS conversion. Reference `wall_seconds` measures warmed CPU pair traversal
and deterministic reduction; `worker_wall_seconds` also includes imports,
JIT, covariance and I/O. Setup/JIT/analysis times are separately recorded.
These scopes differ; neither is silently reported as an end-to-end speedup.
The reference uses the requested number of CPU threads on rank zero,
not `MPI ranks * threads`; no GPU timing is reported.

Run regression and optional upstream integration tests:

```sh
LYA2PCF_SOURCE=/path/to/lya2pcf python3 -m pytest \
  tests/make_tests/test_lya_corr_all_engines.py \
  tests/make_tests/test_lya_analysis.py
```
