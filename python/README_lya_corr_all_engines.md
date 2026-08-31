# Ly-alpha Multi-Engine Driver

Run commands below from the cTreeBalls project root. Use a Python environment
with NumPy, SciPy, Astropy, Matplotlib, and the freshly rebuilt `cyballs`.
MPI additionally requires `mpi4py` linked to the same MPI as cTreeBalls.
`DEFDIMENSION=3` remains required even for radial-only searches.

Build with `LYAFORESTOMPON=1`; additionally enable `LYAFORESTMPION=1` for MPI.
After changing these settings:

```bash
make -j4 cballs cyballs-static-lib
CBALLS_STATIC_LIBRARY_READY=1 python3 setup.py build_ext --inplace --force
# This tree keeps a second import location; synchronize it when present:
cp cyballs*.so python/
python3 python/lya_corr_all_engines.py --list-engines
```

## Download One DESI File

The [official DR1 delta data model](https://data.desi.lbl.gov/doc/releases/dr1/vac/lya-deltas/)
documents these files. A small file in the requested directory is
`delta-1019.fits.gz` (12 forests, 157534 compressed bytes when verified).
Download only this file first, not the entire release:

```bash
bash examples/download_desi_lya_example.sh data/desi-lya-example
```

Equivalent standalone command:

```bash
wget -c -P data/desi-lya-example \
  https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/v1.0/delta-lya-0-0/Delta/delta-1019.fits.gz
```

The driver supports the DR1 image layout: shared `LAMBDA`, per-forest
`METADATA`, `DELTA_BLIND` (or `DELTA`), and `WEIGHT`. It does not pretend that
these files are six-column Cartesian catalogs or legacy one-table-per-forest
FITS. Legacy lya2pcf pickles can first be converted with
`addons/lya_forest_omp/lya2pcf_to_cballs.py` and loaded using `--ascii`.

## Small First Runs

```bash
# Three dedicated 2PCF engines: 3D, radial scan, radial tree.
python3 python/lya_corr_all_engines.py \
  --fits data/desi-lya-example/delta-1019.fits.gz \
  --max-forests 6 --pixel-stride 30 \
  --engine all-omp --threads 2 --rp-bins 8 --rt-bins 8 \
  --save-catalog data/desi-lya-example/pixels.npz \
  --output Output_lya_2pcf

# All seven OpenMP engines, including standalone and combined 2PCF/3PCF.
python3 python/lya_corr_all_engines.py \
  --catalog data/desi-lya-example/pixels.npz \
  --engine all-omp --statistics both --threads 2 --r3-max 60 \
  --output Output_lya_both

# All fourteen engines; OpenMP engines run on rank zero only.
python3 python/lya_corr_all_engines.py \
  --catalog data/desi-lya-example/pixels.npz \
  --engine all --statistics both --mpi-ranks 2 --threads 2 \
  --mpi-extra-arg='--oversubscribe --map-by slot --bind-to none' \
  --output Output_lya_mpi
```

Use a new output directory per run: existing histograms are rejected to prevent
reading stale results after a failed run. `--mpiexec /path/to/mpiexec` selects
the launcher. Alternatively invoke the script directly under `mpiexec -n 2`.
Oversubscription/mapping arguments above are Open MPI-specific; omit or change
them for another implementation. The MPI child must import the same `cyballs`
and use the same filesystem output path on every rank.

A no-download demonstration and notebook are available:

```bash
python3 python/lya_corr_all_engines.py --synthetic \
  --engine all-omp --statistics both --output Output_lya_synthetic
```

Open `examples/lya_corr_all_engines.ipynb` in the same environment.

## Engines and Parameters

Replace `-omp` with `-mpi` for the MPI counterpart of each row.

| Engine | Statistic | Bin controls |
|---|---|---|
| `lya-2pcf-omp` | 3D xi(r_parallel, r_transverse) | `--rp-max`, `--rt-max`, `--rp-bins`, `--rt-bins` |
| `lya-3pcf-omp` | 5D zeta(r1,r2,theta1,theta2,mu) | `--r3-max`, `--r3-bins`, `--theta-bins`, `--mu-bins` |
| `lya-2pcf-3pcf-omp` | Both 3D estimators | Both sets |
| `lya-1d-2pcf-omp` | Radial xi(abs(chi_j-chi_i)) | `--rp-max`, `--rp-bins` |
| `lya-1d-tree-2pcf-omp` | Same radial 2PCF, interval tree | Same |
| `lya-1d-3pcf-omp` | Radial zeta(signed lag1, signed lag2) | `--r3-max`, `--r3-bins` |
| `lya-1d-2pcf-3pcf-omp` | Both radial estimators | Both sets |

Engine groups: `all`, `all-omp`, `all-mpi`, `all-1d`, `all-3d`.
`--statistics 2pcf` (default) selects only dedicated 2PCF engines.
`--statistics 3pcf` selects dedicated 3PCF engines.
`--statistics both` permits dedicated and combined methods.
Explicitly requested unavailable or incompatible engines cause an error.

`octree-3pcf-3d-omp/mpi` are intentionally not included. Their scalar Legendre
estimator differs from the forest 5D estimator and their LOS exclusion does
not require all three forests to be distinct. Angular convergence, shear,
and point-count engines are likewise not interchangeable forest estimators.

## Input and Scientific Contract

Coordinates are observer-centered Cartesian comoving distances in Mpc/h.
FITS absorption redshift is `LAMBDA / 1215.67 - 1`; quasar redshift is not
substituted for pixel redshift. The fiducial conversion is Astropy
`FlatLambdaCDM(H0=100*h, Om0=omega_m, Tcmb0=0)`, default `omega_m=0.315`,
`h=0.674`. These are explicit example choices, not a reproduction of every
DESI analysis setting. `--omega-m`, `--h`, `--z-min`, `--z-max` are configurable;
all distances/range arguments use the same units. ASCII/NPZ coordinates must
already be in those units; no cosmology conversion is applied to them.

The reader preserves the supplied delta and weights, removing pixels with
nonfinite delta/weight or nonpositive weight. `--max-forests` takes an initial
subset, and `--pixel-stride` subsamples valid pixels within each forest; neither
is an unbiased survey selection nor a rebinning procedure. They are intended
for demonstration/testing. Omit both for a full input run.

DESI IDs exceed exact float64 integer precision. They remain int64 in Python
and cached NPZ; the Cython API compacts them to distinct small integer IDs
without changing same-quasar membership. Repeated input files or LOS IDs
across FITS files are rejected, avoiding accidental double counting. Do not
mix different delta regions or releases in one autocatalog.

The downloaded sample uses `DELTA_BLIND` with `BLINDING=desi_y1`. The driver
records and displays this; it does not unblind data, refit continua, recenter
delta, reweight redshift, apply distortion matrices, or estimate covariance.
This is a raw weighted correlation measurement, not the complete published
DESI cosmological analysis.

Pairs exclude the same forest. Triplets require three different forests.
Radial methods ignore transverse separation. Their results cannot be equated
with 3D results from a general sky catalog. All bins contain normalized
`numerator / denominator`, with zero returned for empty bins.

## Reuse in Python

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("python").resolve()))
from lya_corr_all_engines import ForestCatalog, RunConfig, run_engine_suite

catalog = ForestCatalog(positions, delta, weights, forest_ids).normalized()
results = run_engine_suite(
    catalog,
    RunConfig(("lya-2pcf-omp", "lya-1d-tree-2pcf-omp"),
              Path("Output_my_forests"), threads=4),
)
```

Or directly: `balls.set_forest_catalog(positions, delta, weights, forest_ids)`.
Do not also set `infile` or `infileformat`. `struct_cleanup()` releases the C
run but retains registered NumPy arrays; `clear_catalogs()` releases both.
One NumPy catalog is registered per rank, then reused across methods; each
method still constructs its own C-owned bodies/tree/histograms.

## Outputs and Costs

Each engine directory contains its native text histograms (including raw
numerators and denominators), log, and used parameters. `summary.json` records
catalog provenance, cosmology, pixel selection, configuration, timings and
same-statistic comparison metrics. Timing is end-to-end C wall time, including
startup, tree building, search, writing and cleanup, excluding FITS conversion
and Python plotting. It is not isolated kernel CPU time.

Comparison CSVs store bin indices, reference, candidate, absolute difference
and signed relative difference `(candidate-reference)/reference`. Zero
reference values yield `nan`, not misleading zero relative differences.
The reference is the first selected method producing the same statistic.
The summary separately reports numerator/denominator differences.

PNG figures keep the four statistic families separate. The 5D 3PCF image is
a weighted projection: sum all angular numerators and denominators at fixed
(r1,r2), then divide. It is not the unweighted mean of normalized bins.
Empty bins are blank in plots but remain zero in native text output.
Use `--no-plots` on headless systems without Matplotlib.

Radial 3PCF ignores transverse separation and may enumerate many more
neighbors than 3D 3PCF. Start small and increase gradually. Each MPI rank
replicates the catalog, tree and histogram memory. `--max-hist-mib` provides
a configurable per-rank histogram guard; this does not include catalog/tree
storage or guarantee that an arbitrarily large run fits in RAM.

## Regression

```bash
make test-lya-corr-all-engines
```

Tests cover DESI layout, 64-bit IDs, filtering, cached input, invalid domains,
same-forest exclusions, all compiled OpenMP forest engines against independent
2PCF/3PCF oracles, sparse comparisons, and one-time registration. The MPI CI
job also runs all built forest engines together on two ranks.
