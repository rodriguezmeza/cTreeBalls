# In-memory kappa correlations

`kappa_corr_all_engines.py` reads and converts a convergence catalog once,
registers the resulting NumPy arrays with one `cyballs.cballs` object, and runs
any number of scalar angular search methods without rereading the catalog.

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
available engines in a group. Multi-rank execution needs `mpi4py` built against
the same MPI implementation used for cTreeBalls:

```sh
python3 python/kappa_corr_all_engines.py \
  --fits map.fits \
  --engine all-mpi \
  --mpi-ranks 4 \
  --threads 8 \
  --outdir Output_kappa_mpi
```

The script may also be launched directly by MPI:

```sh
mpiexec -n 4 python3 python/kappa_corr_all_engines.py \
  --fits map.fits --engine octree-ggg-mpi --threads 8
```

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

`positions` must have shape `(N, 3)` and contain Cartesian unit vectors.
`kappa`, optional `weights`, and optional `mask` must have shape `(N,)`.

This driver intentionally excludes shear, Ly-alpha, cross-catalog, periodic-box,
and physical-3D methods. Those engines require different catalog fields or
estimator semantics and cannot correctly consume a scalar angular kappa map.

Each engine writes `histograms.npz`, `result.json`, and compatible text files
under `<outdir>/<engine>/python/`. The top-level `summary.json` records timings,
failures, MPI rank count, and the one-read/one-registration contract.
