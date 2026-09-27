# Six-method Ly-alpha optimization benchmark

`benchmark_lya_optimization.py` compares two compiled checkouts using the same
FITS/NPZ catalog, forest prefix, pixel stride, bins and thread count. Both
checkouts must contain a `cyballs` extension compatible with the Python running
the script. Use the same compiler, flags, precision, dependencies and machine
for an algorithm comparison. Preserve the old checkout before rebuilding it.

```bash
python scripts/benchmark_lya_optimization.py \
  --baseline-root /path/to/baseline/cTreeBalls \
  --candidate-root /path/to/updated/cTreeBalls \
  --catalog /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --forest-counts 64 256 1000 --pixel-stride 10 \
  --threads 1 4 16 --radius 200 \
  --pair-bins 8 --radial-bins 4 --polar-bins 4 --mu-bins 8 \
  --warmups 1 --repeats 3 --output results/lya-before-after
```

This can be expensive for 3PCF at 1,000 forests. Start at 64/256, then retain
the same bins and radius for the larger test. Radius units are the input's
Cartesian distance units; the script performs no cosmological conversion.
FITS input needs Astropy. NPZ uses `positions`, `delta`, `weights`, `forest_ids`.
Selection sorts forest IDs and preserves input row order within each forest.

All six original/LOS-tree 3D OpenMP methods run by default. `--methods` selects
a subset. The default kernels are the exact pixel pair walker and segment
triplets (`lya2Kernel=0`, `lya3Kernel=0`), with all slops and pivot smoothing
left at zero. `--baseline-parameters` and `--candidate-parameters` accept JSON
overrides, for example `'{"lya3PivotBlock":8}'` to measure scheduling explicitly.
Keep parameter sets identical when measuring a source-code speedup.

If only the updated build is available, use that same checkout for both roots
and compare its independent direct triplet loop with exact segments:

```bash
python scripts/benchmark_lya_optimization.py \
  --baseline-root . --candidate-root . \
  --catalog /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --forest-counts 64 --threads 1 16 --radius 200 \
  --methods lya-3pcf-omp lya-los-tree-3pcf-omp \
    lya-2pcf-3pcf-omp lya-los-tree-2pcf-3pcf-omp \
  --baseline-parameters '{"lya3Kernel":1}' \
  --candidate-parameters '{"lya3Kernel":0}' \
  --output results/lya-direct-versus-segments
```

For the two pair-only methods, compare `'{"lya2Kernel":0}'` with
`'{"lya2Kernel":1}'` in a separate invocation. This measures existing exact
cell aggregation, not just the latest source changes. Do not apply pair kernel
1 to a 3PCF-only method. A larger kernel number is not necessarily faster.

## Measurement and scientific contract

- Each sample is a fresh subprocess, running the requested warmups and one
  measured native `Run`. Baseline/candidate order alternates per repeat.
  A changed output directory forces recomputation. Processes run sequentially.
- Run wall and total process CPU exclude Python catalog loading and result
  getters. Native timing is retained too. Peak RSS includes the interpreter,
  selected catalog, warmups and native allocations; it is not incremental
  scratch memory. No CPU affinity is imposed. Record external affinity or
  scheduler settings when comparing hosts.
- Every sample is checked against the first baseline for that input/method,
  including later thread counts. Raw numerators, denominators and ratios must
  satisfy per-element `rtol=3e-11`, `atol=1e-10`. Occupancy and finite masks
  must match. This permits floating summation roundoff, not a 5% approximation.
- `benchmark.json` contains every sample, median/range, speedup, memory,
  parameters, build identity, extension hash, catalog hash and error metrics.
  Selected NPZ fixtures, original row indices, raw products and logs remain
  beside it. A failed comparison stops the campaign and marks it failed.
- This tests the supplied workload. Retain the independent release-gate
  oracles as well. Neither this benchmark nor a successful slop calibration
  supplies a universal correlation-error bound.

For approximate geometry or pivot smoothing, use the existing
`tests/python/benchmark_lya_cell_approximation.py` and
`tests/python/benchmark_lya_pivot_frontier.py`. They record errors against exact
output and treat small signals/empty bins explicitly. Legacy `theta`,
`THETA`, `rsmooth`, `nsmooth`, `SMOOTHPIVOT` and `BALLS4SCANLEV` are not synonyms
for the anisotropic forest controls.
