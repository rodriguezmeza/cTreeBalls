# Hierarchical scalar 3PCF reuse

The opt-in `scalar-pivot-reuse` option is shared by `octree-2balls-omp`,
`kdtree-2balls-omp`, `balltree-2balls-omp`, and their `-mpi` counterparts.
It applies to scalar angular 3PCF, including combined 2PCF/3PCF runs.
The independent 2PCF traversal is retained.

## Runtime controls

These are environment variables, not Makefile definitions:

| Control | Default | Range | Meaning |
|---|---:|---:|---|
| `CBALLS_SCALAR_PIVOT_TOL` | 0.1 | [0,3] | Maximum modeled phase perturbation of the highest pair harmonic, radians |
| `CBALLS_SCALAR_BIN_THETA` | 0 | [0,1] | Optional internal bin assignment allowance, in bin widths |

Both require finite values. Set the same values on every MPI rank. Invalid
values and rank mismatches fail collectively before traversal. Existing
`CBALLS_SHEAR_*` controls affect the shear engines only.

Add `scalar-pivot-reuse,no-smooth-pivot` to `options`, or pass
`--more-options scalar-pivot-reuse --no-smooth-pivot` to a compatible Python
driver. `theta` must be positive; the two environment variables govern the
new path's acceptance. `nsmooth=1` often exposes more useful pivot hierarchy
in binary trees, at a higher tree-storage cost; benchmark it against the
usual leaf setting. `THETA` and `BALLS4SCANLEV` retain their existing build
roles. No new build flag is required beyond scalar 3PCF support (`TPCFON=1`).

Zero phase budget, `no-one-ball`, `no-two-balls`, smoothing, periodic runs,
and `only-2pcf` disable the new traversal and retain the existing path.
For an exact scalar reference use `no-two-balls,no-smooth-pivot`.
`legacy-one-ball`, `dual-node-direct-triples`, and unrelated engines reject
the new option. This does not change the legacy smoothing or scan algorithms.

## What is reused

Accepted neighbor Fourier moments, normalization sums, and second moments
are retained down the pivot hierarchy. Only unresolved neighbor nodes are
revisited. Bounding intervals mark individual unresolved radial bins rather
than discarding all work below the furthest unresolved neighbor. A radial
pair is reduced at the first pivot node where both bins are complete; an
ownership mask prevents descendants from repeating that contribution.
Mixed completed/uncompleted pairs wait for both legs.

Four native angular products (`coscos`, `sinsin`, `sincos`, `cossin`) and
window modes retain their conventions. Repeated-neighbor subtraction and
cardinality checks are unchanged. Minimum/maximum radial cuts stay strict;
a positive bin theta can move assignments across internal bin edges.

In 3D, the native basis selects the least-aligned Cartesian reference axis.
Acceptance rejects pivot caps that could cross that axis-selection boundary.
It bounds changes in both the separation vector and the native basis. An
accepted source's phase remains frozen in its original coordinates, with a
bound relative to every descendant pivot; there is no chain of repeated
transport errors. At maximum harmonic order M, each leg receives at most
`phase_budget/(2*max(M,1))` radians. This is a phase bound, **not a relative
coefficient-error guarantee**: signed fields, weak coefficients, bin migration,
and poorly conditioned window correction can amplify observable errors.

The shared radial-product reducer now visits contiguous histogram rows per
harmonic, including a SIMD hint for independent symmetric pairs. This layout
change also applies when hierarchical reuse is disabled.

MPI owns each deterministic frontier task on one rank, and reduces raw
histograms before normalization/window correction. Hierarchy storage is live
only for active workers, not every queued task. Moment and ownership-mask
storage has a 64 MiB per-worker cap; checked neighbor frontiers are separate.
An allocation or storage-limit failure is returned collectively.

`getRunMetadata()['scalar_hierarchical_reuse']` freezes effective controls,
activation, actual/represented radial-pair counts, and parent reductions.
Counters are global on the publishing MPI rank. `dual-node-profile` prints
them along with existing native timing diagnostics. Results remain
numerically unqualified until measured against a reference for the intended
observable and geometry; existing qualification infrastructure is retained.

## Reproducible measurements

From the repository root (Python must import the rebuilt `cyballs`):

```sh
python3 scripts/benchmark_scalar_pivot_reuse.py --engine balltree \
  --geometry clustered --n 262144 --threads 16 --repeat 3 \
  --leaf 16 --output results/scalar-before
python3 scripts/benchmark_scalar_pivot_reuse.py --engine balltree \
  --geometry clustered --n 262144 --threads 16 --repeat 3 --leaf 1 \
  --reuse --phase-tol 0.1 --bin-theta 0 --output results/scalar-reuse
python3 scripts/benchmark_scalar_pivot_reuse.py --engine balltree \
  --geometry clustered --n 262144 --threads 16 --leaf 1 \
  --exact --output results/scalar-exact
```

Change `--engine` to `octree` or `kdtree`; add `--mpi` under `mpiexec -n 2`
for MPI (requires mpi4py), reducing threads per rank to keep the same CPU
budget. `--statistic both` includes 2PCF, and `--no-edge` measures raw 3PCF
without window correction. `--catalog-npz` accepts arrays `positions`,
`kappa`, and optional `weights`; `--min-sep` and `--max-sep` are in catalog
coordinate units. Supply the same data and radial limits to all three runs.
`--max-n` selects the highest signal harmonic (default 3), and `--bins`
sets the radial-bin count (default 20). NPZ subsamples are deterministic
permutation prefixes controlled by `--sampling-seed`. Both native caches
are disabled for these measurements. Outputs are `.json`
(timings/provenance) and `.npz` (all four raw components, complex signal,
corrected signal and W0 when requested). The script never labels a timing
measurement as a passed accuracy qualification.

Measure full-sky, survey, and clustered geometries separately. More reuse or
fewer products does not itself establish a wall-time speedup. Small catalogs,
strict bin containment, basis-boundary caps, and irregular unresolved frontiers
can make this path slower. Defaults remain unchanged.

## Regression checks

```sh
make test-scalar-pivot-reuse
mpiexec -n 2 python3 -O tests/python/test_scalar_pivot_reuse.py --mpi
```

The retained tests exercise independent ordered triangles and all four raw
components, weighted masks and cross-catalogs, linear/logarithmic bins,
partial parent completion, exact/fallback paths, strict radial cutoffs,
thread determinism, MPI invalid/mismatched controls, and randomized native
phase-bound checks. Current end-to-end numerical validation uses the 3D build.
Masked octree cross-catalog dispatch now retains the requested two catalogs
instead of silently invoking an auto-correlation.

## Measured limits (2026-10-02)

The implementation is experimental and remains opt-in. On one macOS host,
three-repeat median cold native times for 1,048,576 clustered points and
16 threads improved by 1.45x (octree), 2.93x (kd-tree), and 2.12x (ball-tree).
The corresponding MPI runs improved by 1.32–2.49x relative to their previous
MPI implementations, with 2x8 or 4x4 ranks/threads on the same host. These are
**loose-profile timings**, at phase budget 3 and bin theta 1, with 20 log bins
and signal harmonics 0–3; they are not matched-accuracy speedups.

That profile is not qualified for scientific results. On an 8,192-point
clustered subsample, relative L2 errors in the four raw components were
3.18–6.50%; window-corrected errors were much larger, and some engines changed
the finite-output support. Even phase 0.1 / bin theta 0 gave 3.04–4.62%
window-corrected error, despite raw errors of roughly 2e-6–2e-5. Strict bin
containment also imposed substantial runtime cost on the large fixture.

On a one-million-object DES geometry with one shear component used solely as
a scalar proxy, the loose reuse profile was slower for all three engines:
5.91/8.42/5.13 seconds versus 4.60/5.12/4.69 seconds previously (octree/kd-tree/
ball-tree). Its full-catalog window-corrected relative L2 error was 7.50–9.48%.
This measurement supports keeping the original path for that workload.

Always compare raw and corrected observables and their valid-output support
against exact results at the intended geometry, binning, and harmonic order.
A raw-moment phase bound does not constrain the conditioning of the window
solve. No multi-node scaling or general numerical qualification is claimed.
