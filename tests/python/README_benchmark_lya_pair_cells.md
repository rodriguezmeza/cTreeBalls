# Ly-alpha 2PCF forest-cell acceleration

`lya2Kernel=1` enables persistent cell-pair traversal for `lya-2pcf-omp`,
`lya-los-tree-2pcf-omp`, and their `2pcf-3pcf-omp` variants. The original
pixel traversal remains the default (`lya2Kernel=0`) and calibration reference.
Radial-only 1D methods do not accept these controls. The original 3D MPI
methods now support pair cells; existing 1D
methods remain available with their original settings.

## Scientific contract

For pixel distances chi_i, chi_j and observer directions n_i, n_j, write
c = clamp(n_i dot n_j, -1, 1). The estimator uses

```
rp = abs(chi_i - chi_j) * sqrt((1+c)/2)
rt = (chi_i + chi_j) * sqrt((1-c)/2)
xi[bin] = sum(w_i delta_i w_j delta_j) / sum(w_i w_j)
```

Only unordered pairs in distinct forests contribute, with `0 <= rp < lya2RpMax`
and `0 <= rt < lya2RtMax`. The reference discovery-sphere cutoff is preserved.
Zero-weight pairs count toward the accepted pair count but contribute no weight.
Empty bins have xi=0. Forest IDs remain full-width integers.

Each node contains pixels from one forest and retains actual radial/Cartesian/LOS
extrema, weight sum W, and weighted-field sum Q. An accepted pair of nodes adds
`N += Qp*Qq`, `D += Wp*Wq`, and `count += np*nq`. Signed fields never define the
geometric center. With zero slop, all enclosed pixel pairs must be certified to
lie in one bin. Unresolved products split either node or evaluate an exact small
tile. Floating-point summation is reassociated, so exact means the same estimator
and bin membership, with roundoff-level differences in raw sums.

The two approximation controls are **independent of every lya3 control**:

| Parameter | Default | Meaning |
|---|---:|---|
| `lya2Kernel` | 0 | 0 reference pixel traversal; 1 persistent cell pairs |
| `lya2RpSlop` | 0 | Allowed parallel bin leakage, fraction of RpMax/RpBins |
| `lya2RtSlop` | 0 | Allowed transverse bin leakage, fraction of RtMax/RtBins |

Slops must be finite in [0,1]. Positive slop requires kernel 1. It relaxes only
bin assignment, never forest exclusions or domain cuts. It is **not a relative
xi error bound**. A small slop can cause large relative errors near a zero crossing.
The native output and `getRunMetadata()['lya_2pcf']` record these settings.

## How the 3PCF improvements transfer

* A hierarchy is built once per search. Combined runs with `lya3Kernel=3/4`
  reuse the same hierarchy for pairs and triples. Other combined kernels run an
  independent pair pass, followed by 3PCF discovery using only its own radius.
* Padded Cartesian bounds reject separated cells before angular calculations;
  bounds on rp and rt reject additional products before descendant work.
* Observer-angle enclosures pass down the recursion when LOS extrema differ only
  at roundoff scale. A child remains inside its inherited enclosure. Actual
  angular changes trigger tighter bounds. This works for bent forests too.
* Narrow forests with sparse pixels and fine bins use conservative radial
  windows in the distance ordering. This avoids subdividing products that
  cannot usefully aggregate; ambiguous pixel bins still use exact arithmetic.
* Both cells can aggregate and split. The split score measures radial and LOS
  uncertainty in units of the two bin widths.
* Small direct tiles reuse certified angular intervals to avoid repeated pixel
  dot products and square roots. Ambiguous pixels use the reference arithmetic.
  SIMD geometry is separated from histogram writes. Forest-pair
  rows use dynamic OpenMP scheduling and fixed ordered reductions.
* Every Cartesian product is visited once. A quadratic pixel cache would have
  no useful reuse here, so the pair path uses stack-local bounds and tiles.
  It does not allocate the 3PCF per-pivot geometry caches or scalar 3PCF tensors.

Tree and histogram allocations use checked dimensions and the existing memory
budget. Combined persistent 3PCF runs still require their usual 3PCF caches.
Non-double compute precision or fast-math profiles descend to reference pixel arithmetic instead
of accepting interval-certified aggregates; scientific calibration targets the
normal double, non-fast-math build.

## Reproducible benchmark

From the repository root, using a Python environment with NumPy (and the existing
FITS reader's dependencies for FITS input):

```bash
python tests/python/benchmark_lya_pair_cells.py \
  --fits /path/to/catalogs/lya_15_xyz_raw_with_losid.fits \
  --max-forests 1000 --pixel-stride 10 \
  --method lya-2pcf-omp --threads 16 \
  --rp-max 160 --rt-max 160 --rp-bins 50 --rt-bins 50 \
  --warmups 1 --repeats 3 \
  --case 0:0 --case .001:.001 --case .01:.01 \
  --outdir results/lya-pair-cells-desi1000-t16
```

Replace `/path/to/catalogs` with your catalog directory. Use the same cuts and
bin counts as your intended workload. The default benchmark
case is exact `0:0`; adding positive-slop cases is an explicit calibration request.
Use `--method lya-los-tree-2pcf-omp` to compare against the legacy LOS-tree engine.
Use a new output directory for every run. NPZ, six-column ASCII and synthetic
catalogs are supported too. `--baseline-cballs /path/to/old/cballs` compares timing
against a previous executable; both executable hashes are retained.

For combined runs:

```bash
python tests/python/benchmark_lya_pair_cells.py \
  --catalog results/lya-pair-cells-desi1000-t16/catalog.npz \
  --method lya-2pcf-3pcf-omp --kernel3 4 --threads 16 \
  --r3-max 160 --r3-bins 4 --theta-bins 4 --mu-bins 4 \
  --case 0:0 --outdir results/lya-pair-cells-combined
```

This times the **entire combined run**, so a dominant 3PCF cost can hide a large
pair speedup. Start combined tests on a smaller catalog. The benchmark verifies
that the 3PCF product is unchanged while varying only the pair controls.

`summary.json` retains commands, input/executable hashes, settings, raw counts,
work counters, whole-process wall/CPU time, peak process RSS and per-case accuracy.
Reference and candidate cases alternate within each repetition to reduce timing
drift. Each native output directory retains its histograms, metadata and log. Per-bin
comparison tables retain reference/candidate xi and raw numerator/denominator.
FITS conversion is excluded from native timing; native input/tree/output are
included. CPU is total process CPU, not elapsed time divided by thread count.
The C log's `build_CPU` and `search_CPU` split forest construction and pair traversal;
work counters (`nodes`, `pruned`, `angle_reuses`, etc.) are counts, not timings.

Accuracy acceptance checks finite output, identical counts and occupancy, at most
5% relative xi error where `abs(xi_ref)>1e-6`, and absolute error at most `5e-8`
in other occupied bins. Thresholds are configurable. Exact cases also require
raw sums to match within roundoff (`rtol=3e-11`, `atol=1e-9`). Every repeat must
produce identical arrays. A failed exact case exits 1. `--require-accepted` exits
2 unless at least one requested positive-slop case passes. A passing calibration
only supports the measured catalog, domain, binning and build.

## Hierarchical radial moments

`lya2Kernel=1` and `lya3Kernel=5` select certified pair-range sums and adaptive
radial/polar moment combinations for the original 3D OpenMP/MPI methods.
Combined runs share the forest tree; zero slop preserves the histogram estimator.
Persistent kernels 3/4/5 also support MPI. The all-engines driver accepts
`--lya2-kernel 1 --lya3-kernel 5 --lya3-pivot-cell-max 8`; these controls apply
only to eligible 3D Ly-alpha statistics. Default kernels remain 0.
See the [implementation and benchmark guide](../../docs/LYA_HIERARCHICAL_REUSE.md).
