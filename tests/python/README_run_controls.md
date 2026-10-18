# Speed, memory, accuracy and plotting controls

These notes describe the current native comparison drivers and CPU benchmarks.
They do not change the engine defaults or qualify a new speed profile.

## Measure the requested error

For a limit on **every nonzero coefficient**, use
`abs(fast - exact) / abs(exact)` and take the maximum over all required bins,
multipoles and components. Use complex magnitudes for complex coefficients;
comparing only their real parts is insufficient. A ratio of `0.10` means 10%.
Exact-zero references must remain zero. Report nonfinite coefficients and
unsupported bins explicitly; do not silently drop them or apply a denominator
floor. Small nonzero coefficients can dominate this maximum.

The reference must use the same selected objects, sampling seed, weights, field
centering/projection, mask, radial edges, geometry, normalization, correction
policy and multipole orders. Test each catalog to which a claimed bound applies.
A result on a smaller sample does not certify a larger one. Agreement between
cTreeBalls and dual node with approximations enabled measures agreement between
those runs, not either run's error relative to exact computation.

The all-engine drivers' `--qualification-rtol` controls a **relative L2** target.
Setting it to `0.10` does not establish a 10% limit on each coefficient. Without a
matching reference the approximation penalty is `UNMEASURED`.

Kappa/Ly-alpha benchmark comparison reports distinguish:

- `median_relative`, `p95_relative`, `max_relative`: historical symmetric
  difference `2*abs(fast-reference)/max(abs(fast)+abs(reference), atol)`.
- `max_reference_relative`: maximum unfloored reference-denominator ratio over
  finite nonzero references, in fractional units. The `reference` column names
  the baseline; this number is an exact error only if that baseline is exact.
- `reference_nonzero_bins`, `zero_reference_mismatches`, `nonfinite_bins`:
  coverage and failure counts. Inspect these alongside the ratio.
- `coefficient_status`: checks the reference-relative ratios against `rtol`,
  requires zero references to remain zero, and fails on any nonfinite input.
  Empty arrays are unmeasured. This is supplemental; it does not change the
  historical `status`/exit-code gate.
- `status`: the existing NumPy `allclose` gate,
  `abs(fast-reference) <= atol + rtol*abs(reference)`. Use `--atol 0` with
  `--rtol 0.10` for a strict 10% comparison against a matching exact reference.
  Equal-error sweeps use their own `--equal-error-rtol` and
  `--equal-error-atol`; set the latter to zero for the same strict requirement.

`--no-exact` in the kappa benchmark does not independently run an exact baseline.
Keep exact and approximate runs in different output directories. The shear addon
CPU benchmark records timings and checksums, not full per-coefficient accuracy;
a matching checksum is not an error bound. Use saved all-engine qualification
packets or complete coefficient arrays for scientific accuracy comparisons.
ENCORE reports retain their existing metric and physical 3D estimator contracts.

## Scalar two-ball speed controls

These controls apply to eligible angular scalar 3PCF paths in
`octree-2balls-*`, `balltree-2balls-*` and `kdtree-2balls-*`:

- `scalar-pivot-reuse,no-smooth-pivot` opts into hierarchical moment reuse.
- `CBALLS_SCALAR_PIVOT_TOL` is a phase budget in radians, in `[0,3]`.
- `CBALLS_SCALAR_BIN_THETA` is an internal bin allowance in bin widths, in `[0,1]`.
- `dual-node-bin-theta` is a separate native options token. Drivers with
  `--dual-node-bin-theta` add that token; benchmark scripts pass it through
  `--ctreeballs-extra-options`.
- Runtime `theta` (`--tree-theta` or `--ctreeballs-theta`) is an independent
  geometric acceptance control. Positive theta is required for scalar reuse.
- `nsmooth` affects leaf capacity/work according to the engine. Disabling pivot
  smoothing does not force singleton leaves. Sweep it at fixed scientific input.

The shell prefix
`env CBALLS_SCALAR_PIVOT_TOL=3 CBALLS_SCALAR_BIN_THETA=0.5`
is an aggressive candidate to measure, **not a 5% or 10% error guarantee**.
Reuse can be slower on some catalogs. Periodic, exact, smoothed or 2PCF-only
paths do not gain scalar 3PCF reuse. Inspect the effective run metadata.

For an unsmoothed scalar two-ball body reference use `theta=0` and native options
`no-one-ball,no-two-balls,no-smooth-pivot`, omitting `scalar-pivot-reuse`.
For spherical shear use `theta=0` and `no-one-ball,no-smooth-pivot`;
`no-two-balls` alone does not make the shear 3PCF neighbor scan exact.
Shear uses its own `CBALLS_SHEAR_*` controls. Ly-alpha uses its own forest/LOS
and kernel controls. Private `*_dual_node.py` peer-backend slops are separate
from native cTreeBalls controls. Exact traversal still retains the specified
bin widths, finite multipole cutoff and sampled catalog.

## Timing and memory

Measure repeated **wall times** on an otherwise idle machine. Record OpenMP
threads, MPI ranks, binding and tree/cache reuse. Report setup/loading, Python
Run and native MainLoop scopes separately. CPU time includes work across
threads; under MPI use maximum participating-rank wall time and summed CPU.
Neither the largest thread count nor enabling every approximation is guaranteed
to be fastest. Independent catalogs can be distributed across processes or
machines, subject to total cores and memory.

Use `--more-options no-out-Hist,dual-node-profile` for a scalar diagnostic run,
or the corresponding `--ctreeballs-extra-options` in the benchmark. Inspect
build, traversal, scratch, multipole products and reductions before optimizing.
Remove `dual-node-profile` for timing trials. Unset `CBALLS_SHEAR_PROFILE` to
remove shear instrumentation: even a value of `0` enables it. `no-out-Hist`
suppresses native histogram files; retain the arrays needed for validation.

HEALPix coordinate chunks, early temporary release, opt-in
`--transfer-catalog-ownership`, and automatically sized compact tree nodes
preserve selected input and approximation controls. They reduce particular
allocations, not necessarily peak RSS or compute time. Native particle/tree
storage and MPI rank replication remain. Transfer is a Python registration
option, not a native `--more-options` token. See `README_memory_controls.md`.
Positive `--max-points` thins a catalog and lower NSIDE changes map resolution;
these are changes to the scientific input, not equivalent memory optimizations.

Use a fresh output directory when changing controls unless the workflow's resume
fingerprint explicitly includes those controls, their environment values and
build identity. A resumed old result is not a newly measured trial.

## Angular edges, centers and plots

Min/max separation controls specify **outer edges**. With 20 logarithmic bins,
8-to-200 arcmin edges give centers approximately 8.67036-to-184.53294 arcmin.
The Sofia figure-1 convention instead targets 20 nominal centers at 8-to-200
arcmin; its spherical chord-bin edges are approximately 7.35028-to-217.68454
arcmin. Read effective saved edges/centers, not just the command line.

Only frontends that advertise a calculation `--binning sofia-fig1` preset can
supply those edges automatically. Their explicit separation limits override
preset defaults; omit the limits to use the preset. Some shear/private drivers
use `--binning`/`--plot-binning` solely to choose plotted pairs. Follow each
parser's help: a plotting-only preset never changes the calculation grid.

In the separate DES flattening workflow, `--selection sofia-fig1` checks the
saved grid. `--selection index-45` retains its pair-index pattern on another
20-bin grid without rebinning or claiming the Sofia angular centers. A failed
plot command can leave older files in an existing directory; use a fresh output
directory and inspect the saved sample/region metadata. Multipole flags can
mean maximum order or a count: native `mChebyshev=7`/driver `--multipoles 7`
retains orders 0 through 7; the DES flattening plotter's `--multipoles 8` selects
those eight orders.

Kappa benchmark plots default on (`--plots`/`--no-plots`). Shear addon benchmark
plots also default on; omit `--no-plots`. Plotting existing coefficients does
not recompute an exact reference or measure an approximation error by itself.
