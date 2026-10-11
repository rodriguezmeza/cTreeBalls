# Fast defaults and measured approximation penalties

Fast traversal is the default. A 5% relative L2 error is a comparison target,
not a universal bound enforced by a geometric tolerance. Runs retain their
chosen controls; supplying an exact reference makes their penalty measurable.

| Control | Default | Meaning |
| --- | ---: | --- |
| Runtime `theta`, Python `--tree-theta` | 1 | Approximate angular tree traversal |
| Make `THETA` | 1 | Separate legacy scan criterion |
| `CBALLS_SHEAR_PIVOT_TOL` | 3 | Shear reuse phase budget, radians |
| `CBALLS_SHEAR_BIN_THETA` | 0.5 | Shear reuse internal bin allowance, bin widths |
| `CBALLS_SCALAR_PIVOT_TOL` | 0.1 | Scalar reuse phase budget, radians |
| `CBALLS_SCALAR_BIN_THETA` | 0 | Scalar reuse requires full radial-bin containment |
| `nsmooth`, Python `--nsmooth` | 8 | Leaf capacity on applicable shear/kappa paths; no automatic smoothing |
| `smooth-pivot` | off | Explicit option changes the unsmoothed estimator |
| `lya2RpSlop`, `lya2RtSlop` | 0 | Exact parallel/transverse bin membership |
| `lya3MuSlop`, `lya3RadialSlop`, `lya3PolarSlop` | 0 | Exact triplet bin membership |
| `lyaPivotRadius` | 0 | No approximate same-forest pivot grouping |

The three angular trees share these controls in OpenMP and MPI. The native
CLI, Cython parser, current benchmark drivers, private dual node comparison
drivers and portable DES launchers use opening 1. Reuse remains explicit:
select `shear-pivot-reuse` or `scalar-pivot-reuse`. A zero phase budget disables
reuse, but ordinary aggregation can still occur at positive opening. Lyα
retains accelerated certified cell aggregation at zero slop; positive slops
remain explicit choices and use the same measured-error reporting.
`nsmooth` remains a leaf-capacity control on applicable paths and does not
silently enable smoothing.

Make settings in `addons/Makefile_addons_settings` use `?=`, so they can be
overridden before compilation. Runtime environment values override their
corresponding compiled reuse controls. They do not overwrite runtime `theta`.
Restart Python after rebuilding the extension. Existing shell overrides and
explicit command arguments take precedence over these defaults.

## Scope of the controls

The two `CBALLS_SHEAR_*` tolerances affect opt-in spherical shear 3PCF reuse
only. They do not change kappa or shear 2PCF acceptance. Kappa retains its
separate `CBALLS_SCALAR_PIVOT_TOL=0.1` and `CBALLS_SCALAR_BIN_THETA=0` defaults.
The active public, private dual node-comparison and portable DES kappa/shear
drivers now use `nsmooth=8`, matching native cballs. Explicit arguments win.

Leaf capacity changes tree size and the balance between node and body work.
It can also change which approximate groups are accepted, so its effect on
speed and error is not monotonic. Some scalar 2PCF policies choose adaptive
capacities; the shear octree retains individual-body leaves. In the exact
unsmoothed limit a leaf-capacity change should preserve the estimator up to
floating-point summation. No percentage error is implied by nsmooth=8.

## Measuring the penalty

Use `relative_L2 = norm(fast - reference) / norm(reference)` separately for
each observable. The drivers print percentages and save
`ENGINE/qualification/approximation-penalty.json`. They report raw numerators,
windows, shear pair weights, normalized and corrected results where the estimator exports them,
and the four shear components separately. A zero reference norm has undefined
relative error; the absolute error is retained. Different finite-output masks,
infinities and failed corrections are visible, rather than counted as zero
error. A small L2 error does not bound every individual bin.

First produce exact results with the same input selection, weights, bins,
statistics and multipoles. For angular kappa and shear, use
`--tree-theta 0 --no-smooth-pivot --more-options no-one-ball,no-two-balls`.
Then run the desired fast settings and append:

```sh
--accuracy-reference results/exact --accuracy-report-only --qualification-rtol 0.05
```

The reference root must contain `ENGINE/qualification` for every selected
native engine. For example, run the same `shear_corr_all_engines.py` input
arguments twice, changing only the following controls and output directory:

```sh
# Exact reference
--tree-theta 0 --no-smooth-pivot --more-options no-one-ball,no-two-balls --outdir results/exact
# Fast reuse with measured errors; retains results above the target
--tree-theta 1 --no-smooth-pivot --more-options shear-pivot-reuse --outdir results/fast \
  --accuracy-reference results/exact --accuracy-report-only
```

Without a reference the penalty is `UNMEASURED`. Computing an exact reference
costs time; it is never silently launched by the fast run. Reusing saved exact
results avoids paying that cost on every tuning trial. A subsample comparison
measures that subsample only, not the error of a subsequent full-catalog run.
For Lyα, use the same direct estimator with zero slops and pivot radius for
the reference. External dual node adapters retain their separate comparisons;
these native result packets do not measure dual node error.

`--accuracy-report-only` retains outputs even if the target is exceeded.
Without it, a failed explicitly requested numerical qualification raises an
error after saving the evidence. The older qualification status uses
`L2(error) <= atol + rtol*L2(reference)` and exact integer counts; the new
`relative_target_met` reports the strict relative target independently of
`atol`. Thus a very weak signal can satisfy the absolute criterion yet exceed
5% relative error. Invalid references or mismatched catalogues still fail.

Single-rank comparisons of the same engine and thread count also report the
native MainLoop wall-time ratio when timings are available. This excludes
validation and I/O. MPI benchmark timing uses maximum rank wall time and summed
rank CPU time; rank-zero timing alone is not presented as job speedup.

## Measured DES example, 2026-10-10

One-million-object DES r002/bin1/region1 subsample, spherical ball-tree shear
3PCF, m=7, 20 Sofia-figure-1 bins (nominal centers 8–200 arcmin;
outer angular edges 7.35028–217.68454 arcmin), 16 threads, no smoothing.
These are single MainLoop measurements on this host, not repeated medians.
These historical trials used nsmooth=1; they do not measure the current
nsmooth=8 default for the binary trees.

| Phase / bin allowance | Wall seconds | Speedup vs exact | Raw Upsilon X L2 | Corrected Gamma X L2 |
| --- | ---: | ---: | ---: | ---: |
| Exact | 73.45 | 1.00 | 0 | 0 |
| **3 / 0.05, tighter reuse trial** | **45.73** | **1.61** | **0.764%** | **0.283%** |
| 3 / 0.1 | 43.32 | 1.70 | 1.003% | 0.313% |
| 3 / 0.5 | 23.98 | 3.06 | 2.816% | 2.172% |
| 1 / 0.05 | 62.68 | 1.17 | 0.499% | 0.049% |

For the tighter 3 / 0.05 trial, the worst individual shear-component errors
were 0.829% raw and 0.362% corrected. These are measurements on this catalogue,
not guaranteed errors on other catalogues. On a 100,000-object sample, the
aggressive 3 / 0.5 profile instead produced raw errors of 17.4–19.3% and
corrected errors of 5.1–10.7% across the three trees. Small samples can also
lose the speed benefit because setup and histogram work dominate.

Scalar signed fields and nearly singular window corrections can amplify tiny
raw errors into large corrected errors. A previous clustered scalar test at
opening 0.025, phase 0.1 and zero bin allowance had about 0.065% raw error but
104% corrected error with ordinary traversal; reuse made the corrected error
even larger. These settings are not a general error bound. Inspect corrected
results and finite-output support as well as raw moments.

The tighter shear reuse profile (phase 3, bin allowance 0.05) on that
100,000-object sample had pooled raw errors 3.33%, 3.52%, and 3.59% for
octree, KD-tree, and ball-tree respectively; corrected errors were 0.70%,
1.02%, and 0.93%. The largest individual raw-component error was 4.25%.
At this smaller size reuse gave no consistent speedup.

For an actual KAPPA field on the same 100,000-object DES selection and binning,
with nsmooth=1 (a historical trial, not the current leaf default),
ordinary opening-1 scalar traversal gave raw 3PCF errors 23.41%, 19.69%, and
27.53%, with corrected errors 19.18%, 15.39%, and 10.81% (octree, KD-tree,
ball-tree). Shipped scalar reuse controls reduced raw errors to
0.01343%, 0.00746%, and 0.00950%, and corrected errors below 0.00053%, but
were slower than exact traversal on this sample. These are distinct from
the earlier shear-component-as-scalar proxy measurements.

Shear 2PCF needs a separate check: at opening 1 on the 100,000-object sample,
the octree had 1.56% xi-plus and 3.16% xi-minus error, while its pair-weight
histogram differed by 12.91%. KD-tree and ball-tree xi-minus errors were about
0.144%, with unchanged pair weights. The result packet includes `xi_weight`
so normalization changes are visible. Older 2PCF shear reference packets that
lack this array must be regenerated.

## What this comparison covers

The error is measured against the same retained estimator, bins and multipoles.
It does not measure finite-multipole truncation, sample selection, bin width,
or observational/statistical errors. The Lyα multipole reconstruction requires
its own convergence check or comparison with direct angular bins; its finite
series is not accepted as an exact hard-bin reference. Angular kappa and shear
reconstructions likewise require multipole-convergence checks.

`tests/python/test_search_defaults.py` checks implicit versus explicit fast
controls for the three kappa and shear trees and preserves the zero-opening
reference check. `test_approximation_penalty.py` checks true relative errors,
weak components, undefined errors, and report-only output retention.

## Required exact-reference science campaign

Fast defaults are unchanged. The active gate now runs
`scripts/survey_science_acceptance.py`: signed, clustered, masked and weak-signal
catalogs with exact controls, every returned observable, and selected radial
bins under a pointwise science metric. Execution PASS does not promote rejected
default candidates. The separate `science_acceptance` field remains REJECTED
when any required candidate fails. Qualify the actual production packet with
explicit tolerances and selection; see [survey semantics](physical_survey_semantics.rst).
