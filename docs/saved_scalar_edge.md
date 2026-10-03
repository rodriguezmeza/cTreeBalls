# Saved scalar angular-window correction (format v1)

`edge-corrections-from-files` solves the same complex, truncated scalar Fourier
window system as the modern in-run two-ball correction. It takes exactly two
comma-separated **prefixes**, in this order: signal, window. Prefixes may be the
same, reside in different directories, or contain spaces. Quote the entire
`in=...` argument in a shell. Commas, newlines and leading/trailing whitespace
are not supported inside a prefix. No current-directory `rbins` file is used.

## Generate compatible exports

Run a supported scalar two-ball engine with
`options=KKKCorrelation,edge-corrections,no-normalize-HistZeta,...` and histogram
output enabled. The KD-tree, PCA ball-tree and octree two-ball implementations
share this exporter, including their MPI root publication. Use `no-smooth-pivot` and `theta=0` when exact raw-catalog moments are
required. Approximate traversal settings remain approximate after correction.
Legacy compatibility traversals do not necessarily emit this format; use the
modern traversal when a version-1 manifest is not produced.

For the default histogram name, the prefix is `OUTPUT_DIRECTORY/histZetaM`.
The exporter writes additional `edge_*` signal matrices at 17-digit precision,
before any subsequent normalization/output processing. It retains the existing
ordinary signal, window and corrected products.

For example, after producing compatible signal and random-window exports:

```sh
./cballs \
  'in=/path/to/signal/histZetaM,/path/to/window/histZetaM' \
  rootDir=/path/to/corrected \
  options=edge-corrections-from-files numberThreads=1
```

The two prefixes must contain version-1 manifests. `infileformat`, `iCatalogs`,
and catalog-loading options are irrelevant to this preprocessing task.
`sizeHistN` and `mChebyshev` are read from the signal manifest, not inferred by
probing filenames. The window must supply every mode through twice the signal
maximum. Extra window modes are allowed and are not needed by the smaller solve.
`no-out-Hist` and `full-sky` are rejected for saved correction: this task produces
files and always solves the supplied window system.

Python may use `cballs.set(infile='signal_prefix,window_prefix',
options='edge-corrections-from-files', rootDir='corrected')` followed by `Run()`.
Errors raise normally and the same object can be reconfigured and reused.
Successful preprocessing writes files and cleans its native state; it does not
publish catalog-based `run_settings` or expose native histogram getters.

## Scientific contract

For modes `-M <= ell,n <= M`, solve

```text
sum_n [W_(ell-n) / W_0] zeta_n = S_ell / W_0
S_m = cos_m + sin_m + i (sincos_m - cossin_m)
S_-m = conjugate(S_m); W_-m = conjugate(W_m)
```

The raw signal and window use ordered distinct-neighbor triples and the same
phase convention. The window is the complex weighted count moment, with orders
`0..2M`. Its zeroth mode is real. The signal's zeroth mode is also real.
All four signal products are required; no imaginary term or high window mode is
assumed zero when its file is missing. This is a finite-mode deconvolution, not
a guarantee against leakage from true signal modes beyond the chosen M.

The shared solver scales by positive W0 and uses complex LU with partial
pivoting and the existing rejection threshold `128 * DBL_EPSILON * (2M+1)`.
This preserves the in-run numerical acceptance policy. It does not regularize
ill-conditioned systems. The pivot ratio is a diagnostic proxy, **not** a
condition number, covariance estimate, or accuracy guarantee.

Independent signal/window catalogs may differ in their objects, weights and
sizes. Their geometry, periodic box where applicable, bin edges, lower cutoff,
logarithmic convention, normalization and phase convention must agree. The
loader checks these properties. It cannot establish that an independently
chosen random catalog models the intended survey selection, or supply any
missing relative count/density normalization: the supplied raw S and W must
already represent the intended equation. No automatic random-catalog rescaling
is performed. These files do not contain hashes of the originating catalogs.

## Files and manifest

For prefix `P`, each numeric file has exactly B*B whitespace-separated values
in row-major order. Exporters write B rows of B values. Indices on filenames
are **mode+1**:

| Required payload | Files |
| --- | --- |
| Signal modes 0..M | `P_edge_cos_1.txt` ... `P_edge_cos_(M+1).txt` |
| Other signal components | Corresponding `edge_sin`, `edge_sincos`, `edge_cossin` files |
| Window modes 0..2M | `P_window_Re_1.txt` ... `P_window_Re_(2M+1).txt` |
| Imaginary window | Corresponding `window_Im` files |

`P_edge_manifest.txt` has these tokens in this exact order (example B=2, M=1):

```text
CBALLS_SCALAR_EDGE 1
bins 2
signal_mmax 1
window_mmax 2
dimensions 3
periodic 0
box 2 2 2
geometry observer-tangent-fourier
normalization raw-ordered-distinct
phase cc+ss+i(sc-cs)
logarithmic 0
lower_cutoff 0.02
edges 0.02 0.76 1.5
complete
```

For two dimensions, `box` has two values and geometry is `planar-fourier`.
Grid comparisons are exact after reading the 17-digit decimal representation.
Nonperiodic box extents may differ. The logarithmic zero-cutoff first-bin
extension is represented by the native export's actual lower edge and separate
zero lower cutoff. All edges must be finite and strictly increasing.

The first prefix supplies the signal files; the second supplies the window
files. The first prefix's window matrices are not silently substituted. A
manifest on each prefix defines its own conventions and available orders.
Malformed numbers, missing required files, extra/missing values, unsupported
versions, incompatible grids or nonreal finite monopoles fail explicitly.
Explicit `nan`/`inf` numerical payloads produce nonfinite **bin** statuses;
representable subnormals are retained, while textual overflow or underflow
to zero is rejected as malformed input.

## Results, invalid bins, and write failures

Outputs are `histZetaM_EE_(m+1).txt` (real) and
`histZetaM_EE_Im_(m+1).txt` (imaginary), respecting a configured histogram prefix.
`histZetaM_window_diagnostics.txt` has one row per radial-bin pair:
`bin1 bin2 status window_monopole pivot_ratio`.

| Status | Meaning | Corrected values |
| --- | --- | --- |
| 1 | Valid solve, including a measured zero | Finite |
| 2 | Empty/nonpositive window monopole | NaN |
| 3 | Singular or numerically rejected pivot | NaN |
| 4 | Nonfinite input or solve | NaN |

Nonfinite inputs take precedence over empty support. The completion file
`histZetaM_edge_result.txt` is written only after all corrected matrices and
diagnostics succeed. An export similarly writes its `edge_manifest.txt` after
its raw matrices and diagnostics. Existing completion markers are invalidated
before payload overwrites. Input validation completes before result writing.
A failed write can leave partial payload files; require the corresponding
completion marker and a successful process exit before using a newly produced
set. Markers are not cryptographic integrity checks or directory-wide atomic
transactions. Use a fresh output directory for each scientific run.

## Migration and validation

Old unversioned `cos/sin/cos_N/sin_N` sets are rejected. Regenerate exports from
the catalogs with the updated modern scalar engine. Simply adding a manifest
to an old set is unsafe: cross terms, high window modes, normalization and
precision may be missing or incompatible. The historical `matrixClm` helper
remains solely for older addon source compatibility; this preprocessing route
does not call it.

Run `python -m pytest -q tests/python/test_saved_scalar_edge.py` after building
both the native program and Python extension. The regression suite includes
the W2=0.4 analytic case (zeta1=1/7), independent-prefix rescaling, complex
NumPy solves, support statuses, malformed/incomplete inputs, writer failure,
Python recovery, and native round trips for three engines and both linear and
logarithmic bins. Existing in-run brute-force oracles independently check the
shared numerical path.
