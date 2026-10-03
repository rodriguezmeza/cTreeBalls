# Exact angular bins alongside anisotropic Ly-alpha multipoles

`lya-anisotropic-multipole-3pcf-omp` supports `lya3MuMode=1` to compute exact
opening-angle bins alongside the raw Legendre moments. This eliminates angular
truncation error in the **additional hard-bin product**, independently of Lmax.
The default `lya3MuMode=0` retains the existing cost and approximate reconstruction.

## Why a finite series cannot guarantee exact hard bins

The raw moments describe the numerator and denominator measures separately:
`T_l = sum_triangles weight * P_l(mu)`, with the numerator also including the
three fields. Reconstructing a top-hat bin uses
`sum_(l=0..L) [(2l+1)/2 * integral_bin P_l(mu) dmu] T_l`.
This is a polynomial approximation to a discontinuous window. Angular atoms,
bin edges, sparse support, and cancellations in the numerator can produce
leakage, signed reconstructed denominators, and large ratio errors.

There is no universal finite-L postprocessing fix. For example, the positive
five-node and six-node Gauss-Legendre measures have identical moments through
L=8, but different masses in the hard bin [0,0.5). A function of those moments
alone cannot return both exact answers. This counterexample follows from
Gaussian quadrature's polynomial exactness through degree 2n-1 (see
[NIST DLMF 3.5](https://dlmf.nist.gov/3.5)); the regression constructs both
measures independently with NumPy. Increasing L or damping ringing may help a
particular catalog, but neither establishes a general per-bin accuracy bound.
Never clip signed reconstructed denominators and call the result exact.

## Usage and output contracts

Add these parameters to the usual native catalog command:

```text
search=lya-anisotropic-multipole-3pcf-omp
lya3MuMode=1
lya3Kernel=0
lya3LMax=8
```

This is a runtime parameter, also accepted by `cyballs.set()`. Rebuild both the
native binary and the Python extension after changing command fields.

| Product with the default histogram basename | Meaning |
|---|---|
| `histZetaM_lya_multipoles.txt` | Raw signed Legendre numerator/denominator moments, exact up to floating arithmetic |
| `histZetaM_lya5d_multipole.txt` | Existing approximate finite-L reconstruction; unchanged meaning, retained as a diagnostic |
| `histZetaM_lya5d.txt` | Additional exact hard-bin numerator, denominator and ratio, only in mode 1 |

The exact file uses the ordinary five-dimensional 13-column layout and its
existing zero-denominator convention (`zeta=0`). Empty bins are omitted unless
`options=lya-output-empty-bins` is included. The approximate file retains all
bins, signed raw sums, and NaN ratios for nonpositive denominators. Its finite-L
diagnostic can still fail accuracy tests when the exact companion passes.
Different product names prevent existing readers from mistaking it for an
exact histogram. Output directory contents from earlier runs are not deleted;
use a fresh directory and the current run metadata to identify available products.

With `options=no-out-Hist`, `getForestResults()['arrays']` returns independent
copies of all four raw arrays in mode 1:

- `moments_numerator`, `moments_denominator`: shape `(R,R,T,T,L+1)`.
- `triple_numerator`, `triple_denominator`: shape `(R,R,T,T,MuBins)`.

`run-metadata.json: lya_mu_output` records mode, exact-bin availability, product
suffixes and shared discovery. The older `lya_3pcf.mu_reconstruction_approximate`
flag stays true: it correctly describes the still-approximate finite-L file.
The exact algorithm does not by itself certify sampling error, catalog
systematics, or observational accuracy; general qualification metadata remains
UNMEASURED until an actual reference comparison is supplied.

## Computation and cost

Each pivot discovers and sorts its neighbors once. The existing harmonic
hierarchy computes moments. The same neighbor array then supplies the certified
hard-bin hierarchy: accept products only when conservative opening-cosine bounds
fit inside one bin; otherwise refine and evaluate individual pixel pairs.
Mixed-forest radial moments are reused on suitable sparse frontiers. Radial and
polar assignments, distinct-forest exclusions and original pixel geometry are
preserved. Both leg orders are included, and geometric triangle counts are
published once, including zero-weight triangles. No pivot smoothing or geometry
slop is allowed by the multipole engine.

There is an additional histogram and exact angular computation; this does not
promise the speed of moment-only output. The extra histogram plan is
`R^2*T^2*MuBins * [2*sizeof(REAL) + threads*(2*sizeof(REAL)+sizeof(size_t))]`
bytes, plus segment/mixed-forest scratch. Checked allocation preflight includes
both grids and their live scratch. Workers release both on success or failure.
Fixed publication blocks preserve thread-count reproducibility within a build.
If only hard bins are needed, use an exact `lya-3pcf-omp` or
`lya-los-tree-3pcf-omp` engine directly; computing unused multipoles adds work.
Mode 1 is available only for the anisotropic multipole OpenMP engine; this change
does not introduce an anisotropic multipole MPI engine.

## Calibration and regression

```sh
python tests/python/benchmark_lya_triplet_kernels.py \
  --synthetic --threads 2 --lmax 4 8 16 32 --mu-mode exact \
  --warmups 1 --repeats 3 --require-accepted-multipole \
  --outdir results/lya-multipole-exact-mu
make test-lya-multipole-hierarchy
```

Replace `--synthetic` by `--ascii pixels.txt`, `--catalog saved.npz`, or
`--fits /path/to/forest.fits` for representative data. The calibration driver
compares the selected mu product to the direct reference and separately records
`reconstruction_comparisons` in exact mode. Acceptance of the exact companion
is never reported as acceptance of its approximate reconstruction. Full process
timing includes both outputs and is not a moment-only kernel benchmark.

Regression coverage includes an independent ordered-triangle oracle, raw moments,
orders 0/1/4/12/32, 1/3/4/20/65 mu bins, boundary and clustered geometry, bent
forests, zero/tiny/dominant weights, signed fields, one/two-forest exclusions,
OpenMP determinism, unchanged discovery counts, independence of exact bins from
Lmax, ownership after cleanup, mode changes, memory failure/recovery, invalid
controls, and the synthetic case where finite-L reconstruction fails.

## Retained local measurements (2026-10-03)

Double precision, 16 OpenMP threads, Lmax=8, 10 polar bins and 20 mu bins,
Rmax=160; 8 radial bins for the long-forest case and 20 for the other cases.
Medians of three runs, each in a fresh Python 3.11 process with one warmup.
Timing covers native MainLoop; catalog registration, histogram serialization and
Python result copies are excluded. These are bounded fixtures, not full-survey
scaling or guarantees for another catalog.

| Catalog | Pixels | Moments only (s) | Moments + exact bins (s) | Exact bins only (s) | Combined / moments |
|---|---:|---:|---:|---:|---:|
| long-2048 | 2048 | 0.3077 | 2.0650 | 1.9311 | 6.71x |
| cluster-4096 | 4096 | 0.3276 | 0.3727 | 0.3768 | 1.14x |
| survey-6991 | 6991 | 0.7339 | 1.4199 | 0.6672 | 1.93x |

Exact-only timing uses `lya-los-tree-3pcf-omp lya3Kernel=5`, whose automatic
hierarchy choices can differ from the combined path. Combined raw moments were
bitwise identical to moment-only output. Exact raw products agreed with the
exact-only engine: worst relative L2 discrepancy 8.44e-16; occupied-bin support
was identical. The clustered case benefits from bulk certified products; the
long-forest case needs substantially more angular work. Sharing discovery does
not guarantee that the combined run is faster than the sum of two separate runs.

On the 96-pixel synthetic calibration with four bins per axis, the exact
companion's maximum relative correlation error was 1.10e-13 over 132 eligible
bins, with zero invalid eligible bins and zero empty-bin leakage, at each of
Lmax=0,4,8,16,32. The finite-L diagnostic failed the 5% criterion at every tested
order. Its L=4 maximum error was 11.44 (about 1144%), with four invalid eligible
bins. This establishes the correction on that fixture, not a universal relative
error bound near vanishing signals.
