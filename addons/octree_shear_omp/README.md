# octree-shear-omp

`octree-shear-omp` computes flat-sky weak-lensing shear two- and three-point
correlation functions in the Porth x-projection (`Gamma^x`). It uses the
production cTreeBalls octree for exact radial pruning and descends accepted
cells to bodies so that spin-2 phase factors are never approximated by a cell
average.

## Estimators

The coordinate system is a right-handed Cartesian tangent plane. Angles are
measured counterclockwise from the positive x axis, and input shear follows
`gamma = gamma1 + i gamma2`. For pair direction `phi`, the weighted two-point
outputs are

```text
xi+ = sum(wa wb gamma_a conj(gamma_b)) / sum(wa wb)
xi- = sum(wa wb gamma_a gamma_b exp(-4 i phi_ab)) / sum(wa wb)
```

The 3PCF implementation evaluates the four natural components
`Gamma_0^x` to `Gamma_3^x` through ring multipoles `G_n` and survey-window
multipoles `W_n`. This is the Porth x-projection with
`zeta1=(phi1+phi2)/2`, `zeta2=phi1`, and `zeta3=phi2`. It is not the centroid
projection; comparisons with centroid-projected measurements must first apply
the projection conversion of Lucas-Porth et al. Terms with the same neighbor
in both rings are removed when the two legs use the same catalog, so all
auto-correlation triplets contain three distinct objects. Finite survey
geometry is corrected by solving the mode-coupling system

```text
C[ell,n] = N[ell-n] / N[0]
C Gamma = Upsilon / N[0]
```

for each radial-bin pair. Singular or empty systems use the project-wide zero
denominator policy and return zero rather than NaN or infinity.

## Build

```bash
make OCTREESHEAROMPON=1 clean
make OCTREESHEAROMPON=1 PYTHON=python3 all
```

For a source installation:

```bash
OCTREESHEAROMPON=1 python3 -m pip install . --no-build-isolation
```

The compiled profile contains `OCTREESHEAROMP`, `THREEPCFSHEAR`, and `TWOPCF`.

## Inputs

The method requires at least two Cartesian coordinates and non-periodic,
flat-sky geometry. A 3D build is supported when the third coordinate is
constant for every object. Every selected catalog must lie on the same
tangent plane; internally flat catalogs at different constant z coordinates
are rejected. Catalog weights must be finite and nonnegative.

One, two, or at least three input catalogs select the tomography tuple as
`Z1,Z1,Z1`, `Z1,Z2,Z2`, or `Z1,Z2,Z3`, respectively, using the first three
entries of `iCatalogs`. The 2PCF always uses the `Z1,Z2` pair. For three or
more input catalogs, `iCatalogs` must still contain one entry per input file,
as required by the core parameter parser; entries after the third do not
participate in this estimator. Startup applies one common translation to all
catalogs before tree construction, preserving every cross-catalog separation.

Python can pass shear and weights without an intermediate file:

```python
model.set_catalog(positions, gamma1=g1, gamma2=g2, weights=weights)
```

The existing `multi-columns-ascii` reader accepts `x y z gamma1 gamma2` in a
3D build with `options=pos-and-shear` and `columns=1,2,3,4,5`; that reader sets
weights to one. FITS/HEALPix readers that populate `Gamma1` and `Gamma2` are
also compatible only after projecting the positions to one flat patch.

Periodic geometry is rejected because wrapping a pair changes its survey
direction and is not the flat-sky window model used by this estimator.

## Outputs

Without `options=no-out-Hist`, the method writes:

- `histShearXi*.txt`: radial center, complex `xi+`, complex `xi-`, pair weight.
- `histShearGammaMultipoles*.txt`: corrected `Gamma_mu,n` and raw `Upsilon_mu,n`.
- `histShearGamma*.txt`: angularly reconstructed `Gamma_0` to `Gamma_3`.

The preferred projection-explicit Cython getters are
`getShearUpsilonXMultipoles()`, `getShearGammaXMultipoles()`, and
`getShearGammaX()`. The older generic names remain compatibility aliases.
The other getters are `getShearXiPlus()`, `getShearXiMinus()`,
`getShearXiWeight()`, and `getShearWindowMultipoles()`.

## Regression test

```bash
make OCTREESHEAROMPON=1 PYTHON=python3 cyballs
tests/make_tests/run_test_shear_octree_omp
```

The test compares every estimator stage with an independent NumPy oracle,
checks factorized rings against explicit distinct-object triplets, verifies
active spin-2 rotation invariance, tests `Z1,Z2,Z3` tomography and common-plane
rejection, checks one-thread versus multi-thread bitwise equality, and
exercises recoverable invalid-geometry paths.
