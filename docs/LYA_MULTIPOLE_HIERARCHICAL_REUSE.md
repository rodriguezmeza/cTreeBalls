# Anisotropic Legendre 3PCF hierarchy

The optimization applies to `lya-anisotropic-multipole-3pcf-omp`.
It retains the existing octree neighbor discovery and exact pixel geometry.
It does not translate moments between distinct pivots or add an MPI engine.

## Algebra and exclusions

For each pivot and radial/polar bin, one forest supplies two harmonic vectors:
`sum(w * delta * H_lm)` and `sum(w * H_lm)`. Their dot products give
Legendre moments through the addition theorem. Two different forests are
contracted with the pivot's field/weight; both ordered leg permutations are
published, including two contributions when their output bin is the same.
The pivot forest is removed by discovery. All geometric counts include
zero-weight pixels, as before.

The hierarchy first contracts forests within a block. Their summed harmonic
vectors then contract once against previous complete blocks. Every pair of
different forests meets exactly once, at its within-block or between-block
boundary. No total-minus-self subtraction is used, including with a dominant
forest. Radial and polar classifications and the physical cutoff remain exact.

An integer occupancy pass tries blocks of 1, 4, 16 and 64 complete forests.
The default uses the hierarchy only if its harmonic product estimate is below
85% of the ordinary prefix estimate. This is a workload choice, not a numerical
tolerance. Sparse clearing, cached recurrence square roots, precomputed output
strides, sparse output publication and the existing typed neighbor sort also accelerate the prefix path.

## Controls and outputs

For this engine only:

- `lya3Kernel=0`: automatic choice (default).
- `lya3Kernel=1`: original prefix reference implementation.
- `lya3Kernel=2`: force four-forest groups for testing or workload calibration.
- `lya3LMax=0..32`: retained Legendre orders, unchanged.

`histZetaM_lya_multipoles.txt` contains signed raw moments. The companion
`histZetaM_lya5d_multipole.txt` is an **approximate** finite-Lmax top-hat
reconstruction; a nonpositive reconstructed denominator retains a NaN ratio.
Faster accumulation does not establish agreement with the hard-bin estimator.
Use `tests/python/benchmark_lya_triplet_kernels.py` for that qualification.

`run-metadata.json: lya_multipole_reuse` records hierarchical/prefix pivot
counts, completed forest blocks, actual products and reference-prefix cost.
Kernel 1 retains the original implementation and leaves these reuse counters
zero. Different kernels can differ by floating-point rounding; a fixed kernel,
build, input order and pivot block size is deterministic across OpenMP threads.

## Resources and regression

Let B = radial bins times polar bins and H = (Lmax+1)^2. Moment scratch is
`(6*B*H + 3*H)*sizeof(REAL) + B*(5*sizeof(size_t)+2) + B*B` per worker, excluding
histograms and neighbor storage. The original moment workspace used about
`4*B*H*sizeof(REAL)`. The new storage is bounded independently of forest count.
Checked memory preflight includes live neighbors, all workers and the histogram
plan. Worker cleanup releases this storage on success or failure.

Run:

```sh
make test-lya-multipole-hierarchy
```

Tests compare raw signed moments to an independent ordered-triangle Legendre
recurrence, cover zero/dominant/tiny weights, signed fields, one/two forests,
radial/cutoff boundaries, bent and clustered geometries, rotations, orders
through 32, automatic/forced/reference modes and thread determinism. Existing
window reconstruction and resource failure tests remain required.

A reproducible binary-version comparison uses the existing native MainLoop
benchmark (copy the old Python extension into a separate directory first):

```sh
python scripts/benchmark_lya_hierarchy.py \
  --reference-module /path/to/before --candidate-module . \
  --methods lya-anisotropic-multipole-3pcf-omp \
  --geometry forests --forests 128 --pixels 32 --threads 16 \
  --reference-parameters '{"lya2Kernel":0,"lya3Kernel":0,"lya3LMax":8}' \
  --candidate-parameters '{"lya2Kernel":0,"lya3Kernel":0,"lya3LMax":8}' \
  --warmups 1 --repeats 3 --outdir results/lya-multipole-before-after
```

The benchmark alternates versions, retains raw moment arrays and provenance,
and reports MainLoop wall time separately from process CPU time. It excludes
catalog loading/registration and output-file serialization. Dataset, cutoff,
bin occupancy, Lmax and thread count affect the gain; qualify representative
catalogs before selecting nondefault controls.


## Exact opening-angle bins with multipoles

Set `lya3MuMode=1` on `lya-anisotropic-multipole-3pcf-omp` to compute an
additional exact `_lya5d` product using the same neighbor discovery and sort.
`_lya_multipoles` retains raw moments and `_lya5d_multipole` remains an explicitly
approximate diagnostic. Exact bins are independent of Lmax and require extra
histogram memory and angular computation. Default mode 0 is unchanged.
The Python getter returns both moment and exact triple arrays even with
`no-out-Hist`; metadata records their separate meanings. See
[the reconstruction guide](LYA_MULTIPOLE_RECONSTRUCTION.md) for use, resource costs and validation.
