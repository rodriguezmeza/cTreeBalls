# Runtime stabilization regressions

Run with the active 3D, double-precision, CLASSLIB, scalar-tree, KD-box,
neighbor-box, IOLIB and CFITSIO profile:

```sh
make test-runtime-stabilization
make test-runtime-stabilization-cython PYTHON=python3
```

The second target builds the local extension in place; it does not install it
into Python's global environment. Its Python environment needs NumPy and Cython.
Both targets run the native analytic and constructor tests. The first runs the
CLI cases and skips the five extension-dependent test groups.

`counts_65536.json` records the small deterministic reproduction used during the
active-addon audit. Generate the catalog with NumPy's default PCG64 generator,
seed 811, 65,536 points uniform in [-1, 1)^3. The recorded unordered pair counts
are [41037, 224782, 1222140, 6501418] for four logarithmic bins from 0.02 to 0.2.
The old optional count correlation returned infinity in all four bins because
the signed 32-bit product N*N wrapped to zero on the tested build. The regression
checks those same pair counts and the independently calculated finite result
2*DD*8/(N^2*shell_volume)-1. The catalog is generated at test time, not stored as a
large binary fixture.

Other permanent cases in `test_runtime_stabilization.c` and
`make_tests/test_runtime_stabilization.py` cover:

- Analytic normalization at N=46,340, 46,341, 65,536 and above INT_MAX, without
  allocating large catalogs. The last count tests the normalization interface,
  not end-to-end tree capacity.
- Exact linear and logarithmic shells, positive and zero lower cutoffs, empty
  bins, invalid domains/counts and count-only operation. The zero-cutoff log
  first bin intentionally reflects the existing producers' truncation toward
  zero: its lower edge is one log step below the nominal grid.
- Independent all-pairs counts and shell-volume oracles for default octree, KD,
  ball and KD-box methods, and legacy KD/ball paths. Legacy octree uses the B4
  pivot approximation even with `no-one-ball`; its shell normalization is
  checked against the DD it actually publishes. This is not a claim that its
  approximate pair counts equal brute force. Linear/logarithmic runs with and
  without three-point multipoles must publish identical count results.
- KD context, pointer array, node array and body-order allocation failures;
  partial-tree destruction; repeated construction; recoverable Cython errors.
- Literal spaces, quotes and command-substitution characters in directory
  names; repeated/deep directory creation; existing regular files at leaf,
  intermediate and log-directory positions.
- Failed publication in all seven active catalog formats and in the owning
  legacy octree and neighbor-box histogram drivers, followed by Cython reuse.

The earlier Gadget/FITS fixtures remain in `../io_stabilization/` and are tested
by `make test-io-stabilization-cython`. MPI multi-rank behavior and other build
profiles require their own validation.
Public startup currently rejects logarithmic `rminHist=0`; that historical
normalizer branch is covered only through the direct native interface.
