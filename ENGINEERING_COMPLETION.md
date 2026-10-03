# Input, MPI, numerical qualification and scaling contracts

This document describes the repaired active Makefile profile. It does not enable additional addons. Reproduce validation with `python scripts/active_release_gate.py --output <new-directory> --jobs 4 --mpi-command "mpiexec"`; provide launcher options appropriate to the host. The final machine-readable gate and measurement evidence is kept separately from the source tree.

## Catalog validation

All active input routes converge on the shared catalog contract before geometry/tree construction: nonempty allocated data, finite coordinates/scalars/weights, binary masks, finite shear where consumed, and finite box dimensions. Raw selected fields are checked before constant-field overrides, coordinate transformations and map selection. Header-only and mask-only operations have explicit separate paths.

Takahashi input now checks every header/record/map read, count and allocation, validates the HEALPix count, rejects every nonfinite map component, and propagates empty/invalid selection errors. Its historical native-endian/native-long layout remains unchanged. FITS fixed/vector shapes are checked before indexing; radial samples use the correct angle row; companion masks are binary. Optional XYZ LOS identifiers require a scalar numeric FITS column without nulls. Integer storage preserves values above 2**53; floating storage must contain finite, integral values within the signed 64-bit range and is never silently rounded.

`test_input_failure_contracts.py` exercises every truncation byte of the small binary/Takahashi layouts, each field's NaN/infinity cases, malformed headers, ASCII tokens, FITS null/shape/LOS cases and same-object recovery. This is systematic finite fixture coverage, not proof against every possible malformed file or a general fuzz campaign.

## Collective MPI application-error contract

Participating ranks must call an MPI method's `Run` in matching order on a healthy communicator. Errors in CLI/Python parsing, startup setup, catalog loading, geometry/bin preparation, tree/histogram/scratch construction, workers, reductions, output and result publication propagate at common phase boundaries. Control signatures reject incompatible reductions, including different histogram shapes. Replicated catalog fingerprints reject divergent input content. Cached Python calls also participate in entry/state checks. Destructors/cleanup are noncollective; serial methods may run on root alone after MPI work. Neither cleanup nor an individual run finalizes borrowed MPI.

`test_mpi_boundary_contracts.py` enumerates every boundary reached by all active MPI methods in the fixture matrix, injects one failure on either rank, checks identical errors and freed catalogs, then successfully reruns the same object. Ordinary, scalar compatibility and edge-correction routes run with both library-owned and externally initialized MPI. Real parameter, unread-key, conflicting input, different-control, different-catalog, invalid-stage and divergent cached-state cases supplement injection. Existing native CLI and malformed-reader MPI tests remain release gates.

“Exhaustive” here means every **reached named application boundary × either rank × both MPI ownership modes** for the active route matrix. It excludes killed processes, damaged MPI communicators/transports, a rank that never calls the collective API, mixed serial/MPI calls between ranks, failures inside MPI itself and arbitrary hardware/resource exhaustion. No finite regression establishes those as recoverable contracts.

## Core and compatibility numerical corrections

Core `octree-sincos-omp` retains its historical estimator: an unweighted pivot field average of weighted neighbor means, including neighbor self-products. It is not the modern raw distinct-neighbor sum. Smoothing now preserves pivot and neighbor group contributions without multiplying neighbor counts by pivot multiplicity or reapplying a global skipped-pivot factor. Linear-bin and accepted-cell pair denominators include the missing contributions.

The native octree `legacy-one-ball` path uses body pivots, conservative radial-bin acceptance for neighbor cells, weighted pair denominators, and sum-of-squares cell moments for diagonal/self subtraction. A cell's squared aggregate is not the sum of its members' squares. Body-pivot scans remove a geometry/frame/binning ambiguity in the old cell-pivot approximation. This correctness choice can increase compatibility runtime; modern two-ball traversal is preserved.

Independent oracle tests cover unequal tightly clustered groups, signed rapidly varying fields, zero weights, linear/log bins, several multipole orders, exact controls and smoothing that actually combines points. Core and raw-GGG references remain separate. Accuracy acceptance now gates both repaired routes rather than classifying them as diagnostic-only. Approximation tolerances remain catalog- and observable-specific.

## User and driver qualification

`model.getResults()` returns owned arrays plus provenance and a qualification record. The default is **UNMEASURED**, including native CLI metadata: a successful run is not evidence of approximation accuracy. `model.qualifyAgainst(exact_model_or_packet, rtol=.02, atol=1e-10)` attaches evidence to subsequent metadata/results.

Qualification requires matching catalog-content fingerprints, estimator, geometry, bins, weights, masks, coordinate conventions and semantic options. References must disable aggregation and smoothing, or use an exact supported forest method. Each actual returned observable is checked separately with `L2(error) <= atol + rtol*L2(reference)`, matching finite-bin masks and no infinities; integer count arrays must match exactly. Reports expose weak-reference bins, maximum absolute bin error and pointwise tolerance exceedances. Global L2 acceptance does not promise every bin is within a relative tolerance. Empty/undefined floating bins may match as NaNs; an entirely undefined observable cannot qualify.

The convergence, shear and Ly-alpha all-engines drivers write `ENGINE/qualification/observables.npz` and `qualification-result.json`, with an array-file digest and the precise provenance. They display the status and accept `--qualification-reference <exact-driver-output-root>`, `--qualification-rtol` and `--qualification-atol`. A requested rejected qualification retains its evidence and fails the driver run. The exact reference must use the same input and requested products; incompatible estimators are rejected. References are explicit because recomputing large exact catalogs can be expensive.

For scalar/shear reference runs use theta zero, `no-one-ball`, `no-two-balls` where supported, and `no-smooth-pivot`. Forest references must disable all approximation options; approximate anisotropic-multipole references are rejected. The active physical Legendre engines traverse original bodies and can serve as exact references for their own estimator. Keep a different output directory for each run. The generic saved-result workflow is:

```sh
python scripts/qualify_results.py --candidate candidate/ENGINE/qualification \
  --reference exact/ENGINE/qualification --output qualification.json
```

The utility exits nonzero on rejected qualification. `save_result_packet` and `load_result_packet` preserve/verify copied arrays; no qualification transfers to a different catalog, observable, smoothing radius or opening tolerance.

## Representative scalar scaling

```sh
python scripts/benchmark_representative.py --output measurements \
  --counts 2048 8192 32768 --thread-counts 1 2 4 --mpi-command "mpiexec"
```

This is a reproducible scalar suite over nested fixed-seed full-sky and clustered geometry, signed fields, nonuniform weights and masks. It retains full-catalog exact KD references at every size, actual observable packets for every measured configuration, cold/warm calls, three modern scalar trees, 1/2/4 OpenMP threads and two MPI ranks. Each candidate is qualified against the full matching catalog rather than a small subsample. Timing is maximum per-rank wall time; CPU is summed across ranks. RSS is per-process high-water memory including Python/imports, reported as the maximum rank value, not summed node memory. Runs are sequential fresh workers on one host; cache reuse occurs only between cold/warm calls.

The JSON distinguishes successful measurement from numerical qualification: rejected approximations remain visible evidence. These are single cold/warm samples without confidence intervals, dedicated multinode/interconnect testing or a universal speedup claim. They extend representative scalar scaling, not the scaling coverage of every forest/shear estimator. Use the supplied suite with larger sizes and deployment hardware as needed.
