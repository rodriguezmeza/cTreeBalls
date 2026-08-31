# balltree-2balls-omp_3pcf

This 3PCF-only addon implements TreeCorr's production `LogMultipole` idea for
scalar fields. For every accepted pivot ball it scans a second ball tree and
accumulates radial complex moments of the neighboring field. Products of
those moments generate all `(r1,r2,m)` bins, changing the expensive part of
the calculation from explicit triples to pivot-neighbor pairs. Same-bin
second moments remove `q == r` exactly.

The production tree scan uses TreeCorr's `Log` slop contract: the sum of pivot
and neighbor radii must fit within `theta` times the local radial-bin width and
satisfy the angular phase bound `theta*pi/(2*mChebyshev+1)`. Controlled leakage
at radial-bin boundaries is part of this approximation. Large-radius bins are
completed on coarse pivots; their moments are inherited while only unresolved
small-radius bins descend the pivot tree. `no-two-balls` disables node
acceptance and gives exact body-level pivot-neighbor moments.

In 3D, the angular bound uses projected tangent bearings, including pivot
motion. Inherited moments and their second moments are transported into
each descendant pivot's tangent basis before mixing with newly scanned bins.
Keep observer-centered positions and chord bins. See
[the scalar contract](../../docs/3pcf.rst) for mode handedness, weighting,
undefined bearings, and distinct-neighbor normalization.

`nsmooth` is the default leaf capacity. Add `treecorr-singleton-leaves` for
the most aggressive node acceptance at the cost of a larger tree. The former
genuine triple-node traversal remains available for validation with
`treecorr-direct-triples`; it is not the performance path.

Use `theta=1` for the TreeCorr-style production setting. Reduce `theta` for a
tighter approximation; the `theta -> 0` limit converges to the exact
`no-two-balls` result.

Build with `BALLTREE2BALLSOMP3PCFON=1 TPCFON=1` and select it with:

```text
search=balltree-2balls-omp_3pcf
options=pos-and-convergence,KKKCorrelation,only-3pcf
```

OpenMP works on a fixed pivot-node frontier. Every task owns private moment
and histogram arrays, and task results are reduced in fixed order so changing
the worker count does not change the output ordering.

Histogram normalization is selected at runtime. By default each multipole is
divided by the distinct-triplet count, or by the triplet-weight sum when
`weights-norm` is active. Add `no-normalize-HistZeta` to retain raw triplet
sums. This runtime choice applies even when the build has `NONORMHISTON=1`.
