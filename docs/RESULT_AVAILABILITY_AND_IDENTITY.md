# Result availability and qualification identity

Scalar pair getters require both successful MainLoop publication and an explicit
computed-product flag. Allocating or clearing an array does not make it a result.
`getHistCF()` therefore raises when the selected computation did not produce CF;
for the modern scalar two-ball engines, request `compute-HistN,and-CF` with a pair
computation. `getHistNN()` and scalar auto/cross getters likewise reject unavailable
products. Shear, physical and forest products belong to their dedicated getters.

Changing settings/catalogs, cleanup and failed runs invalidate publication. MPI
scalar products are available only on the publishing rank. A copied NumPy result
remains independent of later cleanup or reconfiguration. A computed zero remains
valid; it is not used as an availability sentinel. Neighbor-box empty bins now
publish the same zero convention already used by that engine's file output.

`getResults()` metadata now includes `inputs.effective_catalogs`. Immediately
before MainLoop, after input interpretation and before tree traversal or smoothing,
the wrapper hashes native catalog rows in bounded chunks. Floating fields use
little-endian float64 with signed zeros canonicalized. Mask/body/forest/LOS
identities use little-endian int64, without conversion through floating point.
The descriptor records the row count, dimension, encoding and schema version.

The identity includes coordinates, scalar or shear fields, weights, masks, and
applicable forest distances, sightlines and integer membership IDs. Unused
family-specific fields are zero. Catalog and row order are significant. In-memory
forest IDs are already compacted by the existing loader: renaming IDs without
changing their dense membership representation is equivalent, while changing
which rows belong to a forest is not. All loaded catalogs are fingerprinted,
including auxiliary mask/random catalogs.

Qualification compares these effective identities before applying numerical error
tolerances. Equivalent loaders/column selections can compare when they produce
identical canonical native rows and all other estimator controls match. File paths,
raw file SHA-256 hashes and selected columns/formats remain provenance; raw-byte
identity is not a substitute for the interpreted catalog. Changes to unused file
columns may therefore be compatible, while changing the selected field, weights,
masks or membership is incompatible even if numerical arrays remain close.

This policy is deliberately exact about input values: alternate coordinate
conversions with roundoff differences are not assumed equivalent. It does not
establish row-permutation invariance or equivalence between different estimators.
Existing result packets without `effective_catalogs` must be recomputed; numerical
agreement alone cannot reconstruct their missing input identity. Saved packets
with the new identity retain qualification support after reload.

The hashing step is excluded from MainLoop timing. A local identity-generation
failure participates in `MPI Python catalog identity` consensus before any rank
enters the collective calculation. This preserves coherent failure and recovery.

Focused regression suites:

```bash
PYTHONPATH="$PWD" python3 -m pytest -q \
  tests/python/test_histogram_availability.py \
  tests/python/test_effective_catalog_identity.py \
  tests/python/test_cython_in_memory_catalog.py
```

On macOS, run the same tests with `MallocScribble=1` to verify that unavailable
products cannot expose allocator bytes. The count-CF test uses an independent
pair enumeration and shell-volume formula, including invalidation and recovery.
