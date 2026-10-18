# Catalog memory and timing controls

The scalar compact octree, PCA ball-tree and KD-tree builders allocate node
capacity from the required tree shape. This is automatic in rebuilt binaries;
there is no smaller-tree option. It preserves leaf policy, point order and
approximation settings. The KD parallel builder and shear PCA builder already
used exact capacity. Other engine trees are not all affected by this change.
Fewer reserved nodes do not imply an equal reduction in peak RSS, particularly
with lazily committed calloc pages. There is no guaranteed speedup.

The kappa and shear drivers and HEALPix CPU benchmarks accept
`--fits-chunk-pixels 1048576` (default; kappa minimum 1024, shear any positive size). Smaller chunks bound
coordinate and field-gather temporary arrays; they preserve the retained pixels
and their order. Native-resolution kappa also scans the FITS map in chunks.
Downsampling can still require full input/output maps. Shear DES table loading
retains `--fits-chunk-rows` as an alias for the same chunk-size setting. These are Python arguments, not
native `options=` tokens or environment variables.

Kappa, shear and Ly-alpha drivers, including the private `*_dual_node.py`
counterparts, accept `--transfer-catalog-ownership`. The native Python benchmark
paths expose the same opt-in flag. Copy-on-registration remains the default.
Transfer avoids another NumPy copy of owning, aligned, C-contiguous native
float64 input buffers (uint8 for masks). Frontends materialize views or converted
fields during setup when needed; these conversions still allocate memory.
The strict API `set_catalog(..., take_ownership=True)` rejects conversions and
views. Transferred arrays become read-only; all existing writable aliases must
also stop being mutated or resized. The model retains strong references through
`struct_cleanup` for reuse; replacement, `clear_catalogs`, `clean_all` or
destruction releases them. Input arrays remain read-only afterwards.

`set_forest_catalog(..., take_ownership=True)` transfers positions, delta and
weights under the same contract. Integer forest IDs are still compacted into
private storage. Apply continuum projection before transfer; the API rejects
combining `project_delta=True` with ownership transfer. Ly-alpha drivers already
project before registration. Native `Run` still copies registered arrays to
particle storage. MPI maintains a catalog and tree in each rank. External peer
backends retain their own input behavior. ENCORE executable benchmarks use file
input and report that explicitly; a NumPy transfer flag would not apply there.

Memory policy is saved in driver/benchmark metadata. Input preparation and
registration remain outside correlation timing. Compare both setup and compute
costs when assessing transfer. Wall time is elapsed time; process CPU includes
all worker threads. MPI wall uses the maximum participating rank and CPU sums
participating ranks. Never divide CPU time by threads and call it wall time.
Repeated benchmarks retain catalogs/results and can reuse caches; process
high-water RSS from one long sweep is not an independent per-engine peak.
Use fresh processes with identical catalogs, threads, accuracy controls and
cache policy for memory or performance comparisons.

`CBALLS_MEMORY_BUDGET_MB` limits checked native allocations/plans per rank; it is
not a process RSS cap. Raising it does not add RAM. Fewer threads can reduce
private histogram scratch. `no-balltree-tree-cache` trades cache retention for
rebuilding on repeated runs in supporting engines. None of these memory changes
smooths cells or approximates the field. `--max-points` selects a different
sample; lowering NSIDE changes map resolution. Those scientific changes need
their own accuracy comparison. Neither theta nor reuse tolerances alone specify
a measured relative error.

Example additions to an existing HEALPix kappa or shear command:

```text
--fits-chunk-pixels 262144 --transfer-catalog-ownership
```

Keep the same input, selection, bins and approximation settings when comparing
against the original run. Use `options=make-info`, `print-options`,
`print-search-methods`, `print-environment` and `IA-tutor` for the corresponding
native build/help reports.
