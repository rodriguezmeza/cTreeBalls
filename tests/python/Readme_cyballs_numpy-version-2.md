# Building cyballs with NumPy 2

Build the extension against the NumPy version in the target environment. If
an isolated pip build selects an incompatible NumPy ABI, rebuild from the
checkout with:

```bash
python3 -m pip install . --no-build-isolation
```

The maintained profile discovers CFITSIO as an external dependency through
`pkg-config cfitsio`. Set `PKG_CONFIG_PATH` when the installation is outside
the system search path. Build `cballs`, `libcballs.a`, and `cyballs` with the
same Makefile profile; mixing archives or generated PXD files from different
profiles is unsupported.
