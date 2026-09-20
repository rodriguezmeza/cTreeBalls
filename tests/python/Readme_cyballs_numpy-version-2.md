# Building cyballs with NumPy 2

Build against the NumPy version in the target environment. For a checkout
build, `make -j4 all` updates the native targets and in-place Cython extension;
it does not install a package into the environment. To install using the
current environment's NumPy and Cython:

```sh
python3 -m pip install . --no-build-isolation
```

Install the build requirements first when disabling isolation. Restart Python
after rebuilding, and check `cyballs.__file__` and `cyballs.build_info()` to
identify the extension actually loaded. Native C and Cython must use identical
feature flags; do not mix a generated PXD or archive from another profile.

The public profile uses external GSL and CFITSIO. Configure `gsl-config`,
`pkg-config cfitsio`, or explicit include/library paths for that environment.
