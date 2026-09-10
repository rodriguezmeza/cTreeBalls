Checked. The `CFITSIOLIBON=1` fix works for CFITSIO now.

What happened:

`make all` built and installed `cyballs`, but import initially failed with a NumPy ABI error because `pyproject.toml` says:

```toml
numpy>=1.22,<2
```

while your runtime Python has:

```text
NumPy 2.1.3
```

So pip’s isolated build compiled `cyballs` against NumPy 1.x, then Python tried to import it with NumPy 2.1.3.

I rebuilt with:

```bash
python -m pip install . --no-build-isolation
```

and now import passes:

```text
/opt/homebrew/anaconda3/lib/python3.13/site-packages/cyballs.cpython-313-darwin.so
2.1.3
```

Also confirmed the vendored CFITSIO path is working: `otool -L` no longer shows external `libcfitsio`, and `nm -u` no longer shows unresolved FITS symbols like `_ffclos`.

To make `make all` reliable, change the `cyballs` target in `Makefile` from:

```make
$(PYTHON) -m pip install .
```

to:

```make
$(PYTHON) -m pip install . --no-build-isolation
```

or update `pyproject.toml` to allow NumPy 2:

```toml
numpy>=2.0
```

For your current environment, `--no-build-isolation` is the quickest correct fix.
