Benchmark workflows

Start with docs/benchmarks.rst and docs/search_methods.rst. Select an equivalent
field, geometry, binning, estimator, and normalization before comparing times.

Current optional CPU suite:
    addons/python_env/cputime_comparison/benchmark_kappa_corr.py

Its Readme.txt describes the Conda/Jupyter environment, external TreeCorr,
Corrfunc, FCFC, lya2pcf and ENCORE sources, MPI, timing scopes, and CSV outputs.
addons/python_env is a local optional workspace and may not be distributed.
Do not confuse its general suite with the older tests/python script of the
same name, which has a different CLI and a narrower benchmark purpose.

Current catalog drivers:
    python/kappa_corr_all_engines.py
    python/lya_corr_all_engines.py

Run each with --help and --list-engines. Read their READMEs for angle units,
masking, normalization, and rank ownership. They retain in-memory input across
engines rather than rereading files. Examples and notebooks are in examples/.

Before timing, from the root with matching C/Cython builds:
    python3 -m pytest -q tests/make_tests/test_scalar_numerical_contract.py

Use no-normalize-HistZeta,weights-norm for shared raw scalar 3PCF comparisons.
Native kernels now exclude repeated neighbors; do not subtract them again.
SMOOTHPIVOTON=1 does not turn smoothing on without options=smooth-pivot.
Rebuild both C and Cython after changing flags, then restart the notebook kernel.

A compiled 3PCF capability does not mean every pair-only timing performs the
same work. Use supported runtime order selectors and inspect work_contract.
Keep metadata, used-values files, counts, and numerical differences with timings.
Only benchmark differences within the same estimator; forest radial/3D and
angular/physical-3D outputs are not interchangeable.

Historical Takahashi download/plot workflows remain in tests/handson material.
Use bounded small-catalog checks before downloading or processing full-sky maps.
