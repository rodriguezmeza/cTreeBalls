kdtree-omp

Enable KDTREEOMPON=1 in addons/Makefile_addons_settings.
TWOPCFON and TPCFON select compiled 2PCF and scalar angular 3PCF support.
nsmooth sets leaf capacity. Exact traversal is the default; behavior-ball
enables cell aggregation with logarithmic bins. no-one-ball disables it.

In 3D, retain the observer origin and use chord separations on unit-sphere
catalogs. Higher scalar modes use projected tangent bearings, not chord
interior angles. See docs/3pcf.rst for the shared angular convention.

no-normalize-HistZeta selects raw weighted distinct-triplet multipoles.
Cell moments retain both sum(weight*kappa) and sum((weight*kappa)^2);
the latter removes repeated neighbors correctly. Use weights-norm alongside
the raw option for portable comparisons across engine families.
The legacy normalized path is a different contract.

SMOOTHPIVOTON=1 does not activate smoothing. An explicit smooth-pivot option
is required on supported paths; it is rejected with raw KD-tree multipoles.
Setting rsmooth without that option leaves smoothing inactive.

Validate with tests/make_tests/test_scalar_numerical_contract.py using matching
C and Cython builds. See docs/search_methods.rst and docs/benchmarks.rst.
There is no kdtree-mpi sibling in the current maintained scalar engine list.


