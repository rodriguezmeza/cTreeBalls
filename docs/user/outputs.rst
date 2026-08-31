Outputs and File Formats
========================

Run products are written beneath ``rootDir``.  Exact files depend on enabled
statistics, search method, and ``options``.

Common Outputs
--------------

``rBins.txt``
    Radial-bin values used by the histogram outputs.

``histNN*.txt``
    Neighbor-count histogram when ``compute-HistN`` is requested.

``histCF*.txt``
    Count-derived correlation output when ``and-CF`` is requested.

``histXi2pcf*.txt``
    Two-point correlation function.

``histZetaM*.txt``
    Three-point multipole matrices.  Component and multipole suffixes depend
    on the active 3PCF implementation.

``mHistZetaM*.txt`` and ``histZetaG*.txt``
    Additional multipole/angle representations produced by selected options
    and add-ons.

Provenance and Logs
-------------------

Each run writes a ``*-usedvalues`` parameter snapshot beneath ``rootDir``.
When ``verbose_log`` is positive, log output is written under the run's
temporary/log directory.  Treat the used-values file as part of every analysis
artifact.

Naming Controls
---------------

``histNNFileName``, ``histXi2pcfFileName``, and ``histZetaFileName`` change the
base names.  ``suffixOutFiles`` appends a run-specific suffix and is useful
when comparing methods without overwriting previous outputs.

Recommended Output Practice
---------------------------

* use a unique ``rootDir`` or suffix for each parameter set;
* keep the parameter file and generated used-values file together;
* record the commit hash and Makefile settings;
* verify bin units before plotting or comparing runs;
* do not interpret the reduced quickstart settings as science-ready values.

Scalar Complex Modes
--------------------

The four ``histZetaM_cos_N``, ``sin_N``, ``sincos_N`` and ``cossin_N``
matrices reconstruct ``coscos+sinsin + 1j*(sincos-cossin)``.
``N=1`` denotes order zero. Compare reconstructed modes, not individual
basis-dependent component matrices, under sky rotations.

``no-normalize-HistZeta,weights-norm`` selects the shared weighted raw
distinct-triplet contract. KD-tree, legacy balltree, and BALLS4 remove repeated
neighbors natively; no Python subtraction is required. See :doc:`../3pcf`.

For complex edge correction, ``histZetaM_EE_N.txt`` contains the real part and
``histZetaM_EE_Im_N.txt`` the imaginary part. The kappa driver saves
``zeta_edge_complex_N`` in its NPZ products. Read
``getHistZetaM_EE_complex(N)`` after a successful compatible run.
Empty/singular windows yield zero; physical-3D survey validity/conditioning
columns are separate products. MPI output is owned by rank 0.

Use fresh directories after changing a build or estimator. Older affected
raw/multipole files cannot be made current simply by renaming their headers.
