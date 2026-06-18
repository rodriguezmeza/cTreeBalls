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
