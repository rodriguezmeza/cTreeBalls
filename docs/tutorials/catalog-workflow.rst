Catalog and Map Workflow
========================

This page outlines a production-style catalog workflow without requiring a
specific external data download.

Choose and Validate the Format
------------------------------

Select an ``infileformat`` supported by the active build.  For a HEALPix map:

.. code-block:: text

   infile = map.fits
   infileformat = fits-healpix
   iCatalogs = 1

For a plain three-dimensional scalar catalog, use ``columns-ascii`` with rows
containing ``x y z value``.

Convert Without Searching
-------------------------

Use ``outfile`` with ``options=stop`` to test reading or convert formats:

.. code-block:: bash

   ./cballs infile=catalog.txt infileformat=columns-ascii \
      outfile=checked_catalog outfileformat=columns-ascii options=stop

Inspect the converted catalog before launching a long search.

Configure Angular Bins
----------------------

For unit-sphere scalar angular data, define the radial range in chord distance:

.. code-block:: text

   rminHist = 0.00213811
   rangeN = 0.0633205
   sizeHistN = 20

Convert intended angles with ``r=2*sin(alpha/2)``, where ``alpha`` is in
radians; these small chord values are approximately 7.35--217.68 arcmin.
The distinction matters at wide angles. Choose bins for the science case,
preserve the observer origin, and verify the generated used-values file.

Run and Plot
------------

.. code-block:: bash

   ./cballs parameters_map.txt

The current multi-engine driver is ``python/kappa_corr_all_engines.py``;
see :doc:`../user/python` and its README for masks, edge corrections, input
reuse, and angle controls. ``tests/python`` also contains historical plotting
examples. Use each script's ``--help`` before applying it to new data.

Large Takahashi Data
--------------------

Continue with :doc:`../handson` for the original Takahashi download example.
That realization is several gigabytes and is not required to verify cTreeBalls.
