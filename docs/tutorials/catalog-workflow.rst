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

For unit-sphere data, define the radial range in radians:

.. code-block:: text

   rminHist = 0.00213811
   rangeN = 0.0633205
   sizeHistN = 20

These values correspond approximately to 7.35--217.68 arcmin.  Choose bins for
the science case and verify them in the generated used-values file.

Run and Plot
------------

.. code-block:: bash

   ./cballs parameters_map.txt

The scripts under ``tests/python`` provide examples for loading HEALPix data,
running ``cyballs``, comparing outputs, and plotting 2PCF/3PCF products.  Use
``python script_name.py --help`` before applying a script to new data.

Large Takahashi Data
--------------------

Continue with :doc:`../handson` for the original Takahashi download example.
That realization is several gigabytes and is not required to verify cTreeBalls.
