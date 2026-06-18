Inputs and Catalog Formats
==========================

cTreeBalls accepts generated test models or one or more input catalogs.
Runtime parameters may be supplied on the command line or in a parameter file.

Generated Test Models
---------------------

The compact validation path uses ``testmodel`` and ``nbody`` rather than an
external file.  Common generated geometries include simple cubic and random
unit-sphere catalogs.  ``seed`` controls reproducibility and ``lengthBox`` sets
the Cartesian test volume when applicable.

Catalog Parameters
------------------

``infile``
    One filename or a comma-separated list.

``infileformat``
    One format or a comma-separated list aligned with ``infile``.

``iCatalogs``
    Integer labels selecting catalogs for auto- or cross-correlation paths.

``outfile`` and ``outfileformat``
    Optional conversion output.  Pair with ``options=stop`` when only checking
    or converting an input catalog.

Standard Formats
----------------

``columns-ascii``
    Text rows containing positions and a scalar value.  The conventional
    three-dimensional layout is ``x y z value``.

``binary``
    Native binary catalog format used by several bundled tests.

``takahashi``
    Takahashi all-sky ray-tracing map input.

Additional Formats
------------------

With the matching add-ons enabled, the code includes ``multi-columns-ascii``,
Gadget, and ``fits-healpix`` readers.  The detailed format list is retained in
:doc:`../catalog_files`.

Multiple Catalogs
-----------------

A parameter-file pattern for multiple maps is:

.. code-block:: text

   infile = map_a.fits,map_b.fits,map_c.fits
   infileformat = fits-healpix,fits-healpix,fits-healpix
   iCatalogs = 1,2,3

Confirm that each list has the expected length and that the selected search
method supports the intended cross-correlation.

Units and Geometry
------------------

For unit-sphere maps, ``rminHist`` and ``rangeN`` are angular separations in
radians.  For Cartesian catalogs they follow the position units in the input
file.  Record this convention with every output because the histogram files do
not independently encode the input-unit metadata.
