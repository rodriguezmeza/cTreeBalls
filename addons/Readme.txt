cTreeBalls addons

Enable features in addons/Makefile_addons_settings. Core and machine settings
are in Makefile_settings and Makefile_machine. Rebuild C and Cython together.

Discover the actual executable:
    ./cballs options=make-info
    ./cballs options=print-options
    ./cballs options=print-search-methods

The maintained method guide is docs/search_methods.rst. The active profile has:

- Scalar angular KD-tree, PCA-ball-tree, and native-octree 2-ball methods,
  each with OpenMP and MPI entry points.
- Full-sky spin-2 octree, KD-tree, and PCA-ball-tree 2-ball methods.
- Forest-aware 3D and radial estimators with OpenMP and MPI entry points.
- Physical 3D scalar/survey octree estimators with OpenMP and MPI entry points.
- KD-tree and neighbor-box methods for periodic Cartesian 2PCF work.

The scalar angular contract is in docs/3pcf.rst: observer-centered coordinates,
chord distances, distinct weighted raw triplets, and runtime-only smoothing.
Different field/geometry estimators are not interchangeable.

Utilities include CLASS parsing, IOLIB, Gadget I/O, CFITSIO/HEALPix, and PXD
hooks. GSL and CFITSIO are external libraries in this profile. Some active
engines privately link shared or compatibility components; those components do
not register additional public search names.

