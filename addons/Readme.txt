cTreeBalls addons

Enable features in addons/Makefile_addons_settings. Core and machine settings
are in Makefile_settings and Makefile_machine. Rebuild C and Cython together.

Discover the actual executable:
    ./cballs options=make-info
    ./cballs options=print-options
    ./cballs options=print-search-methods

The maintained method guide is docs/search_methods.rst. Current families:

- Scalar angular: octree-ggg-omp/mpi, kdtree-omp, balltree-omp/mpi,
  octree-balls4-omp/mpi, balltree-2balls-omp/mpi, octree-2balls-omp/mpi,
  balltree-2balls-omp_3pcf and balltree-2balls-mpi_3pcf (3PCF only).
- Spin-2 flat-sky shear: octree-shear-omp.
- Forest-aware 3D/radial estimators: seven lya-...-omp and seven lya-...-mpi
  searches; radial modes still require DEFDIMENSION=3.
- Physical 3D scalar/survey: octree-3pcf-3d-omp/mpi.
- Other box, historical, and experimental engines remain profile-dependent.

The scalar angular contract is in docs/3pcf.rst: observer-centered coordinates,
chord distances, distinct weighted raw triplets, and runtime-only smoothing.
Different field/geometry estimators are not interchangeable.

Utilities include CLASS parsing, IOLIB, Gadget I/O, CFITSIO/HEALPix, PXD hooks,
and bundled GSL/CFITSIO. Each active addon directory documents its own flags.

Archived docs, backups, development queues, and vendored third-party READMEs
are not descriptions of the active build. addons/python_env is an optional
local benchmark workspace, not a required installed-package dependency.


