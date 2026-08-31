Two-point Correlation Functions
===============================

The selected field and estimator determine the meaning of a pair histogram.
See :doc:`search_methods` before comparing outputs across addons.

Scalar and Counts
-----------------

``histNN`` counts accepted pairs; it is not by itself a random-catalog-corrected
galaxy correlation estimator. ``compute-HistN`` requests this product in the
applicable engines. ``and-CF`` selects additional count-derived output where
supported.

``histXi2pcf`` is the scalar pair-correlation product in the corresponding
engine. Check whether normalization uses counts or statistical weights
and whether 2PCF was compiled/enabled. In the angular scalar engines,
unit-sphere separation bins are Euclidean chords, not exactly radians.
The corrected angular-bearing policy excludes undefined legs from 3PCF but
does not remove otherwise valid radial pairs.

Other Fields
------------

* :doc:`shear`: complex spin-2 two-point functions with their pair weights.
* :doc:`lyman_alpha`: anisotropic or radial-only forest products with same-quasar exclusions.
* :doc:`scalar_3d`: physical-3D scalar and data-minus-random survey estimators.

Where available, ``only-2pcf`` avoids unrelated scalar 3PCF work.
The physical-3D family instead uses ``only-2pcf-3d``; forest methods encode
their selected orders in the search name. The two ``balltree-2balls-..._3pcf``
addons have no 2PCF mode. See :doc:`benchmarks` for comparable timing contracts.
