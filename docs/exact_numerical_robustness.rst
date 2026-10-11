Exact-mode numerical range and distinct products
================================================

This contract covers the production scalar ``octree-2balls``, ``kdtree-2balls``
and ``balltree-2balls`` engines and the physical ``octree-3pcf-3d`` engine,
including their OpenMP and MPI variants. Use unsmoothed body-exact settings
(``theta=0,no-one-ball,no-two-balls,no-smooth-pivot``) for numerical validation.
Historical compatibility estimators selected by ``legacy-one-ball`` and the
self-inclusive ``octree-sincos`` estimator have separate definitions; this
contract does not change those definitions. Survey normalization has its own
acceptance tests and is not certified by these raw-estimator checks.

Distinct-product accumulation
-----------------------------

For disjoint sets A and B, with sums S and ordered distinct products D,

.. math::

   D(A\cup B) = D(A)+D(B)+S(A)\overline{S(B)}+S(B)\overline{S(A)}.

For real values the cross term is twice the product. Singleton sets have
D=0. Scalar tree nodes and per-pivot angular moments merge these quantities
directly. This avoids subtracting a dominant self square from an almost equal
squared sum. Same-bin complex diagonals are real by definition; exchanging
radial bins conjugates the complex moment. Physical harmonics use the same
prefix cross-products. With ``exclude-all-same-los``, complete line-of-sight
groups form the disjoint sets, so only allowed cross-group terms accumulate.

This is binary64 arithmetic, not arbitrary precision. Signed cancellation,
angular conditioning and reduction order can still affect the last bits.
Agreement is tested against independently enumerated Decimal sums, with
numerators and denominators checked separately. Ordinary exact-test tolerances
are unchanged. Subnormal checks also allow a small absolute number of minimum
binary64 quanta, because a fixed relative error is unattainable near zero.

Derived-arithmetic policy
-------------------------

Finite catalogs do not guarantee representable products. Raw scalar 3PCF,
window and physical numerator/denominator products must be representable in
the compiled accumulator. Nonfinite intermediate markers, task totals,
MPI totals and defined normalized results are checked before histogram
publication. Nonzero checked products that underflow to zero produce a
``numerical range`` error. Representable subnormals are retained. A request
may fail conservatively if an intermediate angular component is outside the
range, even when a particular final component could be representable.
No clipping or replacement of invalid arithmetic by a zero result occurs.

Scalar pair outputs expose a normalized field correlation and pair counts.
For extreme fields or weights (binary exponent outside [-200,200]), an
unsmoothed fallback accumulates products as compensated mantissa sums with
separate binary exponents. It forms the final ratio after combining these
scales. Thus constant-field correlations remain one in occupied bins even
for weights of 1e-300 or 1e300. Empty bins and all-zero-weight denominators
retain the documented zero representation. A genuinely unrepresentable
normalized field result fails. Smoothing at extreme range fails explicitly.

The fallback uses O(N^2) time and O(radial bins) memory and replicates its
work on each MPI rank. It is intended for numerical range protection, not
large-catalog acceleration. Ordinary-range catalogs retain the tree path.
Physical results expose raw products, so they reject unrepresentable raw
numerators or denominators even when their ratio alone would be finite.

MPI arithmetic failures reach a collective consensus. A failed Python run
has no valid results; replacing its catalog and rerunning the same object
must recover. Singular or empty edge-deconvolution cells retain their existing
status and NaN representation; those intentional statuses are separate from
an invalid raw calculation.

Acceptance and release evidence
-------------------------------

``tests/python/test_exact_numerical_robustness.py`` checks the analytic
four-point target 18H+6 for H from 1 through 1e16, signed multibin Decimal
references, permutations, threads, distinct LOS groups, empty and zero-weight
bins, subnormal products, overflow boundaries and recovery. Its MPI entry
point is run with one, two and four ranks by the active release gate.
``test_exact_numerical_cli.py`` checks raw results and rejects overflow and
underflow before files are written; the sanitizer component runs it too.
The clean installed-wheel oracle repeats the dynamic-range and recovery
contracts. The source-distribution policy requires these tests and headers.
