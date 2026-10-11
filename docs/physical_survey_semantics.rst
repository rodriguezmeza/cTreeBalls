Physical survey results and acceptance
======================================

The physical Legendre survey engines compute correction after globally summing
the D-alpha*R and alpha*R measurements and **before** optional file output.
``getPhysicalResults()`` and ``getResults()`` return independent NumPy copies
of both raw and corrected products, including with ``no-out-Hist``. Saved
result packets retain these products and their validity after model cleanup.
This contract applies to ``octree-3pcf-3d-omp`` and its MPI implementation.

Orders and normalization
------------------------

``surveySignalLmax`` is the highest returned corrected order. Its default -1
preserves the historical output extent: ``mChebyshev-1``, or ``mChebyshev``
with ``survey-keep-top-multipole``. Explicit nonnegative values override that
mapping. Non-survey physical mode retains its existing ``mChebyshev`` contract.

``surveySolvePadding`` (default 1) adds nuisance signal modes to the solve.
Let S = signal order + padding. The numerator is measured through S, the solve
uses orders 0..S, and only the requested signal orders are returned.
``surveyWindowLmax`` (default -1) measures the random window through 2S.
An explicit window order must be at least 2S and at most 64; S is at most 32.
Insufficient bandwidth is rejected before traversal. Higher window modes have
zero coupling inside the chosen finite signal space. Increase padding to test
signal truncation; merely raising an already sufficient window order cannot
remove power from omitted signal modes.

For standard coefficients multiplying Legendre polynomials, the system is

.. math::

   M_{ij}=(2i+1)\sum_{k=0}^{2S}\frac{R_k}{R_0}
             \begin{pmatrix}i&j&k\\0&0&0\end{pmatrix}^{2},
   \qquad M\zeta=N/R_0.

Every random moment uses distinct catalog members, including same-shell pairs.
Only R0 is necessarily positive; higher random multipoles can have either sign.
The ENCORE output coefficient is
``4*pi*(-1)**ell*zeta_legendre/sqrt(2*ell+1)``.
The pair correction is ``(D-alpha*R)^2/(alpha*R)^2`` with the same normalization
of random weights. Random density and angular bandwidth are independent
convergence controls. A small solve residual alone cannot certify an adequate
selection model, signal bandwidth, or random density.

Owned products and validity
----------------------------

Raw ``data_pair_*`` and ``data_triple_*`` contain D-alpha*R products; raw
``random_pair_*`` and ``random_triple_*`` contain alpha*R products. The
``*_numerator`` and ``*_denominator`` names remain for compatibility.
These raw triples have different leading extents: S+1 for data and W+1 for
randoms. They are not corrected observables and their ratio is not a general
3PCF correction.

* ``zeta_corrected_legendre`` and ``zeta_corrected_encore`` have shape
  ``(signal_lmax+1, nbins, nbins)``.
* ``zeta_valid``, ``zeta_status``, ``zeta_rcond``, and
  ``zeta_relative_residual`` have shape ``(nbins, nbins)`` and apply to every
  coefficient of the corresponding radial pair.
* ``xi_corrected``, ``xi_valid`` and ``xi_status`` have shape ``(nbins,)``.

Unsupported coefficients are NaN, never a valid measured zero. Status codes
are 0 empty support, 1 valid, 2 numerically singular, 3 poorly conditioned,
4 excessive residual, and 5 nonfinite or unrepresentable derived arithmetic.
An invalid bin does not invalidate an otherwise completed run; callers must
select supported bins. A valid zero signal has valid=1 and finite coefficients.

``zeta_rcond`` is ``1/(norm_inf(M)*norm_inf(inv(M)))`` in the standard
Legendre basis. ``surveyRcondMin`` (default 1e-10, valid range (0,1]) is the
minimum accepted value. ``zeta_relative_residual`` is
``norm_inf(M*x-b)/(norm_inf(M)*norm_inf(x)+norm_inf(b))``; a zero system RHS
has residual zero. The acceptance limit is ``128*(S+1)*DBL_EPSILON``.
Empty or singular support has rcond=0; a residual that cannot be evaluated is
NaN. A well-conditioned solve does not bound the astrophysical noise.

Text files retain their first 14 triple columns, replacing the old pivot proxy
with the infinity-norm condition at column 13 (one based), and append rcond,
relative residual, and status. The validity column remains 14. Pair files
append status after valid. Files emit only bin1<bin2 for compatibility; API
arrays also expose same-shell and transposed pairs. The condition in a text
row can be infinity for invalid support. All coefficients use 17-digit output.

``solve_physical_survey(numerator, random_window, signal_lmax=1, padding=1)``
provides the same pure native solve for measured moments. The first axis is
ell and trailing axes must agree. It returns corrected standard coefficients
and diagnostics without modifying inputs. This is useful for injected-signal
and bandwidth tests; callers retain responsibility for moment normalization.

Qualification and science acceptance
-------------------------------------

Fast search defaults are unchanged. Execution success, exact tuple agreement,
and science acceptance are distinct outcomes. ``qualification_report`` and
``model.qualifyAgainst`` require a matching explicit exact reference and compare
corrected physical products whenever survey raw products are present. They
reject incomplete legacy raw-only survey packets. Survey orders and condition
policy are part of the comparison scope; validity codes must match exactly.
Rcond/residual arrays are checked as support diagnostics, separately from
scientific error norms. Entirely unsupported corrected products cannot qualify.

The default metric is per-observable L2 error. Set ``metric='pointwise'`` to
require ``abs(error)<=atol+rtol*abs(reference)`` in every selected bin.
``selection={name: boolean_array}`` supplies explicit masks with the exact
observable shapes; omitted masks select all bins. Masks must select something.
Their shapes, counts and hashes are retained in the qualification scope.
Shear components and forest normalized products are additionally compared
separately, inheriting their parent selection. Empty/nonfinite masks must
match. Absolute tolerances give weak or zero signals a stated error budget.

Run the mandatory campaign independently with::

   python scripts/survey_science_acceptance.py --output /new/evidence/directory

It retains exact, repeated exact, and default packets on signed, clustered,
masked and weak-signal fixtures for scalar, shear and forest engines. All
returned observables, normalized forest products and individual shear
components are evaluated, with full-bin L2 and predeclared interior-radial-bin
pointwise metrics. Campaign ``status=PASS`` requires complete comparisons and
passing exact controls. ``science_acceptance`` independently records whether
the defaults meet the stated science threshold; rejected candidates remain
rejected in the release evidence. This is bounded evidence for these fixtures,
not automatic acceptance of the user's production survey.

The survey regression also retains a positive-window analytic counterexample,
Gauss quadrature polynomial recovery, and a catalog tuple injection through
the actual measured random window at 24, 48 and 96 random points and window
orders 4 and 6. The latter isolates correction semantics from stochastic
D-R signal estimation; it is not an unbiased-cosmology experiment. Additional
tests exercise solve padding, high window order 64, singular/empty support,
condition rejection, serialization, copied arrays, reused objects and MPI.
