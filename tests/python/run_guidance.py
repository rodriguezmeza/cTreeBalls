"""Shared frontend guidance and explicitly reference-relative diagnostics.

These diagnostics do not declare a reference exact or replace existing gates.
"""
import numpy as np

REFERENCE_COLUMNS = [
    "reference_nonzero_bins", "zero_reference_mismatches", "nonfinite_bins",
    "max_reference_relative", "coefficient_status",
]


def reference_relative_summary(candidate, reference, rtol):
    """Compare every coefficient, including complex coefficients, without a floor.

    Ratios are fractions (0.10 = 10%). Exact-zero references must remain zero.
    Nonfinite data cannot pass; empty arrays are unmeasured.
    The caller is responsible for matching catalogs, bins and estimators.
    """
    candidate = np.asarray(candidate)
    reference = np.asarray(reference)
    if candidate.shape != reference.shape:
        return {"coefficient_status": "shape-mismatch"}
    if not np.isfinite(rtol) or rtol < 0:
        raise ValueError("rtol must be nonnegative and finite")
    finite = np.isfinite(candidate) & np.isfinite(reference)
    nonzero = finite & (reference != 0)
    zero = finite & (reference == 0)
    nonfinite_count = int(np.count_nonzero(~finite))
    zero_mismatch = int(np.count_nonzero(zero & (candidate != 0)))
    # Conversion avoids unsigned subtraction and preserves complex magnitudes.
    dtype = np.result_type(candidate.dtype, reference.dtype, np.float64)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        error = np.abs(candidate[nonzero].astype(dtype) - reference[nonzero].astype(dtype))
        relative = error / np.abs(reference[nonzero].astype(dtype))
    maximum = float(np.max(relative)) if relative.size else 0.0
    status = "pass"
    if nonfinite_count or zero_mismatch or not np.isfinite(maximum) or maximum > rtol:
        status = "fail"
    elif candidate.size == 0:
        status = "unmeasured"
    return {
        "reference_nonzero_bins": int(np.count_nonzero(nonzero)),
        "zero_reference_mismatches": zero_mismatch,
        "nonfinite_bins": nonfinite_count,
        "max_reference_relative": maximum if candidate.size else float("nan"),
        "coefficient_status": status,
    }


def add_run_guidance(parser, *, field, benchmark=False, private=False):
    """Append advice without changing scientific defaults or parser behavior."""
    parts = [
        "Accuracy: a per-coefficient limit uses abs(fast-reference)/abs(reference) "
        "for every nonzero coefficient (complex magnitude for complex multipoles); "
        "a zero reference must remain zero. A relative L2 target does not bound "
        "every coefficient. Require the same catalog selection, weights, bins, "
        "geometry, normalization and multipole orders in the exact reference. "
        "Agreement between two approximate engines is not an exact-error measurement.",
        "Timing: benchmark on an otherwise idle machine, record cache/warmup state, "
        "and compare repeated wall times at matched measured accuracy. CPU time "
        "sums worker work and can exceed elapsed time. Loading, Python Run, native "
        "MainLoop and plotting are different timing scopes. Chunking, explicit "
        "catalog ownership transfer and automatic compact node sizing preserve "
        "the sample; thinning or downsampling changes it. See README_run_controls.md.",
    ]
    if field == "kappa":
        parts.append(
            "Scalar two-ball speed controls: scalar-pivot-reuse,no-smooth-pivot; "
            "dual-node-bin-theta is a separate native options token. "
            "CBALLS_SCALAR_PIVOT_TOL is a phase budget in radians [0,3] and "
            "CBALLS_SCALAR_BIN_THETA is a bin-width allowance [0,1]. Values 3 "
            "and 0.5 do not guarantee 5% or 10% error. Reuse requires positive "
            "theta and eligible nonperiodic 3PCF. An unsmoothed scalar two-ball "
            "body reference uses theta=0 and no-one-ball,no-two-balls,no-smooth-pivot, "
            "without scalar-pivot-reuse. Profile separately with native options "
            "no-out-Hist,dual-node-profile; remove profiling for timing trials. "
            "nsmooth controls leaf work as well as engine-specific smoothing; "
            "no-smooth-pivot does not force singleton leaves."
        )
    elif field == "shear":
        parts.append(
            "Shear controls use CBALLS_SHEAR_* with shear-pivot-reuse, not "
            "CBALLS_SCALAR_*. Phase/bin tolerances and theta are independent, "
            "and none is a percentage error guarantee. For a spherical two-ball "
            "body reference use theta=0 and no-one-ball,no-smooth-pivot; "
            "no-two-balls alone is insufficient for 3PCF. Unset CBALLS_SHEAR_PROFILE "
            "for timing trials: even a value of 0 enables instrumentation."
        )
    elif field == "lya":
        parts.append(
            "Ly-alpha uses its own forest, LOS, kernel and bin-slop controls. "
            "Angular scalar/shear pivot-reuse controls do not configure forest "
            "estimators. Compare only matching estimator families and exclusions."
        )
    else:
        parts.append(
            "ENCORE comparisons use physical 3D estimators and file input; angular "
            "scalar/shear reuse and Python catalog-transfer controls do not apply."
        )
    if field in ("kappa", "shear"):
        parts.append(
            "Angular min/max separation parameters are outer bin edges, not "
            "first/last bin centers. Read the saved effective bins when plotting. "
            "The Sofia figure-1 layout targets 20 centers from 8 to 200 arcmin; "
            "8-to-200 arcmin edges define a different grid. A plot selection "
            "cannot rebin existing coefficients."
        )
    if benchmark and field in ("kappa", "lya"):
        parts.append(
            "Comparison reports: median_relative/p95_relative/max_relative are "
            "symmetric differences with an atol floor. max_reference_relative "
            "is the unfloored ratio to the named reference (multiply by 100 for "
            "percent); coefficient_status checks all coefficients at rtol with "
            "exact-zero matching. The existing status/exit gate uses allclose: "
            "abs(delta)<=atol+rtol*abs(reference). Use atol=0 for a strict "
            "fractional gate. Approximate runs alone supply no exact baseline."
        )
    if benchmark and field == "shear":
        parts.append(
            "This shear benchmark reports timings and checksums, not a complete "
            "coefficient-error comparison or exact-reference qualification. "
            "Plots are enabled by default; --no-plots disables them. Named "
            "calculation binning presets fill omitted limits; explicit limits override them."
        )
    if private:
        parts.append(
            "Private dual-node comparison backend controls remain separate from "
            "native cTreeBalls options and environment variables. Native controls "
            "do not configure the peer backend. Consult this driver's binning "
            "help: a plotting-only preset does not change computation bins."
        )
    parser.epilog = " ".join(filter(None, [parser.epilog, *parts]))
