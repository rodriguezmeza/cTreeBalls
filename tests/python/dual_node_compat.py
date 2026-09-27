"""Optional external dual-node backend with project-facing tolerance names.

Only this boundary knows the dependency's import and keyword spellings.
Importing it does not load or require the external backend. Returned catalogs
and correlations retain the dependency's native numerical behavior.
"""
from functools import lru_cache
from importlib import import_module
import math


class BackendUnavailable(ImportError):
    """The optional external dual-node package is not installed."""


@lru_cache(maxsize=1)
def _backend():
    try:
        return import_module("treecorr")
    except ModuleNotFoundError as exc:
        if exc.name != "treecorr":
            raise
        raise BackendUnavailable("Install the optional external dual-node backend") from exc


def version():
    return _backend().__version__


def catalog(**settings):
    """Construct a catalog; units, fields, weights and masks remain explicit."""
    return _backend().Catalog(**settings)


def correlation(kind, *, bin_theta=0.0, angle_theta=0.0, **settings):
    """Construct an NN, KK, GG, NNN, KKK or GGG correlation.

    Tolerances pass through unchanged; they are not a guaranteed relative-error
    bound. Zero disables these geometric tolerances. Binning, metric and
    estimator matching remain the caller's responsibility.
    """
    constructors = {"NN": "NNCorrelation", "KK": "KKCorrelation", "GG": "GGCorrelation",
                    "NNN": "NNNCorrelation", "KKK": "KKKCorrelation", "GGG": "GGGCorrelation"}
    if kind not in constructors:
        raise ValueError(f"Unsupported external dual-node correlation: {kind}")
    for name, value in (("bin_theta", bin_theta), ("angle_theta", angle_theta)):
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if "bin_slop" in settings or "angle_slop" in settings:
        raise TypeError("Use bin_theta and angle_theta in the project interface")
    # Translation must precede construction; upstream validates its own keywords.
    kwargs = dict(settings, bin_slop=bin_theta, angle_slop=angle_theta)
    return getattr(_backend(), constructors[kind])(**kwargs)
