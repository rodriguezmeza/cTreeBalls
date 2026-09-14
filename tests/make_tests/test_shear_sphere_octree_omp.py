#!/usr/bin/env python3
"""Full-sky spin-2 oracle and determinism tests."""

from __future__ import annotations

import ctypes
import os
import re
import tempfile

import numpy as np

from cyballs import CosmoComputationError, cballs, search_method_id


NMAX = 2
BINS = 3
PHI_BINS = 16
RMIN = 0.03
RMAX = 1.75
OPTIONS = "no-out-Hist,no-smooth-pivot,no-one-ball"
FAST_OPTIONS = "no-out-Hist,no-smooth-pivot"
SMOOTH_OPTIONS = "no-out-Hist,smooth-pivot,no-one-ball"
SMOOTH_RMIN_ARCMIN = 40.0
ENGINE = os.environ.get(
    "CBALLS_SHEAR_SPHERE_ENGINE", "octree-shear-sphere-2balls-omp",
)
DUAL_NODE_TWO_BALLS = ENGINE in {
    "octree-shear-sphere-2balls-omp",
    "kdtree-shear-sphere-2balls-omp",
    "balltree-shear-sphere-2balls-omp",
}


def fixture(nbody: int = 45):
    rng = np.random.default_rng(7182818)
    positions = rng.normal(size=(nbody, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    longitude = np.arctan2(positions[:, 1], positions[:, 0])
    latitude = np.arcsin(positions[:, 2])
    gamma = (0.018 + 0.006*np.cos(3.0*latitude)) \
        * np.exp(2j*(longitude + 0.21*np.sin(latitude)))
    weights = 0.65 + 0.7*rng.random(nbody)
    return positions, gamma, weights


def octant_fixture(nbody: int = 320):
    rng = np.random.default_rng(314159)
    positions = rng.uniform(0.05, 1.0, size=(nbody, 3))
    positions /= np.linalg.norm(positions, axis=1)[:, None]
    gamma = rng.normal(scale=0.02, size=nbody) \
        + 1j*rng.normal(scale=0.02, size=nbody)
    weights = 0.5 + rng.random(nbody)
    return positions, gamma, weights


def smooth_fixture(nbase: int = 18):
    positions, gamma, weights = fixture(nbase)
    separation = np.deg2rad(30.0/60.0)
    grouped_positions = []
    grouped_gamma = []
    grouped_weights = []
    for index, position in enumerate(positions):
        east, _ = tangent_basis(position)
        close_neighbor = np.cos(separation)*position + np.sin(separation)*east
        close_neighbor /= np.linalg.norm(close_neighbor)
        grouped_positions.extend((position, close_neighbor))
        grouped_gamma.extend((gamma[index], gamma[index]*(0.83 + 0.11j)))
        grouped_weights.extend((weights[index], 0.7*weights[index] + 0.2))
    return (
        np.asarray(grouped_positions),
        np.asarray(grouped_gamma),
        np.asarray(grouped_weights),
    )


def tangent_basis(unit):
    equatorial_norm = np.hypot(unit[0], unit[1])
    if equatorial_norm > 64.0*np.finfo(float).eps:
        east = np.array([-unit[1], unit[0], 0.0])/equatorial_norm
    else:
        east = np.array([1.0, 0.0, 0.0])
    north = np.cross(unit, east)
    return east, north


def geometry_and_transport(pivot, neighbor, gamma):
    chord = np.linalg.norm(neighbor - pivot)
    east, north = tangent_basis(pivot)
    tangent = neighbor - np.dot(pivot, neighbor)*pivot
    tangent /= np.linalg.norm(tangent)
    phase = complex(np.dot(tangent, east), np.dot(tangent, north))

    neighbor_east, _ = tangent_basis(neighbor)
    denominator = 1.0 + np.dot(pivot, neighbor)
    transported_east = east - np.dot(east, neighbor)/denominator \
        * (pivot + neighbor)
    transported_north = north - np.dot(north, neighbor)/denominator \
        * (pivot + neighbor)
    c = np.dot(neighbor_east, transported_east)
    s = np.dot(neighbor_east, transported_north)
    orientation_norm = np.hypot(c, s)
    rotation = complex(c/orientation_norm, s/orientation_norm)**2
    return chord, phase, gamma*rotation


def radial_bin(distance):
    if not RMIN < distance < RMAX:
        return None
    result = int((distance - RMIN)/(RMAX - RMIN)*BINS)
    return result if 0 <= result < BINS else None


def solve_coupling(upsilon, window):
    multipoles = 2*NMAX + 1
    corrected = np.zeros_like(upsilon)
    for first in range(BINS):
        for second in range(BINS):
            n_zero = window[2*NMAX, first, second]
            if abs(n_zero) == 0.0:
                continue
            coupling = np.empty((multipoles, multipoles), dtype=np.complex128)
            for row, ell in enumerate(range(-NMAX, NMAX + 1)):
                for column, order in enumerate(range(-NMAX, NMAX + 1)):
                    coupling[row, column] = (
                        window[ell - order + 2*NMAX, first, second]/n_zero
                    )
            if np.linalg.cond(coupling) <= 1.0e10:
                corrected[:, :, first, second] = np.linalg.solve(
                    coupling, (upsilon[:, :, first, second]/n_zero).T,
                ).T
    return corrected


def smooth_pivot_state(positions, gamma, weights, smooth_radius):
    active = np.ones(positions.shape[0], dtype=bool)
    pivot_gamma = weights*gamma
    pivot_weight = weights.copy()
    for pivot in range(positions.shape[0]):
        if not active[pivot]:
            continue
        for claimed in range(positions.shape[0]):
            if pivot == claimed or not active[claimed]:
                continue
            chord, _, transported_gamma = geometry_and_transport(
                positions[pivot], positions[claimed], gamma[claimed],
            )
            if chord < RMAX and chord <= smooth_radius:
                active[claimed] = False
                pivot_gamma[pivot] += weights[claimed]*transported_gamma
                pivot_weight[pivot] += weights[claimed]
    return active, pivot_gamma, pivot_weight


def oracle(positions, gamma, weights, smooth_radius=None):
    multipoles = 2*NMAX + 1
    denominator_modes = 4*NMAX + 1
    ring_max = max(2*NMAX, NMAX + 3)
    xi_plus = np.zeros(BINS, dtype=np.complex128)
    xi_minus = np.zeros(BINS, dtype=np.complex128)
    xi_weight = np.zeros(BINS)
    upsilon = np.zeros((4, multipoles, BINS, BINS), dtype=np.complex128)
    window = np.zeros((denominator_modes, BINS, BINS), dtype=np.complex128)
    active = np.ones(positions.shape[0], dtype=bool)
    pivot_gamma = weights*gamma
    pivot_weight = weights.copy()
    if smooth_radius is not None:
        active, pivot_gamma, pivot_weight = smooth_pivot_state(
            positions, gamma, weights, smooth_radius,
        )

    for pivot in range(positions.shape[0]):
        if not active[pivot]:
            continue
        rings = np.zeros((2*ring_max + 1, BINS), dtype=np.complex128)
        weight_rings = np.zeros_like(rings)
        diagonal_g6 = np.zeros(BINS, dtype=np.complex128)
        diagonal_g2 = np.zeros(BINS, dtype=np.complex128)
        diagonal_abs2 = np.zeros(BINS, dtype=np.complex128)
        diagonal_w2 = np.zeros(BINS)
        weighted_pivot_gamma = pivot_gamma[pivot]

        for neighbor in range(positions.shape[0]):
            if neighbor == pivot:
                continue
            distance, phase, transported_gamma = geometry_and_transport(
                positions[pivot], positions[neighbor], gamma[neighbor],
            )
            bin_index = radial_bin(distance)
            if bin_index is None:
                continue
            weighted_gamma = weights[neighbor]*transported_gamma
            for order in range(-ring_max, ring_max + 1):
                rings[order + ring_max, bin_index] += weighted_gamma*phase**order
                weight_rings[order + ring_max, bin_index] += \
                    weights[neighbor]*phase**order
            diagonal_g6[bin_index] += weighted_gamma**2*phase**-6
            diagonal_g2[bin_index] += weighted_gamma**2*phase**-2
            diagonal_abs2[bin_index] += \
                weights[neighbor]**2*abs(transported_gamma)**2*phase**-2
            diagonal_w2[bin_index] += weights[neighbor]**2
            xi_plus[bin_index] += weighted_pivot_gamma*np.conj(weighted_gamma)
            xi_minus[bin_index] += \
                weighted_pivot_gamma*weighted_gamma*phase**-4
            xi_weight[bin_index] += pivot_weight[pivot]*weights[neighbor]

        def ring(order, bin_index):
            return rings[order + ring_max, bin_index]

        def weight_ring(order, bin_index):
            return weight_rings[order + ring_max, bin_index]

        for first in range(BINS):
            for second in range(BINS):
                same_bin = first == second
                for order in range(-2*NMAX, 2*NMAX + 1):
                    value = weight_ring(order, first) \
                        * np.conj(weight_ring(order, second))
                    if same_bin:
                        value -= diagonal_w2[first]
                    window[order + 2*NMAX, first, second] += \
                        pivot_weight[pivot]*value
                for order in range(-NMAX, NMAX + 1):
                    products = np.array([
                        ring(order - 3, first)*ring(-order - 3, second),
                        ring(order - 1, first)*ring(-order - 1, second),
                        np.conj(ring(-order - 1, first))*ring(-order - 3, second),
                        ring(order - 3, first)*np.conj(ring(order - 1, second)),
                    ])
                    if same_bin:
                        products -= np.array([
                            diagonal_g6[first], diagonal_g2[first],
                            diagonal_abs2[first], diagonal_abs2[first],
                        ])
                    factors = np.array([
                        weighted_pivot_gamma,
                        np.conj(weighted_pivot_gamma),
                        weighted_pivot_gamma,
                        weighted_pivot_gamma,
                    ])
                    upsilon[:, order + NMAX, first, second] -= factors*products

    xi_plus = np.divide(xi_plus, xi_weight, out=np.zeros_like(xi_plus),
                        where=xi_weight != 0.0)
    xi_minus = np.divide(xi_minus, xi_weight, out=np.zeros_like(xi_minus),
                         where=xi_weight != 0.0)
    return {
        "xi_plus": xi_plus,
        "xi_minus": xi_minus,
        "xi_weight": xi_weight,
        "upsilon": upsilon,
        "window": window,
        "multipoles": solve_coupling(upsilon, window),
    }


def run_native_catalogs(catalogs, threads, options=OPTIONS, theta=1.0,
                        rsmooth=None, masks=None, engine=None):
    model = cballs()
    parameters = {
        "searchMethod": ENGINE if engine is None else engine,
        "iCatalogs": ",".join(str(index + 1)
                                for index in range(len(catalogs))),
        "usePeriodic": "false",
        "useLogHist": "false",
        "rangeN": RMAX,
        "rminHist": RMIN,
        "sizeHistN": BINS,
        "sizeHistPhi": PHI_BINS,
        "mChebyshev": NMAX,
        "theta": theta,
        "lengthBox": 2.2,
        "numberThreads": threads,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-sphere-"),
        "options": options,
    }
    if rsmooth is not None:
        parameters["rsmooth"] = str(rsmooth)
    model.set(parameters)
    for catalog, (positions, gamma, weights) in enumerate(catalogs):
        keyword_arguments = {
            "weights": weights,
            "gamma1": gamma.real,
            "gamma2": gamma.imag,
            "catalog": catalog,
        }
        if masks is not None:
            keyword_arguments["mask"] = masks[catalog]
        model.set_catalog(positions, **keyword_arguments)
    try:
        model.Run(level=["MainLoop"])
        option_tokens = set(options.split(","))
        result = {}
        if "only-3pcf" not in option_tokens:
            result.update({
                "xi_plus": model.getShearXiPlus().copy(),
                "xi_minus": model.getShearXiMinus().copy(),
                "xi_weight": model.getShearXiWeight().copy(),
            })
        if "only-2pcf" not in option_tokens:
            result.update({
                "upsilon": model.getShearUpsilonXMultipoles().copy(),
                "window": model.getShearWindowMultipoles().copy(),
                "multipoles": model.getShearGammaXMultipoles().copy(),
            })
        return result
    finally:
        model.struct_cleanup()


def run_native(positions, gamma, weights, threads, options=OPTIONS, theta=1.0,
               rsmooth=None, mask=None, engine=None):
    return run_native_catalogs(
        [(positions, gamma, weights)], threads, options, theta, rsmooth,
        None if mask is None else [mask], engine,
    )


def run_native_with_profile(positions, gamma, weights, threads, options):
    descriptor, path = tempfile.mkstemp(prefix="ctreeballs-shear-profile-")
    saved_stderr = os.dup(2)
    previous = os.environ.get("CBALLS_SHEAR_PROFILE")
    try:
        os.environ["CBALLS_SHEAR_PROFILE"] = "1"
        os.dup2(descriptor, 2)
        os.close(descriptor)
        result = run_native(
            positions, gamma, weights, threads, options=options,
        )
        ctypes.CDLL(None).fflush(None)
    finally:
        os.dup2(saved_stderr, 2)
        os.close(saved_stderr)
        if previous is None:
            os.environ.pop("CBALLS_SHEAR_PROFILE", None)
        else:
            os.environ["CBALLS_SHEAR_PROFILE"] = previous
    try:
        with open(path, encoding="utf-8") as stream:
            profile = stream.read()
    finally:
        os.unlink(path)
    return result, profile


def assert_results_identical(reference, compatibility, label):
    if set(reference) != set(compatibility):
        raise AssertionError(f"{label}: result fields differ")
    for name, values in reference.items():
        np.testing.assert_array_equal(
            compatibility[name], values,
            err_msg=f"{label}: compatibility changed {name}",
        )


def test_spherical_oracle_and_determinism():
    positions, gamma, weights = fixture()
    expected = oracle(positions, gamma, weights)
    one_thread = run_native(positions, gamma, weights, 1)
    for name, values in expected.items():
        np.testing.assert_allclose(
            one_thread[name], values, rtol=4.0e-12, atol=4.0e-13,
            err_msg=name,
        )
    many_threads = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
    )
    for name, values in one_thread.items():
        if not np.array_equal(values, many_threads[name]):
            raise AssertionError(f"spherical OpenMP result changed: {name}")
    narrow_theta = run_native(positions, gamma, weights, 1, theta=0.05)
    for name, values in one_thread.items():
        np.testing.assert_allclose(
            narrow_theta[name], values, rtol=4.0e-12, atol=4.0e-13,
            err_msg=f"theta changed exact spherical result: {name}",
        )

    fast_one = run_native(
        positions, gamma, weights, 1, options=FAST_OPTIONS, theta=0.05,
    )
    fast_many = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
        options=FAST_OPTIONS, theta=0.05,
    )
    for name, values in one_thread.items():
        np.testing.assert_allclose(
            fast_one[name], values, rtol=4.0e-12, atol=4.0e-13,
            err_msg=f"accepted-node limit differs from exact: {name}",
        )
        if not np.array_equal(fast_one[name], fast_many[name]):
            raise AssertionError(
                f"accepted-node result changed with OpenMP threads: {name}"
            )

    exact_pair = run_native(
        positions, gamma, weights, 1,
        options=f"{OPTIONS},only-2pcf", theta=0.05,
    )
    if DUAL_NODE_TWO_BALLS:
        no_two_balls_pair = run_native(
            positions, gamma, weights, 1,
            options=f"{FAST_OPTIONS},no-two-balls,only-2pcf", theta=1.0,
        )
        for name, values in exact_pair.items():
            np.testing.assert_array_equal(
                no_two_balls_pair[name], values,
                err_msg=f"no-two-balls did not select exact pairs: {name}",
            )
    fast_pair_one = run_native(
        positions, gamma, weights, 1,
        options=f"{FAST_OPTIONS},only-2pcf", theta=0.05,
    )
    fast_pair_many = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
        options=f"{FAST_OPTIONS},only-2pcf", theta=0.05,
    )
    for name, values in exact_pair.items():
        np.testing.assert_allclose(
            fast_pair_one[name], values, rtol=4.0e-12, atol=4.0e-13,
            err_msg=f"dual-tree limit differs from exact: {name}",
        )
        if not np.array_equal(fast_pair_one[name], fast_pair_many[name]):
            raise AssertionError(
                f"dual-tree result changed with OpenMP threads: {name}"
            )

    large_theta = run_native(positions, gamma, weights, 1, theta=2.0)
    for name, values in one_thread.items():
        np.testing.assert_allclose(
            large_theta[name], values, rtol=4.0e-12, atol=4.0e-13,
            err_msg=f"theta > 1 changed exact spherical result: {name}",
        )


def test_spherical_accepted_cell_transport():
    positions, gamma, weights = fixture(600)
    exact = run_native(
        positions, gamma, weights, 1, options=OPTIONS, theta=1.0,
    )
    fast_one = run_native(
        positions, gamma, weights, 1, options=FAST_OPTIONS, theta=1.0,
    )
    fast_many = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
        options=FAST_OPTIONS, theta=1.0,
    )
    limits = {
        "xi_plus": 5.0e-3 if DUAL_NODE_TWO_BALLS else 5.0e-4,
        "xi_minus": 2.5e-1 if DUAL_NODE_TWO_BALLS else 1.0e-1,
        "xi_weight": 5.0e-2 if DUAL_NODE_TWO_BALLS else 1.0e-12,
        "upsilon": 5.0e-2,
        "window": 1.0e-3,
        "multipoles": 6.0e-2,
    }
    accepted_a_cell = False

    for name, values in exact.items():
        scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
        relative_error = float(
            np.max(np.abs(fast_one[name] - values))/scale
        )
        if relative_error > limits[name]:
            raise AssertionError(
                f"accepted spherical cells exceed {name} error contract: "
                f"{relative_error} > {limits[name]}"
            )
        if not np.array_equal(fast_one[name], values):
            accepted_a_cell = True
        if not np.array_equal(fast_one[name], fast_many[name]):
            raise AssertionError(
                f"accepted-cell result changed with OpenMP threads: {name}"
            )
    if not accepted_a_cell:
        raise AssertionError("accepted-cell fixture exercised only body nodes")


def test_spherical_octant_frontier_parallelism():
    if ENGINE != "octree-shear-sphere-2balls-omp":
        return

    positions, gamma, weights = octant_fixture()
    options = f"{FAST_OPTIONS},only-3pcf"
    serial = run_native(
        positions, gamma, weights, 1, options=options,
    )
    threaded, profile = run_native_with_profile(
        positions, gamma, weights, min(4, os.cpu_count() or 1), options,
    )
    match = re.search(r"frontier_tasks=(\d+)", profile)
    if match is None or int(match.group(1)) <= 1:
        raise AssertionError(
            f"single-octant frontier did not branch:\n{profile}"
        )
    worker_pivots = [
        int(value) for thread, value in re.findall(
            r"thread=(\d+) pivots=(\d+)", profile,
        ) if int(thread) > 0
    ]
    if not any(value > 0 for value in worker_pivots):
        raise AssertionError(
            f"single-octant frontier left every worker idle:\n{profile}"
        )
    assert_results_identical(
        serial, threaded, "single-octant frontier determinism",
    )


def test_spherical_runtime_order_switches():
    positions, gamma, weights = fixture(39)
    combined = run_native(positions, gamma, weights, 1, theta=0.05)
    triple_only = run_native(
        positions, gamma, weights, 1,
        options=f"{OPTIONS},only-3pcf", theta=0.05,
    )

    assert set(triple_only) == {"upsilon", "window", "multipoles"}
    for name, values in triple_only.items():
        np.testing.assert_allclose(
            values, combined[name], rtol=4.0e-12, atol=4.0e-13,
            err_msg=f"only-3pcf changed {name}",
        )


def test_spherical_smooth_pivot_transport():
    positions, gamma, weights = smooth_fixture()
    oracle_radius = 2.0*np.sin(
        0.5*np.deg2rad(SMOOTH_RMIN_ARCMIN/60.0)
    )
    expected = oracle(positions, gamma, weights, smooth_radius=oracle_radius)
    raw = oracle(positions, gamma, weights)
    if np.array_equal(expected["xi_weight"], raw["xi_weight"]):
        raise AssertionError("smooth-pivot fixture did not group any bodies")

    one_thread = run_native(
        positions, gamma, weights, 1, options=SMOOTH_OPTIONS,
        theta=0.05, rsmooth=SMOOTH_RMIN_ARCMIN,
    )
    many_threads = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
        options=SMOOTH_OPTIONS, theta=0.05,
        rsmooth=SMOOTH_RMIN_ARCMIN,
    )
    for name, values in expected.items():
        np.testing.assert_allclose(
            one_thread[name], values, rtol=5.0e-12, atol=5.0e-13,
            err_msg=f"spherical smooth-pivot transport: {name}",
        )
        if not np.array_equal(one_thread[name], many_threads[name]):
            raise AssertionError(
                f"spherical smooth-pivot result changed with threads: {name}"
            )

    pair_options = "no-out-Hist,smooth-pivot,only-2pcf"
    pair_result = run_native(
        positions, gamma, weights, 1, options=pair_options,
        theta=1.0e-4, rsmooth=SMOOTH_RMIN_ARCMIN,
    )
    for name in ("xi_plus", "xi_minus", "xi_weight"):
        scale = max(float(np.max(np.abs(expected[name]))), np.finfo(float).tiny)
        relative_error = float(
            np.max(np.abs(pair_result[name] - expected[name]))/scale
        )
        if relative_error > 2.0e-4:
            raise AssertionError(
                f"smoothed only-2pcf accepted-node error is too large for "
                f"{name}: {relative_error}"
            )
    smooth_error = np.max(np.abs(pair_result["xi_plus"] - expected["xi_plus"]))
    raw_error = np.max(np.abs(pair_result["xi_plus"] - raw["xi_plus"]))
    if not smooth_error < raw_error:
        raise AssertionError("smoothed only-2pcf used raw-pivot dual-tree fields")

    exact_pair_options = f"{SMOOTH_OPTIONS},only-2pcf"
    single_catalog_pair = run_native(
        positions, gamma, weights, 1, options=exact_pair_options,
        theta=0.05, rsmooth=SMOOTH_RMIN_ARCMIN,
    )
    two_catalog_pair = run_native_catalogs(
        [(positions, gamma, weights),
         (positions.copy(), gamma.copy(), weights.copy())],
        1, options=exact_pair_options, theta=0.05,
        rsmooth=SMOOTH_RMIN_ARCMIN,
    )
    for name, values in single_catalog_pair.items():
        np.testing.assert_allclose(
            two_catalog_pair[name], values, rtol=5.0e-12, atol=5.0e-13,
            err_msg=f"repeated tree build changed explicit rsmooth: {name}",
        )


def test_spherical_smooth_radius_contract():
    positions, gamma, weights = fixture(18)
    parameters = {
        "searchMethod": ENGINE,
        "iCatalogs": "1",
        "usePeriodic": "false",
        "useLogHist": "false",
        "rangeN": RMAX,
        "rminHist": RMIN,
        "sizeHistN": BINS,
        "sizeHistPhi": PHI_BINS,
        "mChebyshev": NMAX,
        "lengthBox": 2.2,
        "numberThreads": 1,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-radius-"),
        "options": SMOOTH_OPTIONS,
        "rsmooth": str(SMOOTH_RMIN_ARCMIN),
    }
    model = cballs()
    model.set(parameters)
    model.set_catalog(
        positions, weights=weights, gamma1=gamma.real, gamma2=gamma.imag,
    )
    try:
        model.Run(level=["MainLoop"])
        expected = 2.0*np.sin(
            0.5*np.deg2rad(SMOOTH_RMIN_ARCMIN/60.0)
        )
        np.testing.assert_allclose(
            model.getrsmooth(), expected, rtol=2.0e-12, atol=5.0e-15,
            err_msg="explicit rsmooth was rescaled by an unrelated tree theta",
        )
    finally:
        model.struct_cleanup()

    unsafe = cballs()
    parameters["rootDir"] = tempfile.mkdtemp(
        prefix="ctreeballs-shear-radius-negative-"
    )
    parameters["rsmooth"] = "60"
    unsafe.set(parameters)
    unsafe.set_catalog(
        positions, weights=weights, gamma1=gamma.real, gamma2=gamma.imag,
    )
    try:
        unsafe.Run(level=["MainLoop"])
    except CosmoComputationError as error:
        if "safe limit 0.5*rminHist" not in str(error):
            raise
    else:
        raise AssertionError("a smooth group spanning measured pairs was accepted")
    finally:
        unsafe.struct_cleanup()


def test_spherical_mask_equivalence():
    positions, gamma, weights = fixture(57)
    mask = np.ones(positions.shape[0], dtype=np.uint8)
    mask[::4] = 0
    masked_options = f"{OPTIONS},read-mask"
    masked = run_native(
        positions, gamma, weights, 1, options=masked_options,
        theta=0.05, mask=mask,
    )
    filtered = run_native(
        positions[mask != 0], gamma[mask != 0], weights[mask != 0], 1,
        options=OPTIONS, theta=0.05,
    )
    many_threads = run_native(
        positions, gamma, weights, min(4, os.cpu_count() or 1),
        options=masked_options, theta=0.05, mask=mask,
    )
    for name, values in filtered.items():
        np.testing.assert_allclose(
            masked[name], values, rtol=5.0e-12, atol=5.0e-13,
            err_msg=f"read-mask differs from filtered catalog: {name}",
        )
        if not np.array_equal(masked[name], many_threads[name]):
            raise AssertionError(
                f"masked result changed with OpenMP threads: {name}"
            )


def test_spherical_contract_failures():
    positions, gamma, weights = fixture(12)
    model = cballs()
    model.set({
        "searchMethod": ENGINE,
        "iCatalogs": "1",
        "usePeriodic": "false",
        "useLogHist": "false",
        "rangeN": 2.1,
        "rminHist": RMIN,
        "sizeHistN": BINS,
        "sizeHistPhi": PHI_BINS,
        "mChebyshev": NMAX,
        "lengthBox": 2.2,
        "numberThreads": 1,
        "verbose": 0,
        "verbose_log": 0,
        "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-sphere-negative-"),
        "options": OPTIONS,
    })
    model.set_catalog(
        positions, weights=weights, gamma1=gamma.real, gamma2=gamma.imag,
    )
    try:
        model.Run(level=["MainLoop"])
    except CosmoComputationError as error:
        if "cannot exceed 2" not in str(error):
            raise
    else:
        raise AssertionError("spherical rangeN > 2 was accepted")
    finally:
        model.struct_cleanup()


if __name__ == "__main__":
    test_spherical_oracle_and_determinism()
    test_spherical_accepted_cell_transport()
    test_spherical_octant_frontier_parallelism()
    test_spherical_runtime_order_switches()
    test_spherical_smooth_pivot_transport()
    test_spherical_smooth_radius_contract()
    test_spherical_mask_equivalence()
    test_spherical_contract_failures()
    print(
        f"PASS: {ENGINE} full-sky spin-2 oracle, smoothing, "
        "contracts, and determinism"
    )
