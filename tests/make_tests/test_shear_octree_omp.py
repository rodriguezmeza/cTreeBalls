#!/usr/bin/env python3
"""Independent oracle and OpenMP regression tests for octree-shear-omp."""

from __future__ import annotations

import math
import os
import tempfile

import numpy as np

from cyballs import CosmoComputationError, CosmoSevereError, cballs


NMAX = 2
BINS = 2
PHI_BINS = 16
RMIN = 0.01
RMAX = 2.0


def fixture():
    index = np.arange(48, dtype=np.float64)
    angle = 2.0*np.pi*index/index.size + 0.013*np.sin(3.0*index)
    radius = 0.28 + 0.58*((index.astype(np.int64) % 5)/4.0)
    positions = np.zeros((index.size, 3), dtype=np.float64)
    positions[:, 0] = radius*np.cos(angle)
    positions[:, 1] = radius*np.sin(angle)
    gamma = (0.025 + 0.004*np.cos(2.0*angle)) \
        * np.exp(2j*(angle + 0.17*np.sin(angle)))
    weights = 0.7 + 0.6*((index.astype(np.int64) % 7)/6.0)
    return positions, gamma, weights


def radial_bin(distance):
    if not RMIN < distance < RMAX:
        return None
    result = int((distance - RMIN)/(RMAX - RMIN)*BINS)
    return result if 0 <= result < BINS else None


def oracle(positions, gamma, weights):
    multipoles = 2*NMAX + 1
    denominator_modes = 4*NMAX + 1
    ring_max = max(2*NMAX, NMAX + 3)
    xi_plus = np.zeros(BINS, dtype=np.complex128)
    xi_minus = np.zeros(BINS, dtype=np.complex128)
    xi_weight = np.zeros(BINS, dtype=np.float64)
    upsilon = np.zeros((4, multipoles, BINS, BINS), dtype=np.complex128)
    window = np.zeros((denominator_modes, BINS, BINS), dtype=np.complex128)

    for pivot in range(positions.shape[0]):
        rings = np.zeros((2*ring_max + 1, BINS), dtype=np.complex128)
        weight_rings = np.zeros_like(rings)
        diagonal_g6 = np.zeros(BINS, dtype=np.complex128)
        diagonal_g2 = np.zeros(BINS, dtype=np.complex128)
        diagonal_abs2 = np.zeros(BINS, dtype=np.complex128)
        diagonal_w2 = np.zeros(BINS, dtype=np.float64)

        for neighbor in range(positions.shape[0]):
            if neighbor == pivot:
                continue
            separation = positions[neighbor, :2] - positions[pivot, :2]
            distance = np.linalg.norm(separation)
            bin_index = radial_bin(distance)
            if bin_index is None:
                continue
            phase = complex(*(separation/distance))
            weighted_gamma = weights[neighbor]*gamma[neighbor]
            for order in range(-ring_max, ring_max + 1):
                rings[order + ring_max, bin_index] += weighted_gamma*phase**order
                weight_rings[order + ring_max, bin_index] += \
                    weights[neighbor]*phase**order

            diagonal_g6[bin_index] += weighted_gamma**2*phase**-6
            diagonal_g2[bin_index] += weighted_gamma**2*phase**-2
            diagonal_abs2[bin_index] += \
                weights[neighbor]**2*abs(gamma[neighbor])**2*phase**-2
            diagonal_w2[bin_index] += weights[neighbor]**2
            pivot_gamma = weights[pivot]*gamma[pivot]
            xi_plus[bin_index] += pivot_gamma*np.conj(weighted_gamma)
            xi_minus[bin_index] += pivot_gamma*weighted_gamma*phase**-4
            xi_weight[bin_index] += weights[pivot]*weights[neighbor]

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
                    window[order + 2*NMAX, first, second] += weights[pivot]*value

                for order in range(-NMAX, NMAX + 1):
                    products = np.array(
                        [
                            ring(order - 3, first)*ring(-order - 3, second),
                            ring(order - 1, first)*ring(-order - 1, second),
                            np.conj(ring(-order - 1, first))
                            * ring(-order - 3, second),
                            ring(order - 3, first)
                            * np.conj(ring(order - 1, second)),
                        ],
                        dtype=np.complex128,
                    )
                    if same_bin:
                        products -= np.array(
                            [
                                diagonal_g6[first],
                                diagonal_g2[first],
                                diagonal_abs2[first],
                                diagonal_abs2[first],
                            ]
                        )
                    pivot_factors = weights[pivot]*np.array(
                        [gamma[pivot], np.conj(gamma[pivot]),
                         gamma[pivot], gamma[pivot]]
                    )
                    upsilon[:, order + NMAX, first, second] -= \
                        pivot_factors*products

    xi_plus = np.divide(xi_plus, xi_weight, out=np.zeros_like(xi_plus),
                        where=xi_weight != 0.0)
    xi_minus = np.divide(xi_minus, xi_weight, out=np.zeros_like(xi_minus),
                         where=xi_weight != 0.0)

    corrected = np.zeros_like(upsilon)
    solved_bins = 0
    for first in range(BINS):
        for second in range(BINS):
            n_zero = window[2*NMAX, first, second]
            if abs(n_zero) == 0.0:
                continue
            coupling = np.empty((multipoles, multipoles), dtype=np.complex128)
            for row, ell in enumerate(range(-NMAX, NMAX + 1)):
                for column, order in enumerate(range(-NMAX, NMAX + 1)):
                    coupling[row, column] = \
                        window[ell - order + 2*NMAX, first, second]/n_zero
            if np.linalg.cond(coupling) > 1.0e10:
                continue
            corrected[:, :, first, second] = np.linalg.solve(
                coupling, (upsilon[:, :, first, second]/n_zero).T
            ).T
            solved_bins += 1
    if solved_bins == 0:
        raise AssertionError("oracle fixture produced no invertible mode system")

    delta_phi = 2.0*np.pi/PHI_BINS
    phi = -np.pi + (np.arange(PHI_BINS) + 0.5)*delta_phi
    angular = np.zeros((4, PHI_BINS, BINS, BINS), dtype=np.complex128)
    for order in range(-NMAX, NMAX + 1):
        window_factor = np.sinc(order*delta_phi/(2.0*np.pi))
        angular += corrected[:, order + NMAX, None, :, :] \
            * np.exp(1j*order*phi)[None, :, None, None] \
            * window_factor/(2.0*np.pi)

    return {
        "xi_plus": xi_plus,
        "xi_minus": xi_minus,
        "xi_weight": xi_weight,
        "upsilon": upsilon,
        "window": window,
        "multipoles": corrected,
        "angular": angular,
    }


def direct_triplet_oracle(catalogs):
    """Enumerate distinct objects without using the factorized ring identity."""
    pivot_catalog = catalogs[0]
    first_catalog = catalogs[1] if len(catalogs) >= 2 else pivot_catalog
    second_catalog = catalogs[2] if len(catalogs) >= 3 else first_catalog
    first_is_pivot = len(catalogs) < 2
    second_is_pivot = len(catalogs) < 2
    same_neighbor_catalog = len(catalogs) < 3
    multipoles = 2*NMAX + 1
    denominator_modes = 4*NMAX + 1
    xi_plus = np.zeros(BINS, dtype=np.complex128)
    xi_minus = np.zeros(BINS, dtype=np.complex128)
    xi_weight = np.zeros(BINS, dtype=np.float64)
    upsilon = np.zeros((4, multipoles, BINS, BINS), dtype=np.complex128)
    window = np.zeros((denominator_modes, BINS, BINS), dtype=np.complex128)
    p_pos, p_gamma, p_weight = pivot_catalog
    f_pos, f_gamma, f_weight = first_catalog
    s_pos, s_gamma, s_weight = second_catalog

    for pivot in range(p_pos.shape[0]):
        for first in range(f_pos.shape[0]):
            if first_is_pivot and first == pivot:
                continue
            first_separation = f_pos[first, :2] - p_pos[pivot, :2]
            first_distance = np.linalg.norm(first_separation)
            first_bin = radial_bin(first_distance)
            if first_bin is None:
                continue
            first_phase = complex(*(first_separation/first_distance))
            weighted_pivot = p_weight[pivot]*p_gamma[pivot]
            weighted_first = f_weight[first]*f_gamma[first]
            pair_weight = p_weight[pivot]*f_weight[first]
            xi_plus[first_bin] += weighted_pivot*np.conj(weighted_first)
            xi_minus[first_bin] += (
                weighted_pivot*weighted_first*first_phase**-4
            )
            xi_weight[first_bin] += pair_weight

            for second in range(s_pos.shape[0]):
                if second_is_pivot and second == pivot:
                    continue
                if same_neighbor_catalog and second == first:
                    continue
                second_separation = s_pos[second, :2] - p_pos[pivot, :2]
                second_distance = np.linalg.norm(second_separation)
                second_bin = radial_bin(second_distance)
                if second_bin is None:
                    continue
                second_phase = complex(*(second_separation/second_distance))
                weighted_second = s_weight[second]*s_gamma[second]
                triplet_weight = p_weight[pivot]*f_weight[first]*s_weight[second]

                for order in range(-2*NMAX, 2*NMAX + 1):
                    window[order + 2*NMAX, first_bin, second_bin] += (
                        triplet_weight*first_phase**order
                        * np.conj(second_phase**order)
                    )

                for order in range(-NMAX, NMAX + 1):
                    products = np.array(
                        [
                            weighted_first*first_phase**(order - 3)
                            * weighted_second*second_phase**(-order - 3),
                            weighted_first*first_phase**(order - 1)
                            * weighted_second*second_phase**(-order - 1),
                            np.conj(
                                weighted_first*first_phase**(-order - 1)
                            )
                            * weighted_second*second_phase**(-order - 3),
                            weighted_first*first_phase**(order - 3)
                            * np.conj(
                                weighted_second*second_phase**(order - 1)
                            ),
                        ],
                        dtype=np.complex128,
                    )
                    pivot_factors = p_weight[pivot]*np.array(
                        [p_gamma[pivot], np.conj(p_gamma[pivot]),
                         p_gamma[pivot], p_gamma[pivot]],
                        dtype=np.complex128,
                    )
                    upsilon[:, order + NMAX, first_bin, second_bin] -= (
                        pivot_factors*products
                    )

    xi_plus = np.divide(xi_plus, xi_weight, out=np.zeros_like(xi_plus),
                        where=xi_weight != 0.0)
    xi_minus = np.divide(xi_minus, xi_weight, out=np.zeros_like(xi_minus),
                         where=xi_weight != 0.0)
    corrected = np.zeros_like(upsilon)
    for first_bin in range(BINS):
        for second_bin in range(BINS):
            n_zero = window[2*NMAX, first_bin, second_bin]
            if abs(n_zero) == 0.0:
                continue
            coupling = np.empty((multipoles, multipoles), dtype=np.complex128)
            for row, ell in enumerate(range(-NMAX, NMAX + 1)):
                for column, order in enumerate(range(-NMAX, NMAX + 1)):
                    coupling[row, column] = window[
                        ell - order + 2*NMAX, first_bin, second_bin
                    ]/n_zero
            if np.linalg.cond(coupling) > 1.0e10:
                continue
            corrected[:, :, first_bin, second_bin] = np.linalg.solve(
                coupling,
                (upsilon[:, :, first_bin, second_bin]/n_zero).T,
            ).T

    delta_phi = 2.0*np.pi/PHI_BINS
    phi = -np.pi + (np.arange(PHI_BINS) + 0.5)*delta_phi
    angular = np.zeros((4, PHI_BINS, BINS, BINS), dtype=np.complex128)
    for order in range(-NMAX, NMAX + 1):
        window_factor = np.sinc(order*delta_phi/(2.0*np.pi))
        angular += corrected[:, order + NMAX, None, :, :] \
            * np.exp(1j*order*phi)[None, :, None, None] \
            * window_factor/(2.0*np.pi)
    return {
        "xi_plus": xi_plus,
        "xi_minus": xi_minus,
        "xi_weight": xi_weight,
        "upsilon": upsilon,
        "window": window,
        "multipoles": corrected,
        "angular": angular,
    }


def configured_model(catalogs, threads, *, periodic=False):
    balls = cballs()
    catalog_selection = ",".join(str(index + 1)
                                 for index in range(len(catalogs)))
    balls.set(
        {
            "searchMethod": "octree-shear-omp",
            "iCatalogs": catalog_selection,
            "usePeriodic": str(periodic).lower(),
            "useLogHist": "false",
            "rangeN": RMAX,
            "rminHist": RMIN,
            "sizeHistN": BINS,
            "sizeHistPhi": PHI_BINS,
            "mChebyshev": NMAX,
            "lengthBox": 2.2,
            "numberThreads": threads,
            "verbose": 0,
            "verbose_log": 0,
            "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-output-"),
            "options": "no-out-Hist",
        }
    )
    for catalog_index, (positions, gamma, weights) in enumerate(catalogs):
        balls.set_catalog(
            positions, weights=weights,
            gamma1=gamma.real, gamma2=gamma.imag,
            catalog=catalog_index,
        )
    return balls


def run_catalogs(catalogs, threads):
    balls = configured_model(catalogs, threads)
    try:
        balls.Run(level=["MainLoop"])
        upsilon_x = balls.getShearUpsilonXMultipoles().copy()
        multipoles_x = balls.getShearGammaXMultipoles().copy()
        angular_x = balls.getShearGammaX().copy()
        if not np.array_equal(upsilon_x,
                              balls.getShearUpsilonMultipoles()):
            raise AssertionError("Upsilon x-projection alias disagrees")
        if not np.array_equal(multipoles_x,
                              balls.getShearGammaMultipoles()):
            raise AssertionError("Gamma multipole x-projection alias disagrees")
        if not np.array_equal(angular_x, balls.getShearGamma()):
            raise AssertionError("Gamma angular x-projection alias disagrees")
        return {
            "xi_plus": balls.getShearXiPlus().copy(),
            "xi_minus": balls.getShearXiMinus().copy(),
            "xi_weight": balls.getShearXiWeight().copy(),
            "upsilon": upsilon_x,
            "window": balls.getShearWindowMultipoles().copy(),
            "multipoles": multipoles_x,
            "angular": angular_x,
        }
    finally:
        balls.struct_cleanup()


def run_once(positions, gamma, weights, threads):
    return run_catalogs([(positions, gamma, weights)], threads)


def assert_close(actual, expected):
    for name in expected:
        np.testing.assert_allclose(
            actual[name], expected[name], rtol=3.0e-12, atol=3.0e-13,
            err_msg=name,
        )


def test_oracle_and_openmp_determinism():
    positions, gamma, weights = fixture()
    expected = oracle(positions, gamma, weights)
    one_thread = run_once(positions, gamma, weights, 1)
    assert_close(one_thread, expected)
    direct = direct_triplet_oracle([(positions, gamma, weights)])
    assert_close({name: one_thread[name] for name in direct}, direct)

    many_threads = run_once(
        positions, gamma, weights, min(4, os.cpu_count() or 1)
    )
    for name in one_thread:
        if not np.array_equal(one_thread[name], many_threads[name]):
            raise AssertionError(f"OpenMP result is not deterministic: {name}")


def test_active_spin2_rotation_invariance():
    positions, gamma, weights = fixture()
    reference = run_once(positions, gamma, weights, 1)
    angle = 0.731
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle)],
         [np.sin(angle), np.cos(angle)]],
        dtype=np.float64,
    )
    rotated_positions = positions.copy()
    rotated_positions[:, :2] = positions[:, :2] @ rotation.T
    rotated_gamma = gamma*np.exp(2j*angle)
    rotated = run_once(rotated_positions, rotated_gamma, weights, 1)
    assert_close(rotated, reference)


def tomography_fixture():
    positions, gamma, weights = fixture()
    catalogs = []
    for catalog_index in range(3):
        selection = np.arange(catalog_index, positions.shape[0], 3)
        catalog_positions = positions[selection].copy()
        catalog_positions[:, 0] += 0.017*catalog_index
        catalog_positions[:, 1] -= 0.011*catalog_index
        catalog_gamma = gamma[selection]*np.exp(0.07j*catalog_index)
        catalog_weights = weights[selection]*(1.0 + 0.05*catalog_index)
        catalogs.append((catalog_positions, catalog_gamma, catalog_weights))
    return catalogs


def test_three_catalog_tomography_matches_direct_triplets():
    catalogs = tomography_fixture()
    two_catalog_expected = direct_triplet_oracle(catalogs[:2])
    two_catalog_actual = run_catalogs(catalogs[:2], 1)
    assert_close(two_catalog_actual, two_catalog_expected)

    expected = direct_triplet_oracle(catalogs)
    actual = run_catalogs(catalogs, min(3, os.cpu_count() or 1))
    assert_close(actual, expected)


def test_getter_and_geometry_failures_are_recoverable():
    positions, gamma, weights = fixture()
    balls = cballs()
    try:
        balls.getShearGamma()
    except CosmoSevereError:
        pass
    else:
        raise AssertionError("shear getter accepted an uncomputed object")

    balls.set(
        {
            "searchMethod": "octree-shear-omp",
            "iCatalogs": "1",
            "usePeriodic": "true",
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
            "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-negative-"),
            "options": "no-out-Hist",
        }
    )
    balls.set_catalog(
        positions, weights=weights,
        gamma1=gamma.real, gamma2=gamma.imag,
    )
    try:
        balls.Run(level=["MainLoop"])
    except CosmoComputationError as error:
        if "periodic geometry" not in str(error):
            raise
    else:
        raise AssertionError("periodic flat-sky shear run was accepted")
    finally:
        balls.struct_cleanup()

    first_catalog = (positions[:16].copy(), gamma[:16].copy(),
                     weights[:16].copy())
    shifted_positions = positions[16:32].copy()
    shifted_positions[:, 2] = 0.125
    second_catalog = (shifted_positions, gamma[16:32].copy(),
                      weights[16:32].copy())
    balls = configured_model([first_catalog, second_catalog], 1)
    try:
        balls.Run(level=["MainLoop"])
    except CosmoComputationError as error:
        if "do not share a tangent plane" not in str(error):
            raise
    else:
        raise AssertionError("cross-catalog tangent-plane mismatch was accepted")
    finally:
        balls.struct_cleanup()


if __name__ == "__main__":
    test_oracle_and_openmp_determinism()
    test_active_spin2_rotation_invariance()
    test_three_catalog_tomography_matches_direct_triplets()
    test_getter_and_geometry_failures_are_recoverable()
    print("PASS: shear Gamma^x oracle, spin-2 rotation, tomography, recovery, "
          "and OpenMP determinism")
