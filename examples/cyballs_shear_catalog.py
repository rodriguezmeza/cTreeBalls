#!/usr/bin/env python3
"""Compute flat-sky shear 2PCFs and 3PCFs from an in-memory NumPy catalog."""

import tempfile

import numpy as np

from cyballs import cballs


def main():
    count = 64
    angle = np.linspace(0.0, 2.0*np.pi, count, endpoint=False)
    radius = 0.35 + 0.5*((np.arange(count) % 5)/4.0)
    positions = np.zeros((count, 3), dtype=np.float64)
    positions[:, 0] = radius*np.cos(angle)
    positions[:, 1] = radius*np.sin(angle)
    gamma = (0.025 + 0.005*np.cos(2.0*angle))*np.exp(2j*angle)
    weights = 0.8 + 0.4*((np.arange(count) % 7)/6.0)

    model = cballs()
    model.set(
        {
            "searchMethod": "octree-shear-omp",
            "iCatalogs": "1",
            "usePeriodic": "false",
            "useLogHist": "false",
            "rminHist": 0.01,
            "rangeN": 2.0,
            "sizeHistN": 8,
            "sizeHistPhi": 16,
            "mChebyshev": 4,
            "lengthBox": 2.0,
            "numberThreads": 4,
            "rootDir": tempfile.mkdtemp(prefix="ctreeballs-shear-"),
            "verbose": 0,
            "verbose_log": 0,
            "options": "no-out-Hist",
        }
    )
    model.set_catalog(
        positions, weights=weights,
        gamma1=gamma.real, gamma2=gamma.imag,
    )

    try:
        model.Run(level=["MainLoop"])
        xi_plus = model.getShearXiPlus()
        xi_minus = model.getShearXiMinus()
        gamma_multipoles = model.getShearGammaXMultipoles()
        gamma_angular = model.getShearGammaX()
        print("xi+:", xi_plus)
        print("xi-:", xi_minus)
        print("Gamma^x multipoles shape:", gamma_multipoles.shape)
        print("Gamma^x angular shape:", gamma_angular.shape)
    finally:
        model.clean_all()


if __name__ == "__main__":
    main()
