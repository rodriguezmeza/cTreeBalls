#!/usr/bin/env python3
"""Public convergence benchmark entry point.

Uses the CLI and timing reports of kappa_corr_all_engines.py. Private
multi-backend benchmark environments are not dependencies of this checkout.
Original convergence benchmark by Axel Romero Tisnado.
"""
from kappa_corr_all_engines import main


if __name__ == "__main__":
    raise SystemExit(main())
