"""Optional, isolated adapter to an unmodified external lya2pcf CPU checkout.

The 3D product uses the upstream estimator. Radial and same-LOS products are
explicit kernel adaptations, not features claimed for the upstream pipeline.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np


REFERENCE_NAMES = {
    "3d_2pcf": "lya2pcf-cpu",
    "1d_2pcf": "lya2pcf-kernel-radial",
    "1d_same_los_2pcf": "lya2pcf-kernel-same-los",
}
CONTRACTS = {
    "3d_2pcf": "upstream CPU pair kernel; distinct unordered forest pairs; actual sky geometry",
    "1d_2pcf": "adaptation: upstream pair kernel on collinear sightlines; distinct forest IDs",
    "1d_same_los_2pcf": (
        "adaptation: upstream collinear pair kernel on disjoint within-forest pixel slices; "
        "unordered pairs, no self pairs; equal mean of occupied forest/bin correlations"
    ),
}


def product_table(numerator, denominator, config, family):
    xi = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    index = np.indices(numerator.shape).reshape(numerator.ndim, -1).T
    centers = (index + .5) * np.array(
        [config["rp_max"]/config["rp_bins"]] +
        ([config["rt_max"]/config["rt_bins"]] if family == "3d_2pcf" else []))
    return np.column_stack((index, centers, xi.ravel(), numerator.ravel(), denominator.ravel()))


def run_references(catalog, config, families):
    """One fresh process per binning: Numba captures upstream module globals."""
    source = Path(config.lya2pcf_source).expanduser().resolve()
    for name in ("parameters.py", "correlation_procedures_cpu.py", "post_processing.py"):
        if not (source / name).is_file():
            raise ValueError(f"--lya2pcf-source: missing {source/name}")
    output = config.output_dir / "lya2pcf_reference"
    output.mkdir(exist_ok=False)
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in vars(config).items()}
    settings_path = output / "settings.json"
    settings_path.write_text(json.dumps(settings, indent=2) + "\n")
    result = {}
    with tempfile.TemporaryDirectory(prefix="lya-reference-") as temporary:
        input_path = Path(temporary) / "pixels.npz"
        np.savez(input_path, positions=catalog.positions, delta=catalog.delta,
                 weights=catalog.weights, forest_ids=catalog.forest_ids)
        for family in sorted(set(families) & REFERENCE_NAMES.keys()):
            print(f"Running {REFERENCE_NAMES[family]}: {CONTRACTS[family]}", flush=True)
            started = time.perf_counter()
            target = output / family
            command = [sys.executable, str(Path(__file__).resolve()),
                       "--source", str(source), "--catalog", str(input_path),
                       "--settings", str(settings_path), "--family", family,
                       "--output", str(target)]
            # MPI rank zero launches a non-MPI worker; no global import mutations.
            with target.with_suffix(".log").open("w") as log:
                status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                        timeout=config.reference_timeout or None, check=False)
            if status.returncode:
                raise RuntimeError(f"lya2pcf failed; see {target.with_suffix('.log')}:\n" +
                                   target.with_suffix(".log").read_text()[-3000:])
            metadata = json.loads(target.with_suffix(".json").read_text())
            with np.load(target.with_suffix(".npz"), allow_pickle=False) as data:
                table = data["table"]
            elapsed = time.perf_counter() - started
            result[REFERENCE_NAMES[family]] = dict(
                wall_seconds=metadata["compute_seconds"], worker_wall_seconds=elapsed,
                products={family: table}, ranks=1, smooth_pivot="unsupported",
                provenance=metadata, reference_archive=str(target.with_suffix(".npz")))
            print(f"  completed: pair traversal {metadata['compute_seconds']:.6f}s; "
                  f"worker including imports/JIT/analysis {elapsed:.6f}s", flush=True)
    return result


def compute_reference(source, catalog, config, family):
    """Called only in the disposable worker, or an isolated integration test."""
    started = time.perf_counter()
    sys.path.insert(0, str(source))
    cpu = importlib.import_module("correlation_procedures_cpu")
    if Path(cpu.__file__).resolve().parent != source.resolve():
        raise RuntimeError("reference module was not loaded from the requested checkout")
    radial = family != "3d_2pcf"
    shape = (config["rp_bins"], 1 if radial else config["rt_bins"])
    cpu.rpmax, cpu.rtmax = config["rp_max"], config["rt_max"]
    cpu.numpix_rp, cpu.numpix_rt = shape
    cpu.shape_hist = shape
    positions, delta, weights, identifiers = (
        catalog[k] for k in ("positions", "delta", "weights", "forest_ids"))
    distances = np.hypot(np.hypot(positions[:, 0], positions[:, 1]), positions[:, 2])
    order = np.argsort(identifiers, kind="stable")
    splits = np.flatnonzero(np.diff(identifiers[order]) != 0) + 1
    forests = []
    for indices in np.split(order, splits):
        indices = indices[weights[indices] > 0]
        if not len(indices):
            continue
        unit = positions[indices[0]] / distances[indices[0]]
        if not radial and not np.allclose(positions[indices]/distances[indices, None],
                                         unit, rtol=0, atol=1e-10):
            raise ValueError("lya2pcf 3D requires each forest ID to lie on one observer sightline")
        ra = float(np.arctan2(unit[1], unit[0]) % (2*np.pi))
        dec = float(np.arcsin(np.clip(unit[2], -1, 1)))
        forests.append((ra, dec, np.ascontiguousarray(weights[indices]),
                        np.ascontiguousarray(weights[indices]*delta[indices]),
                        np.ascontiguousarray(distances[indices]), int(identifiers[indices[0]])))
    forests.sort(key=lambda f: (f[0], f[5]))
    if not forests:
        raise ValueError("no positive-weight forests")
    angmax = np.nextafter(np.pi, np.inf)

    def pair(left, right):
        ra1, dec1, w1, dw1, dc1, _ = left
        ra2, dec2, w2, dw2, dc2, _ = right
        if radial:
            ra1 = dec1 = ra2 = dec2 = 0.0
        elif abs(ra1-ra2) >= cpu.chiquito or abs(dec1-dec2) >= cpu.chiquito:
            cosine = (math.sin(dec1)*math.sin(dec2) +
                      math.cos(dec1)*math.cos(dec2)*math.cos(ra1-ra2))
            if not -1 <= cosine <= 1:
                raise ValueError("upstream angular arccos is out of range for a forest pair; "
                                 "refusing to silently drop its NaN-angle histogram")
        return cpu.pair_correlation(angmax, ra1, dec1, w1, dw1, 0, dc1, 0,
                                    ra2, dec2, w2, dw2, 0, dc2, 0)

    # Warm only the exact signature that the threaded traversal will use.
    empty = np.empty(0, dtype=np.float64)
    dummy = (0.0, 0.0, empty, empty, empty, 0)
    jit_started = time.perf_counter()
    pair(dummy, dummy)
    jit_seconds = time.perf_counter() - jit_started
    tree = None
    if not radial:
        from scipy.spatial import cKDTree
        ra = np.array([f[0] for f in forests])
        dec = np.array([f[1] for f in forests])
        units = np.column_stack((np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)))
        tree = cKDTree(units)
        radius = min(2.0, config["rt_max"]/float(distances[weights > 0].min()))
        radius += 1e-9  # Conservative guard for upstream small-angle geometry.
    covariance_requested = config["reference_covariance"] and not radial
    regions = np.zeros(len(forests), dtype=np.int64)
    if covariance_requested:
        import healpy as hp
        regions = hp.ang2pix(config["reference_nside"], np.pi/2-dec, ra, nest=False)
    region_ids, region_index = np.unique(regions, return_inverse=True)
    cells = math.prod(shape)
    arrays_bytes = 16*len(region_ids)*cells + 16*(2*config["threads"]+2)*cells
    if covariance_requested:
        arrays_bytes += 32*cells*cells
    if arrays_bytes > config["max_hist_mib"]*2**20:
        raise ValueError("reference histogram/covariance workspace exceeds --max-hist-mib")
    partial_num = np.zeros((len(region_ids), *shape))
    partial_den = np.zeros_like(partial_num)

    def visit(i):
        left = forests[i]
        num, den = np.zeros(shape), np.zeros(shape)
        if family == "1d_same_los_2pcf":
            # Disjoint slices avoid both diagonal self pairs and cancellation
            # from subtracting self weights out of the zero-lag bin.
            for p in range(len(left[2])-1):
                one = (*left[:2], *(a[p:p+1] for a in left[2:5]), left[5])
                rest = (*left[:2], *(a[p+1:] for a in left[2:5]), left[5])
                w, dw = pair(one, rest)
                den += w
                num += dw
            num = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
            den = (den > 0).astype(float)
        else:
            neighbors = (range(i+1, len(forests)) if radial else
                         sorted(j for j in tree.query_ball_point(units[i], radius) if j > i))
            for j in neighbors:
                w, dw = pair(left, forests[j])
                den += w
                num += dw
        return i, num, den

    setup_seconds = time.perf_counter() - started
    compute_started = time.perf_counter()
    last_progress = compute_started
    with ThreadPoolExecutor(max_workers=config["threads"]) as pool:
        # Bounded batches prevent completed dense histograms accumulating for
        # the entire survey. Reduction order stays independent of thread count.
        step = max(1, 2*config["threads"])
        for first in range(0, len(forests), step):
            for i, num, den in pool.map(visit, range(first, min(first+step, len(forests)))):
                partial_num[region_index[i]] += num
                partial_den[region_index[i]] += den
            if time.perf_counter() - last_progress > 5:
                print(f"forests completed {min(first+step, len(forests))}/{len(forests)}", flush=True)
                last_progress = time.perf_counter()
    compute_seconds = time.perf_counter() - compute_started
    numerator, denominator = partial_num.sum(axis=0), partial_den.sum(axis=0)
    payload = dict(table=product_table(numerator[:, 0] if radial else numerator,
                                      denominator[:, 0] if radial else denominator, config, family))
    covariance_status = "not requested for this product"
    analysis_started = time.perf_counter()
    if covariance_requested:
        payload.update(region_ids=region_ids, partial_numerator=partial_num,
                       partial_denominator=partial_den)
        occupied_regions = np.any(partial_den > 0, axis=(1, 2))
        if occupied_regions.sum() < 2:
            covariance_status = "unavailable: fewer than two occupied sky regions"
        else:
            post = importlib.import_module("post_processing")
            occupied = denominator.ravel() > 0
            w = np.ascontiguousarray(partial_den.reshape(len(region_ids), -1)[:, occupied])
            n = partial_num.reshape(len(region_ids), -1)[:, occupied]
            xi = np.divide(n, w, out=np.zeros_like(n), where=w > 0)
            covariance = np.zeros((cells, cells))
            if occupied.any():
                occupied_covariance, means = post.cov(xi, w)
                covariance[np.ix_(occupied, occupied)] = occupied_covariance
                np.testing.assert_allclose(means, numerator.ravel()[occupied]/denominator.ravel()[occupied],
                                           rtol=1e-12, atol=1e-14)
            payload["covariance"] = covariance
            covariance_status = "upstream post_processing.cov; unsmoothed weighted HEALPix subsampling"
    metadata = dict(contract=CONTRACTS[family], source=str(source), threads=config["threads"],
                    kernel_sha256=hashlib.sha256(Path(cpu.__file__).read_bytes()).hexdigest(),
                    scheduling="RA then int64-ID ordering; each unordered forest pair once; cKDTree sky pruning",
                    setup_seconds=setup_seconds, jit_seconds=jit_seconds,
                    compute_seconds=compute_seconds,
                    analysis_seconds=time.perf_counter()-analysis_started,
                    covariance=covariance_status, regions=len(region_ids),
                    occupied_regions=int(np.any(partial_den > 0, axis=(1, 2)).sum()),
                    timing="CPU kernel traversal and reduction; excludes setup, JIT, covariance and I/O",
                    geometry="upstream small-angle approximation retained; boundary-bin differences can occur")
    return payload, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "catalog", "settings", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--family", choices=REFERENCE_NAMES, required=True)
    args = parser.parse_args()
    config = json.loads(args.settings.read_text())
    with np.load(args.catalog, allow_pickle=False) as catalog:
        payload, metadata = compute_reference(args.source, catalog, config, args.family)
    np.savez_compressed(args.output.with_suffix(".npz"), **payload)
    args.output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
