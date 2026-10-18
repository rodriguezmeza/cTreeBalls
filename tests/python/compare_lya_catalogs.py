#!/usr/bin/env python3
"""Compare saved cTreeBalls Ly-alpha 2PCF/3PCF results from several catalogs.

Example (run from run2):
  python3 python/compare_lya_catalogs.py Output_DESI Output_mock \
      --engine lya-2pcf-3pcf-omp --labels DESI Mock --reference DESI \
      --output Output_catalog_comparison

Inputs are driver output roots, individual engine directories, or native
histogram files. Only matching estimator families and physical grids are
compared. This script needs NumPy and Matplotlib, but does not import cyballs
or read the original FITS catalogs. See README_compare_lya_catalogs.txt.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys

import numpy as np


@dataclass(frozen=True)
class Product:
    filename: str
    indices: tuple[str, ...]
    coordinates: tuple[str, ...]
    coordinate_axes: tuple[int, ...]
    grid: tuple[str, ...]
    order: int

    @property
    def columns(self):
        return len(self.indices) + len(self.coordinates) + 3


PRODUCTS = {
    "3d_2pcf": Product("histXi2pcf_lya.txt", ("bp", "bt"),
                       ("rp_center", "rt_center"), (0, 1),
                       ("rp_max", "rt_max", "rp_bins", "rt_bins"), 2),
    "3d_3pcf": Product("histZetaM_lya5d.txt", ("b1", "b2", "t1", "t2", "bmu"),
                       ("r1", "r2", "theta1", "theta2", "mu"), (0, 1, 2, 3, 4),
                       ("r3_max", "r3_bins", "theta_bins", "mu_bins"), 3),
    "1d_2pcf": Product("histXi2pcf_lya1d.txt", ("bin",),
                       ("radial_separation",), (0,), ("rp_max", "rp_bins"), 2),
    "1d_3pcf": Product("histZetaM_lya1d.txt", ("b1", "b2"),
                       ("lag1", "lag2"), (0, 1), ("r3_max", "r3_bins"), 3),
    "1d_same_los_2pcf": Product("histXi2pcf_lya1d_same_los.txt", ("bin",),
                                ("radial_separation",), (0,),
                                ("rp_max", "rp_bins"), 2),
    "multipole_3pcf": Product("histZetaM_3d.txt", ("ell", "bin1", "bin2"),
                              ("r1_center", "r2_center"), (1, 2),
                              ("r3_max", "r3_bins", "multipole_rmin"), 3),
}


@dataclass
class Input:
    label: str
    directory: Path
    files: dict[str, Path]
    config: dict
    catalog: dict
    summary: Path | None


@dataclass
class Aligned:
    indices: np.ndarray
    coordinates: np.ndarray
    numerator: np.ndarray
    denominator: np.ndarray
    labels: list[str]
    index_names: tuple[str, ...]
    coordinate_names: tuple[str, ...]
    coordinate_axes: tuple[int, ...]

    @property
    def occupied(self):
        return self.denominator > 0

    @property
    def values(self):
        return np.divide(self.numerator, self.denominator,
                         out=np.full_like(self.numerator, np.nan),
                         where=self.occupied)


def resolve_input(path, label, engine=None):
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise ValueError(f"input does not exist: {path}")
    if path.is_file():
        files = {key: path for key, spec in PRODUCTS.items() if path.name == spec.filename}
        if not files:
            raise ValueError(f"{path}: expected a native histogram filename")
        directory = path.parent
    else:
        directory = path
        if engine and (path / engine).is_dir():
            directory = path / engine
        files = {key: directory / spec.filename for key, spec in PRODUCTS.items()
                 if (directory / spec.filename).is_file()}
        if not files:
            children = [p for p in directory.iterdir() if p.is_dir() and
                        any((p / s.filename).is_file() for s in PRODUCTS.values())]
            if engine:
                raise ValueError(f"{path}: no histogram for engine {engine!r}")
            if len(children) != 1:
                names = ", ".join(sorted(p.name for p in children)) or "none"
                raise ValueError(f"{path}: select --engine or pass an engine directory; found {names}")
            directory = children[0]
            files = {key: directory / spec.filename for key, spec in PRODUCTS.items()
                     if (directory / spec.filename).is_file()}
    summary = next((p / "summary.json" for p in (directory, directory.parent)
                    if (p / "summary.json").is_file()), None)
    metadata = {} if summary is None else json.loads(summary.read_text())
    if not isinstance(metadata, dict):
        raise ValueError(f"{summary}: expected a JSON object")
    config, catalog = metadata.get("config", {}), metadata.get("catalog", {})
    if not isinstance(config, dict) or not isinstance(catalog, dict):
        raise ValueError(f"{summary}: config and catalog must be JSON objects")
    return Input(label, directory, files, config, catalog, summary)


def read_table(path, spec):
    with path.open() as stream:
        has_rows = any(line.strip() and not line.lstrip().startswith("#") for line in stream)
    table = np.loadtxt(path, ndmin=2) if has_rows else np.empty((0, spec.columns))
    if table.shape[1] != spec.columns or not np.isfinite(table).all():
        raise ValueError(f"{path}: expected {spec.columns} finite columns")
    ndim = len(spec.indices)
    indices = table[:, :ndim]
    if np.any(indices < 0) or np.any(indices != np.floor(indices)) or np.any(indices > 2**53-1):
        raise ValueError(f"{path}: histogram indices must be nonnegative exact integers")
    if len(np.unique(indices, axis=0)) != len(table):
        raise ValueError(f"{path}: duplicate histogram bin indices")
    if np.any(table[:, -1] < 0):
        raise ValueError(f"{path}: negative denominator")
    ratio = np.divide(table[:, -2], table[:, -1], out=np.zeros(len(table)),
                      where=table[:, -1] > 0)
    # Multipole output is printed with 12 digits; forest output with 17.
    if not np.allclose(table[:, -3], ratio, rtol=1e-10, atol=1e-14):
        raise ValueError(f"{path}: correlation is inconsistent with numerator/denominator")
    if np.any((table[:, -1] == 0) & (table[:, -2] != 0)):
        raise ValueError(f"{path}: nonzero numerator in an empty bin")
    return table


def check_metadata(inputs, spec):
    notes = []
    for field in spec.grid:
        known = [(item.label, item.config[field]) for item in inputs if field in item.config]
        if known and any(not np.isclose(v, known[0][1], rtol=1e-11, atol=1e-12)
                         for _, v in known[1:]):
            raise ValueError(f"incompatible bin grids: {field}: {known}")
    units = {item.catalog["distance_unit"] for item in inputs if "distance_unit" in item.catalog}
    if len(units) > 1:
        raise ValueError(f"incompatible distance units: {sorted(units)}")
    if any(item.summary is None or not all(k in item.config for k in spec.grid) for item in inputs):
        notes.append("Some bin metadata are absent: coordinate/index mappings were checked, "
                     "but complete bin edges and units must agree in the original runs.")
    for field in ("project_delta", "redshift_weight_exponent", "weight_z_ref", "h", "omega_m"):
        known = [(item.label, item.catalog[field]) for item in inputs if field in item.catalog]
        if known and any(value != known[0][1] for _, value in known[1:]):
            notes.append(f"Catalog metadata differ for {field}: {known}")
    return notes


def align_tables(tables, spec, labels):
    """Align sparse native bins without interpolating or replacing missing data by zero."""
    ndim, ncoord = len(spec.indices), len(spec.coordinates)
    indices = np.concatenate([t[:, :ndim].astype(np.int64) for t in tables])
    keys, inverse = np.unique(indices, axis=0, return_inverse=True)
    coordinates = np.full((len(keys), ncoord), np.nan)
    num = np.zeros((len(keys), len(tables)))
    den = np.zeros_like(num)
    # Validate coordinate mappings even for bins which are disjoint in other axes.
    for c, axis in enumerate(spec.coordinate_axes):
        bins = indices[:, axis]
        values = np.concatenate([t[:, ndim+c] for t in tables])
        order = np.argsort(bins, kind="stable")
        ordered_bins, ordered_values = bins[order], values[order]
        if len(order):
            starts = np.r_[0, np.flatnonzero(np.diff(ordered_bins)) + 1]
            group = np.cumsum(np.r_[True, np.diff(ordered_bins) != 0]) - 1
            if not np.allclose(ordered_values, ordered_values[starts][group], rtol=1e-10, atol=1e-11):
                raise ValueError(f"incompatible coordinate/index mapping for {spec.coordinates[c]}")
            centers = ordered_values[starts]
            if np.any(np.diff(centers) <= 0):
                raise ValueError(f"non-increasing bin coordinates for {spec.coordinates[c]}")
            # Forest grids have known origins in these axes. This also detects
            # different bin widths when the two sparse catalogs occupy disjoint bins.
            name = spec.coordinates[c]
            if name not in ("lag1", "lag2", "r1_center", "r2_center"):
                origin = -1.0 if name == "mu" else 0.0
                widths = (centers-origin)/(ordered_bins[starts]+.5)
                if np.any(widths <= 0) or not np.allclose(widths, widths[0], rtol=1e-10, atol=1e-11):
                    raise ValueError(f"incompatible bin widths for {name}")
    offset = 0
    for j, table in enumerate(tables):
        slots = inverse[offset:offset+len(table)]
        coordinates[slots] = table[:, ndim:ndim+ncoord]
        num[slots, j], den[slots, j] = table[:, -2], table[:, -1]
        offset += len(table)
    return Aligned(keys, coordinates, num, den, labels,
                   spec.indices, spec.coordinates, spec.coordinate_axes)


def differences(data, reference, floor, relative_floor=0.0):
    values = data.values
    ref = values[:, reference]
    threshold = max(floor, relative_floor * float(np.max(np.abs(ref[np.isfinite(ref)]), initial=0)))
    shared = data.occupied & data.occupied[:, reference, None]
    delta = np.where(shared, values - ref[:, None], np.nan)
    relative = np.divide(delta, ref[:, None], out=np.full_like(delta, np.nan),
                         where=shared & (np.abs(ref[:, None]) > threshold))
    return delta, relative, threshold


def metrics(data, reference, floor, relative_floor):
    values = data.values
    delta, relative, threshold = differences(data, reference, floor, relative_floor)
    result = {}
    for j, label in enumerate(data.labels):
        if j == reference:
            continue
        common = data.occupied[:, reference] & data.occupied[:, j]
        valid = np.isfinite(relative[:, j])
        a, d = values[common, reference], delta[common, j]
        # Scale the norm before squaring to avoid overflow/underflow in ordinary norms.
        scale = float(np.max(np.abs(a), initial=0))
        norm = float(np.linalg.norm(a/scale)) if scale else 0.0
        rel_l2 = float(100*np.linalg.norm(d/scale)/norm) if norm else None
        result[label] = dict(
            common_occupied_bins=int(common.sum()),
            reference_only_bins=int((data.occupied[:, reference] & ~data.occupied[:, j]).sum()),
            candidate_only_bins=int((~data.occupied[:, reference] & data.occupied[:, j]).sum()),
            relative_defined_bins=int(valid.sum()),
            near_zero_reference_bins=int((common & ~valid).sum()),
            max_abs_difference=float(np.max(np.abs(d))) if len(d) else None,
            rms_difference=float(np.sqrt(np.mean(d*d))) if len(d) else None,
            median_abs_relative_percent=float(100*np.median(np.abs(relative[valid, j]))) if valid.any() else None,
            max_abs_relative_percent=float(100*np.max(np.abs(relative[valid, j]))) if valid.any() else None,
            relative_l2_percent=rel_l2,
        )
    return dict(relative_floor_effective=threshold, union_bins=len(data.indices),
                all_catalogs_occupied_bins=int(np.all(data.occupied, axis=1).sum()), comparisons=result)


def aggregate(data, axes, mask):
    """Pool raw sums on one shared native-bin support, never average xi values."""
    axes = tuple(axes)
    keys, first, inverse = np.unique(data.indices[mask][:, axes], axis=0,
                                    return_index=True, return_inverse=True)
    num, den = (np.zeros((len(keys), len(data.labels))) for _ in range(2))
    np.add.at(num, inverse, data.numerator[mask])
    np.add.at(den, inverse, data.denominator[mask])
    coord_columns = [data.coordinate_axes.index(a) for a in axes]
    coordinates = data.coordinates[mask][first][:, coord_columns]
    return Aligned(keys, coordinates, num, den, data.labels,
                   tuple(data.index_names[a] for a in axes),
                   tuple(data.coordinate_names[c] for c in coord_columns), tuple(range(len(axes))))


def slug(label):
    return re.sub(r"[^a-zA-Z0-9_.-]+", "_", label).strip("._")[:60] or "catalog"


def write_csvs(data, reference, floor, relative_floor, output, tag):
    delta, relative, _ = differences(data, reference, floor, relative_floor)
    values = data.values
    files = []
    for j, label in enumerate(data.labels):
        if j == reference:
            continue
        path = output / f"{tag}__{j+1}_{slug(label)}_vs_{reference+1}_{slug(data.labels[reference])}.csv"
        header = [*data.index_names, *data.coordinate_names, "reference", "candidate",
                  "difference", "relative_difference", "relative_percent", "reference_denominator",
                  "candidate_denominator", "common_occupied", "relative_defined"]
        with path.open("w") as stream:
            stream.write(",".join(header) + "\n")
            # Bound temporary export memory for large five-dimensional histograms.
            for start in range(0, len(values), 10000):
                s = slice(start, start+10000)
                rows = np.column_stack((data.indices[s], data.coordinates[s], values[s, reference],
                    values[s, j], delta[s, j], relative[s, j], 100*relative[s, j],
                    data.denominator[s, reference], data.denominator[s, j],
                    data.occupied[s, reference] & data.occupied[s, j], np.isfinite(relative[s, j])))
                np.savetxt(stream, rows, delimiter=",", fmt="%.17g")
        files.append(str(path))
    return files


def pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_figure(fig, output, tag, formats, dpi):
    paths = []
    for fmt in formats:
        path = output / f"{tag}.{fmt}"
        fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight", pad_inches=.15)
        paths.append(str(path))
    pyplot().close(fig)
    return paths


def curve_plot(data, reference, args, tag, title, xlabel, *, physical=False):
    plt = pyplot()
    values = data.values
    delta, relative, threshold = differences(data, reference, args.relative_floor, args.relative_scale_floor)
    if len(values) > args.max_plot_bins:
        sample = np.unique(np.linspace(0, len(values)-1, args.max_plot_bins).astype(int))
        title += f"\nPlot samples {len(sample):,} of {len(values):,} bins; CSV retains all bins"
    else:
        sample = np.arange(len(values))
    x = data.coordinates[:, 0] if physical else np.arange(len(values))
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                              gridspec_kw={"height_ratios": [2, 1, 1]})
    fig.subplots_adjust(left=.12, right=.97, bottom=.12, top=.89, hspace=.12)
    colors = plt.get_cmap("tab10")
    for j, label in enumerate(data.labels):
        style = dict(color=colors(j % 10), lw=1.15, marker="." if len(sample) < 200 else None,
                     ms=3, label=label + (" (reference)" if j == reference else ""))
        axes[0].plot(x[sample], values[sample, j], **style)
        if j != reference:
            axes[1].plot(x[sample], delta[sample, j], **style)
            axes[2].plot(x[sample], 100*relative[sample, j], **style)
    for ax in axes:
        ax.axhline(0, color=".35", lw=.65, zorder=0)
        ax.grid(alpha=.18)
        ax.tick_params(labelsize=9)
    axes[0].set_ylabel(r"$\xi$" if "2pcf" in tag else r"$\zeta$")
    axes[1].set_ylabel(r"$C-C_{\rm ref}$")
    axes[2].set_ylabel(r"$100(C-C_{\rm ref})/C_{\rm ref}$ [%]", fontsize=10)
    if args.relative_plot_scale == "symlog":
        axes[2].set_yscale("symlog", linthresh=args.relative_linthresh)
    axes[2].set_xlabel(xlabel)
    axes[0].legend(loc="best", fontsize=9, ncol=min(3, len(data.labels)))
    candidates = [j for j in range(len(data.labels)) if j != reference]
    if not np.any(np.isfinite(delta[:, candidates])):
        axes[1].text(.5, .5, "No shared occupied bins", transform=axes[1].transAxes, ha="center")
    if not np.any(np.isfinite(relative[:, candidates])):
        axes[2].text(.5, .5, "Relative differences undefined", transform=axes[2].transAxes, ha="center")
    fig.suptitle(title, fontsize=12, y=.97)
    scale_note = f" Percent axis: {args.relative_plot_scale}."
    fig.text(.12, .035, f"Reference: {data.labels[reference]}. Percent differences omit |reference| <= {threshold:.3g}.{scale_note}\n"
             "Missing/empty bins are undefined. Differences do not estimate statistical uncertainty.", fontsize=9)
    return save_figure(fig, args.output, tag, args.formats, args.dpi)


AXIS_LABELS = {
    "rp_center": r"$r_\parallel$ [Mpc/h]", "rt_center": r"$r_\perp$ [Mpc/h]",
    "r1": r"$r_1$ [Mpc/h]", "r2": r"$r_2$ [Mpc/h]",
    "lag1": r"$\Delta\chi_1$ [Mpc/h]", "lag2": r"$\Delta\chi_2$ [Mpc/h]",
    "r1_center": r"$r_1$ [Mpc/h]", "r2_center": r"$r_2$ [Mpc/h]",
    "radial_separation": r"$|\Delta\chi|$ [Mpc/h]",
}


def center_edges(centers, singleton_width=None):
    if len(centers) == 1:
        width = singleton_width or abs(centers[0]) or 1.
        return np.array([centers[0]-width/2, centers[0]+width/2])
    mid = (centers[1:] + centers[:-1])/2
    return np.r_[2*centers[0]-mid[0], mid, 2*centers[-1]-mid[-1]]


def map_plot(data, reference, args, tag, title, widths=(None, None)):
    """Plot sparse cells by their index lattice; absent rows remain masked."""
    if not len(data.indices):
        return []
    spans = [int(np.max(data.indices[:, a]))-int(np.min(data.indices[:, a]))+1 for a in range(2)]
    if spans[0]*spans[1] > 4_000_000:
        return []  # Check before constructing potentially very large index arrays.
    # Native grids are uniform in these coordinates. Fill index gaps with NaNs;
    # do not stretch an occupied sparse cell across missing neighbors.
    grids, row_indices = [], []
    for axis in range(2):
        index, first = np.unique(data.indices[:, axis], return_index=True)
        coords = data.coordinates[first, axis]
        grid = np.arange(index.min(), index.max()+1)
        if len(index) > 1:
            slope = np.diff(coords)/np.diff(index)
            if not np.allclose(slope, slope[0], rtol=1e-9, atol=1e-10):
                return []  # Nonuniform detached grids still have exact CSV/native-bin curves.
            centers = coords[0] + (grid-index[0])*slope[0]
        else:
            centers = coords
        grids.append(centers)
        row_indices.append(data.indices[:, axis]-index[0])
    plt = pyplot()
    values = data.values
    delta, relative, threshold = differences(data, reference, args.relative_floor, args.relative_scale_floor)
    fields = [values, delta, 100*relative]
    limits = [max(float(np.max(np.abs(f[np.isfinite(f)]), initial=0)), 1e-15) for f in fields]
    n = len(data.labels)
    fig, axes = plt.subplots(n, 3, figsize=(13.5, 3.1*n+1.3), squeeze=False)
    fig.subplots_adjust(left=.08, right=.98, bottom=.1, top=1-.7/(3.1*n+1.3), hspace=.45, wspace=.4)
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad="#e5e7eb")
    xedges, yedges = center_edges(grids[1], widths[1]), center_edges(grids[0], widths[0])
    for j, label in enumerate(data.labels):
        for k, (field, limit) in enumerate(zip(fields, limits)):
            ax = axes[j, k]
            if j == reference and k:
                ax.set_axis_off()
                ax.text(.5, .5, "Reference catalog" if k == 1 else
                        f"Gray cells: undefined\nPercent floor: {threshold:.3g}",
                        transform=ax.transAxes, ha="center", va="center", color=".35")
                continue
            array = np.full((len(grids[0]), len(grids[1])), np.nan)
            array[tuple(row_indices)] = field[:, j]
            if k == 2 and args.relative_plot_scale == "symlog":
                from matplotlib.colors import SymLogNorm
                normalization = dict(norm=SymLogNorm(args.relative_linthresh, vmin=-limit, vmax=limit))
            else:
                normalization = dict(vmin=-limit, vmax=limit)
            plot = ax.pcolormesh(xedges, yedges, np.ma.masked_invalid(array),
                                 cmap=cmap, **normalization, shading="flat", rasterized=True)
            ax.set(xlabel=AXIS_LABELS.get(data.coordinate_names[1], data.coordinate_names[1]),
                   ylabel=AXIS_LABELS.get(data.coordinate_names[0], data.coordinate_names[0]))
            ax.set_title(label if k == 0 else ("Difference" if k == 1 else
                         f"Relative difference [%], {args.relative_plot_scale}"), fontsize=10)
            fig.colorbar(plot, ax=ax, fraction=.045, pad=.03)
    fig.suptitle(title, fontsize=12)
    fig.text(.08, .025, "Differences use shared occupied bins. These are catalog differences, not statistical error bars.", fontsize=9)
    return save_figure(fig, args.output, tag, args.formats, args.dpi)


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", type=Path, help="two or more output roots, engine folders or histogram files")
    parser.add_argument("--labels", nargs="+", help="one distinct label per input")
    parser.add_argument("--engine", help="engine subdirectory when inputs are driver result roots")
    parser.add_argument("--reference", help="reference label (default: first input)")
    parser.add_argument("--statistics", choices=("2pcf", "3pcf", "both"), default="both",
                        help="compare matching saved products; both includes all shared families")
    parser.add_argument("--product", action="append", choices=tuple(PRODUCTS),
                        help="restrict to a family; may be repeated (also useful with --select)")
    parser.add_argument("--select", action="append", default=[], metavar="INDEX=VALUE",
                        help="native-bin slice, e.g. --product 3d_3pcf --select t1=1 --select bmu=2")
    parser.add_argument("--relative-floor", type=float, default=1e-12,
                        help="absolute reference threshold for percentage differences")
    parser.add_argument("--relative-scale-floor", type=float, default=0.,
                        help="also mask references <= this fraction of the largest |reference|")
    parser.add_argument("--relative-plot-scale", choices=("linear", "symlog"), default="symlog",
                        help="percent plot scale (symlog retains zero/sign and shows large outliers)")
    parser.add_argument("--relative-linthresh", type=float, default=1.,
                        help="symlog linear region in percentage points (default +/-1 percent)")
    parser.add_argument("--output", type=Path, default=Path("Output_lya_catalog_comparison"))
    parser.add_argument("--formats", nargs="+", choices=("png", "pdf", "svg"), default=["png"])
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--max-plot-bins", type=int, default=6000, help="native plot sample cap; CSV/metrics keep every bin")
    parser.add_argument("--no-plots", action="store_true", help="write CSV and JSON only; Matplotlib is then unnecessary")
    args = parser.parse_args(argv)
    if len(args.inputs) < 2:
        parser.error("provide at least two saved catalog results")
    if args.labels is None:
        args.labels = [p.parent.name if p.is_file() else p.name for p in args.inputs]
    if len(args.labels) != len(args.inputs) or len(set(args.labels)) != len(args.labels):
        parser.error("--labels must contain one distinct label per input")
    if args.reference is None:
        args.reference = args.labels[0]
    if args.reference not in args.labels:
        parser.error("--reference must be one of the catalog labels")
    for name in ("relative_floor", "relative_scale_floor"):
        if not np.isfinite(getattr(args, name)) or getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and nonnegative")
    if args.dpi <= 0 or args.max_plot_bins < 2:
        parser.error("require --dpi > 0 and --max-plot-bins >= 2")
    if not np.isfinite(args.relative_linthresh) or args.relative_linthresh <= 0:
        parser.error("--relative-linthresh must be finite and positive")
    args.selections = {}
    for selection in args.select:
        try:
            name, number = selection.split("=")
            value = int(number)
            if value < 0 or name in args.selections:
                raise ValueError
        except ValueError:
            parser.error("--select needs a nonnegative integer INDEX=VALUE; specify each index once")
        args.selections[name] = value
    args.output = args.output.expanduser().resolve()
    return args


def run(args):
    inputs = [resolve_input(p, label, args.engine) for p, label in zip(args.inputs, args.labels)]
    common = set.intersection(*(set(item.files) for item in inputs))
    requested = set(args.product) if args.product else common
    if args.product and requested-common:
        raise ValueError(f"products absent from one or more inputs: {sorted(requested-common)}")
    keys = [key for key, spec in PRODUCTS.items() if key in requested & common and
            (args.statistics == "both" or args.statistics == f"{spec.order}pcf")]
    if not keys:
        raise ValueError("no shared products for the requested statistic; estimator families are kept separate")
    reference = args.labels.index(args.reference)
    # Read and validate all products before writing output, so incompatible grids fail early.
    loaded, notes = {}, []
    available = set.union(*(set(item.files) for item in inputs))
    for key in sorted(available-common):
        if (args.statistics == "both" or args.statistics == f"{PRODUCTS[key].order}pcf") and (not args.product or key in args.product):
            missing = ", ".join(item.label for item in inputs if key not in item.files)
            notes.append(f"Omitted {key}: no saved histogram in {missing}")
    for key in keys:
        spec = PRODUCTS[key]
        notes.extend(check_metadata(inputs, spec))
        tables = [read_table(item.files[key], spec) for item in inputs]
        for name, value in args.selections.items():
            if name not in spec.indices:
                raise ValueError(f"{key} has no index {name!r}; use --product to select the intended family")
            axis = spec.indices.index(name)
            tables = [table[table[:, axis] == value] for table in tables]
        loaded[key] = align_tables(tables, spec, args.labels)
    args.output.mkdir(parents=True, exist_ok=True)
    report = dict(schema_version=1, reference=args.reference,
                  definition="difference = candidate - reference; relative_percent = 100 * difference / reference",
                  interpretation="Catalog differences include sampling, weighting and selection effects. "
                                 "They are not statistical uncertainties or significance estimates.",
                  projection="Sum numerator / sum denominator on native bins occupied in EVERY catalog; "
                             "each catalog retains its own pair/triplet weights. No interpolation.",
                  relative_floor=args.relative_floor, relative_scale_floor=args.relative_scale_floor,
                  plotting=dict(relative_scale=args.relative_plot_scale,
                                relative_linthresh_percent=args.relative_linthresh,
                                max_plot_bins=args.max_plot_bins, formats=args.formats),
                  selections=args.selections, notes=list(dict.fromkeys(notes)),
                  inputs=[dict(label=i.label, directory=str(i.directory),
                               summary=str(i.summary) if i.summary else None,
                               histograms={k: str(v) for k, v in i.files.items()},
                               catalog=i.catalog, binning={k: i.config[k] for k in
                               set().union(*(set(PRODUCTS[p].grid) for p in keys)) if k in i.config})
                          for i in inputs], products={})
    for key, data in loaded.items():
        spec = PRODUCTS[key]
        item = dict(native=metrics(data, reference, args.relative_floor, args.relative_scale_floor), files=[])
        item["files"].extend(write_csvs(data, reference, args.relative_floor, args.relative_scale_floor, args.output, key))
        title = key.replace("_", " ")
        if args.selections:
            title += " | " + ", ".join(f"{k}={v}" for k, v in args.selections.items())
        if not args.no_plots:
            physical = len(spec.indices) == 1
            xlabel = AXIS_LABELS[spec.coordinates[0]] if physical else "Native bin index (lexicographic order; coordinates in CSV)"
            item["files"].extend(curve_plot(data, reference, args, key+"_bins", title, xlabel, physical=physical))
        shared = np.all(data.occupied, axis=1)
        if key == "multipole_3pcf":
            shared &= data.indices[:, 0] == 0
            axes, projection_name = (1, 2), "ell=0 monopole"
        elif key == "3d_3pcf":
            axes, projection_name = (0, 1), "angular-pooled triplet sums"
        elif key == "3d_2pcf":
            axes, projection_name = (0, 1), "native pair bins"
        elif key == "1d_3pcf":
            axes, projection_name = (0, 1), "native signed-lag bins"
        else:
            axes, projection_name = (0,), "native radial bins"
        if len(axes) == 2 and np.any(shared):
            projected = aggregate(data, axes, shared)
            item["map_projection"] = metrics(projected, reference, args.relative_floor, args.relative_scale_floor)
            item["map_projection"]["description"] = projection_name
            item["files"].extend(write_csvs(projected, reference, args.relative_floor, args.relative_scale_floor,
                                            args.output, key+"_map"))
            if not args.no_plots:
                widths = []
                for axis, name in enumerate(projected.coordinate_names):
                    if name in ("lag1", "lag2"):
                        cfg = inputs[reference].config
                        width = cfg["r3_max"]/cfg["r3_bins"] if all(k in cfg for k in ("r3_max", "r3_bins")) else None
                    elif name in ("r1_center", "r2_center"):
                        cfg = inputs[reference].config
                        width = (cfg["r3_max"]-cfg.get("multipole_rmin", 0))/cfg["r3_bins"] if all(k in cfg for k in ("r3_max", "r3_bins")) else None
                    else:
                        width = projected.coordinates[0, axis]/(projected.indices[0, axis]+.5)
                    widths.append(width)
                item["files"].extend(map_plot(projected, reference, args, key+"_maps",
                    title+" | "+projection_name+"\nCommon occupied native bins in all catalogs", widths))
            curve = aggregate(data, (axes[0],), shared)
            item["radial_projection"] = metrics(curve, reference, args.relative_floor, args.relative_scale_floor)
            item["radial_projection"]["description"] = "Raw sums pooled over all other retained axes on common occupied native bins"
            item["files"].extend(write_csvs(curve, reference, args.relative_floor, args.relative_scale_floor,
                                            args.output, key+"_radial"))
            if not args.no_plots:
                item["files"].extend(curve_plot(curve, reference, args, key+"_radial",
                    title+" | weighted radial projection\nSum numerator / sum denominator; common occupied native bins",
                    AXIS_LABELS[curve.coordinate_names[0]], physical=True))
        report["products"][key] = item
        for label, result in item["native"]["comparisons"].items():
            percent = result["median_abs_relative_percent"]
            display = "undefined" if percent is None else f"{percent:.6g}%"
            print(f"{key}: {label} vs {args.reference}: {result['common_occupied_bins']} shared bins; "
                  f"median |relative difference| = {display}; "
                  f"{result['near_zero_reference_bins']} near-zero reference bins")
    report_path = args.output / "comparison_summary.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for note in report["notes"]:
        print("NOTE:", note, file=sys.stderr)
    print(f"Saved plots, CSV comparisons and summary to {args.output}")
    return report


def main(argv=None):
    try:
        run(parse_arguments(argv))
    except (ValueError, OSError, ImportError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
