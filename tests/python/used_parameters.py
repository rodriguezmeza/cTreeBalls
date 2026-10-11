"""Persistent, readable parameter records for native and dual node drivers."""
from collections.abc import Mapping
from dataclasses import asdict
import json
import os
from pathlib import Path
import tempfile

from benchmark_timing import runtime_environment


def _plain(value):
    if isinstance(value, Mapping):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def _value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _record(prefix, values, *, comment=True):
    """Extra information uses native '%' comments, not unknown C parameters."""
    lines = []
    for key, value in values.items():
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, Mapping):
            lines.extend(_record(name, value, comment=comment))
        else:
            encoded = json.dumps(_plain(value), ensure_ascii=False, allow_nan=False)
            lines.append(f"{'% ' if comment else ''}{name:<35} = {encoded}\n")
    return lines


def _write(directory, lines):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "parameters-usedvalues"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=directory,
                prefix=".parameters-usedvalues-", delete=False) as stream:
            temporary = Path(stream.name)
            stream.writelines(lines)
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return str(target)


def _driver_details(config, catalog):
    lines = ["\n% Python driver arguments; backend-effective settings are recorded separately.\n"]
    lines.extend(_record("driver", asdict(config)))
    if catalog is not None:
        lines.extend(_record("catalog", dict(nbody=int(catalog.nbody),
            metadata=getattr(catalog, "metadata", {}))))
    lines.extend(_record("environment", runtime_environment()))
    return lines


def save_native_parameters(model, config, catalog=None):
    """Retain the native parameter file after MainLoop, before scratch cleanup.

    Only global rank zero publishes. Values exposed by run_settings replace
    startup values at full Python float precision. Other native parameters
    remain in the original file, including addon-specific parameters.
    """
    metadata = model.getRunMetadata()
    if metadata.get("parallel", {}).get("rank", 0) != 0:
        return None
    effective = model.run_settings["effective"]
    source = Path(effective["rootDir"]) / "cyballs_param.txt-usedvalues"
    lines = ["% Native cTreeBalls parameters after MainLoop.\n",
             "% run_stage = \"MainLoop complete\"\n",
             "% Input labels such as m0 denote data supplied from Python, not a disk catalog.\n",
             "% rootDir is the native I/O directory used during this run; it may be temporary.\n"]
    for line in source.read_text(encoding="utf-8").splitlines(keepends=True):
        key, separator, value = line.partition("=")
        if separator and key.strip() in effective:
            value = effective[key.strip()]
            line = f"{key.strip():<35} = {_value(value)}\n"
        lines.append(line)
    lines.extend(_driver_details(config, catalog))
    lines.extend(_record("native", dict(catalog_sizes=effective["catalog_sizes"])))
    for section in ("effective_smoothing", "opening_tolerance", "bin_edges",
                    "parallel", "shear_hierarchical_reuse", "scalar_hierarchical_reuse"):
        if section in metadata:
            lines.extend(_record(section, metadata[section]))
    for section, family, bin_key in (
            ("shear_hierarchical_reuse", "SHEAR", "radial_bin_slop"),
            ("scalar_hierarchical_reuse", "SCALAR", "radial_bin_theta")):
        controls = metadata.get(section, {})
        if controls.get("enabled"):
            lines.extend(_record("used", {
                f"CBALLS_{family}_PIVOT_TOL": controls["phase_budget_radians"],
                f"CBALLS_{family}_BIN_THETA": controls[bin_key]}))
    build = metadata.get("build", {})
    compiled = build.get("resolved_settings", {})
    lines.extend(_record("build", dict(id=build.get("id"), **{
        name: value for name, value in compiled.items()
        if name == "THETA" or name.startswith(("CBALLS_SHEAR_", "CBALLS_SCALAR_"))})))
    return _write(Path(config.output_dir)/metadata["engine"], lines)


def save_native_startup(native_root, config, engine, catalog=None):
    """Publish the native startup record before counting, also on failed runs.

    The native writer creates the source only on global rank zero. Shear gives
    each rank its own scratch directory; kappa calls this only on rank zero.
    """
    source = Path(native_root)/"cyballs_param.txt-usedvalues"
    if not source.is_file():
        return None
    lines = ["% cTreeBalls startup parameters; counting has not completed.\n",
             "% run_stage = \"initialized\"\n",
             "% Effective backend controls are added after MainLoop completes.\n",
             source.read_text(encoding="utf-8")]
    lines.extend(_driver_details(config, catalog))
    return _write(Path(config.output_dir)/engine, lines)


def save_dual_node_parameters(config, catalog, result, correlations):
    """Record actual dual node properties, including defaults resolved by it."""
    lines = ["% dual node run settings; this is not a native cTreeBalls parameter file.\n",
             "% Native theta/nsmooth/reuse options in driver arguments do not control dual node.\n"]
    lines.extend(_record("", {name: result[name] for name in (
        "engine", "version", "dual_node_version", "geometry", "nbody",
        "threads_per_rank", "requested_threads_per_rank", "window_correction_applied",
        "edge_corrections", "three_pcf_normalization", "scalar_weights") if name in result},
        comment=False))
    for label, correlation in correlations.items():
        if correlation is None:
            continue
        fields = ("min_sep", "max_sep", "nbins", "bin_size", "bin_type", "sep_units",
                  "bin_slop", "angle_slop", "max_n", "brute", "metric", "coords")
        values = {name: getattr(correlation, name) for name in fields
                  if hasattr(correlation, name)}
        lines.extend(_record("dual_node."+label, values, comment=False))
    lines.extend(_driver_details(config, catalog))
    return _write(Path(config.output_dir)/result["engine"], lines)
