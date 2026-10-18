"""Scientific contracts for the saved-catalog plotting utility (no cyballs required)."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_lya_catalogs as compare


def pair_table(xi=(.1, -.2, .0), den=(1., 2., 3.)):
    den = np.asarray(den)
    return np.column_stack((np.arange(3), np.zeros(3), [2., 6., 10.], np.full(3, 2.),
                            xi, np.asarray(xi)*den, den))


def save_run(path, table, *, key="3d_2pcf", config=None):
    path.mkdir(parents=True, exist_ok=True)
    np.savetxt(path/compare.PRODUCTS[key].filename, table, fmt="%.17g")
    if config is not None:
        (path/"summary.json").write_text(json.dumps(dict(config=config, catalog={"distance_unit": "Mpc/h"})))
    return path


def test_relative_sign_zero_floor_and_non_first_reference():
    a = pair_table()
    b = a.copy()
    b[:, -3:-1] *= 1.21  # Scaling every delta by 1.1 scales 2PCF by 1.21.
    d = compare.align_tables([b[::-1], a], compare.PRODUCTS["3d_2pcf"], ["mock", "data"])
    delta, rel, _ = compare.differences(d, 1, 1e-12)
    np.testing.assert_allclose(delta[:2, 0], [.021, -.042])
    np.testing.assert_allclose(rel[:2, 0], [.21, .21])
    assert np.isnan(rel[2, 0])
    report = compare.metrics(d, 1, 1e-12, 0.)["comparisons"]["mock"]
    assert report["relative_l2_percent"] == pytest.approx(21)
    assert report["near_zero_reference_bins"] == 1


def test_sparse_bins_are_undefined_and_reordered_bins_align():
    a = pair_table()
    b = a[[2, 1]].copy()
    b[0, -3:-1] = [3., 9.]
    d = compare.align_tables([a, b], compare.PRODUCTS["3d_2pcf"], ["A", "B"])
    assert np.isnan(d.values[0, 1])
    delta, rel, _ = compare.differences(d, 0, 1e-12)
    assert np.isnan(delta[0, 1])
    assert np.isnan(rel[2, 1])  # A zero correlation is valid data, but an invalid relative denominator.
    assert delta[2, 1] == 3
    stats = compare.metrics(d, 0, 1e-12, 0.)["comparisons"]["B"]
    assert stats["reference_only_bins"] == 1
    assert stats["common_occupied_bins"] == 2


def test_zero_denominator_is_not_a_zero_correlation():
    a = pair_table(xi=(.1, -.2, 0.), den=(1., 2., 0.))
    d = compare.align_tables([a, a], compare.PRODUCTS["3d_2pcf"], ["A", "B"])
    assert not d.occupied[2].any()
    assert np.isnan(d.values[2]).all()


def test_projection_uses_raw_weights_and_one_common_support():
    a = pair_table(xi=(1., 3., 500.), den=(1., 3., 100.))
    a[:, 0], a[:, 1] = 0, np.arange(3)
    a[:, 2], a[:, 3] = 2., [2., 6., 10.]
    b = a[:2].copy()
    b[:, -3:-1] *= 2
    d = compare.align_tables([a, b], compare.PRODUCTS["3d_2pcf"], ["A", "B"])
    projected = compare.aggregate(d, (0,), np.all(d.occupied, axis=1))
    # (1*1 + 3*3)/(1+3) = 2.5, not the arithmetic mean 2 or a pool including 500.
    np.testing.assert_allclose(projected.values, [[2.5, 5.]])
    np.testing.assert_allclose(projected.denominator, [[4., 4.]])


def test_three_catalog_projection_uses_intersection_of_all_three():
    a = pair_table(xi=(1., 2., 3.))
    d = compare.align_tables([a, a[:2], a[1:]], compare.PRODUCTS["3d_2pcf"], ["A", "B", "C"])
    np.testing.assert_array_equal(np.flatnonzero(np.all(d.occupied, axis=1)), [1])


@pytest.mark.parametrize("change,match", [
    (lambda t: np.vstack([t, t[:1]]), "duplicate"),
    (lambda t: np.column_stack([t[:, :6], -np.ones(len(t))]), "negative denominator"),
    (lambda t: t + np.array([.1, 0, 0, 0, 0, 0, 0]), "exact integers"),
    (lambda t: t + np.array([0, 0, 0, 0, 1, 0, 0]), "inconsistent"),
    (lambda t: np.column_stack([t[:, :6], [1., 2., np.nan]]), "finite columns"),
])
def test_reader_rejects_malformed_histograms(tmp_path, change, match):
    path = tmp_path/"hist.txt"
    np.savetxt(path, change(pair_table()))
    with pytest.raises(ValueError, match=match):
        compare.read_table(path, compare.PRODUCTS["3d_2pcf"])


def test_coordinate_mapping_rejects_same_indices_different_physical_bins():
    a, b = pair_table(), pair_table()
    b[:, 2] *= 2
    with pytest.raises(ValueError, match="coordinate/index"):
        compare.align_tables([a, b], compare.PRODUCTS["3d_2pcf"], ["A", "B"])


def test_sparse_disjoint_bin_widths_rejected_without_metadata():
    a, b = pair_table()[:1], pair_table()[1:2]
    b[:, 2] *= 2
    with pytest.raises(ValueError, match="bin widths"):
        compare.align_tables([a, b], compare.PRODUCTS["3d_2pcf"], ["A", "B"])


def test_changed_grid_metadata_fail_before_output_is_written(tmp_path):
    a = save_run(tmp_path/"a", pair_table(), config=dict(rp_max=12., rp_bins=3, rt_max=4., rt_bins=1))
    b = save_run(tmp_path/"b", pair_table(), config=dict(rp_max=24., rp_bins=3, rt_max=4., rt_bins=1))
    output = tmp_path/"bad_output"
    args = compare.parse_arguments([str(a), str(b), "--no-plots", "--output", str(output)])
    with pytest.raises(ValueError, match="incompatible bin grids"):
        compare.run(args)
    assert not output.exists()


def test_engine_resolution_never_silently_picks_one_of_many(tmp_path):
    save_run(tmp_path/"run"/"engineA", pair_table())
    save_run(tmp_path/"run"/"engineB", pair_table())
    with pytest.raises(ValueError, match="select --engine"):
        compare.resolve_input(tmp_path/"run", "data")
    resolved = compare.resolve_input(tmp_path/"run", "data", "engineB")
    assert resolved.directory.name == "engineB"


def test_no_shared_estimator_family_is_an_error(tmp_path):
    a = save_run(tmp_path/"A", pair_table())
    b = save_run(tmp_path/"B", np.array([[0, 2., .1, .1, 1.]]), key="1d_2pcf")
    args = compare.parse_arguments([str(a), str(b), "--no-plots"])
    with pytest.raises(ValueError, match="families are kept separate"):
        compare.run(args)


def test_5d_triplet_scaling_and_angular_pooling(tmp_path):
    # Two angular bins for one radial pair: weighted zeta=(1+9)/4=2.5.
    a = np.array([[0, 1, 0, 0, 0, 2, 6, np.pi/4, np.pi/4, -.5, 1, 1, 1],
                  [0, 1, 0, 0, 1, 2, 6, np.pi/4, np.pi/4,  .5, 3, 9, 3.]])
    b = a.copy()
    b[:, -3:-1] *= 1.331  # delta multiplied by 1.1 -> 3PCF multiplied by 1.1^3.
    paths = [save_run(tmp_path/name, table, key="3d_3pcf") for name, table in [("A", a), ("B", b)]]
    args = compare.parse_arguments([*map(str, paths), "--no-plots", "--output", str(tmp_path/"out")])
    result = compare.run(args)
    native = result["products"]["3d_3pcf"]["native"]["comparisons"]["B"]
    assert native["median_abs_relative_percent"] == pytest.approx(33.1)
    csv = next(args.output.glob("3d_3pcf_map__*.csv"))
    row = np.genfromtxt(csv, delimiter=",", names=True)
    assert row["reference"] == 2.5
    assert row["candidate"] == pytest.approx(2.5*1.331)
    assert row["relative_percent"] == pytest.approx(33.1)


def test_relative_scale_floor_and_all_zero_reference():
    a = pair_table(xi=(1., 1e-8, 0.))
    d = compare.align_tables([a, a], compare.PRODUCTS["3d_2pcf"], ["A", "B"])
    _, relative, threshold = compare.differences(d, 0, 0., 1e-6)
    assert threshold == 1e-6
    assert np.isnan(relative[1:]).all()
    a[:, -3:-1] = 0
    zero = compare.align_tables([a, a], compare.PRODUCTS["3d_2pcf"], ["A", "B"])
    assert compare.metrics(zero, 0, 0., 0.)["comparisons"]["B"]["relative_l2_percent"] is None


def test_header_only_empty_catalogs_remain_valid(tmp_path):
    spec = compare.PRODUCTS["3d_3pcf"]
    path = tmp_path/spec.filename
    path.write_text("# no accepted triplets\n")
    table = compare.read_table(path, spec)
    assert table.shape == (0, 13)
    d = compare.align_tables([table, table], spec, ["A", "B"])
    result = compare.metrics(d, 0, 1e-12, 0.)["comparisons"]["B"]
    assert result["common_occupied_bins"] == 0
    assert result["median_abs_relative_percent"] is None


def test_selection_and_three_catalog_cli_exports_full_csv(tmp_path):
    paths = [save_run(tmp_path/name, pair_table()) for name in ("A", "B", "C")]
    command = [sys.executable, str(Path(compare.__file__)), *map(str, paths), "--reference", "B",
               "--product", "3d_2pcf", "--select", "bp=1", "--output", str(tmp_path/"out"), "--no-plots"]
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads((tmp_path/"out"/"comparison_summary.json").read_text())
    assert summary["reference"] == "B"
    assert set(summary["products"]["3d_2pcf"]["native"]["comparisons"]) == {"A", "C"}
    rows = np.genfromtxt(next((tmp_path/"out").glob("3d_2pcf__*.csv")), delimiter=",", names=True)
    assert rows["bp"] == 1
    assert rows["reference"] == -.2


def test_plot_outputs_are_created_without_cyballs(tmp_path):
    paths = [save_run(tmp_path/name, pair_table()) for name in ("A", "B")]
    args = compare.parse_arguments([*map(str, paths), "--output", str(tmp_path/"out"), "--formats", "png", "pdf"])
    result = compare.run(args)
    files = result["products"]["3d_2pcf"]["files"]
    for name in ("3d_2pcf_bins", "3d_2pcf_radial", "3d_2pcf_maps"):
        for suffix in ("png", "pdf"):
            path = args.output/(name+"."+suffix)
            assert str(path) in files
            assert path.stat().st_size > 1000


@pytest.mark.parametrize("key,table", [
    ("1d_2pcf", [[0, 2., .1, 1., 10.], [1, 6., -.2, -2., 10.]]),
    ("1d_same_los_2pcf", [[0, 2., .1, 1., 10.], [1, 6., -.2, -2., 10.]]),
    ("1d_3pcf", [[0, 1, -6., -2., .1, 1., 10.], [1, 2, -2., 2., -.2, -2., 10.]]),
    ("multipole_3pcf", [[0, 1, 2, 2., 6., .1, 1., 10.], [1, 1, 2, 2., 6., -.2, -2., 10.]])
])
def test_other_native_estimator_families(tmp_path, key, table):
    a = np.asarray(table)
    b = a.copy()
    b[:, -3:-1] *= 1.25
    paths = [save_run(tmp_path/name, t, key=key) for name, t in [("A", a), ("B", b)]]
    args = compare.parse_arguments([*map(str, paths), "--output", str(tmp_path/"out"), "--no-plots"])
    result = compare.run(args)
    stats = result["products"][key]["native"]["comparisons"]["B"]
    assert stats["median_abs_relative_percent"] == pytest.approx(25)
    assert stats["relative_l2_percent"] == pytest.approx(25)


def test_no_shared_support_and_all_zero_reference_can_be_plotted(tmp_path):
    a = pair_table()[:1]
    b = pair_table()[1:2]
    paths = [save_run(tmp_path/name, t) for name, t in [("A", a), ("B", b)]]
    args = compare.parse_arguments([*map(str, paths), "--output", str(tmp_path/"out")])
    report = compare.run(args)
    assert "radial_projection" not in report["products"]["3d_2pcf"]
    a[:, -3:-1] = 0
    save_run(paths[0], a)
    save_run(paths[1], a)
    compare.run(args)
    assert (args.output/"3d_2pcf_bins.png").stat().st_size > 1000
