"""Public benchmark timing contracts without an MPI runtime."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from kappa_corr_all_engines import aggregate_rank_timings, _timing_metadata, parse_arguments
from shear_corr_all_engines import RunConfig
from kappa_corr_all_engines import write_timing_report


def test_rank_timings_sum_cpu_and_preserve_critical_path():
    rows = [_timing_metadata(5, 2, 1, 3, "test"),
            _timing_metadata(1, 4, 4, 7, "test")]
    for rank, row in enumerate(rows):
        row.update(rank=rank, native_reported_cpu_time=rank + 0.5)
    result = aggregate_rank_timings(rows)
    assert result["setup_wall_time"] == 5
    assert result["compute_wall_time"] == 4
    assert result["total_wall_time"] == 6
    assert result["compute_cpu_time"] == 10
    assert result["total_cpu_time"] == 16
    assert result["native_reported_cpu_time"] == 2
    assert result["ranks"] == 2
    assert result["rank_timings"] == rows


def test_no_participating_ranks_is_an_error():
    with pytest.raises(ValueError):
        aggregate_rank_timings([])


def test_documented_scalar_exact_switches_parse():
    args = parse_arguments([
        "--catalog-npz", "catalog.npz", "--engine", "all",
        "--more-options", "no-one-ball,no-two-balls", "--no-smooth-pivot",
    ])
    assert args.more_options == ["no-one-ball,no-two-balls"]
    assert args.no_smooth_pivot


def test_timing_report_describes_rank_aggregation(tmp_path):
    output = tmp_path / "timing_report.txt"
    write_timing_report(output, {})
    text = output.read_text()
    assert "rank 0 only" not in text
    assert "CPU is summed across participating ranks" in text


def test_shear_leaf_and_opening_parameters():
    assert RunConfig(nsmooth=8, tree_theta=0).normalized().nsmooth == 8
    for value in (float("nan"), float("inf"), -1):
        with pytest.raises(ValueError, match="tree_theta"):
            RunConfig(tree_theta=value).normalized()
    with pytest.raises(ValueError, match="nsmooth"):
        RunConfig(nsmooth=0).normalized()
