"""Legacy text plot limits use the units already stored in the input columns."""
from pathlib import Path
import importlib.util
import sys
from unittest.mock import patch

import matplotlib.figure
import numpy as np
import pytest


PLOTTERS = [p for p in ((Path(__file__).resolve().parents[1] / "python/compare_xi2pcf_curves.py"),
                        (Path(__file__).resolve().parents[1] / "python/compare_xi2pcf_curves_3files.py"))
            if p.exists()]


@pytest.fixture(params=PLOTTERS, ids=lambda path: path.stem)
def plotter(request):
    spec = importlib.util.spec_from_file_location("legacy_plot_under_test", request.param)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def command(plotter, path, *args):
    result = ["plot", "--file-a", str(path), "--file-b", str(path)]
    if "_3files" in plotter.__file__:
        result += ["--file-c", str(path)]
    return result + list(args)


def test_explicit_limits_set_both_axes_without_relabeling_input(plotter, tmp_path):
    path = tmp_path / "xi.txt"
    np.savetxt(path, np.column_stack(([3., 10., 80., 200.], [4., 3., 2., 1.])))
    captured = []

    def capture(fig, *args, **kwargs):
        captured.extend(fig.axes)

    argv = command(plotter, path, "--min-sep=2.5", "--max-sep=250", "--xscale=arcmin",
                   "--scale=semilogx", "--outdir", str(tmp_path))
    with patch.object(sys, "argv", argv), patch.object(matplotlib.figure.Figure, "savefig", capture):
        plotter.main()
    assert len(captured) == 2
    for ax in captured:
        np.testing.assert_allclose(ax.get_xlim(), [2.5, 250])
    np.testing.assert_array_equal(captured[0].lines[0].get_xdata(), [3., 10., 80., 200.])
    assert "arcmin" in captured[-1].get_xlabel()


def test_empty_cut_has_clear_error(plotter, tmp_path):
    path = tmp_path / "xi.txt"
    np.savetxt(path, [[1, 2], [3, 4]])
    argv = command(plotter, path, "--theta-min=8", "--theta-max=200")
    with patch.object(sys, "argv", argv), pytest.raises(ValueError, match="empty angular selection"):
        plotter.main()


@pytest.mark.parametrize("limits", [
    ["--min-sep=nan"], ["--max-sep=inf"], ["--min-sep=-1"],
    ["--theta-min=250", "--theta-max=2.5"],
])
def test_invalid_limits_fail_before_loading_data(plotter, tmp_path, limits):
    argv = command(plotter, tmp_path / "not-read.txt", *limits)
    with patch.object(sys, "argv", argv), pytest.raises(SystemExit) as exc:
        plotter.main()
    assert exc.value.code == 2
