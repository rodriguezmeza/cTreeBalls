"""Saved-window F2/F4 regressions: independent NumPy systems and native exports."""
from pathlib import Path
import os
import subprocess

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
EXE = Path(os.environ.get("CBALLS", ROOT / "cballs")).resolve()


def write_export(prefix, signal, window, edges=None, logarithmic=0):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    mmax, bins = len(signal)-1, signal.shape[1]
    if edges is None:
        edges = np.linspace(.02, 1.5, bins+1)
    prefix.with_name(prefix.name + "_edge_manifest.txt").write_text(
        f"CBALLS_SCALAR_EDGE 1\nbins {bins}\nsignal_mmax {mmax}\n"
        f"window_mmax {len(window)-1}\ndimensions 3\nperiodic 0\nbox 2 2 2\n"
        "geometry observer-tangent-fourier\nnormalization raw-ordered-distinct\n"
        f"phase cc+ss+i(sc-cs)\nlogarithmic {logarithmic}\nlower_cutoff .02\nedges "
        + " ".join(format(x, '.17g') for x in edges) + "\ncomplete\n")
    # All four components matter; splitting cross terms catches swapped signs.
    cross = np.full_like(signal.real, .03125)
    cross[0] = 0
    for suffix, array in (("edge_cos", .4*signal.real), ("edge_sin", .6*signal.real),
                          ("edge_sincos", signal.imag+cross), ("edge_cossin", cross),
                          ("window_Re", window.real), ("window_Im", window.imag)):
        for m, matrix in enumerate(array):
            np.savetxt(f"{prefix}_{suffix}_{m+1}.txt", matrix, fmt="%.17g")
    return prefix


def arrays(mmax=1, bins=2):
    s = np.zeros((mmax+1, bins, bins), dtype=complex)
    w = np.zeros((2*mmax+1, bins, bins), dtype=complex)
    s[0] = 1; s[1:] = .2
    w[0] = 1
    return s, w


def saved(tmp_path, signal, window=None, expect=None, extra=()):
    output = tmp_path / "corrected"
    cwd = tmp_path / "unrelated cwd"
    cwd.mkdir(exist_ok=True)
    cmd = [str(EXE), "searchMethod=kdtree-2balls-omp",
           "in=" + str(signal) + ("," + str(window) if window is not None else ""),
           f"rootDir={output}", "numberThreads=1", "verbose=0", "verbose_log=0",
           "options=edge-corrections-from-files", *extra]
    proc = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=30)
    text = proc.stdout + proc.stderr
    if expect is not None:
        assert proc.returncode != 0, text
        assert expect in text, text
        return output
    assert proc.returncode == 0, text
    assert (output / "histZetaM_edge_result.txt").read_text().endswith("complete\n")
    return output


def result(output, mmax):
    return np.array([np.loadtxt(output / f"histZetaM_EE_{m+1}.txt", ndmin=2)
                     + 1j*np.loadtxt(output / f"histZetaM_EE_Im_{m+1}.txt", ndmin=2)
                     for m in range(mmax+1)])


def oracle(signal, window):
    mmax = len(signal)-1
    modes = np.arange(-mmax, mmax+1)
    difference = modes[:, None]-modes[None, :]
    out = np.empty_like(signal)
    for i in range(signal.shape[1]):
        for j in range(signal.shape[2]):
            a = window[np.abs(difference), i, j]
            a = np.where(difference < 0, a.conj(), a)
            rhs = signal[np.abs(modes), i, j]
            rhs = np.where(modes < 0, rhs.conj(), rhs)
            out[:, i, j] = np.linalg.solve(a, rhs)[mmax:]
    return out


def test_f2_keeps_window_order_twice_signal_order(tmp_path):
    s, w = arrays(); w[2] = .4
    prefix = write_export(tmp_path / "data" / "hist", s, w)
    corrected = result(saved(tmp_path, prefix, prefix), 1)
    np.testing.assert_allclose(corrected[1], 1/7, rtol=2e-15)
    np.testing.assert_allclose(corrected, oracle(s, w), rtol=2e-15)


def test_f4_uses_second_prefix_and_ignores_cwd(tmp_path):
    s, w = arrays(); w[0] = 2
    first = write_export(tmp_path / "signal directory" / "data", s, w)
    second = write_export(tmp_path / "window directory" / "random", s, w)
    a = result(saved(tmp_path, first, second), 1)
    w[0] = 100
    write_export(second, s, w)
    b = result(saved(tmp_path, first, second), 1)
    np.testing.assert_allclose(a[1], .1, rtol=2e-15)
    np.testing.assert_allclose(b[1], .002, rtol=2e-15)
    saved(tmp_path, first, tmp_path / "missing", expect="cannot open v1 manifest")


@pytest.mark.parametrize("mmax", [0, 1, 2, 4])
def test_complex_window_matches_independent_solve(tmp_path, mmax):
    s, w = arrays(mmax)
    rng = np.random.default_rng(104+mmax)
    s[0] = .5
    s[1:] = rng.normal(size=s[1:].shape) + 1j*rng.normal(size=s[1:].shape)
    w[0] = 4
    w[1:] = .03*(rng.normal(size=w[1:].shape) + 1j*rng.normal(size=w[1:].shape))
    prefix = write_export(tmp_path / "fixture", s, w)
    actual = result(saved(tmp_path, prefix, prefix), mmax)
    np.testing.assert_allclose(actual, oracle(s, w), atol=2e-15, rtol=3e-14)
    if mmax:
        assert np.max(np.abs(actual[1:].imag)) > .01


def test_statuses_preserve_valid_zero_and_rejected_bins(tmp_path):
    s, w = arrays(bins=3)
    s[:, 0, 0] = 0  # Valid measured zero.
    w[:, 0, 1] = 0  # Empty.
    w[:, 0, 2] = 1  # Rank one.
    w[2, 1, 0] = np.nan
    s[1, 1, 1] = complex(0, np.inf)
    w[0, 1, 2] = -1  # Nonpositive support.
    w[:, 2, 0] = 1; w[1:, 2, 0] -= 1e-15  # Nearly singular.
    prefix = write_export(tmp_path / "status", s, w)
    output = saved(tmp_path, prefix, prefix)
    actual = result(output, 1)
    diagnostics = np.loadtxt(output / "histZetaM_window_diagnostics.txt")
    states = diagnostics[:, 2].reshape(3, 3)
    np.testing.assert_array_equal(states, [[1, 2, 3], [4, 4, 2], [3, 1, 1]])
    assert np.all(actual[:, 0, 0] == 0)
    assert np.all(np.isnan(actual[:, states != 1]))
    assert np.all(np.isfinite(actual[:, states == 1]))


@pytest.mark.parametrize("damage,expected", [
    ("manifest", "cannot open v1 manifest"), ("version", "unsupported manifest"),
    ("normalization", "unsupported manifest"), ("phase", "unsupported manifest"),
    ("mode", "cannot open required mode"), ("cross", "cannot open required mode"),
    ("short", "wrong value count"), ("surplus", "wrong value count"),
    ("text", "malformed matrix"), ("overflow", "malformed matrix"),
    ("grid", "incompatible radial edges"), ("monopole", "monopoles must be real"),
    ("trailing", "unsupported manifest"), ("dimensions", "unsupported manifest"),
    ("large", "resource"),
])
def test_invalid_files_fail_before_publication(tmp_path, damage, expected):
    s, w = arrays()
    first = write_export(tmp_path / "first", s, w)
    second = write_export(tmp_path / "second", s, w)
    manifest = Path(f"{second}_edge_manifest.txt")
    mode = Path(f"{second}_window_Re_3.txt")
    if damage == "manifest": manifest.unlink()
    elif damage == "version": manifest.write_text(manifest.read_text().replace("EDGE 1", "EDGE 2"))
    elif damage == "normalization": manifest.write_text(manifest.read_text().replace("raw-ordered-distinct", "normalized"))
    elif damage == "phase": manifest.write_text(manifest.read_text().replace("cc+ss+i(sc-cs)", "cc+ss"))
    elif damage == "mode": mode.unlink()
    elif damage == "cross": Path(f"{first}_edge_cossin_2.txt").unlink()
    elif damage == "short": mode.write_text("0 0 0\n")
    elif damage == "surplus": mode.write_text("0 0 0 0 0\n")
    elif damage == "text": mode.write_text("0 not-a-number 0 0\n")
    elif damage == "overflow": mode.write_text("0 1e9999 0 0\n")
    elif damage == "grid": write_export(second, s, w, edges=[.02, .7, 1.5])
    elif damage == "monopole": Path(f"{second}_window_Im_1.txt").write_text("1 0 0 0\n")
    elif damage == "trailing": manifest.write_text(manifest.read_text()+"unexpected\n")
    elif damage == "dimensions": manifest.write_text(manifest.read_text().replace("dimensions 3", "dimensions 4"))
    elif damage == "large": manifest.write_text(manifest.read_text().replace("bins 2", "bins 2147483646"))
    output = saved(tmp_path, first, second, expect=expected)
    assert not (output / "histZetaM_edge_result.txt").exists()
    assert not (output / "histZetaM_EE_1.txt").exists()


def test_two_inputs_required_and_writer_failure(tmp_path):
    s, w = arrays(); prefix = write_export(tmp_path / "data", s, w)
    saved(tmp_path, prefix, expect="exactly two prefixes")
    output = tmp_path / "corrected"
    (output / "histZetaM_EE_2.txt").mkdir()
    saved(tmp_path, prefix, prefix, expect="cannot open output")
    assert not (output / "histZetaM_edge_result.txt").exists()
    (output / "histZetaM_EE_2.txt").rmdir()
    saved(tmp_path, prefix, prefix)


@pytest.mark.parametrize("engine", ["kdtree-2balls-omp", "balltree-2balls-omp", "octree-2balls-omp"])
@pytest.mark.parametrize("logarithmic", [False, True])
def test_native_export_round_trip(tmp_path, engine, logarithmic):
    from test_two_ball_edge_corrections import catalog, write_catalog
    write_catalog(tmp_path / "input.txt", catalog())
    output = tmp_path / "native"
    command = [str(EXE), f"searchMethod={engine}", f"in={tmp_path / 'input.txt'}",
        "infmt=columns-ascii-all", f"rootDir={output}", "numberThreads=2",
        f"useLogHist={str(logarithmic).lower()}", "usePeriodic=false", "rangeN=1.5",
        "rminHist=.02", "sizeHistN=4", "mChebyshev=2", "sizeHistPhi=8", "nsmooth=2",
        "theta=0", "verbose=0", "verbose_log=0",
        "options=KKKCorrelation,only-3pcf,weights-norm,read-mask,no-smooth-pivot,edge-corrections,no-normalize-HistZeta"]
    proc = subprocess.run(command, text=True, capture_output=True, timeout=60)
    assert proc.returncode == 0, proc.stdout+proc.stderr
    prefix = output / "histZetaM"
    corrected = saved(tmp_path, prefix, prefix)
    np.testing.assert_allclose(result(corrected, 2), result(output, 2), rtol=2e-13, atol=2e-14)
    a = np.loadtxt(corrected / "histZetaM_window_diagnostics.txt")
    b = np.loadtxt(output / "histZetaM_window_diagnostics.txt")
    np.testing.assert_array_equal(a[:, :4], b[:, :4])
    np.testing.assert_allclose(a[:, 4], b[:, 4], rtol=2e-13, atol=2e-14)


def test_python_failure_and_same_object_recovery(tmp_path):
    from cyballs import cballs
    s, w = arrays(); prefix = write_export(tmp_path / "data", s, w)
    obj = cballs()
    try:
        params = dict(searchMethod="kdtree-2balls-omp", infile=f"{prefix},{tmp_path/'missing'}",
            rootDir=str(tmp_path / "python"), numberThreads=1, verbose=0, verbose_log=0,
            options="edge-corrections-from-files")
        obj.set(params)
        with pytest.raises(Exception, match="cannot open v1 manifest"): obj.Run()
        obj.set(infile=f"{prefix},{prefix}")
        obj.Run()
        assert obj.Run() == 0.0  # Completed preprocessing can be reused.
        assert (tmp_path / "python" / "histZetaM_edge_result.txt").exists()
        assert obj.run_settings is None
        with pytest.raises(Exception): obj.getHistNN()
        obj.clean_all()
        from test_provenance_window import parameters, register
        from test_two_ball_edge_corrections import catalog
        register(obj, catalog())
        obj.set(parameters(tmp_path / "normal"))
        obj.Run()
        assert np.all(np.isfinite(obj.getHistNN()))
    finally:
        obj.struct_cleanup()


def test_subnormal_signal_survives_text_round_trip(tmp_path):
    s, w = arrays(mmax=0); s[0] = 1e-310
    prefix = write_export(tmp_path / "tiny", s, w)
    actual = result(saved(tmp_path, prefix, prefix), 0)
    np.testing.assert_allclose(actual.real, s.real, rtol=0, atol=4*np.nextafter(0., 1.))


def test_larger_window_order_is_compatible(tmp_path):
    s, w = arrays(); w[2] = .4
    signal = write_export(tmp_path / "signal", s, w)
    ws, ww = arrays(mmax=2); ww[:3] = w
    window = write_export(tmp_path / "window", ws, ww)
    np.testing.assert_allclose(result(saved(tmp_path, signal, window), 1), oracle(s, w), rtol=2e-15)


@pytest.mark.parametrize("suffix", [",", ",extra", ",,extra"])
def test_empty_or_extra_prefix_rejected(tmp_path, suffix):
    s, w = arrays(); prefix = write_export(tmp_path / "data", s, w)
    saved(tmp_path, prefix, str(prefix)+suffix, expect="exactly two prefixes")
