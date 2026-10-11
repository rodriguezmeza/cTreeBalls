"""Keep runtime help in step with the cached option registry."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]


def test_cached_options_have_help_entries():
    cache = (ROOT / "include/options_cache.h").read_text()
    cached = set(re.findall(r'X\([^,]+,[^,]+,\s*"([^"]+)"', cache))
    startup = (ROOT / "source/startrun.c").read_text()
    options = startup.split("local int print_options(")[-1]
    documented = set(re.findall(r'\{"([^"]+)",', options))
    assert cached
    # Keep retired cache enum slots stable for the native ABI, but do not
    # advertise options whose implementations are absent from this profile.
    retired = {'NNLandySzalay1', 'NNLandySzalay2', 'NNStandard', 'edge-effects',
               'pivot-loop', 'set-Nb-noSel', 'set-default-param', 'smooth',
               'smooth-min-cell'}
    assert not cached - documented - retired
    assert not documented & retired


def test_make_info_covers_resolved_build_profiles():
    source = (ROOT / "source/startrun.c").read_text()
    names = set(re.findall(r"MAKE_SETTING\((\w+)\)", source))
    expected = {"BUILD_MPI", "BUILD_LYA", "BUILD_SHEAR", "BUILD_PRECISION", "MPICC",
                "KDTREEMPION", "KDTREE2BALLSOMPON", "KDTREE2BALLSMPION",
                "LYAFORESTMPION", "OCTREE3PCF3DMPION",
                "OCTREEBALLS4MPION",
                "LYA1D_OMP_PIVOT_BLOCK_SIZE", "LYA1D_TREE3_LEAF_SIZE",
                "CB3D_OMP_PIVOT_BLOCK_SIZE", "GGG_OMP_PIVOT_CHUNK_SIZE",
                "GGG_OMP_FRONTIER_TARGET_TASKS", "GGG_OMP_FRONTIER_MIN_ACTIVE",
                "GGG_MPI_PIVOT_CLAIM_SIZE", "SLEEFON", "SLEEF_PKG",
                "SLEEF_AVAILABLE", "SLEEF_ENABLED", "SLEEF_CFLAGS",
                "SLEEF_LDFLAGS", "SLEEF_LIBS", "SLEEF_RPATH", "THETA",
                "CBALLS_SHEAR_PIVOT_TOL", "CBALLS_SHEAR_BIN_THETA"}
    assert expected <= names


def test_runtime_help_documents_current_progress_and_count_policy():
    source = (ROOT / "source/startrun.c").read_text()
    assert "NMultipolesON is a legacy build switch" in source
    assert "CBALLS_SHEAR_TASK_MEMORY_MB" in source
    assert "even if set to 0" in source
    assert "bin slop does not replace theta" in source
    # MPI-only profiles must also advertise the reuse option.
    options = source.split("local int print_options(")[-1]
    before_reuse = options.split('{"shear-pivot-reuse"')[0].rsplit("#if", 1)[-1]
    for family in ("OCTREE", "KDTREE", "BALLTREE"):
        assert family + "SHEARSPHERE2BALLSMPI" in before_reuse
