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
    assert not cached - documented


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
                "SLEEF_LDFLAGS", "SLEEF_LIBS", "SLEEF_RPATH"}
    assert expected <= names
