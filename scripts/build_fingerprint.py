#!/usr/bin/env python3
"""Deterministic source/profile identity, evaluated by the resolved Make graph.

Hash packaged source conservatively (including disabled implementations), so
changes can invalidate a build but cannot silently retain its source identity.
Artifacts and execution timestamps are deliberately not part of the build ID.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source_files(root):
    for directory in ("source", "main", "include", "general_libs", "getparam",
                      "addons", "python", "scripts", "tests", "capabilities"):
        for parent, dirs, files in os.walk(root / directory):
            dirs[:] = sorted(d for d in dirs if not d.startswith((".", "backup"))
                             and d not in {"build", "python_env", "__pycache__"})
            for name in sorted(files):
                path = Path(parent) / name
                if path.is_symlink() or name == "cyballs.c":
                    continue
                if (path.suffix in {".c", ".h", ".cpp", ".py", ".pyx", ".sh"}
                        or (directory == "capabilities" and path.suffix == ".json")
                        or name.endswith(".pxd.in") or name.startswith(("Makefile", "run_test"))
                        or "fixtures" in path.parts):
                    yield path
    for name in ("Makefile", "Makefile_settings", "Makefile_machine", "setup.py",
                 "pyproject.toml", "MANIFEST.in"):
        yield root / name


def query(command):
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    if result.returncode:
        raise RuntimeError(f"toolchain query failed: {command}: {result.stderr}")
    return result.stdout.strip()


def resolve(root):
    names = os.environ["CBALLS_FINGERPRINT_VARIABLES"].split()
    # Paths are retained with a relocatable source-root token in the identity.
    settings = {key: os.environ.get(key, "").replace(str(root), "${SOURCE_ROOT}")
                for key in sorted(set(names))}
    compiler = shlex.split(os.environ["CC"])
    sources = {str(p.relative_to(root)): digest(p.read_bytes())
               for p in sorted(set(source_files(root)))}
    source_digest = digest(json.dumps(sources, sort_keys=True, separators=(",", ":")).encode())
    tools = {"compiler_version": query(compiler + ["--version"]),
             "compiler_target": query(compiler + ["-dumpmachine"]),
             "python_version": platform.python_version(), "python_implementation": platform.python_implementation(),
             "system": platform.system(), "machine": platform.machine(), "byteorder": sys.byteorder}
    tools["build_packages"] = {name: importlib.metadata.version(name) for name in ("Cython", "numpy", "setuptools")}
    for name, command in (("gsl", [os.environ.get("GSL_CONFIG", "gsl-config"), "--version"]),
                          ("cfitsio", [os.environ.get("PKG_CONFIG", "pkg-config"), "--modversion", "cfitsio"])):
        try:
            tools[name+"_version"] = query(command)
        except (OSError, RuntimeError):
            tools[name+"_version"] = "unavailable; resolved link flags retained"
    if settings.get("BUILD_MPI") == "mpi":
        wrapper = shlex.split(os.environ.get("MPICC", os.environ["CC"]))
        tools["mpi_wrapper"] = query(wrapper + ["-show"])
        try:
            tools["mpi_version"] = query(wrapper + ["--showme:version"])
        except RuntimeError:
            tools["mpi_version"] = "wrapper does not expose --showme:version"
    identity = {"schema_version": 1, "resolved_settings": settings,
                "toolchain": tools, "source_sha256": source_digest,
                "source_count": len(sources),
                "source_scope": "conservative packaged C/header, Python, Make and regression inputs"}
    identity["id"] = digest(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode())
    return identity, sources


def write_changed(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(text)
        temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--header", type=Path)
    args = parser.parse_args()
    info, sources = resolve(args.root.resolve())
    if args.header:
        packed = json.dumps(info, sort_keys=True, separators=(",", ":"))
        header = ('#ifndef CBALLS_BUILD_FINGERPRINT_H\n#define CBALLS_BUILD_FINGERPRINT_H\n'
                  '#define CBALLS_BUILD_ID ' + json.dumps(info["id"]) + '\n'
                  '#define CBALLS_BUILD_JSON ' + json.dumps(packed) + '\n#endif\n')
        write_changed(args.header, header)
        write_changed(args.header.parent / "build-fingerprint.json",
                      json.dumps({**info, "source_files": sources}, sort_keys=True, indent=2) + "\n")
    else:
        print(json.dumps(info, sort_keys=True))


if __name__ == "__main__":
    main()
