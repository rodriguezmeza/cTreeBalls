#!/usr/bin/env python3
"""Build/install cyballs only when its inputs or installed artifacts change."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import sysconfig
import tempfile

from build_fingerprint import write_changed


def file_hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() if hasattr(
            hashlib, "file_digest") else hashlib.sha256(stream.read()).hexdigest()


def artifact_state(path):
    path = Path(path)
    stat = path.stat()
    return [str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def installation_probe(root, verify=False):
    # Do not mistake an in-place extension (or PYTHONPATH entry) for an install.
    code = """import importlib.util, json
spec = importlib.util.find_spec('cyballs')
result = {'path': spec.origin if spec else None}
"""
    if verify:
        code += """import cyballs
result['id'] = cyballs.build_info()['id']
"""
    code += "print(json.dumps(result))"
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory(prefix="cyballs-install-check-") as directory:
        result = subprocess.run([sys.executable, "-c", code], cwd=directory,
                                env=env, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError("Cannot verify installed cyballs:\n" + result.stderr)
    return json.loads(result.stdout)


def input_state(root, profile):
    info = json.loads((profile / "build-fingerprint.json").read_text())
    names = ("CC", "CXX", "CFLAGS", "CPPFLAGS", "LDFLAGS", "LDSHARED",
             "SDKROOT", "DEVELOPER_DIR", "MACOSX_DEPLOYMENT_TARGET",
             "ARCHFLAGS", "OMPI_CC", "MPICH_CC", "CBALLS_CPPFLAGS")
    environment = {name: os.environ.get(name, "") for name in names}
    for name in ("CC", "CXX", "CFLAGS", "CPPFLAGS", "LDFLAGS", "LDSHARED",
                 "ARCHFLAGS", "CBALLS_CPPFLAGS"):
        environment[name] = shlex.split(environment[name])
    return {
        "schema": 1, "id": info["id"], "python": sys.executable,
        "prefix": sys.prefix, "version": sys.version,
        "ext_suffix": sysconfig.get_config_var("EXT_SUFFIX"),
        "environment": environment,
        "library": file_hash(root / ("lib" + os.environ["CBALLS_LIB"] + ".a")),
    }


def is_current(saved, inputs, paths):
    try:
        return (saved["inputs"] == inputs and
                saved["artifacts"] == [artifact_state(path) for path in paths])
    except (KeyError, OSError, TypeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--profile", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root, profile = args.root.resolve(), args.profile.resolve()
    key = hashlib.sha256((sys.executable + "\n" + sys.prefix).encode()).hexdigest()[:16]
    stamp = profile / ("cyballs-installed-" + key + ".json")
    inputs = input_state(root, profile)
    local = root / ("cyballs" + inputs["ext_suffix"])
    try:
        saved = json.loads(stamp.read_text())
    except (OSError, ValueError):
        saved = {}
    if not isinstance(saved, dict):
        saved = {}
    installed = installation_probe(root)["path"]
    if (not args.force and installed and
            is_current(saved, inputs, [local, installed])):
        print("cyballs: unchanged; skipping build and installation.", flush=True)
        return

    # Distutils does not notice changes to compiler flags or older cached
    # profile headers. Force the wrapper only when its build identity changes.
    command = [sys.executable, "setup.py", "build_ext", "--inplace"]
    if args.force or saved.get("inputs") != inputs:
        command.append("--force")
    print("cyballs: building the matching Python extension...", flush=True)
    subprocess.check_call(command, cwd=root)
    print("cyballs: installing into " + sys.prefix, flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", ".",
                           "--no-build-isolation", "--no-deps", "--no-input"], cwd=root)
    installed = installation_probe(root, verify=True)
    if installed["id"] != inputs["id"]:
        raise RuntimeError("Installed cyballs has the wrong native build fingerprint")
    write_changed(stamp, json.dumps(
        {"inputs": inputs, "artifacts": [artifact_state(local),
                                        artifact_state(installed["path"])]},
        sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    main()
