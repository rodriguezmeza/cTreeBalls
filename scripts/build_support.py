#!/usr/bin/env python3
"""Content-stable build metadata and native command signatures."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess

from build_fingerprint import write_changed


COMPILE_VARIABLES = (
    "CC", "OPTFLAG", "OMPFLAG", "CCFLAG", "INCLUDES",
    "PROJECT_WARNING_FLAGS", "VENDORED_SOURCE_PATTERNS",
    "SDKROOT", "DEVELOPER_DIR", "MACOSX_DEPLOYMENT_TARGET",
    "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH",
    "COMPILER_PATH", "GCC_EXEC_PREFIX", "OMPI_CC", "MPICH_CC",
    "NATIVE_PAIR_VECTOR_LOG",
)
LINK_VARIABLES = (
    "CC", "OPTFLAG", "OMPFLAG", "LDFLAG", "MLIBS", "FITSIOLIBS",
    "AR", "CBALLS_OBJECTS", "LIBRARY_PATH", "SDKROOT", "DEVELOPER_DIR",
    "MACOSX_DEPLOYMENT_TARGET", "OMPI_CC", "MPICH_CC",
)
COMMAND_VARIABLES = {
    "CC", "OPTFLAG", "OMPFLAG", "CCFLAG", "INCLUDES",
    "PROJECT_WARNING_FLAGS", "VENDORED_SOURCE_PATTERNS",
    "LDFLAG", "MLIBS", "FITSIOLIBS", "AR", "CBALLS_OBJECTS",
}


def command_signature(kind, environ=None):
    environ = os.environ if environ is None else environ
    names = COMPILE_VARIABLES if kind == "compile" else LINK_VARIABLES
    values = {name: environ.get(name, "") for name in names}
    # GNU Make releases differ in whitespace when expanding empty additions.
    # Compare shell words, preserving whitespace inside quoted arguments.
    for name in COMMAND_VARIABLES.intersection(values):
        values[name] = shlex.split(values[name])
    compiler = shlex.split(environ["CC"])
    values["compiler_path"] = shutil.which(compiler[0], path=environ.get("PATH"))
    # Also notices upgrades of a compiler or changes behind an MPI wrapper.
    values["compiler_version"] = subprocess.check_output(
        compiler + ["--version"], text=True, env=environ)
    if environ.get("BUILD_MPI") == "mpi":
        wrapper = shlex.split(environ.get("MPICC") or environ["CC"])
        values["mpi_wrapper_flags"] = shlex.split(subprocess.check_output(
            wrapper + ["-show"], text=True, env=environ))
    return values


def make_info(environ=None):
    environ = os.environ if environ is None else environ
    lines = ["#ifndef CBALLS_MAKE_INFO_H", "#define CBALLS_MAKE_INFO_H", ""]
    for name in environ["CBALLS_MAKE_INFO_VARIABLES"].split():
        lines.append("#define CBALLS_MAKE_" + name + " " +
                     json.dumps(environ.get(name, "")))
    return "\n".join(lines + ["", "#endif", ""])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("compile", "link", "make-info"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    text = (make_info() if args.kind == "make-info" else
            json.dumps(command_signature(args.kind), sort_keys=True, indent=2) + "\n")
    write_changed(args.output, text)


if __name__ == "__main__":
    main()
