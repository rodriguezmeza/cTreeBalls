#!/usr/bin/env python3
"""Build-system regressions; small C fixtures use the real project Makefiles."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import build_support
import incremental_cyballs
from build_fingerprint import write_changed


class MetadataTests(unittest.TestCase):
    def test_write_only_when_changed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state"
            write_changed(path, "same")
            before = path.stat().st_mtime_ns
            write_changed(path, "same")
            self.assertEqual(before, path.stat().st_mtime_ns)
            write_changed(path, "different")
            self.assertEqual(path.read_text(), "different")

    def test_c_string_escaping(self):
        text = build_support.make_info({
            "CBALLS_MAKE_INFO_VARIABLES": "CC FLAGS MISSING",
            "CC": "cc", "FLAGS": 'a"b\\c\nd',
        })
        self.assertIn('#define CBALLS_MAKE_FLAGS "a\\"b\\\\c\\nd"', text)
        self.assertIn('#define CBALLS_MAKE_MISSING ""', text)

    def test_compile_and_link_settings(self):
        env = {"CC": "cc", "CCFLAG": "-DA=1", "CBALLS_OBJECTS": "a.o"}
        with patch.object(subprocess, "check_output", return_value="compiler 1"):
            before = build_support.command_signature("compile", env)
            self.assertEqual(before, build_support.command_signature(
                "compile", dict(env, CCFLAG="   -DA=1  ")))
            self.assertNotEqual(
                build_support.command_signature("compile", dict(env, CCFLAG='-DMSG="a b"')),
                build_support.command_signature("compile", dict(env, CCFLAG='-DMSG="a  b"')))
            self.assertEqual(before, build_support.command_signature(
                "compile", dict(env, CBALLS_OBJECTS="b.o", PYTHON="another-python")))
            self.assertNotEqual(before, build_support.command_signature(
                "compile", dict(env, CCFLAG="-DA=2")))
            self.assertNotEqual(before, build_support.command_signature(
                "compile", dict(env, OMPI_CC="another-compiler")))
            self.assertNotEqual(before, build_support.command_signature(
                "compile", dict(env, NATIVE_PAIR_VECTOR_LOG="compiler")))
            self.assertNotEqual(build_support.command_signature("link", env),
                                build_support.command_signature("link", dict(env, CBALLS_OBJECTS="b.o")))
        mpi = dict(env, BUILD_MPI="mpi", MPICC="mpicc")
        with patch.object(subprocess, "check_output", side_effect=["compiler 1", "cc -DMPI_ABI=1"]):
            before = build_support.command_signature("compile", mpi)
        with patch.object(subprocess, "check_output", side_effect=["compiler 1", "cc -DMPI_ABI=2"]):
            self.assertNotEqual(before, build_support.command_signature("compile", mpi))

    def test_installation_stamp_checks_artifacts_and_environment(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cyballs.so"
            path.write_bytes(b"test")
            inputs = {"id": "abc", "prefix": "env1"}
            saved = {"inputs": inputs, "artifacts": [incremental_cyballs.artifact_state(path)]}
            self.assertTrue(incremental_cyballs.is_current(saved, inputs, [path]))
            self.assertFalse(incremental_cyballs.is_current(saved, dict(inputs, prefix="env2"), [path]))
            path.write_bytes(b"changed")
            self.assertFalse(incremental_cyballs.is_current(saved, inputs, [path]))
            path.unlink()
            self.assertFalse(incremental_cyballs.is_current(saved, inputs, [path]))
            self.assertFalse(incremental_cyballs.is_current({}, inputs, [path]))

    def test_failed_install_does_not_publish_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = root / "build"
            profile.mkdir()
            with patch.object(sys, "argv", ["incremental_cyballs.py", "--root", directory,
                                           "--profile", str(profile)]), \
                 patch.object(incremental_cyballs, "input_state",
                              return_value={"id": "abc", "ext_suffix": ".so"}), \
                 patch.object(incremental_cyballs, "installation_probe",
                              return_value={"path": None}), \
                 patch.object(subprocess, "check_call",
                              side_effect=[0, subprocess.CalledProcessError(1, ["pip"])]):
                with self.assertRaises(subprocess.CalledProcessError):
                    incremental_cyballs.main()
            self.assertEqual(list(profile.glob("cyballs-installed-*.json")), [])


class MakeGraphTests(unittest.TestCase):
    def test_incremental_native_graph(self):
        make = os.environ.get("TEST_MAKE", "make")
        compiler = os.environ.get("TEST_CC", "cc")
        with tempfile.TemporaryDirectory(prefix="cballs-make-test-") as directory:
            root = Path(directory)
            for name in ("Makefile", "Makefile_machine", "Makefile_settings",
                         "setup.py", "pyproject.toml", "MANIFEST.in"):
                shutil.copy2(ROOT / name, root / name)
            (root / "scripts").mkdir()
            for name in ("build_support.py", "build_fingerprint.py",
                         "generate_make_info_header.sh", "generate_capabilities.py"):
                shutil.copy2(ROOT / "scripts" / name, root / "scripts" / name)
            (root / "source").mkdir()
            (root / "include").mkdir()
            (root / "capabilities").mkdir()
            (root / "capabilities" / "engines.json").write_text(json.dumps({
                "schema_version": 1, "engines": [], "regression_groups": {},
            }))
            subprocess.check_call([sys.executable, "scripts/generate_capabilities.py"],
                                  cwd=root)
            header = root / "include" / "leaf.h"
            header.write_text("#define LEAF 1\n")
            names = ("main", "cballsio", "cballs", "startrun", "testdata", "treeload",
                     "cballsutils", "search", "abi_check", "run_metadata",
                     "engine_registry", "runtime_context", "memory_catalog",
                     "common_histogram", "smooth_pivots", "mpi_runtime",
                     "clib", "mathfns", "inout", "mathutil", "numrec", "getparam")
            for name in names:
                source = "int " + name + "_fixture(void) { return 0; }\n"
                if name == "main":
                    source = "int main(void) { return 0; }\n"
                if name == "testdata":
                    source = '#include "leaf.h"\nint testdata_fixture(void) { return LEAF; }\n'
                if name == "run_metadata":
                    source = '#include "cballs_build_fingerprint.h"\nconst char *fixture_id = CBALLS_BUILD_ID;\n'
                if name == "startrun":
                    source = '#include "cballs_make_info.h"\n'
                (root / "source" / (name + ".c")).write_text(source)
            command = [make, "--no-print-directory", "-j4",
                       "ADDONSON=0", "CLASSLIBON=0", "OPENMPMACHINE=0", "OMPFLAG=",
                       "USEGSL=0", "SLEEFON=0", "CC=" + compiler,
                       "PYTHON=" + sys.executable, "cballs", "libcballs.a", "cyballs-static-lib"]

            def run(*options):
                # macOS ships Make 3.81, with whole-second file timestamps.
                time.sleep(1.05)
                result = subprocess.run(command + list(options), cwd=root,
                                        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                self.assertEqual(result.returncode, 0, result.stdout[-12000:])
                return result.stdout

            def objects(profile):
                return {p.name: p.stat().st_mtime_ns for p in profile.glob("*.o")}

            def changed(before, after):
                return {name for name in before if before[name] != after[name]}

            cold = run()
            self.assertEqual(cold.count("ar rv "), 1)
            profile = next((root / "build").glob("*-double"))
            before = objects(profile)
            self.assertEqual(len(before), len(names))
            libtime = (root / "libcballs.a").stat().st_mtime_ns
            exectime = (root / "cballs").stat().st_mtime_ns
            warm = run()
            self.assertNotIn(" -c ", warm)
            self.assertNotIn("ar rv ", warm)
            self.assertEqual(before, objects(profile))
            self.assertEqual(libtime, (root / "libcballs.a").stat().st_mtime_ns)
            self.assertEqual(exectime, (root / "cballs").stat().st_mtime_ns)

            source = root / "source" / "testdata.c"
            time.sleep(1.05)
            source.write_text(source.read_text() + "/* incremental fixture */\n")
            run()
            after = objects(profile)
            self.assertEqual(changed(before, after), {"testdata.o", "run_metadata.o"})
            before = after
            time.sleep(1.05)
            header.write_text("#define LEAF 2\n")
            run()
            after = objects(profile)
            self.assertEqual(changed(before, after), {"testdata.o", "run_metadata.o"})

            # Lost dependency files must regenerate even if source is unchanged.
            before = after
            (profile / "testdata.d").unlink()
            run()
            self.assertEqual(changed(before, objects(profile)), {"testdata.o"})

            before = objects(profile)
            run("LONGINTON=0")
            self.assertEqual(changed(before, objects(profile)), set(before))
            run()  # Restore the original flag set.
            before = objects(profile)

            # An edited Makefile with identical resolved commands is not a
            # reason to recompile every engine (provenance still changes).
            machine = root / "Makefile_machine"
            time.sleep(1.05)
            machine.write_text(machine.read_text() + "\n# metadata-only edit\n")
            run()
            self.assertEqual(changed(before, objects(profile)), {"run_metadata.o"})

            # Switching cached profiles must republish the matching artifacts.
            run("SINGLEPON=1")
            self.assertNotEqual((root / "cballs").read_bytes(), (profile / "cballs").read_bytes())
            before = objects(profile)
            run()
            self.assertEqual(before, objects(profile))
            self.assertEqual((root / "cballs").read_bytes(), (profile / "cballs").read_bytes())
            self.assertEqual((root / "libcballs.a").read_bytes(), (profile / "libcballs.a").read_bytes())

            # Shrinking the member list must remove stale archive members.
            run("EXTERNAL=clib.o")
            members = subprocess.check_output(["ar", "t", str(root / "libcballs.a")], text=True)
            self.assertNotIn("numrec.o", members)
            self.assertIn("clib.o", members)

            before = (profile / "testdata.o").read_bytes()
            source.write_text("this is not valid C;\n")
            time.sleep(1.05)
            result = subprocess.run(command, cwd=root, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(before, (profile / "testdata.o").read_bytes())


if __name__ == "__main__":
    unittest.main()
