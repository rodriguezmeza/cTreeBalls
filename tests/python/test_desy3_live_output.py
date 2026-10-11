"""The batch driver must forward output before its child exits and retain logs."""
import os
from pathlib import Path
import select
import signal
import subprocess
import sys
import tempfile
import unittest

DRIVERS = [Path(__file__).resolve().parents[1] / 'python/run_desy3_shear_catalogs.py']
if os.environ.get('CBALLS_PRIVATE_SHEAR_DRIVER'):
    DRIVERS.append(Path(os.environ['CBALLS_PRIVATE_SHEAR_DRIVER']))

WORKER = '''import importlib.util, pathlib, sys
path = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(path.parent))
spec = importlib.util.spec_from_file_location('batch_driver', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
raise SystemExit(module.run_logged([sys.executable, '-u', '-c', sys.argv[3], sys.argv[4]], pathlib.Path(sys.argv[2])))
'''
CHILD = '''from pathlib import Path
import os, sys, time
os.write(1, b'first\\n')
while not Path(sys.argv[1]).exists(): time.sleep(.01)
os.write(2, b'second\\n')
raise SystemExit(7)
'''


class LiveOutputTests(unittest.TestCase):
    def test_live_forwarding_log_and_failure_status(self):
        for driver in DRIVERS:
            with self.subTest(driver=str(driver)), tempfile.TemporaryDirectory() as tmp:
                log, gate = Path(tmp)/'run.log', Path(tmp)/'release'
                with subprocess.Popen([sys.executable, '-u', '-c', WORKER,
                        str(driver), str(log), CHILD, str(gate)], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE) as process:
                    try:
                        ready, _, _ = select.select([process.stdout], [], [], 10)
                        self.assertTrue(ready, 'first line buffered until child completion')
                        self.assertEqual(process.stdout.readline(), b'first\n')
                        self.assertIsNone(process.poll())
                        self.assertEqual(log.read_bytes(), b'first\n')
                    finally:
                        gate.touch()
                        output, error = process.communicate(timeout=10)
                    self.assertEqual(process.returncode, 7, error)
                    self.assertEqual(output, b'second\n')
                    self.assertEqual(log.read_bytes(), b'first\nsecond\n')

    def test_interrupt_reaps_child(self):
        child = 'import os,time; print(os.getpid(),flush=True); time.sleep(60)'
        for driver in DRIVERS:
            with self.subTest(driver=str(driver)), tempfile.TemporaryDirectory() as tmp:
                log = Path(tmp)/'run.log'
                with subprocess.Popen([sys.executable, '-u', '-c', WORKER,
                        str(driver), str(log), child, 'unused'], stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE) as process:
                    try:
                        ready, _, _ = select.select([process.stdout], [], [], 10)
                        self.assertTrue(ready)
                        child_pid = int(process.stdout.readline())
                        os.kill(process.pid, signal.SIGINT)
                        process.communicate(timeout=10)
                        self.assertNotEqual(process.returncode, 0)
                        with self.assertRaises(ProcessLookupError):
                            os.kill(child_pid, 0)
                    finally:
                        if process.poll() is None:
                            process.terminate()
                            process.communicate(timeout=10)


if __name__ == '__main__':
    unittest.main()
