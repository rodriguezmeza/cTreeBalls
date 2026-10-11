"""Offline job planning, MPI wrapper queries and launcher failure propagation."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'tests/slurm_scripts'
sys.path.insert(0,str(ROOT/'scripts'))
from build_fingerprint import mpi_wrapper_flags
spec=importlib.util.spec_from_file_location('perlmutter_run_job',ROOT/'tests/python/run_job.py')
job=importlib.util.module_from_spec(spec);spec.loader.exec_module(job)

class WrapperTests(unittest.TestCase):
    def wrapper(self,replies):
        def query(argv,**kwargs):
            value=replies.get(argv[-1])
            if value is None:raise subprocess.CalledProcessError(1,argv)
            return value
        with patch.object(subprocess,'check_output',side_effect=query):
            return mpi_wrapper_flags('cc')

    def test_cray(self):
        self.assertEqual(self.wrapper({'--cray-print-opts=cflags':'-I/cray/include',
                                      '--cray-print-opts=libs':'-L/cray/lib -lmpi -Wl,-rpath,/cray/lib'}),
                         (['-I/cray/include'],['-L/cray/lib','-lmpi','-Wl,-rpath,/cray/lib']))

    def test_openmpi(self):
        self.assertEqual(self.wrapper({'--showme:compile':'-I/mpi/include',
                                      '--showme:link':'-L/mpi/lib -lmpi'}),
                         (['-I/mpi/include'],['-L/mpi/lib','-lmpi']))

    def test_mpich(self):
        self.assertEqual(self.wrapper({'-show':'gcc -I/mpi/include -lmpich'}),
                         (['-I/mpi/include','-lmpich'],['-I/mpi/include','-lmpich']))

    def test_missing(self):
        with self.assertRaises(RuntimeError):self.wrapper({})

    def test_empty_compile_flags_allowed(self):
        self.assertEqual(self.wrapper({'--cray-print-opts=cflags':'','--cray-print-opts=libs':'-lmpi'}),([],['-lmpi']))

class JobTests(unittest.TestCase):
    def call(self,argv,env=None):
        out=io.StringIO()
        with patch.dict(os.environ,dict(env or {},CTREEBALLS_ROOT=str(ROOT)),clear=True),contextlib.redirect_stdout(out):
            code=job.main(argv)
        return code,out.getvalue()

    def test_inventory_and_syntax(self):
        registered=set(job.engines(ROOT))
        inventory=json.loads((HERE/'engine_jobs.json').read_text())['engines']
        self.assertEqual({e['engine'] for e in inventory},registered)
        self.assertEqual({p.stem for p in HERE.glob('*.slurm')} - {'validate_mpi', 'benchmark_dedicated'},registered)
        for path in [*HERE.glob('*.slurm'),*HERE.glob('*.sh')]:
            subprocess.run(['bash','-n',str(path)],check=True)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.params';p.write_text('infile = /shared/catalog.txt\nsearchMethod = wrong\nnumberThreads = 999\nrootDir = /old\n')
            for name in registered:
                _,output=self.call([name,'--dry-run','--params',str(p)])
                plan=json.loads(output);shape=plan['layout']
                self.assertEqual(shape['physical_cores'],256 if name.endswith('-mpi') else 128)
                self.assertIn('searchMethod = '+name,plan['effective_parameters'])
                self.assertNotIn('wrong',plan['effective_parameters'])
                self.assertNotIn('/old',plan['effective_parameters'])
                self.assertEqual(plan['command'][0],'srun')

    def test_shear_parent_and_shared_packet(self):
        _,output=self.call(['balltree-shear-sphere-2balls-mpi','--mode','shear','--dry-run','--',
            '--fits','/shared/data with spaces.fits','--more-options','shear-pivot-reuse'],
            {'SCRATCH':'/shared/scratch','CBALLS_SHEAR_PIVOT_TOL':'0.1'})
        plan=json.loads(output);cmd=plan['command']
        self.assertNotEqual(cmd[0],'srun');self.assertEqual(cmd[cmd.index('--mpiexec')+1],'srun')
        self.assertIn('--mpi-extra-arg=--kill-on-bad-exit=1',cmd)
        self.assertIn('/shared/data with spaces.fits',cmd)
        self.assertTrue(plan['environment']['TMPDIR'].startswith('/shared/scratch/'))
        self.assertEqual(plan['environment']['CBALLS_SHEAR_PIVOT_TOL'],'0.1')
        self.assertNotIn('CBALLS_SHEAR_TASKS',plan['environment'])

    def test_python_collective_driver(self):
        for engine,mode in [('octree-2balls-mpi','kappa'),('lya-los-tree-3pcf-mpi','lya'),('octree-3pcf-3d-mpi','lya')]:
            _,output=self.call([engine,'--mode',mode,'--dry-run','--','--catalog-npz' if mode=='kappa' else '--catalog','/shared/in.npz'])
            cmd=json.loads(output)['command'];self.assertEqual(cmd[0],'srun')
            self.assertIn('--ntasks=4',cmd)

    def test_allocations(self):
        with self.assertRaises(ValueError):job.layout('x-mpi',{},ranks_per_node=4,threads=64)
        with self.assertRaises(ValueError):job.layout('x-omp',{},nodes=2)
        env={'SLURM_JOB_ID':'123','SLURM_JOB_NUM_NODES':'2','SLURM_NTASKS':'4',
             'SLURM_NTASKS_PER_NODE':'2','SLURM_CPUS_PER_TASK':'128'}
        self.assertEqual(job.layout('x-mpi',env)['threads_per_rank'],64)
        with self.assertRaises(ValueError):job.layout('x-mpi',env,threads=128)

    def test_memory_catalogue_and_managed_arguments_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.params';p.write_text('infile = m0\n')
            for args in [ ['balltree-shear-sphere-2balls-mpi','--params',str(p),'--dry-run'],
                          ['balltree-shear-sphere-2balls-mpi','--mode','shear','--dry-run','--','--fits','x','--threads','999'] ]:
                with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):self.call(args)

    def test_failures_and_effective_native_file(self):
        with tempfile.TemporaryDirectory() as d:
            temp=Path(d);binary=temp/'cballs';srun=temp/'srun';p=temp/'input.params'
            binary.write_text('#!/bin/sh\necho "- octree-2balls-mpi (id=177)"\nexit 1\n');binary.chmod(0o755)
            srun.write_text('#!/bin/sh\necho simulated-rank-failure\nexit 7\n');srun.chmod(0o755)
            p.write_text('infile = /shared/catalog.txt\nsearchMethod = wrong\nnumberThreads = 1\nrootDir = /old\ntheta = 1\n')
            env={'SLURM_JOB_ID':'123','CBALLS_EXECUTABLE':str(binary),'CBALLS_SRUN':str(srun),'CBALLS_WORKDIR':str(temp)}
            code,_=self.call(['octree-2balls-mpi','--params',str(p),'--outdir',str(temp/'out'),'--','theta=0'],env)
            self.assertEqual(code,7)
            saved=(temp/'out/effective.parameters').read_text()
            self.assertIn('theta = 0',saved);self.assertNotIn('theta = 1',saved)
            self.assertIn('numberThreads = 64',saved)
            self.assertEqual(json.loads((temp/'out/job-status.json').read_text())['exit_code'],7)
            binary.write_text('#!/bin/sh\necho missing\nexit 1\n')
            with contextlib.redirect_stderr(io.StringIO()):
                code,_=self.call(['octree-2balls-mpi','--params',str(p),'--outdir',str(temp/'missing')],env)
            self.assertEqual(code,1)
            self.assertIn('unavailable',json.loads((temp/'missing/job-status.json').read_text())['error'])

    def test_scaling_is_dry_by_default(self):
        result=subprocess.run([sys.executable,str(ROOT/'tests/python/submit_scaling.py'),'octree-2balls-mpi',
            '--account','example','--nodes','1','2','--ranks-per-node','2','4','--repeats','2',
            '--','--params','/shared/input.params'],text=True,capture_output=True,check=True)
        self.assertEqual(len(result.stdout.splitlines()),8)
        self.assertIn('--cpus-per-task=64',result.stdout)

if __name__=='__main__':unittest.main()
