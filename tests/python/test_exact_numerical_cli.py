#!/usr/bin/env python3
"""Native N1 analytic raw products and N2 fail-before-publication checks.

Uses only NumPy, so it also runs against the instrumented executable without
loading a Python extension into the sanitizer process.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile
import numpy as np

ENGINES=('octree-2balls-omp','kdtree-2balls-omp','balltree-2balls-omp','octree-3pcf-3d-omp')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cballs',type=Path,required=True)
    args=parser.parse_args();binary=args.cballs.resolve()
    p=np.vstack((np.eye(3),-np.ones(3)/np.sqrt(3)))
    with tempfile.TemporaryDirectory(prefix='exact-numerical-cli-') as tmp:
        root=Path(tmp)
        for engine in ENGINES:
            for label,w,target in [('dynamic',np.array([1e16,1,1,1]),18e16+6),
                                   ('overflow',np.full(4,1e155),None),
                                   ('underflow',np.full(4,1e-120),None),
                                   ('zero',np.zeros(4),0.)]:
                path=root/(engine+'-'+label);path.mkdir()
                catalog=path/'catalog.txt';np.savetxt(catalog,np.column_stack([p,np.ones(4),w]),fmt='%.17g')
                options='pos-and-convergence-weight,weights-norm,no-smooth-pivot,no-one-ball,no-two-balls,no-normalize-HistZeta'
                physical='3pcf-3d' in engine
                options+=',compute-2pcf-3d,compute-3pcf-3d' if physical else ',edge-corrections'
                command=[str(binary),f'search={engine}',f'in={catalog}',
                         'infmt=multi-columns-ascii','columns=1,2,3,4,5','iCatalogs=1',
                         f'rootDir={path}',f'options={options}','numberThreads=2',
                         'mChebyshev=2','sizeHistN=4','rangeN=1.99','rminHist=.03',
                         'useLogHist=true','usePeriodic=false','theta=0','verbose=0','verbose_log=0']
                result=subprocess.run(command,text=True,capture_output=True,timeout=90)
                if target is None:
                    assert result.returncode!=0,result.stdout+result.stderr
                    assert 'numerical range' in result.stdout+result.stderr
                    assert not list(path.glob('hist*.txt')),list(path.iterdir())
                else:
                    assert result.returncode==0,result.stdout+result.stderr
                    if physical:
                        rows=np.atleast_2d(np.loadtxt(path/'histZetaM_3d.txt'))
                        # Columns: ell, bin1, bin2, r1, r2, estimator, raw numerator, denominator.
                        rows=rows[rows[:,0]==0];values=rows[:,6:8].sum(axis=0)
                    else:
                        values=[np.loadtxt(path/'histZetaM_cos_1.txt').sum(),
                                np.loadtxt(path/'histZetaM_window_Re_1.txt').sum()]
                    np.testing.assert_allclose(values,[target,target],rtol=3e-12,atol=0)
    print('PASS: native exact distinct products, zero weights, overflow/underflow, and no invalid output')

if __name__=='__main__':main()
