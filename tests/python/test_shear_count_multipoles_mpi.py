"""Run with mpiexec -n 2 python this_file.py: MPI N0 and override parity."""
import sys
from pathlib import Path
from unittest.mock import patch
import traceback

sys.path[:0]=[str(Path(__file__).resolve().parents[2]),str(Path(__file__).resolve().parent)]
import numpy as np
from mpi4py import MPI
import test_shear_sphere_octree_omp as ref


def main():
    comm=MPI.COMM_WORLD
    for tree in ('octree','kdtree','balltree'):
        for mode in ('n0-auto','n0-cross','forced-full'):
            catalogs=[ref.fixture(44)]
            if mode=='n0-cross':catalogs.extend([ref.fixture(31),ref.fixture(28)])
            options='no-out-Hist,no-smooth-pivot,no-one-ball,only-3pcf,no-count-multipoles'
            if mode=='forced-full':options+=',edge-corrections'
            with patch.object(ref,'NMAX',7):
                actual=ref.run_native_catalogs(catalogs,2,options=options,
                    engine=f'{tree}-shear-sphere-2balls-mpi',collect_results=comm.rank==0)
                error=None
                if comm.rank==0:
                    try:
                        expected=ref.run_native_catalogs(catalogs,2,options=options,
                            engine=f'{tree}-shear-sphere-2balls-omp')
                        for key in actual:
                            np.testing.assert_allclose(actual[key],expected[key],
                                rtol=1e-8 if key=='multipoles' and mode=='forced-full' else 2e-12,
                                atol=2e-12,err_msg=f'{tree}/{mode}/{key}')
                        assert actual['window'].shape[0]==(29 if mode=='forced-full' else 1)
                    except Exception:error=traceback.format_exc()
                errors=comm.allgather(error)
                if any(errors):raise AssertionError('\n'.join(e for e in errors if e))
            if comm.rank==0:print('MPI_COUNT_MODES_OK',tree,mode,flush=True)


if __name__=='__main__':main()
