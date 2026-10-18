"""Complete-forest kdtree broad phase preserves the exact 2PCF estimator."""
import os
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
import pytest

import test_lya_forest_omp as oracle
import test_lya_forest_los_tree as los


def sample(case, forests=67, pixels=4):
    rng=np.random.default_rng(1904)
    directions=rng.normal(size=(forests,3))
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    # Small groups yield real near-neighbours on a broad, non-axis-aligned sky.
    directions[1:4]=directions[0]+np.arange(1,4)[:,None]*.003
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    chi=rng.uniform(80,260,(forests,pixels,1))
    xyz=directions[:,None,:]*chi
    if case=='bent': xyz+=rng.normal(0,9,xyz.shape)
    elif case=='centres':
        xyz=np.empty_like(xyz)
        xyz[:]=np.linspace(-700,700,pixels)[None,:,None]*directions[:,None,:]
        # Every bounding-box centre is identical; all observer distances > 0.
    elif case=='long':
        xyz[0]=np.column_stack([np.linspace(50,500,pixels),np.zeros((pixels,2))])
        xyz[1]=np.column_stack([np.linspace(499,950,pixels),np.zeros((pixels,2))])
        xyz[2]=np.column_stack([np.linspace(501,900,pixels),np.zeros((pixels,2))])
    elif case=='boundaries':
        xyz[0]=np.array([[100,0,0],[106,0,0],[130,0,0],[160,0,0]])
        xyz[1]=np.array([[100,1e-5,0],[130-1e-5,0,1e-5],[130+1e-5,0,0],[-100,0,0]])
    field=rng.normal(size=forests*pixels)
    weight=rng.uniform(.1,2,forests*pixels)
    if case=='zero': weight[::3]=0
    if case=='cancellation': field[:]=np.where(np.arange(len(field))%2,1,-1)
    points=np.column_stack([xyz.reshape(-1,3),field,weight,
        np.repeat(np.arange(forests)*37-2**40,pixels)])
    rng.shuffle(points)
    return points


def counters(log):
    line=next(x for x in log.splitlines() if 'Ly-alpha 2PCF forest index:' in x)
    return {k:int(v) for k,v in re.findall(r'(enabled|forests|nodes|node_visits|pruned_nodes|candidate_pairs|bitmap_bytes_per_thread)=(\d+)',line)}


@pytest.mark.parametrize('case',['ordinary','bent','long','centres','boundaries','zero','cancellation'])
def test_complete_bounds_exact_oracle_and_thread_order(tmp_path,case):
    points=sample(case)
    cat=los.catalog_file(tmp_path,points)
    expected=oracle.oracle_2pcf(points)
    for method in ('lya-2pcf-omp','lya-los-tree-2pcf-omp'):
        root=tmp_path/method;root.mkdir()
        los.run(cat,root/'reference',method,1,extra=('lya2Kernel=0',))
        log=los.run(cat,root/'one',method,1,extra=('lya2Kernel=1',))
        los.run(cat,root/'four',method,4,extra=('lya2Kernel=1',))
        los.compare_products(root/'reference',root/'one',(2,))
        oracle.assert_histogram_close(oracle.read_2pcf(root/'one'/los.PRODUCTS[2]),expected,'forest-bound independent pairs')
        assert (root/'one'/los.PRODUCTS[2]).read_bytes()==(root/'four'/los.PRODUCTS[2]).read_bytes()
        count=counters(log)
        assert count['enabled']==1 and count['forests']==67 and count['nodes']>1
        assert count['candidate_pairs']<=67*66//2
        if case=='ordinary':
            assert count['pruned_nodes']>0 and count['candidate_pairs']<67*66//2


def test_small_and_single_forest_fallback(tmp_path):
    points=sample('ordinary',forests=17)
    cat=los.catalog_file(tmp_path,points)
    log=los.run(cat,tmp_path/'small','lya-2pcf-omp',extra=('lya2Kernel=1',))
    assert counters(log)['enabled']==0
    points[:,5]=91
    np.savetxt(cat,points,fmt=('%.17g',)*5+('%.0f',))
    log=los.run(cat,tmp_path/'single','lya-2pcf-omp',extra=('lya2Kernel=1',))
    assert counters(log)['enabled']==0


@pytest.mark.parametrize('kernel3',[0,5])
def test_combined_keeps_pair_and_triple_results(tmp_path,kernel3):
    points=sample('ordinary',forests=35,pixels=4)
    cat=los.catalog_file(tmp_path,points)
    for prefix in ('lya','lya-los-tree'):
        root=tmp_path/prefix;root.mkdir()
        method=prefix+'-2pcf-3pcf-omp'
        los.run(cat,root/'old',method,1,extra=(f'lya3Kernel={kernel3}',))
        log=los.run(cat,root/'new',method,4,extra=('lya2Kernel=1',f'lya3Kernel={kernel3}'))
        los.compare_products(root/'old',root/'new')
        assert counters(log)['enabled']==1


def test_index_memory_preflight(tmp_path,monkeypatch):
    points=sample('ordinary',forests=257,pixels=2)
    cat=los.catalog_file(tmp_path,points)
    # End-to-end memory rejection followed by recovery without a cap.
    monkeypatch.setenv('CBALLS_MEMORY_BUDGET_MB','.01')
    with pytest.raises(AssertionError,match='budget|workspace|plan|memory'):
        los.run(cat,tmp_path/'limited','lya-2pcf-omp',4,extra=('lya2Kernel=1',))
    monkeypatch.delenv('CBALLS_MEMORY_BUDGET_MB')
    log=los.run(cat,tmp_path/'recover','lya-2pcf-omp',4,extra=('lya2Kernel=1',))
    assert counters(log)['enabled']==1


def test_candidate_summary_crosses_word_boundaries(tmp_path):
    n=4103
    rows=np.arange(n)
    xyz=np.column_stack([1e6+100*rows,50*rows,30*rows]).astype(float)
    selected=np.array([0,63,64,4095,4096,4102])
    xyz[selected]=np.column_stack([100+np.arange(len(selected)),np.zeros((len(selected),2))])
    points=np.column_stack([xyz,np.sin(rows),np.ones(n),rows*37-2**40])
    cat=los.catalog_file(tmp_path,points)
    log=los.run(cat,tmp_path/'index','lya-2pcf-omp',4,extra=('lya2Kernel=1',))
    expected=oracle.oracle_2pcf(points[selected])
    oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/'index'/los.PRODUCTS[2]),expected,'candidate bitmap boundaries')
    assert counters(log)['bitmap_bytes_per_thread']==8*(65+2)


@pytest.mark.parametrize('slop',[.1,.8])
def test_approximate_bins_keep_exact_pair_cuts_and_total_weights(tmp_path,slop):
    cat=los.catalog_file(tmp_path,sample('ordinary'))
    los.run(cat,tmp_path/'reference','lya-2pcf-omp')
    log=los.run(cat,tmp_path/'index','lya-2pcf-omp',4,
        extra=('lya2Kernel=1',f'lya2RpSlop={slop}',f'lya2RtSlop={slop}'))
    assert counters(log)['enabled']==1
    before=tmp_path/'reference'/los.PRODUCTS[2];after=tmp_path/'index'/los.PRODUCTS[2]
    count=lambda p:re.search(r'pairs: (\d+)',p.read_text())[1]
    assert count(before)==count(after)
    np.testing.assert_allclose(np.loadtxt(before,ndmin=2)[:,-2:].sum(axis=0),
        np.loadtxt(after,ndmin=2)[:,-2:].sum(axis=0),rtol=3e-12,atol=3e-12)


@pytest.mark.skipif(not os.environ.get('CTREEBALLS_MPI_COMMAND'),reason='set CTREEBALLS_MPI_COMMAND for explicit local MPI allocation')
def test_mpi_forest_rows_match_openmp(tmp_path):
    cat=los.catalog_file(tmp_path,sample('long'))
    for prefix in ('lya','lya-los-tree'):
        root=tmp_path/prefix;root.mkdir()
        los.run(cat,root/'omp',prefix+'-2pcf-omp',2,extra=('lya2Kernel=1',))
        params=dict(infile=cat,infileformat='lya-ascii',iCatalogs=1,rootDir=root/'mpi',
            searchMethod=prefix+'-2pcf-mpi',numberThreads=2,usePeriodic='false',useLogHist='false',
            rangeN=30,rminHist=.1,sizeHistN=4,lya2RpMax=30,lya2RtMax=30,
            lya2RpBins=oracle.RP_BINS,lya2RtBins=oracle.RT_BINS,lya2Kernel=1,
            verbose=2,verbose_log=2,options='no-smooth-pivot')
        cmd=shlex.split(os.environ['CTREEBALLS_MPI_COMMAND'])+['-n','2',str(los.BINARY)]
        cmd += [f'{k}={v}' for k,v in params.items()]
        result=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
        assert result.returncode==0,result.stdout+result.stderr
        los.compare_products(root/'omp',root/'mpi',(2,))
