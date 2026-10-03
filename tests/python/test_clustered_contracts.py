"""Independent core/compatibility estimators and clustered approximations."""
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'scripts'),str(ROOT/'tests/python')]
from cyballs import cballs
from test_two_ball_edge_corrections import polar
CHECK=unittest.TestCase()


def fixture(seed,groups=8):
    rng=np.random.default_rng(seed)
    centers=rng.normal(size=(groups,3));centers/=np.linalg.norm(centers,axis=1)[:,None]
    sizes=np.arange(groups)%5+2
    pos=np.repeat(centers,sizes,axis=0)+rng.normal(scale=1e-6,size=(sum(sizes),3))
    pos/=np.linalg.norm(pos,axis=1)[:,None]
    field=.3+pos[:,0]-.2*pos[:,1] if seed%2 else rng.normal(size=len(pos))
    weights=rng.uniform(.2,2,len(pos));weights[::9]=0
    return pos,field,weights


def oracle(data,edges,orders,core):
    p,k,w=data;n=len(p);bins=len(edges)-1
    num=np.zeros(bins);den=np.zeros(bins);z=np.zeros((orders,bins,bins),complex)
    for i in range(n):
        distance,phase=polar(p[i],p)
        valid=(distance>edges[0])&(distance<edges[-1])
        indices=np.flatnonzero(valid);b=np.searchsorted(edges,distance[valid],side='right')-1
        np.add.at(num,b,w[i]*k[i]*w[indices]*k[indices])
        np.add.at(den,b,1. if core else w[i]*w[indices])
        counts=np.bincount(b,minlength=bins)
        for order in range(orders):
            moment=np.zeros(bins,complex);selfterm=np.zeros(bins)
            angular=np.isfinite(phase[indices]);jj=indices[angular];bb=b[angular]
            np.add.at(moment,bb,w[jj]*k[jj]*phase[jj]**order)
            if core:
                moment/=np.maximum(counts,1)
                z[order]+=k[i]/n*np.outer(moment,moment.conj())
            else:
                np.add.at(selfterm,bb,(w[jj]*k[jj])**2)
                z[order]+=w[i]*k[i]*(np.outer(moment,moment.conj())-np.diag(selfterm))
    return np.divide(num,den,out=np.zeros_like(num),where=den!=0),z


def run(data,root,core,logarithmic,orders,mode):
    edges=np.geomspace(.02,1.8,5) if logarithmic else np.linspace(.02,1.8,5)
    options=['KKKCorrelation','no-out-Hist','weights-norm','no-normalize-HistZeta']
    if not core:options+=['legacy-one-ball']
    if mode=='exact':options+=['no-one-ball','no-smooth-pivot']
    else:options+=['smooth-pivot' if mode=='smooth' else 'no-smooth-pivot']
    model=cballs()
    try:
        model.set(searchMethod='octree-sincos-omp' if core else 'octree-2balls-omp',
            rootDir=str(root),numberThreads=2,sizeHistN=4,mChebyshev=orders-1,rangeN=1.8,rminHist=.02,
            useLogHist=logarithmic,usePeriodic=False,theta=0 if mode=='exact' else .01,
            nsmooth=4,rsmooth='1',verbose=0,verbose_log=0,options=','.join(options))
        model.set_catalog(data[0],kappa=data[1],weights=data[2]);model.Run()
        return model.getResults()
    finally:model.struct_cleanup()


def test_exact_and_clustered_approximations(tmp_path):
    for seed in (41,42,43):
        data=fixture(seed)
        for core in (False,True):
            for logarithmic in (False,True):
                for orders in (3,4,9):
                    reference=run(data,tmp_path,core,logarithmic,orders,'exact')
                    edges=np.asarray(reference['metadata']['bin_edges']['radial'])
                    pair,triple=oracle(data,edges,orders,core)
                    np.testing.assert_allclose(reference['arrays']['xi'],pair,rtol=2e-12,atol=1e-12)
                    np.testing.assert_allclose(reference['arrays']['zeta_raw'],triple,rtol=3e-10,atol=2e-10)
                    # Signed, rapidly varying fields must also preserve group sums.
                    for mode in ('approx','smooth'):
                        result=run(data,tmp_path,core,logarithmic,orders,mode)
                        for key,target in reference['arrays'].items():
                            error=np.linalg.norm(result['arrays'][key]-target);norm=np.linalg.norm(target)
                            CHECK.assertLessEqual(error, .02*norm+1e-10,(seed,core,logarithmic,orders,mode,key,error,norm))


def test_qualification_is_observable_and_catalog_specific(tmp_path):
    from cyballs import qualification_report
    data=fixture(41)
    exact=run(data,tmp_path,False,True,4,'exact')
    candidate=run(data,tmp_path,False,True,4,'approx')
    report=qualification_report(candidate,exact)
    CHECK.assertEqual(report['status'],'QUALIFIED')
    damaged=dict(candidate,arrays={k:v.copy() for k,v in candidate['arrays'].items()})
    damaged['arrays']['xi'][0]=np.nan
    report=qualification_report(damaged,exact)
    CHECK.assertEqual(report['status'],'REJECTED')
    CHECK.assertEqual(report['observables']['zeta_raw']['status'],'QUALIFIED')
    with CHECK.assertRaises(ValueError):qualification_report(candidate,candidate)
    other=run(fixture(42),tmp_path,False,True,4,'exact')
    with CHECK.assertRaises(ValueError):qualification_report(candidate,other)
    incompatible=run(data,tmp_path,True,True,4,'exact')
    with CHECK.assertRaises(ValueError):qualification_report(candidate,incompatible)


def test_result_packets_and_integer_counts(tmp_path):
    from cyballs import qualification_report,save_result_packet,load_result_packet
    packet=run(fixture(41),tmp_path/'run',False,True,4,'exact')
    save_result_packet(packet,tmp_path/'packet')
    copied=load_result_packet(tmp_path/'packet')
    CHECK.assertEqual(qualification_report(copied,packet)['status'],'QUALIFIED')
    # Count equality is exact even where conversion to float loses a unit.
    left=dict(packet,arrays={'count':np.array([2**63+1],dtype=np.uint64)})
    right=dict(packet,arrays={'count':np.array([2**63+2],dtype=np.uint64)})
    record=qualification_report(left,right)['observables']['count']
    CHECK.assertEqual(record['status'],'REJECTED')
    CHECK.assertEqual(record['maximum_absolute_bin_error'],1.)
    (tmp_path/'packet/observables.npz').write_bytes(b'corrupt')
    with CHECK.assertRaises(ValueError):load_result_packet(tmp_path/'packet')


def test_driver_reference_packets(tmp_path):
    import sys,json
    from dataclasses import replace
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'python'))
    import kappa_corr_all_engines as k
    import shear_corr_all_engines as s
    import lya_corr_all_engines as f
    data=fixture(41)
    cases=[]
    engine='octree-2balls-omp'
    config=k.RunConfig((engine,),tmp_path/'k-exact',theta_min=2.,theta_max=100.,bins=4,multipoles=2,
        threads=1,tree_theta=0.,options=('only-2pcf','no-one-ball','no-two-balls','no-smooth-pivot'),
        verbose=0,verbose_log=0,plots=False)
    cat=k.KappaCatalog(data[0],data[1],weights=data[2])
    cases.append((engine,config,lambda c:k.run_engine_suite(cat,c)))
    engine='octree-shear-sphere-2balls-omp'
    config=s.RunConfig(statistics='2pcf',min_sep=2.,max_sep=100.,bins=4,multipoles=2,
        threads=1,tree_theta=0.,options=('no-one-ball','no-two-balls','no-smooth-pivot'),
        output_dir=tmp_path/'s-exact',verbose=0,verbose_log=0,plots=False)
    shear=s.synthetic_spherical_shear_catalog(40).normalized()
    cases.append((engine,config,lambda c:s.run_ctreeballs(shear,c,'octree-shear-sphere-2balls-omp')))
    engine='lya-2pcf-omp'
    config=f.RunConfig((engine,),tmp_path/'f-exact',threads=1,rp_max=200.,rt_max=200.,rp_bins=4,rt_bins=4,
                       plots=False,analysis=False)
    forest=f.synthetic_catalog(3,3)
    cases.append((engine,config,lambda c:f.run_engine_suite(forest,c)))
    engine='octree-3pcf-3d-omp'
    config=f.RunConfig((engine,),tmp_path/'physical-exact',threads=1,r3_max=300.,r3_bins=4,
                       multipole_lmax=2,plots=False,analysis=False)
    cases.append((engine,config,lambda c:f.run_engine_suite(forest,c)))
    for engine,config,run_driver in cases:
        run_driver(config)
        reference=config.output_dir
        candidate=reference.with_name(reference.name+'-qualified')
        run_driver(replace(config,output_dir=candidate,qualification_reference=reference))
        record=json.loads((candidate/engine/'qualification/qualification-result.json').read_text())
        CHECK.assertEqual(record['metadata']['qualification']['status'],'QUALIFIED')


def test_native_informational_commands():
    import subprocess
    root=Path(__file__).resolve().parents[2]
    for flag in ('--help','-h','-help','--clue','-c','--version'):
        result=subprocess.run([str(root/'cballs'),flag],capture_output=True,text=True,timeout=10)
        CHECK.assertEqual(result.returncode,0,(flag,result.stderr))
        CHECK.assertTrue(result.stdout.strip())


def test_reference_approximation_metadata(tmp_path):
    from cyballs import qualification_report
    from copy import deepcopy
    packet=run(fixture(41),tmp_path,False,True,4,'exact')
    # Schema control tests: every documented approximate forest flag rejects
    # reference evidence, independently of a caller's array equality.
    packet['metadata']['engine']='lya-3pcf-omp'
    for key,flag in (('lya_pivot_frontier','geometry_approximate'),('lya_2pcf','approximate'),
                     ('lya_3pcf','mu_reconstruction_approximate'),('lya_geometry','three_point_approximate')):
        reference=deepcopy(packet);reference['metadata'][key]={flag:True}
        with CHECK.assertRaisesRegex(ValueError,'reference must explicitly'):
            qualification_report(packet,reference)
