"""Qualification must match interpreted catalogs before applying error tolerances."""
from copy import deepcopy
import numpy as np
import pytest
from cyballs import cballs, qualification_report, save_result_packet, load_result_packet

OPTIONS = 'only-2pcf,no-smooth-pivot,no-one-ball,no-two-balls,no-out-Hist,weights-norm,pos-and-convergence-weight'


def fixture():
    rng = np.random.default_rng(2041)
    p = rng.normal(size=(36, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    f = 1+.1*p[:, 0]; w = np.ones(len(p))
    return np.column_stack((p, f, w, 1.001*f, 1.001*w))


def packet(tmp_path, data, columns='1,2,3,4,5', memory=False, mask=None, staged=False):
    m = cballs()
    m.set(searchMethod='octree-2balls-omp', rootDir=str(tmp_path), rangeN=1.8,
          rminHist=.02, sizeHistN=4, mChebyshev=2, numberThreads=1, theta=0,
          lengthBox=2., useLogHist=False, verbose=0, verbose_log=0, options=OPTIONS)
    if memory:
        m.set_catalog(data[:, :3], kappa=data[:, 3], weights=data[:, 4], mask=mask)
    else:
        tmp_path.mkdir(parents=True, exist_ok=True)
        path = tmp_path/'catalog.txt'; np.savetxt(path, data, fmt='%.17g')
        m.set(infile=str(path), infileformat='multi-columns-ascii', columns=columns)
    try:
        if staged: m.Run(level=['SetNumberThreads'])
        m.Run()
        return m.getResults()
    finally:
        m.clean_all()


@pytest.mark.parametrize('columns', ['1,2,3,6,5', '1,2,3,4,7', '2,1,3,4,5'])
def test_same_file_different_selected_fields_cannot_qualify(tmp_path, columns):
    data = fixture()
    ref = packet(tmp_path/'reference', data)
    candidate = packet(tmp_path/'candidate', data, columns=columns)
    assert ref['metadata']['inputs']['file_fingerprints'] == candidate['metadata']['inputs']['file_fingerprints']
    with pytest.raises(ValueError, match='matching catalog'):
        qualification_report(candidate, ref, rtol=100.)


def test_equivalent_field_selection_and_loader_are_comparable(tmp_path):
    data = fixture(); data[:, 5] = data[:, 3]
    ref = packet(tmp_path/'reference', data)
    selected = packet(tmp_path/'selected', data, columns='1,2,3,6,5')
    memory = packet(tmp_path/'memory', data, memory=True, staged=True)
    for candidate in (selected, memory):
        assert qualification_report(candidate, ref)['status'] == 'QUALIFIED'
    # Unused bytes remain provenance and cannot change semantic catalog identity.
    data[:, 6] = 123.
    unused = packet(tmp_path/'unused', data)
    assert unused['metadata']['inputs']['file_fingerprints'] != ref['metadata']['inputs']['file_fingerprints']
    assert qualification_report(unused, ref)['status'] == 'QUALIFIED'


def test_masks_and_legacy_packets_cannot_bypass_identity(tmp_path):
    data = fixture()
    ref = packet(tmp_path/'reference', data, memory=True)
    mask = np.ones(len(data), dtype=np.uint8); mask[0] = 0
    changed = packet(tmp_path/'mask', data, memory=True, mask=mask)
    with pytest.raises(ValueError, match='matching catalog'):
        qualification_report(changed, ref, rtol=100.)
    old = deepcopy(ref); del old['metadata']['inputs']['effective_catalogs']
    with pytest.raises(ValueError, match='rerun legacy'):
        qualification_report(old, old)
    save_result_packet(ref, tmp_path/'saved')
    assert qualification_report(load_result_packet(tmp_path/'saved'), ref)['status'] == 'QUALIFIED'


def test_forest_integer_ids_are_not_rounded_through_float(tmp_path):
    points=np.array([[10.,0.,0.],[12.,0.,0.],[9.,3.,0.],[11.,4.,0.]])
    def run(ids, root):
        m=cballs();m.set(searchMethod='lya-2pcf-omp',rootDir=str(root),verbose=0,verbose_log=0,
                        numberThreads=1,options='no-out-Hist,no-smooth-pivot',
                        lya2RpMax=30.,lya2RtMax=30.,lya2RpBins=2,lya2RtBins=2)
        try:
            m.set_forest_catalog(points,np.array([.5,-.2,.3,.7]),np.ones(4),ids)
            m.Run();return m.getResults()
        finally:m.clean_all()
    # Dense native IDs preserve forest membership, even beyond float precision.
    ids=np.array([2**53,2**53,2**53+2,2**53+2],dtype=np.int64)
    ref=run(ids,tmp_path/'reference');candidate=run(ids+1,tmp_path/'changed')
    for name in ref['arrays']:np.testing.assert_array_equal(ref['arrays'][name],candidate['arrays'][name])
    assert qualification_report(candidate,ref)['status'] == 'QUALIFIED'
    altered = ids.copy(); altered[1] += 1
    changed = run(altered,tmp_path/'changed-partition')
    with pytest.raises(ValueError,match='matching catalog'):
        qualification_report(changed,ref,rtol=100.)
