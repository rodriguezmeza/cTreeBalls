"""Exact separation-plane quadtree engine regressions."""
import json
import math
import re

import numpy as np
import pytest

import test_lya_forest_omp as oracle
import test_lya_forest_los_tree as los

METHOD = "lya-rprt-quadtree-2pcf-omp"


def stats(log):
    line = next(x for x in log.splitlines() if "Ly-alpha rp-rt quadtree:" in x)
    return {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)", line)}


def independent_grid(points, rp_bins, rt_bins, maximum=30.):
    """Small O(N^2) Cartesian oracle, independent of both native trees."""
    result = np.zeros((rp_bins, rt_bins, 2))
    for i, a in enumerate(points):
        for b in points[i+1:]:
            if a[5] == b[5] or np.linalg.norm(a[:3]-b[:3]) >= math.hypot(maximum, maximum):
                continue
            chi1, chi2 = np.linalg.norm(a[:3]), np.linalg.norm(b[:3])
            cosine = np.clip(np.dot(a[:3]/chi1, b[:3]/chi2), -1., 1.)
            rp = abs(chi1-chi2)*math.sqrt(max(0., (1+cosine)/2))
            rt = (chi1+chi2)*math.sqrt(max(0., (1-cosine)/2))
            if not (0 <= rp < maximum and 0 <= rt < maximum):
                continue
            p, t = int(rp/maximum*rp_bins), int(rt/maximum*rt_bins)
            result[p, t] += [a[4]*b[4]*a[3]*b[3], a[4]*b[4]]
    return result.reshape(-1, 2)


@pytest.mark.parametrize("case", ["ordinary", "bent", "long", "antipodal", "zero", "one", "cancellation", "boundaries"])
def test_reference_oracle_and_thread_order(tmp_path, case):
    points = los.make_forests(7, 23, seed=427)
    if case == "bent": points[::3, 1:3] += 19
    elif case == "long": points[::2, :3] *= 20
    elif case == "antipodal": points[::2, :3] *= -1
    elif case == "zero": points[::3, 4] = 0
    elif case == "one": points[:, 5] = 1
    elif case == "cancellation": points[:, 3] = np.where(np.arange(len(points)) % 2, 1., -1.)
    elif case == "boundaries":
        points = np.array([(100,0,0,.2,1,0), (106,0,0,-.1,2,1),
            (130,0,0,.3,0,2), (130-1e-5,0,1e-5,.4,1,3),
            (130+1e-5,0,0,.1,1,4), (100,1e-5,0,.1,1,5),
            (-100,0,0,.2,1,6), (105,7,0,.6,1,7), (95,-7,0,.1,1,8)])
    catalog = los.catalog_file(tmp_path, points)
    los.run(catalog, tmp_path/"reference", "lya-2pcf-omp")
    log = los.run(catalog, tmp_path/"one", METHOD)
    los.run(catalog, tmp_path/"four", METHOD, 4)
    los.compare_products(tmp_path/"reference", tmp_path/"one", (2,))
    assert (tmp_path/"one"/los.PRODUCTS[2]).read_bytes() == (tmp_path/"four"/los.PRODUCTS[2]).read_bytes()
    if case != "boundaries":
        oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/"one"/los.PRODUCTS[2]),
                                     oracle.oracle_2pcf(points), "independent quadtree pairs")
    assert stats(log)["global_pixel_octree"] == 0


@pytest.mark.parametrize("grid", [(1,1), (1,7), (7,1), (3,5), (8,16), (50,50), (11,6)])
def test_rectangular_odd_and_singleton_grids(tmp_path, grid):
    points = los.make_forests(4, 7, seed=93)
    catalog = los.catalog_file(tmp_path, points)
    rp, rt = grid
    log = los.run(catalog, tmp_path/"quad", METHOD, 3,
                  extra=(f"lya2RpBins={rp}", f"lya2RtBins={rt}"))
    table = np.loadtxt(tmp_path/"quad"/los.PRODUCTS[2], ndmin=2)
    np.testing.assert_allclose(table[:, -2:], independent_grid(points, rp, rt), rtol=2e-12, atol=2e-12)
    assert rp*rt <= stats(log)["bin_nodes"] <= 2*rp*rt-1


def test_whole_cell_bin_acceptance_avoids_pixel_work(tmp_path):
    rows = []
    for f in range(5):
        for i, chi in enumerate(np.linspace(100+7.3*f, 100.01+7.3*f, 64)):
            rows.append((chi,0,0,math.sin(i)*.2,1+(i%5)*.1,f))
    points = np.array(rows)
    catalog = los.catalog_file(tmp_path, points)
    log = los.run(catalog, tmp_path/"quad", METHOD, 4)
    los.run(catalog, tmp_path/"ref", "lya-2pcf-omp")
    los.compare_products(tmp_path/"quad", tmp_path/"ref", (2,))
    count = stats(log)
    assert count["single_bin_products"] > 0
    assert count["represented_pairs"] > 64*64
    assert count["bin_descents"] > 0
    assert count["pixel_tests"] < 5*4//2*64*64


def test_coincident_pixels_do_not_need_cartesian_octree(tmp_path):
    points = np.array([(100,0,0,.5,2,f) for f in range(3) for _ in range(12)])
    catalog = los.catalog_file(tmp_path, points)
    log = los.run(catalog, tmp_path/"quad", METHOD, 4)
    oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/"quad"/los.PRODUCTS[2]),
                                 oracle.oracle_2pcf(points), "coincident forest products")
    assert stats(log)["single_bin_products"] > 0
    count = re.search(r"# distinct-forest unordered pairs: (\d+)", (tmp_path/"quad"/los.PRODUCTS[2]).read_text())
    assert int(count[1]) == 3*12*12


@pytest.mark.parametrize("options", [("lyaPivotRadius=.1",), ("lyaScanLevel=1",)])
def test_frontier_controls_are_rejected(tmp_path, options):
    catalog = los.catalog_file(tmp_path, oracle.POINTS)
    with pytest.raises(AssertionError, match="pivot frontier"):
        los.run(catalog, tmp_path/"bad", METHOD, extra=options)


def test_memory_rejection_then_recovery(tmp_path, monkeypatch):
    catalog = los.catalog_file(tmp_path, los.make_forests(8, 16))
    monkeypatch.setenv("CBALLS_MEMORY_BUDGET_MB", ".1")
    with pytest.raises(AssertionError, match="workspace|memory|budget|plan"):
        los.run(catalog, tmp_path/"limited", METHOD)
    monkeypatch.delenv("CBALLS_MEMORY_BUDGET_MB")
    assert stats(los.run(catalog, tmp_path/"recovered", METHOD))["bin_nodes"] > 0


def test_cython_registry_metadata_and_repeated_calls(tmp_path):
    import cyballs
    import lya_corr_all_engines as driver
    assert driver.LYA_ENGINES[METHOD].orders == (2,)
    points = los.make_forests(5, 12)
    b = cyballs.cballs()
    values = []
    try:
        for i, method in enumerate((METHOD, "lya-2pcf-omp", METHOD)):
            b.set(dict(searchMethod=method, rootDir=str(tmp_path/str(i)), numberThreads=3,
                usePeriodic=False, useLogHist=False, rangeN=30., rminHist=.1, sizeHistN=4,
                lya2RpMax=30., lya2RtMax=30., lya2RpBins=5, lya2RtBins=6, lya2Kernel=0,
                options="no-smooth-pivot", verbose=0, verbose_log=0))
            b.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            b.Run(level=["MainLoop"])
            metadata = b.getRunMetadata()
            if method == METHOD:
                assert metadata["engine_id"] == 213
                assert metadata["lya_rprt_quadtree"]["exact_bin_containment"]
                assert not metadata["lya_rprt_quadtree"]["global_pixel_octree"]
                assert metadata["lya_rprt_quadtree"]["two_ball_traversal"]
                assert metadata["lya_rprt_quadtree"]["two_ball_leaf_block_max"] == 2
                assert metadata["lya_2pcf"]["kernel"] == 1
                assert not metadata["lya_2pcf"]["approximate"]
            values.append(np.loadtxt(tmp_path/str(i)/los.PRODUCTS[2]))
            b.struct_cleanup()
        np.testing.assert_allclose(values[0], values[1], rtol=2e-12, atol=2e-12)
        np.testing.assert_array_equal(values[0], values[2])
    finally:
        b.struct_cleanup()


@pytest.mark.parametrize("case", ["straight", "perturbed-los", "zero", "shuffled"])
def test_compact_sparse_tiles_match_independent_pairs(tmp_path, case):
    rng = np.random.default_rng(1706)
    rows = []
    for f in range(4):
        direction = np.array([1., .002*f, .001*f])
        direction /= np.linalg.norm(direction)
        for i, chi in enumerate(np.linspace(100., 270., 43)):
            xyz = chi*direction
            if case == "perturbed-los":
                xyz[2] += 1e-8*(i % 3)
            weight = 0. if case == "zero" and i % 3 == 0 else 1.+.1*(i % 7)
            rows.append((*xyz, np.sin(i+f)*.2, weight, f))
    points = np.array(rows)
    if case == "shuffled": rng.shuffle(points)
    catalog = los.catalog_file(tmp_path, points)
    log = los.run(catalog, tmp_path/"quad", METHOD)
    los.run(catalog, tmp_path/"many", METHOD, 4)
    los.run(catalog, tmp_path/"reference", "lya-2pcf-omp")
    los.compare_products(tmp_path/"reference", tmp_path/"quad", (2,))
    oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/"quad"/los.PRODUCTS[2]),
                                 oracle.oracle_2pcf(points), "sparse window oracle")
    assert (tmp_path/"quad"/los.PRODUCTS[2]).read_bytes() == (tmp_path/"many"/los.PRODUCTS[2]).read_bytes()
    counts = stats(log)
    assert counts["compact_tiles"] > 0
    assert counts["cached_pixel_pairs"] > 0
    assert counts["pixel_tests"] < 6*43*43


@pytest.mark.parametrize("grid", [(3, 5), (8, 16), (1, 7)])
def test_dense_subdivision_keeps_boundary_bins_and_signed_sums(tmp_path, grid):
    # Fine sampling selects the cell-product path rather than sparse windows.
    # Almost collinear forests and grid edges stress inherited bin regions.
    rng = np.random.default_rng(143)
    points = []
    for f in range(3):
        for i, chi in enumerate(np.linspace(98.+3*f, 101.+3*f, 65)):
            points.append((chi, chi*f*.00001, 0., (-1.)**i*.2, 1.+.1*(i%5), f))
    points = np.array(points)
    rng.shuffle(points)
    cat = los.catalog_file(tmp_path, points)
    args = (f"lya2RpBins={grid[0]}", f"lya2RtBins={grid[1]}")
    log = los.run(cat, tmp_path/"quad", METHOD, extra=args)
    los.run(cat, tmp_path/"reference", "lya-2pcf-omp", extra=args)
    los.compare_products(tmp_path/"reference", tmp_path/"quad", (2,))
    expected = independent_grid(points, *grid)
    np.testing.assert_allclose(np.loadtxt(tmp_path/"quad"/los.PRODUCTS[2])[:, -2:],
                               expected, rtol=2e-12, atol=2e-12)
    assert stats(log)["single_bin_products"] > 0


@pytest.mark.parametrize("axis", ["rp", "rt"])
def test_compact_tile_boundary_fallback_matches_reference(tmp_path, axis):
    # Put distinct pixels within roundoff of every internal bin edge and the
    # outer cutoff; interval uncertainty must retain the reference pixel path.
    rows = []
    edges = np.arange(6., 31., 6.) if axis == "rp" else np.arange(5., 31., 5.)
    for f, edge in enumerate(edges, start=1):
        for i, offset in enumerate((-1e-11, 0., 1e-11)):
            if axis == "rp":
                angle = .01*(i+1)
                chi = 100.+(edge+offset)/math.cos(angle/2.)
            else:
                chi = 100.+4.*i
                angle = 2.*math.asin((edge+offset)/(100.+chi))
            xyz = (chi*math.cos(angle),chi*math.sin(angle),0.)
            rows.append((*xyz, .1*(-1.)**i, 1.+.2*i, f))
    rows.append((100.,0.,0.,-.3,1.,0))
    points = np.array(rows)
    cat = los.catalog_file(tmp_path, points)
    log = los.run(cat,tmp_path/"quad",METHOD)
    los.run(cat,tmp_path/"reference","lya-2pcf-omp")
    los.compare_products(tmp_path/"reference",tmp_path/"quad",(2,))
    counts = stats(log)
    assert counts["compact_tiles"] > 0
    assert counts["cached_pixel_pairs"] < counts["pixel_tests"]


@pytest.mark.parametrize("axis", ["rp", "rt", "both"])
@pytest.mark.parametrize("grid", [(5,6), (3,5), (1,7)])
def test_slop_bounds_for_every_encoded_pair(tmp_path, axis, grid):
    # Powers of two encode each cross-forest pair in one unique denominator
    # bit. Decode the output to check placement, exclusions and duplication
    # independently of the tree's counters or approximate representative.
    npix, nq = 4, 8
    dp, dt = 30./grid[0], 30./grid[1]
    rp_slop = .1 if axis in ("rp", "both") else 0.
    rt_slop = .1 if axis in ("rt", "both") else 0.
    parallel = dp if rp_slop and grid[0]>1 else (1.4*dp if grid[0]>1 else 10.)
    pchi = np.linspace(100.,100.2,npix)
    qchi = np.linspace(100.1+parallel-.1,100.1+parallel+.1,nq)
    transverse = dt if rt_slop else .4*dt
    angle = 2.*math.asin(transverse/(pchi.mean()+qchi.mean()))
    rows = [(chi,0.,0.,1.,2.**i,0) for i,chi in enumerate(pchi)]
    rows += [(chi*math.cos(angle),chi*math.sin(angle),0.,1.,2.**(npix*j),1)
             for j,chi in enumerate(qchi)]
    points = np.array(rows)
    cat = los.catalog_file(tmp_path, points)
    options = (f"lya2RpBins={grid[0]}",f"lya2RtBins={grid[1]}",
               f"lya2RpSlop={rp_slop}",f"lya2RtSlop={rt_slop}")
    log = los.run(cat,tmp_path/"slop",METHOD,extra=options)
    los.run(cat,tmp_path/"many",METHOD,4,extra=options)
    assert (tmp_path/"slop"/los.PRODUCTS[2]).read_bytes() == (tmp_path/"many"/los.PRODUCTS[2]).read_bytes()
    seen = set()
    for row in np.loadtxt(tmp_path/"slop"/los.PRODUCTS[2],ndmin=2):
        bp,bt = map(int,row[:2])
        code = int(row[-1]);assert float(code)==row[-1]
        for bit in range(npix*nq):
            if not (code >> bit)&1:continue
            i,j = bit%npix,bit//npix
            assert (i,j) not in seen
            seen.add((i,j))
            a,b = points[i,:3],points[npix+j,:3]
            ca,cb = np.linalg.norm(a),np.linalg.norm(b)
            c = np.clip(np.dot(a/ca,b/cb),-1.,1.)
            rp = abs(ca-cb)*math.sqrt((1+c)/2)
            rt = (ca+cb)*math.sqrt((1-c)/2)
            assert (bp-rp_slop)*dp-1e-9 <= rp <= (bp+1+rp_slop)*dp+1e-9
            assert (bt-rt_slop)*dt-1e-9 <= rt <= (bt+1+rt_slop)*dt+1e-9
            if rp_slop==0:assert bp==int(rp/30.*grid[0])
            if rt_slop==0:assert bt==int(rt/30.*grid[1])
    assert len(seen)==npix*nq
    if rt_slop or (rp_slop and grid[0]>1):
        assert stats(log)["slop_products"]>0
        assert stats(log)["approximate_pairs"]>0


@pytest.mark.parametrize("slops",[(.1,0.),(0.,.1),(.25,.25),(1.,1.)])
def test_slop_conserves_global_sums_and_exact_outer_cuts(tmp_path, slops):
    points = los.make_forests(6,67,seed=179)
    points[::4,4] = 0.
    points[::11,:3] *= 4.
    cat = los.catalog_file(tmp_path,points)
    los.run(cat,tmp_path/"exact",METHOD)
    los.run(cat,tmp_path/"slop",METHOD,4,extra=(f"lya2RpSlop={slops[0]}",f"lya2RtSlop={slops[1]}"))
    los.run(cat,tmp_path/"slop-one",METHOD,1,extra=(f"lya2RpSlop={slops[0]}",f"lya2RtSlop={slops[1]}"))
    assert (tmp_path/"slop"/los.PRODUCTS[2]).read_bytes()==(tmp_path/"slop-one"/los.PRODUCTS[2]).read_bytes()
    products = [tmp_path/d/los.PRODUCTS[2] for d in ("exact","slop")]
    counts = [int(re.search(r"unordered pairs: (\d+)",p.read_text())[1]) for p in products]
    assert counts[0]==counts[1]
    np.testing.assert_allclose(np.loadtxt(products[0])[:,-2:].sum(axis=0),
        np.loadtxt(products[1])[:,-2:].sum(axis=0),rtol=3e-12,atol=3e-10)
    grids=[]
    for path in products:
        table=np.loadtxt(path);grid=np.zeros((5,6,2))
        grid[table[:,0].astype(int),table[:,1].astype(int)]=table[:,-2:]
        grids.append(grid)
    for dimension,slop in enumerate(slops):
        if slop==0:
            np.testing.assert_allclose(grids[0].sum(axis=1-dimension),grids[1].sum(axis=1-dimension),
                                       rtol=3e-12,atol=3e-10)


@pytest.mark.parametrize("value",[-.1,1.01,float("nan"),float("inf")])
def test_invalid_slop_rejected(tmp_path,value):
    cat = los.catalog_file(tmp_path,oracle.POINTS)
    with pytest.raises(AssertionError,match="pair slops finite|invalid numeric value for lya2RpSlop"):
        los.run(cat,tmp_path/"bad",METHOD,extra=(f"lya2RpSlop={value}",))


def test_slop_metadata_then_exact_recovery(tmp_path):
    import cyballs
    points = los.make_forests(4,31)
    b = cyballs.cballs();values=[]
    try:
        for i,slop in enumerate((0.,.1,0.)):
            b.set(dict(searchMethod=METHOD,rootDir=str(tmp_path/str(i)),numberThreads=3,
                usePeriodic=False,useLogHist=False,rangeN=30.,rminHist=.1,sizeHistN=4,
                lya2RpMax=30.,lya2RtMax=30.,lya2RpBins=5,lya2RtBins=6,
                lya2RpSlop=slop,lya2RtSlop=slop,options="no-smooth-pivot",verbose=0,verbose_log=0))
            b.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            b.Run(level=["MainLoop"])
            meta=b.getRunMetadata()
            assert meta["lya_2pcf"]["approximate"]==bool(slop)
            assert meta["lya_rprt_quadtree"]["exact_bin_containment"]==(not slop)
            assert meta["lya_rprt_quadtree"]["outer_cuts_exact"]
            values.append(b.getForestResults()["arrays"]["pair_numerator"])
            b.struct_cleanup()
        np.testing.assert_array_equal(values[0],values[2])
    finally:b.struct_cleanup()


@pytest.mark.parametrize("grid", [(5, 6), (7, 3), (3, 7)])
def test_two_ball_blocks_decode_every_pair(tmp_path, grid):
    # A 3x7 product cannot fit one parallel bin, but several 2x2 radial
    # balls do fit. Unique denominator bits expose omission/duplication and
    # any wrong bin, independently of field cancellation or work counters.
    rows = [(chi, 0., 0., (-1.)**i, 2.**i, 0)
            for i, chi in enumerate(100.+1.2*np.arange(3))]
    angle = .035
    rows += [(chi*math.cos(angle), chi*math.sin(angle), 0., (-1.)**j,
              2.**(3*j), 1) for j, chi in enumerate(104.3+1.2*np.arange(7))]
    points = np.array(rows)
    cat = los.catalog_file(tmp_path, points)
    args = (f"lya2RpBins={grid[0]}", f"lya2RtBins={grid[1]}")
    log = los.run(cat, tmp_path/"quad", METHOD, extra=args)
    los.run(cat, tmp_path/"four", METHOD, 4, extra=args)
    table = np.loadtxt(tmp_path/"quad"/los.PRODUCTS[2], ndmin=2)
    expected = independent_grid(points, *grid)
    np.testing.assert_array_equal(table[:, -2:], expected)
    assert int(table[:, -1].sum()) == 2**21-1
    assert stats(log)["two_ball_products"] > 0
    assert stats(log)["two_ball_pairs"] > 0
    assert (tmp_path/"quad"/los.PRODUCTS[2]).read_bytes() == (tmp_path/"four"/los.PRODUCTS[2]).read_bytes()


def test_two_ball_split_handles_unequal_forest_extents(tmp_path):
    rows = []
    for f, (n, lo, hi) in enumerate([(65, 100., 115.), (97, 95., 165.), (25, 97., 142.)]):
        angle = .025*f
        for i, chi in enumerate(np.linspace(lo, hi, n)):
            rows.append((chi*math.cos(angle), chi*math.sin(angle), 0.,
                         .2*math.sin(i+f), 0. if i%5==0 else 1.+.1*(i%3), f))
    points = np.array(rows)
    np.random.default_rng(1710).shuffle(points)
    cat = los.catalog_file(tmp_path, points)
    log = los.run(cat, tmp_path/"quad", METHOD, 4)
    los.run(cat, tmp_path/"ref", "lya-2pcf-omp")
    los.compare_products(tmp_path/"ref", tmp_path/"quad", (2,))
    oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/"quad"/los.PRODUCTS[2]),
                                 oracle.oracle_2pcf(points), "two-ball unequal extents")
    counts = stats(log)
    assert counts["split_one"] > 0 and counts["split_both"] > 0
    assert counts["product_splits"] == counts["split_one"]+counts["split_both"]
    assert counts["two_ball_pairs"] > 0


def test_two_ball_sparse_gate_retains_exact_pixel_tiles(tmp_path):
    rows = []
    for f in range(3):
        angle = .002*f
        for i, chi in enumerate(np.linspace(100., 900., 43)):
            rows.append((chi*math.cos(angle), chi*math.sin(angle), 0., .2*math.sin(i), 1., f))
    points = np.array(rows)
    cat = los.catalog_file(tmp_path, points)
    log = los.run(cat, tmp_path/"quad", METHOD)
    oracle.assert_histogram_close(oracle.read_2pcf(tmp_path/"quad"/los.PRODUCTS[2]),
                                 oracle.oracle_2pcf(points), "sparse two-ball scheduling")
    assert stats(log)["two_ball_tests"] == 0
    assert stats(log)["pixel_tests"] > 0


def construction_stats(log):
    result = {}
    for line in log.splitlines():
        if "Ly-alpha 2PCF forest build:" in line or "Ly-alpha rp-rt construction:" in line:
            for key, value in re.findall(r"(\w+)=([\d.eE+-]+)", line):
                result[key] = float(value)
    return result


@pytest.mark.parametrize("large_pixels,slop", [(True,0.), (False,0.), (False,.1)])
def test_parallel_construction_matches_serial_bits(tmp_path, large_pixels, slop):
    import cyballs
    # Large pixel input exercises chunk sort, co-ranked merge levels and
    # independently built forest subtrees. Large odd grids exercise tasks,
    # reserved compact node slots and the optional slop leaf directory.
    nf, npix = (512,137) if large_pixels else (8,19)
    index = np.arange(nf*npix)
    forest = index//npix
    radial_index = (index%npix)//3 if large_pixels else index%npix
    chi = 100.+1000.*(forest//4)+4.13*(forest%4)+.1*radial_index/npix
    points = np.column_stack((chi,np.zeros((len(index),2)),.2*np.sin(index),
                              1.+.1*(index%7),forest*37-2**35))
    np.random.default_rng(1716).shuffle(points)
    bins = (5,6) if large_pixels else (129,131)
    values=[];logs=[]
    for label,nt,serial in [("serial",4,True),("parallel",4,False),("one",1,False)]:
        b=cyballs.cballs()
        try:
            out=tmp_path/label
            b.set(dict(searchMethod=METHOD,rootDir=str(out),numberThreads=nt,
                usePeriodic=False,useLogHist=False,rangeN=30.,rminHist=.1,sizeHistN=4,
                lya2RpMax=30.,lya2RtMax=30.,lya2RpBins=bins[0],lya2RtBins=bins[1],
                lya2RpSlop=slop,lya2RtSlop=slop,lya2Kernel=1,verbose=0,verbose_log=2,
                options="no-out-Hist,no-smooth-pivot"+(",lya-serial-build" if serial else "")))
            b.set_forest_catalog(points[:,:3],points[:,3],points[:,4],points[:,5].astype(np.int64))
            b.Run(level=["MainLoop"])
            values.append(b.getForestResults()["arrays"])
            meta=b.getRunMetadata()["lya_rprt_quadtree"]
            assert meta["implicit_rp_rt_projection"]
            assert meta["construction_mode"] == ("serial override" if serial else "automatic")
        finally:b.struct_cleanup()
        files=list(out.rglob("*.log"))
        assert files
        logs.append(construction_stats("\n".join(p.read_text() for p in files)))
    for a in values[1:]:
        for key in a:np.testing.assert_array_equal(a[key],values[0][key])
    assert logs[0]["sort_threads"]==logs[0]["node_threads"]==1
    assert logs[0]["pixel_threads"]==logs[0]["bin_threads"]==1
    assert logs[0]["projection_threads"]==4
    assert logs[2]["sort_threads"]==logs[2]["node_threads"]==1
    assert logs[2]["pixel_threads"]==logs[2]["bin_threads"]==1
    if large_pixels:
        assert logs[1]["sort_threads"]==logs[1]["node_threads"]==4
        assert logs[1]["pixel_threads"]==4
        assert logs[1]["sort_scratch_bytes"]==len(points)*np.dtype(np.intp).itemsize
        assert logs[1]["bin_threads"]==1
    else:
        assert logs[1]["sort_threads"]==logs[1]["node_threads"]==1
        assert logs[1]["pixel_threads"]==1
        assert logs[1]["bin_threads"]==4
