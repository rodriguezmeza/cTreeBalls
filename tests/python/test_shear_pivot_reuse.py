"""Regression/calibration fixtures for opt-in native octree pivot reuse."""
import os
import re
from unittest.mock import patch

import numpy as np

import test_shear_sphere_octree_omp as reference

ENGINE = os.environ.get("CBALLS_SHEAR_SPHERE_ENGINE", "octree-shear-sphere-2balls-omp")
BASE = "no-out-Hist,no-smooth-pivot,only-3pcf"
REUSE = BASE + ",shear-pivot-reuse"


def compare(a, b, rtol=2.e-11, atol=3.e-12):
    for key in a:
        np.testing.assert_allclose(b[key], a[key], rtol=rtol, atol=atol,
                                   err_msg=key)


def test_exact_limit_and_fallbacks():
    catalog = reference.fixture(80)
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "1e-12"}):
        for use_log in (False, True):
            exact = reference.run_native(*catalog, 1, options=BASE+",no-one-ball",
                                         use_log=use_log, engine=ENGINE)
            actual = reference.run_native(*catalog, 4, options=REUSE,
                                          use_log=use_log, engine=ENGINE)
            compare(exact, actual)
            direct = reference.oracle(*catalog, use_log=use_log)
            compare({key: direct[key] for key in actual}, actual,
                    rtol=1.e-10, atol=2.e-11)
    for token in ("no-one-ball", "no-two-balls", "only-2pcf"):
        options = "no-out-Hist,no-smooth-pivot,"+token
        a = reference.run_native(*catalog, 1, options=options, engine=ENGINE)
        b = reference.run_native(*catalog, 1, options=options+",shear-pivot-reuse",
                                 engine=ENGINE)
        for key in a:
            np.testing.assert_array_equal(a[key], b[key])
    a = reference.run_native(*catalog, 1, options=BASE, engine=ENGINE)
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0"}):
        b = reference.run_native(*catalog, 1, options=REUSE, engine=ENGINE)
    compare(a, b, 0, 0)


def clustered_catalog():
    rng = np.random.default_rng(881)
    centers = rng.normal(size=(12, 3))
    centers /= np.linalg.norm(centers, axis=1)[:, None]
    position = np.repeat(centers, 32, axis=0) + rng.normal(size=(384, 3))*1.e-5
    position /= np.linalg.norm(position, axis=1)[:, None]
    gamma = .02*(rng.normal(size=384) + 1j*rng.normal(size=384))
    weight = rng.uniform(.1, 2., 384)
    weight[::31] = 0.0
    return position, gamma, weight


def test_aggregate_pivots_and_determinism():
    catalog = clustered_catalog()
    with patch.object(reference, "ENGINE", ENGINE), patch.dict(
            os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.1"}):
        exact = reference.run_native(*catalog, 1, options=BASE+",no-one-ball")
        result, profile = reference.run_native_with_profile(*catalog, 1, REUSE)
        count = sum(map(int, re.findall(r"aggregate_pivots=(\d+)", profile)))
        if not (count > 0):
            raise AssertionError(profile)
        # This fixture has nearly coincident groups and deliberately weak modes.
        for key in result:
            relative_l2 = np.linalg.norm(result[key]-exact[key])/np.linalg.norm(exact[key])
            if not (relative_l2 < .01):
                raise AssertionError((key, relative_l2))
        threaded = reference.run_native(*catalog, 4, options=REUSE)
        compare(result, threaded, 0, 0)
        combined = reference.run_native(*catalog, 1,
            options="no-out-Hist,no-smooth-pivot,shear-pivot-reuse")
        compare(result, combined, 0, 0)
        pairs = reference.run_native(*catalog, 1,
            options="no-out-Hist,no-smooth-pivot,only-2pcf")
        compare(pairs, combined, 0, 0)


def test_cross_masks_and_polar_geometry():
    first = reference.fixture(42)
    position, gamma, weight = reference.fixture(51)
    position[:4] = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 1.e-6, 0]]
    position /= np.linalg.norm(position, axis=1)[:, None]
    second = position, gamma, weight
    masks = [np.arange(len(c[0])) % 4 != 0 for c in (first, second)]
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "1e-10"}), \
            patch.object(reference, "RMAX", 2.0):
        for catalogs, order, selected_masks in (
                ([first, second], None, masks),
                ([first, second, first], "1,2,1", masks+[masks[0]])):
            for nmax in (2, 5):
                with patch.object(reference, "NMAX", nmax):
                    kwargs = dict(engine=ENGINE, masks=selected_masks,
                                  catalog_order=order, use_log=False)
                    a = reference.run_native_catalogs(catalogs, 1,
                        options=BASE+",read-mask,no-one-ball", **kwargs)
                    b = reference.run_native_catalogs(catalogs, 4,
                        options=REUSE+",read-mask", **kwargs)
                    compare(a, b, rtol=2.e-9, atol=1.e-10)


def test_inherited_rings_and_masked_aggregation():
    rng = np.random.default_rng(773)
    centers, _, _ = reference.fixture(64)
    position = np.repeat(centers, 32, axis=0) + rng.normal(size=(2048, 3))*.004
    position /= np.linalg.norm(position, axis=1)[:, None]
    gamma = .02*(rng.normal(size=2048) + 1j*rng.normal(size=2048))
    weight = rng.uniform(.5, 1.5, 2048)
    with patch.object(reference, "ENGINE", ENGINE), patch.dict(
            os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.5"}):
        actual, profile = reference.run_native_with_profile(
            position, gamma, weight, 4, REUSE)
        if not (sum(map(int, re.findall(r"ancestor_merges=(\d+)", profile))) > 0):
            raise AssertionError(profile)
        exact = reference.run_native(position, gamma, weight, 1,
                                     options=BASE+",no-one-ball")
        for key in actual:
            error = np.linalg.norm(actual[key]-exact[key])/np.linalg.norm(exact[key])
            if not (error < .05):
                raise AssertionError((key, error))
    catalog = clustered_catalog()
    masks = [(np.arange(384) % 7 == 0).astype(np.int32)]
    with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": "0.03"}):
        for catalogs, selected_masks in (([catalog], masks),
                                         ([catalog, catalog], masks*2)):
            exact = reference.run_native_catalogs(catalogs, 1,
                options=BASE+",no-one-ball,read-mask", masks=selected_masks, engine=ENGINE)
            actual = reference.run_native_catalogs(catalogs, 4,
                options=REUSE+",read-mask", masks=selected_masks, engine=ENGINE)
            for key in actual:
                error = np.linalg.norm(actual[key]-exact[key])/np.linalg.norm(exact[key])
                if not (error < .02):
                    raise AssertionError((key, error))


def test_validation_and_smoothing_fallback():
    catalog = reference.smooth_fixture()
    options = "no-out-Hist,smooth-pivot,only-3pcf"
    kwargs = dict(engine=ENGINE, rsmooth=reference.SMOOTH_RMIN_ARCMIN)
    a = reference.run_native(*catalog, 1, options=options, **kwargs)
    b = reference.run_native(*catalog, 1, options=options+",shear-pivot-reuse", **kwargs)
    compare(a, b, 0, 0)
    try:
        reference.run_native(*catalog, 1, options=REUSE+",legacy-one-ball", engine=ENGINE)
    except Exception as error:
        if not ("legacy-one-ball" in str(error)):
            raise AssertionError(error)
    else:
        raise AssertionError("reuse unexpectedly accepted in legacy mode")
    for text in ("nan", "inf", "-1", "3.1", "", "0.1junk"):
        with patch.dict(os.environ, {"CBALLS_SHEAR_PIVOT_TOL": text}):
            try:
                reference.run_native(*catalog, 1, options=REUSE, engine=ENGINE)
            except Exception as error:
                if not ("CBALLS_SHEAR_PIVOT_TOL" in str(error)):
                    raise AssertionError(error)
            else:
                raise AssertionError("invalid tolerance accepted: "+text)


def test_hierarchical_radial_completion():
    """Each pivot/bin pair is owned once, including partial parent results."""
    rng = np.random.default_rng(921)
    centers = np.array([[0,0,1], [.035,0,1], [.12,.01,1],
                        [.25,-.02,1], [.55,.12,1]], dtype=float)
    centers /= np.linalg.norm(centers, axis=1)[:,None]
    position = np.repeat(centers,64,axis=0) + rng.normal(scale=.0003,size=(320,3))
    position /= np.linalg.norm(position,axis=1)[:,None]
    gamma = .02*np.exp(1j*np.arange(320)*.03)
    weight = rng.uniform(.5,1.5,320)
    with patch.multiple(reference, ENGINE=ENGINE, BINS=8, NMAX=3,
                        RMIN=.0001, RMAX=.65), patch.dict(os.environ,
                        CBALLS_SHEAR_PIVOT_TOL='.3', CBALLS_SHEAR_BIN_THETA='.2'):
        actual, profile = reference.run_native_with_profile(position,gamma,weight,1,REUSE)
        counts = {key: sum(map(int, re.findall(r"\b"+key+r"=(\d+)",profile)))
                  for key in ('partial_reductions','radial_pairs','represented_pairs')}
        if not (counts['partial_reductions'] > 0
                and counts['radial_pairs'] < counts['represented_pairs']
                and counts['represented_pairs'] == 320*8*8):
            raise AssertionError(counts)
        threaded = reference.run_native(position,gamma,weight,4,REUSE)
        compare(actual,threaded,0,0)
        exact = reference.run_native(position,gamma,weight,1,BASE+',no-one-ball')
        for key in ('upsilon','window'):
            error = np.linalg.norm(actual[key]-exact[key])/np.linalg.norm(exact[key])
            if not error < .002: raise AssertionError((key,error))
        # Compact groups leave some angular windows almost singular. Qualify
        # the corrected observable only on jointly well-conditioned windows.
        orders = np.arange(-3,4)
        valid = np.zeros((8,8),bool)
        for i in range(8):
            for j in range(8):
                matrices = [x['window'][orders[:,None]-orders[None,:]+6,i,j]
                            for x in (actual,exact)]
                valid[i,j] = all(abs(m[3,3]) > 0 and np.linalg.cond(m) < 1.e4
                                 for m in matrices)
        if valid.sum() < 10: raise AssertionError('insufficient qualified windows')
        a,b = actual['multipoles'][:,:,valid],exact['multipoles'][:,:,valid]
        error = np.linalg.norm(a-b)/np.linalg.norm(b)
        if not error < .003: raise AssertionError(('corrected multipoles',error))


def test_bin_slop_cutoffs_and_validation():
    catalog = reference.fixture(240)
    # Total monopole weight cannot migrate across the radial range limits,
    # even when assignments across internal bin edges are approximate.
    for use_log in (False,True):
        with patch.dict(os.environ, CBALLS_SHEAR_PIVOT_TOL='3',
                        CBALLS_SHEAR_BIN_THETA='.5'):
            exact = reference.run_native(*catalog,1,BASE+',no-one-ball',
                                         engine=ENGINE,use_log=use_log)
            actual = reference.run_native(*catalog,4,REUSE,
                                          engine=ENGINE,use_log=use_log)
            np.testing.assert_allclose(actual['window'][2*reference.NMAX].sum(),
                                       exact['window'][2*reference.NMAX].sum(),
                                       rtol=3.e-13,atol=1.e-10)
    for text in ('nan','inf','-1','1.1','','0.1junk'):
        with patch.dict(os.environ, CBALLS_SHEAR_BIN_THETA=text):
            try:
                reference.run_native(*catalog,1,REUSE,engine=ENGINE)
            except Exception as error:
                if 'CBALLS_SHEAR_BIN_THETA' not in str(error): raise
            else:
                raise AssertionError('invalid bin slop accepted: '+text)


def test_reuse_provenance_is_frozen():
    import tempfile
    position,gamma,weight = reference.fixture(40)
    model = reference.cballs()
    try:
        model.set(dict(searchMethod=ENGINE,rootDir=tempfile.mkdtemp(),
            usePeriodic=False,useLogHist=False,rminHist=.03,rangeN=1.75,
            sizeHistN=3,sizeHistPhi=16,mChebyshev=2,numberThreads=1,
            theta=1,verbose=0,verbose_log=0,options=REUSE))
        model.set_catalog(position,gamma1=gamma.real,gamma2=gamma.imag,weights=weight)
        with patch.dict(os.environ,CBALLS_SHEAR_PIVOT_TOL='.3',CBALLS_SHEAR_BIN_THETA='.2'):
            model.Run(level=['MainLoop'])
        with patch.dict(os.environ,CBALLS_SHEAR_PIVOT_TOL='2',CBALLS_SHEAR_BIN_THETA='.8'):
            info=model.getRunMetadata()['shear_hierarchical_reuse']
        if not (info['enabled'] and info['phase_budget_radians']==.3
                and info['radial_bin_slop']==.2 and info['radial_cutoffs']=='strict'):
            raise AssertionError(info)
        model.set(options=BASE)
        model.Run(level=['MainLoop'])
        if model.getRunMetadata()['shear_hierarchical_reuse']['enabled']:
            raise AssertionError('reuse provenance leaked into the following run')
    finally:
        model.struct_cleanup()


def test_phase_bound_encloses_geodesic_formula():
    """The fast chord algebra must never relax the original angular bound."""
    import pathlib, shlex, subprocess, tempfile
    root=pathlib.Path(__file__).resolve().parents[2]
    text=(root/'addons/octree_shear_sphere_2balls_omp/shear_pivot_reuse.h').read_text()
    begin=text.index('static bool shear_reuse_phase_bound(')
    end=text.index('/* Radial slop',begin)
    source=r"""
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#define TRUE true
#define FALSE false
#define MIN(a,b) ((a)<(b)?(a):(b))
#define rsqrt sqrt
typedef double real;
typedef double shear_reuse_ref;
typedef struct { real tolerance; } shear_reuse_context;
typedef struct { real radius, pivot_error; struct { int ring_max; } work; } shear_reuse_level;
static real shear_reuse_error(real q) { return q; }
"""+text[begin:end]+r"""
static uint64_t state=1729;
static double uniform(void) {
    state=state*6364136223846793005ULL+1442695040888963407ULL;
    return (state>>11)*0x1p-53;
}
int main(void) {
    for(int i=0;i<200000;i++) {
        double a=3.141592653589793*uniform(),b=3.141592653589793*uniform();
        double d=3.141592653589793*uniform(),qerror=.02*uniform();
        if(i%2==0) { a*=.003; b*=.003; }
        shear_reuse_context c={3*uniform()};
        shear_reuse_level l={2*sin(a/2),.01*uniform(),{1+i%33}};
        bool accepted=shear_reuse_phase_bound(&c,&l,qerror,2*sin(d/2),2*sin(b/2));
        if(!accepted) continue;
        if(!(d>a+b && d+a+b<3.141592653589793)) return 1;
        double triangle=tan((a+b)/2)*tan((d+a+b)/2);
        if(!(triangle<1)) return 2;
        double bound=l.work.ring_max*(a+b)/fmin(sin(d-a-b),sin(d+a+b))
                     +8*asin(triangle)+qerror;
        if(bound>c.tolerance/3+1e-12 || l.pivot_error>c.tolerance/3) return 3;
    }
    return 0;
}
"""
    with tempfile.TemporaryDirectory(prefix='shear-phase-bound-') as directory:
        path=pathlib.Path(directory)
        (path/'bound.c').write_text(source)
        subprocess.run(shlex.split(os.environ.get('CC','cc'))+
                       ['-O3','-std=c99',str(path/'bound.c'),'-lm','-o',str(path/'bound')],check=True)
        subprocess.run([str(path/'bound')],check=True)


if __name__ == "__main__":
    for test in (test_exact_limit_and_fallbacks, test_aggregate_pivots_and_determinism,
                 test_cross_masks_and_polar_geometry, test_inherited_rings_and_masked_aggregation,
                 test_validation_and_smoothing_fallback,
                 test_hierarchical_radial_completion, test_bin_slop_cutoffs_and_validation,
                 test_reuse_provenance_is_frozen,
                 test_phase_bound_encloses_geodesic_formula):
        test()
        print("PASS", test.__name__, flush=True)
