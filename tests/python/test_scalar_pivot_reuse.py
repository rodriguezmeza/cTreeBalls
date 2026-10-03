"""Hierarchical scalar moments: independent triangles, ownership, OMP and MPI.

Run directly; --mpi uses the already launched communicator. Checks remain
active with Python -O. No pytest dependency is required by MPI workers.
"""

from pathlib import Path
import os
import sys
import tempfile
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if "--mpi" in sys.argv:
    from mpi4py import MPI

    COMM = MPI.COMM_WORLD
else:
    COMM = None
from cyballs import cballs

BASE = "KKKCorrelation,no-smooth-pivot,weights-norm,no-normalize-HistZeta,read-mask,no-out-Hist"


def check(callback):
    error = None
    if COMM is None or COMM.rank == 0:
        try:
            callback()
        except Exception as exc:
            error = repr(exc)
    if COMM:
        error = COMM.bcast(error, root=0)
    if error:
        raise AssertionError(error)


def fixture(n=36, cluster=False):
    rng = np.random.default_rng(1423)
    if cluster:
        centers = rng.normal(size=(12, 3))
        centers /= np.linalg.norm(centers, axis=1)[:, None]
        pos = np.repeat(centers, n // 12, axis=0) + rng.normal(size=(n, 3)) * 1e-7
    else:
        pos = rng.normal(size=(n, 3))
        pos /= np.linalg.norm(pos, axis=1)[:, None]
        pos[:6] = [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]
    field = rng.normal(size=n)
    weights = rng.uniform(0.2, 1.7, n)
    weights[::17] = 0
    mask = np.ones(n, dtype=np.uint8)
    mask[::19] = 0
    return pos, field, weights, mask


def run(
    engine,
    data,
    *,
    reuse=False,
    exact=False,
    threads=1,
    log=False,
    leaf=1,
    extra="",
    neighbor=None,
    minimum=0.02,
    maximum=1.9,
    bins=4,
    periodic=False,
    weighted=True,
    smooth=False
):
    with tempfile.TemporaryDirectory(prefix="scalar-reuse-test-") as out:
        out = COMM.bcast(out if COMM.rank == 0 else None, root=0) if COMM else out
        m = cballs()
        try:
            options = BASE + ("," + extra if extra else "")
            if not weighted:
                options = options.replace("weights-norm,", "")
            if smooth:
                options = options.replace("no-smooth-pivot", "smooth-pivot")
            if reuse:
                options += ",scalar-pivot-reuse"
            if exact:
                options += ",no-two-balls"
            params = dict(
                searchMethod=engine,
                usePeriodic=periodic,
                useLogHist=log,
                rminHist=minimum,
                rangeN=maximum,
                sizeHistN=bins,
                sizeHistPhi=16,
                mChebyshev=2,
                theta=1,
                nsmooth=leaf,
                numberThreads=threads,
                lengthBox=4,
                verbose=0,
                verbose_log=0,
                rootDir=out,
                options=options,
            )
            if neighbor is not None:
                params["iCatalogs"] = "1,2"
            m.set(params)
            for i, (p, k, w, mask) in enumerate(
                [data] if neighbor is None else [data, neighbor]
            ):
                m.set_catalog(p, kappa=k, weights=w, mask=mask, catalog=i)
            m.Run(level=["MainLoop"])
            if COMM and COMM.rank != 0:
                return None
            result = {"metadata": m.getRunMetadata()}
            # Provenance describes the completed run, even if the environment
            # changes before the caller reads it again.
            with patch.dict(
                os.environ, CBALLS_SCALAR_PIVOT_TOL="2.9", CBALLS_SCALAR_BIN_THETA=".9"
            ):
                np.testing.assert_equal(
                    m.getRunMetadata()["scalar_hierarchical_reuse"],
                    result["metadata"]["scalar_hierarchical_reuse"],
                )
            if "only-2pcf" not in options:
                result["components"] = np.array(
                    [
                        [m.getHistZetaMsincos(order, c).copy() for c in range(1, 5)]
                        for order in range(1, 4)
                    ]
                )
                if "edge-corrections" in options:
                    result["corrected"] = np.array(
                        [m.getHistZetaM_EE_complex(o).copy() for o in range(1, 4)]
                    )
                    result["w0"] = m.getScalarWindowDiagnostics()["window_monopole"]
            if "only-3pcf" not in options:
                result["pairs"] = m.getHistNN().copy()
                result["xi"] = m.getHistXi2pcf().copy()
            return result
        finally:
            m.struct_cleanup()


def direct(data, neighbor=None, log=False, minimum=0.02, maximum=1.9, bins=4):
    """Explicit distinct j,k enumeration of all four native basis products."""
    p, k, w, mask = data
    q, l, v, qmask = data if neighbor is None else neighbor
    result = np.zeros((3, 4, bins, bins))
    same = neighbor is None
    edges = (
        np.geomspace(minimum, maximum, bins + 1)
        if log
        else np.linspace(minimum, maximum, bins + 1)
    )
    for i in np.flatnonzero(mask):
        norm = np.linalg.norm(p[i])
        if not norm:
            continue
        normal = p[i] / norm
        ref = np.eye(3)[np.argmin(np.abs(normal))]
        a = ref - normal * np.dot(ref, normal)
        a /= np.linalg.norm(a)
        b = np.cross(a, normal)
        legs = q - p[i]
        distance = np.linalg.norm(legs, axis=1)
        x = legs @ a
        y = legs @ b
        transverse = np.hypot(x, y)
        good = (
            qmask.astype(bool)
            & (distance > minimum)
            & (distance < maximum)
            & (transverse > 32 * np.finfo(float).eps * distance)
        )
        if same:
            good[i] = False
        indices = np.flatnonzero(good)
        radial = np.searchsorted(edges, distance, side="right") - 1
        angle = np.arctan2(y, x)
        for j in indices:
            for h in indices:
                if j == h:
                    continue
                factor = w[i] * v[j] * v[h] * k[i] * l[j] * l[h]
                for order in range(3):
                    c1, s1 = np.cos(order * angle[j]), np.sin(order * angle[j])
                    c2, s2 = np.cos(order * angle[h]), np.sin(order * angle[h])
                    result[order, :, radial[j], radial[h]] += factor * np.array(
                        [c1 * c2, s1 * s2, s1 * c2, c1 * s2]
                    )
    return result


def compare(a, b, tolerance=2e-10):
    for key in a:
        if key != "metadata":
            np.testing.assert_allclose(
                b[key],
                a[key],
                rtol=tolerance,
                atol=tolerance,
                equal_nan=True,
                err_msg=key,
            )


def suite(engine):
    data = fixture()
    neighbor = fixture(48)
    with patch.dict(
        os.environ, CBALLS_SCALAR_PIVOT_TOL="1e-12", CBALLS_SCALAR_BIN_THETA="0"
    ):
        for log in (False, True):
            for cross in (False, True):
                other = neighbor if cross else None
                exact = run(engine, data, exact=True, log=log, neighbor=other)
                actual = run(
                    engine, data, reuse=True, threads=3, log=log, neighbor=other
                )
                check(lambda: compare(exact, actual))
                expected = direct(data, other, log)
                check(
                    lambda: np.testing.assert_allclose(
                        actual["components"], expected, rtol=2e-10, atol=2e-10
                    )
                )
    # Undefined directions, axis-switch caps, unweighted catalogs, and
    # a wide dynamic range in observer distance still converge to body sums.
    awkward = list(fixture(48))
    awkward[0][:8] = [
        [0, 0, 0],
        [0, 0, 0.4],
        [0, 0, -0.4],
        [1, 1, 2],
        [1 + 1e-8, 1 - 1e-8, 2],
        [1 - 1e-8, 1 + 1e-8, 2],
        [-1, -1, -2],
        [1e-9, 0, 1e-9],
    ]
    with patch.dict(
        os.environ, CBALLS_SCALAR_PIVOT_TOL="1e-12", CBALLS_SCALAR_BIN_THETA="0"
    ):
        exact = run(engine, awkward, exact=True, weighted=False, maximum=5)
        actual = run(engine, awkward, reuse=True, weighted=False, maximum=5, leaf=16)
        check(lambda: compare(exact, actual))
        for minimum in (0.02,):
            exact = run(
                engine,
                data,
                exact=True,
                log=True,
                minimum=minimum,
                extra="only-3pcf,edge-corrections",
            )
            actual = run(
                engine,
                data,
                reuse=True,
                log=True,
                minimum=minimum,
                extra="only-3pcf,edge-corrections",
            )
            check(lambda: compare(exact, actual))
    for periodic, smooth in ((True, False), (False, True)):
        if smooth and engine.startswith("octree"):
            continue
        a = run(engine, data, periodic=periodic, smooth=smooth)
        b = run(engine, data, reuse=True, periodic=periodic, smooth=smooth)
        check(lambda: compare(a, b, 0))
        check(
            lambda: np.testing.assert_equal(
                b["metadata"]["scalar_hierarchical_reuse"]["enabled"], False
            )
        )
    clustered = fixture(768, cluster=True)
    with patch.dict(
        os.environ, CBALLS_SCALAR_PIVOT_TOL=".01", CBALLS_SCALAR_BIN_THETA="0"
    ):
        for log in (False, True):
            exact = run(engine, clustered, exact=True, log=log, extra="only-3pcf")
            one = run(engine, clustered, reuse=True, log=log, extra="only-3pcf")
            many = run(
                engine, clustered, reuse=True, log=log, extra="only-3pcf", threads=3
            )

            def aggregation():
                compare(one, many, 0)
                info = one["metadata"]["scalar_hierarchical_reuse"]
                if not (
                    info["enabled"]
                    and info["parent_reductions"] > 0
                    and info["represented_pairs"] > info["radial_pairs"]
                ):
                    raise AssertionError(info)
                if (
                    np.linalg.norm(one["components"] - exact["components"])
                    / np.linalg.norm(exact["components"])
                    > 2e-4
                ):
                    raise AssertionError(
                        (
                            "clustered raw error",
                            np.linalg.norm(one["components"] - exact["components"])
                            / np.linalg.norm(exact["components"]),
                            info,
                        )
                    )

            check(aggregation)
        # Corrected estimates use the same scalar window operator, including
        # unstable/empty windows. Compare only a fixture with well behaved bins.
        exact = run(engine, data, exact=True, extra="edge-corrections")
        zero = run(engine, data, reuse=True, exact=True, extra="edge-corrections")
        check(lambda: compare(exact, zero))
        combined = run(engine, clustered, reuse=True)
        triple = run(engine, clustered, reuse=True, extra="only-3pcf")
        pairs = run(engine, clustered, reuse=True, extra="only-2pcf")
        baseline = run(engine, clustered, extra="only-2pcf")
        check(
            lambda: np.testing.assert_array_equal(
                combined["components"], triple["components"]
            )
        )
        check(lambda: compare(baseline, pairs, 0))
        check(
            lambda: np.testing.assert_array_equal(combined["pairs"], baseline["pairs"])
        )
        # Bin migration may occur internally, never across the radial cutoffs.
        constant = tuple(
            np.ones_like(x) if i == 1 else x for i, x in enumerate(clustered)
        )
        for log in (False, True):
            exact = run(engine, constant, exact=True, log=log, extra="only-3pcf")
            with patch.dict(os.environ, CBALLS_SCALAR_BIN_THETA=".5"):
                loose = run(
                    engine, constant, reuse=True, log=log, extra="only-3pcf", leaf=16
                )
            check(
                lambda: np.testing.assert_allclose(
                    loose["components"][0, 0].sum(),
                    exact["components"][0, 0].sum(),
                    rtol=2e-13,
                )
            )
        for extra in ("only-2pcf", "no-two-balls", "no-one-ball"):
            a = run(engine, data, extra=extra)
            b = run(engine, data, reuse=True, extra=extra)
            check(lambda: compare(a, b, 0))
            check(
                lambda: np.testing.assert_equal(
                    b["metadata"]["scalar_hierarchical_reuse"]["enabled"], False
                )
            )
    with patch.dict(
        os.environ, CBALLS_SCALAR_PIVOT_TOL="0", CBALLS_SCALAR_BIN_THETA="0"
    ):
        a = run(engine, data)
        b = run(engine, data, reuse=True)
        check(lambda: compare(a, b, 0))
    for name in ("CBALLS_SCALAR_PIVOT_TOL", "CBALLS_SCALAR_BIN_THETA"):
        for value in ("nan", "inf", "-1", "4", "", ".1junk"):
            try:
                with patch.dict(os.environ, {name: value}):
                    run(engine, data, reuse=True)
            except Exception as error:
                if name not in str(error) and "scalar pivot-reuse controls" not in str(
                    error
                ):
                    raise
            else:
                raise AssertionError("invalid control accepted")
    if COMM and COMM.size > 1:
        for value in ("nan", ".8"):
            try:
                with patch.dict(
                    os.environ,
                    CBALLS_SCALAR_BIN_THETA=value if COMM.rank == 1 else ".2",
                ):
                    run(engine, data, reuse=True)
            except Exception:
                pass
            else:
                raise AssertionError("rank-local invalid/mismatched control accepted")
    if COMM is None or COMM.rank == 0:
        print(
            "PASS",
            engine,
            "triangles, hierarchy, threads, cutoffs, fallbacks, controls",
            flush=True,
        )


def test_native_phase_enclosure():
    """Sample descendants, including least-axis switches and near axial legs."""
    import subprocess, shlex

    angular = (ROOT / "include/angular_contracts.h").read_text()
    basis = angular[
        angular.index("static inline bool cballs_angular_basis") : angular.index(
            "\n#endif", angular.index("static inline bool cballs_angular_basis")
        )
    ]
    phase = angular[
        angular.index("static inline bool cballs_angular_phase") : angular.index(
            "/* Conservative variation"
        )
    ]
    reuse = (ROOT / "addons/balltree_2balls_omp/dual_node_pivot_reuse.h").read_text()
    bound = reuse[
        reuse.index("static real dual_node_reuse_basis_error") : reuse.index(
            "static void dual_node_reuse_mark"
        )
    ]
    code = (
        r"""
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#define NDIM 3
#define FALSE false
#define TRUE true
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define CROSSVP(c,a,b) do { (c)[0]=(a)[1]*(b)[2]-(a)[2]*(b)[1]; (c)[1]=(a)[2]*(b)[0]-(a)[0]*(b)[2]; (c)[2]=(a)[0]*(b)[1]-(a)[1]*(b)[0]; } while(0)
typedef double real;
typedef double cballs_storage_real;
typedef double compute_vector[3];
typedef struct { bool angular_basis_valid; real position_norm,radius,normal[3]; } dual_node_multipole_pivot;
"""
        + basis
        + phase
        + bound
        + r"""
static uint64_t state=1273889;
static double u(void) { state^=state<<13;state^=state>>7;state^=state<<17;return (state>>11)*0x1p-53; }
static void ball(double *v,double r) { double norm;do {for(int k=0;k<3;k++)v[k]=2*u()-1;norm=hypot(hypot(v[0],v[1]),v[2]);}while(norm>1||norm==0);for(int k=0;k<3;k++)v[k]*=r; }
int main(void) {
 int accepted=0;
 for(int trial=0;trial<200000;trial++) {
  double p[3],q[3],dr[3],a[3],b[3],dp[3],dq[3],child[3],leg[3],c,s,x=0,y=0;
  ball(p,2);ball(dr,2);for(int k=0;k<3;k++)q[k]=p[k]-dr[k];
  dual_node_multipole_pivot pivot;
  pivot.radius=pow(10.,-8+6*u());double qr=pow(10.,-8+6*u());
  pivot.position_norm=hypot(hypot(p[0],p[1]),p[2]);
  pivot.angular_basis_valid=cballs_angular_basis(p,pivot.normal,a,b);
  double e=dual_node_reuse_basis_error(&pivot);
  if(!isfinite(e)||!cballs_angular_phase(p,dr,&c,&s))continue;
  for(int k=0;k<3;k++){x-=dr[k]*a[k];y-=dr[k]*b[k];}
  double ratio=(pivot.radius+qr+hypot(hypot(dr[0],dr[1]),dr[2])*e)/hypot(x,y);
  if(!(ratio<.5))continue;
  ball(dp,pivot.radius);ball(dq,qr);
  for(int k=0;k<3;k++){child[k]=p[k]+dp[k];leg[k]=child[k]-q[k]-dq[k];}
  double cc,ss;if(!cballs_angular_phase(child,leg,&cc,&ss))return 2;
  if(fabs(atan2(ss*c-cc*s,cc*c+ss*s))>asin(ratio)+2e-12)return 3;
  accepted++;
 }
 if(accepted<10000)return 4;
 printf("PASS phase enclosure: %d descendant samples\n",accepted);return 0;
}
"""
    )
    with tempfile.TemporaryDirectory(prefix="scalar-reuse-phase-") as temp:
        path = Path(temp)
        (path / "phase.c").write_text(code)
        subprocess.run(
            shlex.split(os.environ.get("CC", "cc"))
            + [
                "-O2",
                "-std=c99",
                str(path / "phase.c"),
                "-lm",
                "-o",
                str(path / "phase"),
            ],
            check=True,
        )
        subprocess.run([str(path / "phase")], check=True)


if __name__ == "__main__":
    if COMM is None:
        test_native_phase_enclosure()
    for engine in ("octree", "kdtree", "balltree"):
        suite(engine + "-2balls-" + ("mpi" if COMM else "omp"))
