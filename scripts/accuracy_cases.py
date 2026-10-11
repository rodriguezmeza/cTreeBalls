from acceptance_metadata import ACCURACY_RTOL
"""Bounded regression acceptance cases shared by the runner and documentation.

These settings are not universal error bounds or the fast runtime defaults.
"""
# Max complex Frobenius error, measured separately for 2PCF and 3PCF.
# The octree uses a much smaller opening angle than the binary trees.
SCALAR = [
    dict(engine='kdtree-2balls-omp',theta=.5,rtol=ACCURACY_RTOL,smoothing=True),
    dict(engine='balltree-2balls-omp',theta=.2,rtol=ACCURACY_RTOL,smoothing=True),
    dict(engine='octree-2balls-omp',theta=.025,rtol=ACCURACY_RTOL,smoothing=False),
    dict(engine='octree-sincos-omp',theta=.025,rtol=ACCURACY_RTOL,smoothing=True,estimator_family='core'),
    dict(engine='octree-2balls-omp',theta=.025,rtol=ACCURACY_RTOL,smoothing=True,compatibility=True),
]
SHEAR = [dict(engine=f'{tree}-shear-sphere-2balls-omp',theta=.05,rtol=ACCURACY_RTOL,smoothing=True)
         for tree in ('octree','kdtree','balltree')]

