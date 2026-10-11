"""Shared release obligations and fixed scientific acceptance policies.

A completed software campaign and scientific qualification are separate states.
Changing a policy invalidates its identity; thresholds never depend on observations.
"""
import hashlib
import json
from pathlib import Path
from build_fingerprint import source_files
ROOT = Path(__file__).resolve().parents[1]

REQUIRED = {
    'active': {'active-gate', 'parameter-file-parser', 'configuration-readers', 'incremental-build', 'accuracy-documentation', 'survey-science-campaign'},
    'package': {'sdist-build', 'sdist-metadata', 'sdist-policy', 'clean-install'},
    'sanitizers': {'sanitizer-profile', 'sanitizer-probe-build', 'sanitizer-probe-run',
                   'sanitizer-build', 'parameter-parser', 'parameter-file-parser',
                   'standalone-parser', 'sanitizer-smoke', 'physical-oracle', 'exact-numerical-oracle', 'physical-survey-oracle'},
}


def source_identity(root=ROOT):
    paths = set(source_files(root))
    paths.update(root.glob('*.md'))
    for directory, suffixes in [('docs', {'.md', '.rst', '.py', '.m', '.1', '.html'}),
                                 ('.github', {'.yml', '.yaml'}), ('requirements', {'.txt'}),
                                 ('examples', {'.py', '.ipynb'})]:
        paths.update(p for p in (root/directory).rglob('*') if p.suffix in suffixes)
    digest = hashlib.sha256()
    for path in sorted(paths):
        if path.is_file() and not path.is_symlink():
            digest.update(path.relative_to(root).as_posix().encode() + b'\0')
            digest.update(hashlib.sha256(path.read_bytes()).digest())
    return digest.hexdigest()


ACCURACY_RTOL = .02
ACCURACY_ATOL = 1e-10
EXACT_RTOL = 1e-11
SURVEY_RTOL = .05
SURVEY_ATOL = 1e-14
SURVEY_EXACT_RTOL = 3e-10

def acceptance_policy():
    # Return a fresh record; callers cannot mutate future policy records.
    from accuracy_cases import SCALAR, SHEAR
    record = dict(schema_version=1,
        release_checks={k: sorted(v) for k,v in REQUIRED.items()},
        bounded_accuracy=dict(rtol=ACCURACY_RTOL, atol=ACCURACY_ATOL,
            exact_rtol=EXACT_RTOL, metric='per-product complex Frobenius L2',
            cases=[dict(c) for c in SCALAR+SHEAR]),
        survey_defaults=dict(rtol=SURVEY_RTOL, atol=SURVEY_ATOL,
            exact_rtol=SURVEY_EXACT_RTOL, metric='pointwise',
            selected_radial_bins=[1,2], fixtures=['signed','clustered','masked','weak_signal']))
    record['id'] = hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':'),
                                             allow_nan=False).encode()).hexdigest()
    return record
