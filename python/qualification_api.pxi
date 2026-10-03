def _qualification_scope(packet):
    metadata = packet['metadata']
    inputs = metadata.get('inputs', {})
    catalogs = inputs.get('effective_catalogs')
    if not isinstance(catalogs, dict) or catalogs.get('schema_version') != 1 or not catalogs.get('catalogs'):
        raise ValueError('qualification requires effective catalog fingerprints; rerun legacy result packets with the current build')
    ignored = {'smooth-pivot','no-smooth-pivot','no-one-ball','no-two-balls','behavior-ball',
               'dual-node-bin-theta','dual-node-profile','ggg-profile','no-out-Hist',
               'no-balltree-tree-cache','no-native-tree-cache','no-kdtree-tree-cache',
               'no-balltree-parallel-build','no-balltree-persistent-frontier'}
    options = sorted(set(metadata.get('options','').split(','))-ignored-{''})
    return dict(catalogs=catalogs, estimator=metadata['estimator'], geometry=metadata['geometry'],
                bin_edges=metadata['bin_edges'], weights=metadata['weights'], masks=metadata['masks'],
                box=metadata['box'], multipole_max=metadata.get('multipole_max'), coordinate_convention=metadata.get('coordinate_convention'),
                options=options, catalog_selection=inputs.get('catalog_selection'))


def _qualification_exact(metadata):
    if metadata.get('effective_smoothing',{}).get('enabled',True):
        return False
    method = metadata.get('engine','')
    # These active physical Legendre implementations visit original bodies;
    # theta only rescales a conservative pruning radius, never aggregates cells.
    if method in ('octree-3pcf-3d-omp','octree-3pcf-3d-mpi'):
        return bool(metadata.get('bin_edges'))
    if method.startswith('lya-'):
        if 'multipole' in method:
            return False
        for key in ('lya_pivot_frontier','lya_2pcf','lya_3pcf','lya_geometry'):
            item = metadata.get(key,{})
            if any(v is True for k,v in item.items() if k.endswith('approximate')):
                return False
        return bool(metadata.get('bin_edges'))
    opening=metadata.get('opening_tolerance',{})
    return opening.get('theta') == 0 and (opening.get('no_one_ball') or opening.get('no_two_balls'))


def qualification_report(candidate, reference, *, rtol=.02, atol=1e-10):
    """Compare actual copied observables on the identical catalog and geometry.

    Acceptance uses absolute-plus-relative L2 error, separately per array. Bin
    errors and weak-reference bins are diagnostic and are never hidden by a
    single combined statistic. Catalog identity uses native interpreted fields
    before tree/smoothing mutation, including integer IDs, masks and weights.
    Equivalent loaders can compare when their canonical rows and estimator
    controls match; file bytes and filenames remain provenance only. Legacy
    packets without this identity must be rerun. A result says nothing about another catalog,
    observable, mask, binning, opening angle or smoothing radius.
    """
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol < 0 or atol < 0:
        raise ValueError('qualification tolerances must be finite and nonnegative')
    candidate_scope = _qualification_scope(candidate)
    reference_scope = _qualification_scope(reference)
    if candidate_scope != reference_scope:
        raise ValueError('qualification requires matching catalog content, estimator, geometry, bins, masks and weights')
    if not _qualification_exact(reference['metadata']):
        raise ValueError('reference must explicitly disable aggregation and smoothing; approximate references cannot qualify results')
    arrays = candidate['arrays']; targets = reference['arrays']
    if set(arrays) != set(targets) or not arrays:
        raise ValueError('qualification requires identical nonempty observable sets')
    records={}
    for name, value in arrays.items():
        value=np.asarray(value); target=np.asarray(targets[name])
        if value.shape != target.shape:
            raise ValueError('observable shape mismatch: '+name)
        valid=np.isfinite(target); actual_valid=np.isfinite(value)
        mask_match=bool(np.array_equal(valid,actual_valid))
        has_inf=bool(np.isinf(value).any() or np.isinf(target).any())
        # Subtract integer counts before conversion: unsigned subtraction wraps,
        # and converting large (>2**53) counters first loses unit differences.
        integer = value.dtype.kind in 'iu' and target.dtype.kind in 'iu'
        delta=(value[valid & actual_valid].astype(object)-target[valid & actual_valid].astype(object)) if integer else value[valid & actual_valid]-target[valid & actual_valid]
        if integer: delta=np.asarray(delta,dtype=np.float64)
        error=float(np.linalg.norm(delta)); norm=float(np.linalg.norm(target[valid]))
        maximum=float(np.max(np.abs(delta))) if delta.size else None
        weak=np.abs(target[valid]) <= atol
        per_bin=np.abs(delta)/(atol+rtol*np.abs(target[valid & actual_valid])+np.finfo(float).tiny)
        passed=bool(valid.any() and mask_match and not has_inf and np.isfinite(error)
                    and (not integer or np.array_equal(value,target))
                    and error <= atol+rtol*norm)
        records[name]=dict(status='QUALIFIED' if passed else 'REJECTED', shape=list(value.shape),
            finite_bins=int(valid.sum()), finite_mask_matches=mask_match, contains_infinity=has_inf,
            integer_counts_require_exact_equality=integer,
            absolute_l2=error if np.isfinite(error) else None,
            reference_l2=norm if np.isfinite(norm) else None,
            relative_l2=error/max(norm,atol,np.finfo(float).tiny) if np.isfinite(error) else None,
            maximum_absolute_bin_error=maximum if maximum is None or np.isfinite(maximum) else None, weak_reference_bins=int(weak.sum()),
            bins_exceeding_pointwise_tolerance=int(np.count_nonzero(per_bin>1)),
            maximum_pointwise_tolerance_ratio=float(per_bin.max()) if per_bin.size and np.isfinite(per_bin.max()) else None)
    fingerprint=hashlib.sha256(json.dumps(candidate_scope,sort_keys=True,allow_nan=False).encode()).hexdigest()
    return dict(schema_version=1,status='QUALIFIED' if all(v['status']=='QUALIFIED' for v in records.values()) else 'REJECTED',
        criterion='L2(error) <= atol + rtol*L2(reference), independently per observable; identical finite masks; no infinities; integer counts must match exactly',
        tolerances=dict(rtol=float(rtol),atol=float(atol)), observables=records,
        scope_sha256=fingerprint, scope=candidate_scope,
        candidate=dict(engine=candidate['metadata']['engine'],build=candidate['metadata']['build']['id'],
            opening=candidate['metadata']['opening_tolerance'],smoothing=candidate['metadata']['effective_smoothing']),
        reference=dict(engine=reference['metadata']['engine'],build=reference['metadata']['build']['id'],exact_controls=True),
        limitation='Only these arrays on this catalog and geometry. No qualification is transferred to other data or observables.')

def save_result_packet(packet, directory):
    """Persist the exact arrays being qualified alongside their provenance."""
    from pathlib import Path
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(directory/'observables.npz',**packet['arrays'])
    document=dict(metadata=packet['metadata'],arrays_sha256=hashlib.sha256((directory/'observables.npz').read_bytes()).hexdigest())
    (directory/'qualification-result.json').write_text(json.dumps(document,indent=2,allow_nan=False)+'\n')
    return str(directory/'qualification-result.json')


def load_result_packet(directory):
    """Load and verify a saved result packet; pickle/object arrays are disallowed."""
    from pathlib import Path
    directory=Path(directory)
    document=json.loads((directory/'qualification-result.json').read_text())
    raw=(directory/'observables.npz').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=document['arrays_sha256']:
        raise ValueError('result packet array digest mismatch')
    with np.load(directory/'observables.npz',allow_pickle=False) as data:
        arrays={name:data[name].copy() for name in data.files}
    return dict(metadata=document['metadata'],arrays=arrays)


def publish_result_packet(model, directory, reference=None, *, rtol=.02, atol=1e-10):
    """Driver publication with optional exact evidence and a failing rejection."""
    if reference is not None:
        model.qualifyAgainst(load_result_packet(reference),rtol=rtol,atol=atol)
    packet=model.getResults()
    path=save_result_packet(packet,directory)
    report=packet['metadata']['qualification']
    print('  numerical qualification: '+report['status'],flush=True)
    for name,row in report.get('observables',{}).items():
        print(f"    {name}: {row['status']}; relative L2={row['relative_l2']}; weak bins={row['weak_reference_bins']}; pointwise exceedances={row['bins_exceeding_pointwise_tolerance']}",flush=True)
    if report['status']=='REJECTED':
        raise ValueError('requested numerical qualification failed; retained evidence: '+path)
    return path
