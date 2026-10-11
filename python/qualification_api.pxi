def _qualification_scope(packet):
    metadata = packet['metadata']
    inputs = metadata.get('inputs', {})
    catalogs = inputs.get('effective_catalogs')
    if not isinstance(catalogs, dict) or catalogs.get('schema_version') != 1 or not catalogs.get('catalogs'):
        raise ValueError('qualification requires effective catalog fingerprints; rerun legacy result packets with the current build')
    ignored = {'smooth-pivot','no-smooth-pivot','no-one-ball','no-two-balls','behavior-ball',
               'dual-node-bin-theta','shear-pivot-reuse','scalar-pivot-reuse','dual-node-profile','ggg-profile','no-out-Hist',
               'no-balltree-tree-cache','no-native-tree-cache','no-kdtree-tree-cache',
               'no-balltree-parallel-build','no-balltree-persistent-frontier'}
    options = sorted(set(metadata.get('options','').split(','))-ignored-{''})
    return dict(catalogs=catalogs, estimator=metadata['estimator'], geometry=metadata['geometry'],
                bin_edges=metadata['bin_edges'], weights=metadata['weights'], masks=metadata['masks'],
                box=metadata['box'], multipole_max=metadata.get('multipole_max'), coordinate_convention=metadata.get('coordinate_convention'),
                options=options, catalog_selection=inputs.get('catalog_selection'),
                physical_survey=metadata.get('physical_survey'))


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


def qualification_report(candidate, reference, *, rtol=.02, atol=1e-10, selection=None, metric="l2"):
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
    if metric not in ('l2','pointwise'):
        raise ValueError('metric must be l2 or pointwise')
    selection=dict(selection or {})
    candidate_scope = _qualification_scope(candidate)
    reference_scope = _qualification_scope(reference)
    if candidate_scope != reference_scope:
        raise ValueError('qualification requires matching catalog content, estimator, geometry, bins, masks and weights')
    if not _qualification_exact(reference['metadata']):
        raise ValueError('reference must explicitly disable aggregation and smoothing; approximate references cannot qualify results')
    arrays = dict(candidate['arrays']); targets = dict(reference['arrays'])
    if set(arrays) != set(targets) or not arrays:
        raise ValueError('qualification requires identical nonempty observable sets')
    if set(selection)-set(arrays):
        raise ValueError('selection names must be returned observables')
    survey = candidate['metadata'].get('physical_survey')
    if ('data_triple_numerator' in arrays or 'data_pair_numerator' in arrays) and not survey:
        raise ValueError('legacy physical survey packets require a rerun with corrected products and survey metadata')
    if survey and set(selection)&{'zeta_rcond','zeta_relative_residual'}:
        raise ValueError('selection applies to scientific observables, not support diagnostics')
    if survey:
        for raw, corrected in [('data_triple_numerator','zeta_corrected_legendre'),('data_pair_numerator','xi_corrected')]:
            if raw in arrays and corrected not in arrays:
                raise ValueError('physical survey qualification requires corrected products: '+corrected)
        required=set()
        if 'data_triple_numerator' in arrays:
            required.update(('zeta_corrected_encore','zeta_valid','zeta_status','zeta_rcond','zeta_relative_residual'))
        if 'data_pair_numerator' in arrays:required.update(('xi_valid','xi_status'))
        if not required<=set(arrays):raise ValueError('physical survey qualification requires validity and diagnostics')
    # Report the normalized forest observables as well as their raw sums.
    # Empty bins remain undefined here, so validity changes are visible.
    if candidate['metadata'].get('engine','').startswith('lya-'):
        for prefix in ('pair','triple'):
            numerator, denominator = prefix+'_numerator', prefix+'_denominator'
            if numerator in arrays and denominator in arrays:
                for products in (arrays, targets):
                    num=np.asarray(products[numerator]); den=np.asarray(products[denominator])
                    products[prefix+'_normalized']=np.divide(num,den,
                        out=np.full_like(num,np.nan,dtype=np.float64),where=den>0)
                if numerator in selection:selection[prefix+'_normalized']=selection[numerator]
    # Resolve all four shear components separately; a strong component must
    # not hide an inaccurate weak component in the pooled array norm.
    for name in tuple(arrays):
        if name in ('upsilon_raw','gamma_corrected','gamma_normalized'):
            for i in range(4):
                arrays[name+'_component_'+str(i)] = np.asarray(arrays[name])[i]
                targets[name+'_component_'+str(i)] = np.asarray(targets[name])[i]
                if name in selection:selection[name+'_component_'+str(i)]=np.asarray(selection[name])[i]
    records={}; selected_scope={}; diagnostics={}
    for name, value in arrays.items():
        value=np.asarray(value); target=np.asarray(targets[name])
        if value.shape != target.shape:
            raise ValueError('observable shape mismatch: '+name)
        if name in ('zeta_rcond','zeta_relative_residual') and survey:
            support=np.asarray(targets['zeta_valid'],dtype=bool)
            acceptable=bool(np.array_equal(support,np.asarray(arrays['zeta_valid'],dtype=bool)))
            for packet, diagnostic in ((candidate,value),(reference,target)):
                good=diagnostic[support]
                acceptable=acceptable and bool(np.isfinite(good).all())
                if name=='zeta_rcond':acceptable=acceptable and bool((good>=survey['rcond_min']).all())
                else:acceptable=acceptable and bool(((good>=0)&(good<=128*(survey['solve_lmax']+1)*np.finfo(float).eps)).all())
            diagnostics[name]=dict(status='QUALIFIED' if acceptable else 'REJECTED',
                role='support diagnostic, excluded from scientific error norm',supported_bins=int(support.sum()))
            continue
        mask=np.asarray(selection.get(name,np.ones(value.shape,dtype=bool)))
        if mask.dtype.kind!='b' or mask.shape!=value.shape or not mask.any():
            raise ValueError('selection must be a nonempty Boolean mask with the observable shape: '+name)
        selected_scope[name]=dict(shape=list(mask.shape),selected_bins=int(mask.sum()),
            sha256=hashlib.sha256(mask.tobytes()).hexdigest())
        value=value[mask];target=target[mask]
        valid=np.isfinite(target); actual_valid=np.isfinite(value)
        mask_match=bool(np.array_equal(valid,actual_valid))
        has_inf=bool(np.isinf(value).any() or np.isinf(target).any())
        # Subtract integer counts before conversion: unsigned subtraction wraps,
        # and converting large (>2**53) counters first loses unit differences.
        integer = value.dtype.kind in 'iu' and target.dtype.kind in 'iu'
        discrete=integer or name in ('zeta_valid','zeta_status','xi_valid','xi_status')
        delta=(value[valid & actual_valid].astype(object)-target[valid & actual_valid].astype(object)) if integer else value[valid & actual_valid]-target[valid & actual_valid]
        if integer: delta=np.asarray(delta,dtype=np.float64)
        error=float(np.linalg.norm(delta)); norm=float(np.linalg.norm(target[valid]))
        maximum=float(np.max(np.abs(delta))) if delta.size else None
        weak=np.abs(target[valid]) <= atol
        per_bin=np.abs(delta)/(atol+rtol*np.abs(target[valid & actual_valid])+np.finfo(float).tiny)
        passed=bool(valid.any() and mask_match and not has_inf and np.isfinite(norm) and np.isfinite(error)
                    and (not discrete or np.array_equal(value,target))
                    and (error <= atol+rtol*norm if metric=='l2' else bool((per_bin<=1).all())))
        records[name]=dict(status='QUALIFIED' if passed else 'REJECTED', shape=list(np.asarray(arrays[name]).shape), selected_bins=int(mask.sum()),
            finite_bins=int(valid.sum()), finite_mask_matches=mask_match, contains_infinity=has_inf,
            integer_counts_require_exact_equality=integer,
            absolute_l2=error if np.isfinite(error) else None,
            reference_l2=norm if np.isfinite(norm) else None,
            relative_l2=(error/norm if norm>0 and np.isfinite(norm) and np.isfinite(error)
                         and mask_match and not has_inf else None),
            relative_l2_percent=(100.0*error/norm if norm>0 and np.isfinite(norm)
                         and np.isfinite(error) and mask_match and not has_inf else None),
            relative_target_met=bool(mask_match and not has_inf and valid.any()
                         and np.isfinite(norm) and np.isfinite(error) and error<=rtol*norm),
            zero_reference_norm=bool(norm==0),
            maximum_absolute_bin_error=maximum if maximum is None or np.isfinite(maximum) else None, weak_reference_bins=int(weak.sum()),
            bins_exceeding_pointwise_tolerance=int(np.count_nonzero(per_bin>1)),
            maximum_pointwise_tolerance_ratio=float(per_bin.max()) if per_bin.size and np.isfinite(per_bin.max()) else None)
    evidence_scope=dict(candidate_scope,selection=selected_scope,metric=metric)
    fingerprint=hashlib.sha256(json.dumps(evidence_scope,sort_keys=True,allow_nan=False).encode()).hexdigest()
    return dict(schema_version=1,status='QUALIFIED' if all(v['status']=='QUALIFIED' for v in list(records.values())+list(diagnostics.values())) else 'REJECTED',
        criterion=('L2(error) <= atol + rtol*L2(reference)' if metric=='l2' else 'abs(error) <= atol + rtol*abs(reference), every selected bin')+'; independently per observable, identical finite masks, exact counts and validity codes',
        metric=metric,selection=selected_scope,diagnostics=diagnostics,
        tolerances=dict(rtol=float(rtol),atol=float(atol)), observables=records,
        scope_sha256=fingerprint, scope=evidence_scope,
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


def publish_result_packet(model, directory, reference=None, *, rtol=.02, atol=1e-10, fail_on_rejection=True):
    """Publish measured penalties; optionally retain results beyond the target."""
    from pathlib import Path
    reference_packet = load_result_packet(reference) if reference is not None else None
    if reference_packet is not None:
        model.qualifyAgainst(reference_packet,rtol=rtol,atol=atol)
    packet=model.getResults()
    path=save_result_packet(packet,directory)
    report=packet['metadata']['qualification']
    penalty=dict(schema_version=1,status='MEASURED' if reference_packet is not None else 'UNMEASURED',
        target_relative_l2=float(rtol), target_percent=100.0*rtol,
        criterion='norm(candidate-reference)/norm(reference), separately per observable; no denominator floor',
        observables=report.get('observables',{}), qualification_status=report['status'],
        scope='Only this catalogue, estimator, bins and multipoles. Subsample measurements do not qualify the full catalogue.')
    metadata=packet['metadata']
    penalty['controls']={key:metadata.get(key) for key in ('engine','opening_tolerance',
        'effective_smoothing','shear_hierarchical_reuse','scalar_hierarchical_reuse',
        'lya_pivot_frontier','lya_2pcf','lya_3pcf') if key in metadata}
    if reference_packet is None:
        print('  approximation penalty: UNMEASURED (supply --accuracy-reference for an exact comparison)',flush=True)
    else:
        print(f"  approximation penalty: MEASURED; relative L2 target={100*rtol:g}%",flush=True)
        for name,row in penalty['observables'].items():
            percent=row['relative_l2_percent']
            value=f'{percent:.6g}%' if percent is not None else 'undefined (zero reference or invalid bins)'
            state='within target' if row['relative_target_met'] else 'EXCEEDS TARGET or invalid reference'
            print(f"    {name}: L2={value}; {state}; absolute L2={row['absolute_l2']:.6g}" if row['absolute_l2'] is not None
                  else f"    {name}: L2={value}; {state}",flush=True)
        # MPI rank-zero timing is not a job timing. Existing driver reports
        # perform the required max-wall/sum-CPU reduction separately.
        refmeta=reference_packet['metadata']
        par=metadata.get('parallel',{}); refpar=refmeta.get('parallel',{})
        if (par.get('estimator_ranks')==refpar.get('estimator_ranks')==1
                and par.get('threads_requested')==refpar.get('threads_requested')
                and metadata.get('engine')==refmeta.get('engine')):
            timing=metadata.get('timings',{}); reftime=refmeta.get('timings',{})
            if timing.get('wall_seconds',0)>0 and reftime.get('wall_seconds',0)>0:
                speedup=reftime['wall_seconds']/timing['wall_seconds']
                penalty['native_wall_speedup']=speedup
                penalty['timing_scope']='single-rank MainLoop only; same engine and thread count; cache/host conditions may differ'
                print(f'  measured MainLoop wall speedup: {speedup:.3g}x (single rank; excludes validation and I/O)',flush=True)
    (Path(directory)/'approximation-penalty.json').write_text(json.dumps(penalty,indent=2,allow_nan=False)+'\n')
    print('  numerical qualification: '+report['status'],flush=True)
    if report['status']=='REJECTED' and fail_on_rejection:
        raise ValueError('requested numerical qualification failed; retained evidence: '+path)
    return path
