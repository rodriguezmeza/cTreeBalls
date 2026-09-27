# Public resource planning independent of catalog allocation. Native layout
# sizes, not hard-coded Python guesses, define all item widths.
def resource_policy():
    """Current environment budget and compiled layout; this is not an RSS cap."""
    budget=cballs_memory_budget()
    if not budget: raise ValueError("CBALLS_MEMORY_BUDGET_MB must be a positive integer")
    return {"default_budget_mib": cballs_resource_default_mib(),
            "budget_bytes": budget, "scope": "per-rank checked plans and allocations; not total RSS",
            "type_bytes": {name: cballs_resource_type_size(i) for i,name in enumerate(
                ("body","octree_cell","binary_node","packed_point","real","size_t","long_double"))}}


def resource_plan(engine, pixels, threads=1, parameters=None, retained_cache_bytes=0):
    """Forecast dominant per-rank storage before loading a catalog.

    Known components and backend/tree estimates are separate. The forecast
    excludes the interpreter, file readers, user arrays, MPI, allocator slack,
    and data-dependent frontier/scratch/correction storage. It is not a promise
    that the process will fit RAM. Python integers are checked before returning.
    """
    if search_method_id(engine)<0: raise ValueError("engine is not in this build")
    policy=resource_policy(); sizes=policy["type_bytes"];p=dict(parameters or {})
    def positive(value,name):
        if isinstance(value,bool) or int(value)!=value or int(value)<1: raise ValueError(name+" must be a positive integer")
        return int(value)
    n=positive(pixels,"pixels");t=positive(threads,"threads")
    b=positive(p.get("sizeHistN",20),"sizeHistN");m=positive(p.get("mChebyshev",7)+1,"multipole orders")
    options=set(str(p.get("options","")).split(','));scalar=not any(x in engine for x in ('lya-','shear','box','3pcf-3d','ggg-3d'))
    scalar3=scalar and 'only-2pcf' not in options
    # Exact common NR shapes for the public PXD + smoothing profile.
    settings=build_info()["resolved_settings"];flags=settings.get('CCFLAG','')
    smooth='-DSMOOTHPIVOT' in flags;pxd='-DPXD' in flags
    vectors=10+2*int(smooth)+2*int(pxd)
    vec=lambda a:(a+1)*8
    mat=lambda a,c:(a*c+1)*8+(a+1)*sizes['size_t']
    tensor=lambda a,c,d:(a*c*d+1)*8+(a*c+1)*sizes['size_t']+(a+1)*sizes['size_t']
    common=vectors*vec(b)
    if scalar3:
        common+=2*mat(m,b)+9*tensor(m,b,b)+(mat(b,b) if pxd else 0)
    known={"native_catalogs":n*sizes['body'],"common_histograms":common,"retained_caches":int(retained_cache_bytes)}
    estimates={}
    if engine.startswith('lya-'):
        rb=positive(p.get('lya3RBins',20),'lya3RBins');tb=positive(p.get('lya3ThetaBins',10),'lya3ThetaBins')
        mu=positive(p.get('lya3MuBins',20),'lya3MuBins');rp=positive(p.get('lya2RpBins',50),'lya2RpBins');rt=positive(p.get('lya2RtBins',50),'lya2RtBins')
        radial='1d-' in engine;pair='2pcf' in engine;triple='3pcf' in engine
        cells2=(rp if radial else rp*rt) if pair else 0
        cells3=(4*rb*rb if radial else rb*rb*tb*tb*(positive(p.get('lya3LMax',8)+1,'lya3LMax orders') if 'multipole' in engine else mu)) if triple else 0
        if not radial:
            known['estimator_global_and_worker_histograms']=(cells2+cells3)*(2*sizes['real']+t*(2*sizes['real']+sizes['size_t']))
        else:
            estimates['radial_global_and_worker_histograms']=(cells2+cells3)*(2*sizes['long_double']*(t+1)+t*sizes['size_t'])
        estimates['native_octree_at_one_cell_per_pixel']=n*sizes['octree_cell']
        estimates['forest_indexes_and_neighbors']=n*(96+96*t)
    elif scalar or 'shear' in engine:
        estimates['compact_tree_at_two_nodes_per_pixel']=n*(2*sizes['binary_node']+sizes['packed_point'])
        if scalar3:estimates['scalar_worker_histograms']=t*(4*m*b*b*8+8*b*b*8+2*m*b*8)
        elif 'shear' in engine:estimates['shear_worker_moments']=t*8*(16*m*b*b+16*m*b)
    else:
        estimates['native_octree_at_one_cell_per_pixel']=n*sizes['octree_cell']
        if '3pcf-3d' in engine or 'ggg-3d' in engine:estimates['physical_moments']=2*(t+1)*m*b*b*sizes['real']
    limit=(1 << (8*sizes['size_t']-1))-1
    if any(v<0 or v>limit for v in [*known.values(),*estimates.values(),sum(known.values())+sum(estimates.values())]):
        raise OverflowError('resource dimensions exceed PTRDIFF_MAX')
    return {"engine":engine,"pixels":n,"threads":t,"parameters":p,"policy":policy,
            "known_components_bytes":known,"estimated_components_bytes":estimates,
            "known_total_bytes":sum(known.values()),"estimated_total_bytes":sum(known.values())+sum(estimates.values()),
            "known_exceeds_budget":sum(known.values())>policy['budget_bytes'],
            "estimate_exceeds_budget":sum(known.values())+sum(estimates.values())>policy['budget_bytes'],
            "exclusions":["Python/user arrays","reader/MPI/libraries","allocator slack","data-dependent scratch/frontier/window solve"],
            "tree_estimate_is_upper_bound":False}
