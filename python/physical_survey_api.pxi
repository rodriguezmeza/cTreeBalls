cdef extern from "legendre_window_solver.h":
    int cb_survey_solve(const double *, const double *, int, int, double,
                        double *, double *, double *)


def solve_physical_survey(numerator, random_window, *, signal_lmax, padding=1, rcond_min=1e-10):
    """Correct standard Legendre moments using the same native survey solver.

    Axis zero is ell. Numerator must contain signal_lmax+padding+1 modes;
    random_window must contain at least 2*(signal_lmax+padding)+1 modes.
    Other axes match. Inputs are copied. Returned invalid coefficients are NaN.
    Status: 0 empty, 1 valid, 2 singular, 3 ill-conditioned, 4 residual,
    5 nonfinite/unrepresentable. rcond uses the infinity norm in this basis;
    relative_residual is a backward error, not a truncation-error estimate.
    """
    cdef int solve_lmax, window_lmax, code
    cdef Py_ssize_t i, j, count
    cdef double rc, residual, threshold
    cdef np.ndarray[np.float64_t, ndim=2] num, ran
    cdef np.ndarray[np.float64_t, ndim=1] x, nvec, rvec
    if isinstance(signal_lmax,bool) or int(signal_lmax)!=signal_lmax or not 0<=signal_lmax<=32:
        raise ValueError('signal_lmax must be an integer in 0..32')
    if isinstance(padding,bool) or int(padding)!=padding or not 0<=padding<=32-signal_lmax:
        raise ValueError('padding must be a nonnegative integer; signal_lmax+padding<=32')
    threshold=float(rcond_min)
    if not np.isfinite(threshold) or not 0<threshold<=1:
        raise ValueError('rcond_min must be finite and in (0,1]')
    solve_lmax=int(signal_lmax+padding)
    if np.iscomplexobj(numerator) or np.iscomplexobj(random_window):
        raise ValueError('physical Legendre moments must be real')
    a=np.asarray(numerator,dtype=np.float64); b=np.asarray(random_window,dtype=np.float64)
    if a.ndim<1 or b.ndim!=a.ndim or a.shape[1:]!=b.shape[1:] or a.shape[0]!=solve_lmax+1 or not 2*solve_lmax+1<=b.shape[0]<=65:
        raise ValueError('matching trailing axes and full numerator/window bandwidths are required')
    shape=a.shape[1:];count=int(np.prod(shape,dtype=np.int64)) if shape else 1
    window_lmax=b.shape[0]-1
    num=np.ascontiguousarray(a.reshape((solve_lmax+1,count)))
    ran=np.ascontiguousarray(b.reshape((window_lmax+1,count)))
    nvec=np.empty(solve_lmax+1);rvec=np.empty(window_lmax+1);x=np.empty(solve_lmax+1)
    values=np.full((int(signal_lmax)+1,count),np.nan)
    statuses=np.empty(count,dtype=np.int32);condition=np.empty(count);residuals=np.empty(count)
    for i in range(count):
        for j in range(solve_lmax+1):nvec[j]=num[j,i]
        for j in range(window_lmax+1):rvec[j]=ran[j,i]
        code=cb_survey_solve(&nvec[0],&rvec[0],solve_lmax,window_lmax,threshold,&x[0],&rc,&residual)
        if code==1:values[:,i]=x[:int(signal_lmax)+1]
        statuses[i]=code;condition[i]=rc;residuals[i]=residual
    return dict(zeta_corrected_legendre=values.reshape((int(signal_lmax)+1,)+shape),
        zeta_valid=(statuses==1).reshape(shape),zeta_status=statuses.reshape(shape),
        zeta_rcond=condition.reshape(shape),zeta_relative_residual=residuals.reshape(shape))
