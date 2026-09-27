#ifndef _startrun_lya_forest_omp_07_h
#define _startrun_lya_forest_omp_07_h

/* These controls are deliberately independent of legacy smooth-pivot/rsmooth:
 * enabling a legacy build flag must not silently approximate a forest result. */
if (cmd->lyaScanLevel<0 || cmd->lyaScanLevel>20
    || !isfinite(cmd->lyaPivotRadius) || cmd->lyaPivotRadius<0
    || cmd->lyaPivotMax<1 || cmd->lyaPivotMax>1024)
    cBALLS_FAIL(cmd, "%s: require lyaScanLevel=0..20, finite lyaPivotRadius>=0 and lyaPivotMax=1..1024", routineName);
if (cmd->lyaScanLevel || cmd->lyaPivotRadius>0) {
    const int pivot_kind=lya_forest_method_kind(cmd->searchMethod);
    if (pivot_kind<0 || pivot_kind>2 || lya_forest_is_mpi_method(cmd->searchMethod)
        || lya_forest_is_multipole_method(cmd->searchMethod)
        || (pivot_kind!=0 && cmd->lya3Kernel>=3)
        || (pivot_kind==0 && cmd->lya2Kernel!=0))
        cBALLS_FAIL(cmd, "%s: pivot frontier/smoothing requires a 3D OpenMP hard-bin pixel-pivot path (2PCF kernel 0 or 3PCF kernels 0..2)", routineName);
}

/* Pair-cell controls are independent of all 3PCF controls. */
if (cmd->lya2Kernel < 0 || cmd->lya2Kernel > 1
    || !isfinite(cmd->lya2RpSlop) || cmd->lya2RpSlop < 0 || cmd->lya2RpSlop > 1
    || !isfinite(cmd->lya2RtSlop) || cmd->lya2RtSlop < 0 || cmd->lya2RtSlop > 1)
    cBALLS_FAIL(cmd, "%s: require lya2Kernel=0..1 and pair slops finite in [0,1]", routineName);
if (cmd->lya2Kernel || cmd->lya2RpSlop > 0 || cmd->lya2RtSlop > 0) {
    if (cmd->lya2Kernel != 1)
        cBALLS_FAIL(cmd, "%s: pair slop requires lya2Kernel=1", routineName);
    if ((lya_forest_method_kind(cmd->searchMethod) != 0
         && lya_forest_method_kind(cmd->searchMethod) != 2)
        || lya_forest_is_mpi_method(cmd->searchMethod))
        cBALLS_FAIL(cmd, "%s: pair cells require a 3D OpenMP Ly-alpha 2PCF or combined method", routineName);
}

if (lya_forest_method_kind(cmd->searchMethod) >= 0
    && lya_forest_method_kind(cmd->searchMethod) < 3) {
#if NDIM != 3
    cBALLS_FAIL(cmd, "%s: %s requires DEFDIMENSION=3",
                routineName, cmd->searchMethod);
#endif
    if (cmd->lya3PivotBlock < 0 || cmd->lya3PivotBlock > 4096)
        cBALLS_FAIL(cmd, "%s: require lya3PivotBlock=0..4096", routineName);
    if (cmd->lya3Kernel < 0 || cmd->lya3Kernel > 4 || cmd->lya3LMax < 0 || cmd->lya3LMax > 32)
        cBALLS_FAIL(cmd, "%s: require lya3Kernel=0..4 and lya3LMax=0..32", routineName);
    if (!isfinite(cmd->lya3MuSlop) || cmd->lya3MuSlop < 0 || cmd->lya3MuSlop > 1
        || !isfinite(cmd->lya3RadialSlop) || cmd->lya3RadialSlop < 0 || cmd->lya3RadialSlop > 1
        || !isfinite(cmd->lya3PolarSlop) || cmd->lya3PolarSlop < 0 || cmd->lya3PolarSlop > 1)
        cBALLS_FAIL(cmd, "%s: Ly-alpha geometry slops must be finite in [0,1]", routineName);
    if (cmd->lya3PivotCellMax < 1 || cmd->lya3PivotCellMax > 64)
        cBALLS_FAIL(cmd, "%s: require lya3PivotCellMax=1..64", routineName);
    if ((cmd->lya3RadialSlop > 0 || cmd->lya3PolarSlop > 0) && cmd->lya3Kernel < 3)
        cBALLS_FAIL(cmd, "%s: radial/polar slop requires lya3Kernel=3 or 4", routineName);
    if (cmd->lya3MuSlop > 0 && (cmd->lya3Kernel == 1 || cmd->lya3Kernel == 2))
        cBALLS_FAIL(cmd, "%s: direct kernels 1/2 do not support geometry slop", routineName);
    if ((cmd->lya3MuSlop > 0 || cmd->lya3RadialSlop > 0 || cmd->lya3PolarSlop > 0 || cmd->lya3Kernel >= 3)
        && (lya_forest_method_kind(cmd->searchMethod) == 0 || lya_forest_is_multipole_method(cmd->searchMethod)))
        cBALLS_FAIL(cmd, "%s: cell geometry controls require a hard-bin 3PCF method", routineName);
    if (cmd->lya3Kernel >= 3 && lya_forest_is_mpi_method(cmd->searchMethod))
        cBALLS_FAIL(cmd, "%s: persistent/pivot cell kernels currently require an OpenMP method", routineName);
    if (cmd->usePeriodic)
        cBALLS_FAIL(cmd, "%s: %s uses observer-centered lines of sight and "
                    "does not support periodic boundaries",
                    routineName, cmd->searchMethod);
    if (strcmp(cmd->infilefmt, "lya-ascii") != 0)
        cBALLS_FAIL(cmd, "%s: %s requires infileformat=lya-ascii",
                    routineName, cmd->searchMethod);
    if (gd->ninfiles != 1)
        cBALLS_FAIL(cmd, "%s: %s requires one flattened input catalog",
                    routineName, cmd->searchMethod);

    if (lya_forest_method_kind(cmd->searchMethod) != 1
        && (!isfinite(cmd->lya2RpMax) || cmd->lya2RpMax <= 0.0
            || !isfinite(cmd->lya2RtMax) || cmd->lya2RtMax <= 0.0
            || cmd->lya2RpBins < 1 || cmd->lya2RtBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid Lyman-alpha 2PCF domain or bin count",
                    routineName);

    if (lya_forest_method_kind(cmd->searchMethod) != 0
        && (!isfinite(cmd->lya3RMax) || cmd->lya3RMax <= 0.0
            || cmd->lya3RBins < 1 || cmd->lya3ThetaBins < 1
            || cmd->lya3MuBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid Lyman-alpha 3PCF domain or bin count",
                    routineName);
}

if (lya_forest_method_kind(cmd->searchMethod) >= 3) {
#if NDIM != 3
    cBALLS_FAIL(cmd, "%s: %s requires DEFDIMENSION=3 for catalog input",
                routineName, cmd->searchMethod);
#endif
    if (cmd->usePeriodic)
        cBALLS_FAIL(cmd, "%s: %s uses observer-centered radial distances and "
                    "does not support periodic boundaries",
                    routineName, cmd->searchMethod);
    if (strcmp(cmd->infilefmt, "lya-ascii") != 0)
        cBALLS_FAIL(cmd, "%s: %s requires infileformat=lya-ascii",
                    routineName, cmd->searchMethod);
    if (gd->ninfiles != 1)
        cBALLS_FAIL(cmd, "%s: %s requires one flattened input catalog",
                    routineName, cmd->searchMethod);

    if (lya_forest_method_kind(cmd->searchMethod) != 4
        && lya_forest_method_kind(cmd->searchMethod) != 7
        && (!isfinite(cmd->lya2RpMax) || cmd->lya2RpMax <= 0.0
            || cmd->lya2RpBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid radial Ly-alpha 2PCF domain or bin count",
                    routineName);

    if (lya_forest_method_kind(cmd->searchMethod) != 3
        && lya_forest_method_kind(cmd->searchMethod) != 6
        && lya_forest_method_kind(cmd->searchMethod) != 8
        && (!isfinite(cmd->lya3RMax) || cmd->lya3RMax <= 0.0
            || cmd->lya3RBins < 1 || cmd->lya3RBins > INT_MAX / 2))
        cBALLS_FAIL(cmd, "%s: invalid radial Ly-alpha 3PCF domain or bin count",
                    routineName);
}

#endif
