#ifndef _startrun_lya_forest_omp_07_h
#define _startrun_lya_forest_omp_07_h

if (strcmp(cmd->searchMethod, "lya-2pcf-omp") == 0
    || strcmp(cmd->searchMethod, "lya-3pcf-omp") == 0
    || strcmp(cmd->searchMethod, "lya-2pcf-3pcf-omp") == 0) {
#if NDIM != 3
    cBALLS_FAIL(cmd, "%s: %s requires DEFDIMENSION=3",
                routineName, cmd->searchMethod);
#endif
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

    if ((strcmp(cmd->searchMethod, "lya-3pcf-omp") != 0)
        && (!isfinite(cmd->lya2RpMax) || cmd->lya2RpMax <= 0.0
            || !isfinite(cmd->lya2RtMax) || cmd->lya2RtMax <= 0.0
            || cmd->lya2RpBins < 1 || cmd->lya2RtBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid Lyman-alpha 2PCF domain or bin count",
                    routineName);

    if ((strcmp(cmd->searchMethod, "lya-2pcf-omp") != 0)
        && (!isfinite(cmd->lya3RMax) || cmd->lya3RMax <= 0.0
            || cmd->lya3RBins < 1 || cmd->lya3ThetaBins < 1
            || cmd->lya3MuBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid Lyman-alpha 3PCF domain or bin count",
                    routineName);
}

if (strcmp(cmd->searchMethod, "lya-1d-2pcf-omp") == 0
    || strcmp(cmd->searchMethod, "lya-1d-3pcf-omp") == 0
    || strcmp(cmd->searchMethod, "lya-1d-2pcf-3pcf-omp") == 0
    || strcmp(cmd->searchMethod, "lya-1d-tree-2pcf-omp") == 0) {
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

    if ((strcmp(cmd->searchMethod, "lya-1d-3pcf-omp") != 0)
        && (!isfinite(cmd->lya2RpMax) || cmd->lya2RpMax <= 0.0
            || cmd->lya2RpBins < 1))
        cBALLS_FAIL(cmd, "%s: invalid radial Ly-alpha 2PCF domain or bin count",
                    routineName);

    if ((strcmp(cmd->searchMethod, "lya-1d-2pcf-omp") != 0)
        && (strcmp(cmd->searchMethod, "lya-1d-tree-2pcf-omp") != 0)
        && (!isfinite(cmd->lya3RMax) || cmd->lya3RMax <= 0.0
            || cmd->lya3RBins < 1 || cmd->lya3RBins > INT_MAX / 2))
        cBALLS_FAIL(cmd, "%s: invalid radial Ly-alpha 3PCF domain or bin count",
                    routineName);
}

#endif
