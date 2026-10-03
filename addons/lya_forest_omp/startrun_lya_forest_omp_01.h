#ifndef _startrun_lya_forest_omp_01_h
#define _startrun_lya_forest_omp_01_h

    cmd->lya2RpMax = GetdParam("lya2RpMax");
    cmd->lya2RtMax = GetdParam("lya2RtMax");
    cmd->lya2RpBins = GetiParam("lya2RpBins");
    cmd->lya2RtBins = GetiParam("lya2RtBins");
    cmd->lya3RMax = GetdParam("lya3RMax");
    cmd->lya3RBins = GetiParam("lya3RBins");
    cmd->lya3ThetaBins = GetiParam("lya3ThetaBins");
    cmd->lya3MuBins = GetiParam("lya3MuBins");

    cmd->lya3Kernel = GetiParam("lya3Kernel");

    cmd->lya3LMax = GetiParam("lya3LMax");
    cmd->lya3MuMode = GetiParam("lya3MuMode");

    cmd->lya3PivotBlock = GetiParam("lya3PivotBlock");

    cmd->lya3MuSlop = GetdParam("lya3MuSlop");

    cmd->lya3RadialSlop = GetdParam("lya3RadialSlop");

    cmd->lya3PolarSlop = GetdParam("lya3PolarSlop");

    cmd->lya3PivotCellMax = GetiParam("lya3PivotCellMax");

    cmd->lya2Kernel = GetiParam("lya2Kernel");
    cmd->lya2RpSlop = GetdParam("lya2RpSlop");
    cmd->lya2RtSlop = GetdParam("lya2RtSlop");

    cmd->lyaScanLevel = GetiParam("lyaScanLevel");
    cmd->lyaPivotRadius = GetdParam("lyaPivotRadius");
    cmd->lyaPivotMax = GetiParam("lyaPivotMax");

#endif
