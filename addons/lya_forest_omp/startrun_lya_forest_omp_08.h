#ifndef _startrun_lya_forest_omp_08_h
#define _startrun_lya_forest_omp_08_h

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya2RpMax", cmd->lya2RpMax);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya2RtMax", cmd->lya2RtMax);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya2RpBins", cmd->lya2RpBins);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya2RtBins", cmd->lya2RtBins);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya3RMax", cmd->lya3RMax);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3RBins", cmd->lya3RBins);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3ThetaBins", cmd->lya3ThetaBins);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3MuBins", cmd->lya3MuBins);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3Kernel", cmd->lya3Kernel);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3LMax", cmd->lya3LMax);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3PivotBlock", cmd->lya3PivotBlock);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya3MuSlop", cmd->lya3MuSlop);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya3RadialSlop", cmd->lya3RadialSlop);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya3PolarSlop", cmd->lya3PolarSlop);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya3PivotCellMax", cmd->lya3PivotCellMax);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lya2Kernel", cmd->lya2Kernel);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya2RpSlop", cmd->lya2RpSlop);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lya2RtSlop", cmd->lya2RtSlop);

    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lyaScanLevel", cmd->lyaScanLevel);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTR, "lyaPivotRadius", cmd->lyaPivotRadius);
    WRITE_OUTPUT_OR_FAIL(fdout, buf, FMTI, "lyaPivotMax", cmd->lyaPivotMax);

#endif
