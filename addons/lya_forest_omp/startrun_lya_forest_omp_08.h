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

#endif
