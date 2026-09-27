#ifndef _input_lya_forest_omp_01_h
#define _input_lya_forest_omp_01_h

    PARSER_READ(parser_read_double(pfc, "lya2RpMax", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya2RpMax = param1;
    PARSER_READ(parser_read_double(pfc, "lya2RtMax", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya2RtMax = param1;
    PARSER_READ(parser_read_int(pfc, "lya2RpBins", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya2RpBins = param;
    PARSER_READ(parser_read_int(pfc, "lya2RtBins", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya2RtBins = param;
    PARSER_READ(parser_read_double(pfc, "lya3RMax", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya3RMax = param1;
    PARSER_READ(parser_read_int(pfc, "lya3RBins", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3RBins = param;
    PARSER_READ(parser_read_int(pfc, "lya3ThetaBins", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3ThetaBins = param;
    PARSER_READ(parser_read_int(pfc, "lya3MuBins", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3MuBins = param;

    PARSER_READ(parser_read_int(pfc, "lya3Kernel", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3Kernel = param;

    PARSER_READ(parser_read_int(pfc, "lya3LMax", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3LMax = param;

    PARSER_READ(parser_read_int(pfc, "lya3PivotBlock", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3PivotBlock = param;

    PARSER_READ(parser_read_double(pfc, "lya3MuSlop", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya3MuSlop = param1;

    PARSER_READ(parser_read_double(pfc, "lya3RadialSlop", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya3RadialSlop = param1;

    PARSER_READ(parser_read_double(pfc, "lya3PolarSlop", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya3PolarSlop = param1;

    PARSER_READ(parser_read_int(pfc, "lya3PivotCellMax", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya3PivotCellMax = param;

    PARSER_READ(parser_read_int(pfc, "lya2Kernel", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lya2Kernel = param;
    PARSER_READ(parser_read_double(pfc, "lya2RpSlop", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya2RpSlop = param1;
    PARSER_READ(parser_read_double(pfc, "lya2RtSlop", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lya2RtSlop = param1;

    PARSER_READ(parser_read_int(pfc, "lyaScanLevel", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lyaScanLevel = param;
    PARSER_READ(parser_read_double(pfc, "lyaPivotRadius", &param1, &flag1, errmsg));
    if (flag1 == TRUE) cmd->lyaPivotRadius = param1;
    PARSER_READ(parser_read_int(pfc, "lyaPivotMax", &param, &flag, errmsg));
    if (flag == TRUE) cmd->lyaPivotMax = param;

#endif
