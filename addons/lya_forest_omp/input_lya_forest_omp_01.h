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

#endif
