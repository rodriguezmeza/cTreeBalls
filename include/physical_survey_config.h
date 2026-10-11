#ifndef CBALLS_PHYSICAL_SURVEY_CONFIG_H
#define CBALLS_PHYSICAL_SURVEY_CONFIG_H
/* Include after cmdline_data and scanopt declarations. */
static inline int cb_survey_signal_lmax(const struct cmdline_data *cmd)
{
    return cmd->surveySignalLmax>=0 ? cmd->surveySignalLmax
        : cmd->mChebyshev-(scanopt(cmd->options,"survey-keep-top-multipole") ? 0 : 1);
}
static inline int cb_survey_solve_lmax(const struct cmdline_data *cmd)
{return cb_survey_signal_lmax(cmd)+cmd->surveySolvePadding;}
static inline int cb_survey_window_lmax(const struct cmdline_data *cmd)
{return cmd->surveyWindowLmax>=0 ? cmd->surveyWindowLmax : 2*cb_survey_solve_lmax(cmd);}
#endif
