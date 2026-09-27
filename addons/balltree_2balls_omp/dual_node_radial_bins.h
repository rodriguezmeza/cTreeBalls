/* Shared exact/approximate radial contract; included after the backend context.
 * Independent enumerations remain in tests/make_tests. */
static inline int dual_node_bin_index(const dual_node_search_context *context,
                                     real distance)
{
    const struct cmdline_data *cmd = context->cmd;
    int n;

    if (!(distance > cmd->rminHist && distance < cmd->rangeN))
        return -1;
    if (cmd->useLogHist) {
#ifdef DUAL_NODE_USE_NATURAL_LOG_BINS
        if (cmd->rminHist == 0.0) {
            n = (int)(context->natural_log_scale
                * (rlog(distance) - context->logarithmic_maximum)
                + cmd->sizeHistN) + 1;
        } else {
            n = (int)((rlog(distance)
                - 0.5 * context->logarithmic_minimum2)
                * context->natural_log_scale) + 1;
        }
#else
        if (cmd->rminHist == 0.0) {
            n = (int)(cmd->logHistBinsPD
                * (rlog10(distance) - rlog10(cmd->rangeN))
                + cmd->sizeHistN) + 1;
        } else {
            n = (int)(rlog10(distance / cmd->rminHist)
                * context->gd->i_deltaR) + 1;
        }
#endif
    } else {
        n = (int)((distance - cmd->rminHist)
            * context->gd->i_deltaR) + 1;
    }
    return n >= 1 && n <= cmd->sizeHistN ? n : -1;
}

static inline int dual_node_bin_index_squared(
        const dual_node_search_context *context, real distance2)
{
    const struct cmdline_data *cmd = context->cmd;
    int n;

    if (!(distance2 > context->minimum2
          && distance2 < context->maximum2)) return -1;
    if (!cmd->useLogHist)
        return dual_node_bin_index(context, rsqrt(distance2));
#ifdef DUAL_NODE_USE_NATURAL_LOG_BINS
    if (cmd->rminHist == 0.0) {
        n = (int)((cmd->logHistBinsPD / rlog(10.0))
            * (0.5 * rlog(distance2) - context->logarithmic_maximum)
            + cmd->sizeHistN) + 1;
    } else {
        n = (int)(0.5
            * (rlog(distance2) - context->logarithmic_minimum2)
            * context->natural_log_scale) + 1;
    }
#else
    if (cmd->rminHist == 0.0) {
        n = (int)(cmd->logHistBinsPD
            * (0.5 * rlog10(distance2) - rlog10(cmd->rangeN))
            + cmd->sizeHistN) + 1;
    } else {
        n = (int)(0.5 * rlog10(distance2 / context->minimum2)
            * context->gd->i_deltaR) + 1;
    }
#endif
    return n >= 1 && n <= cmd->sizeHistN ? n : -1;
}

/* dual-node's b = bin_theta * bin_size, expressed in distance units. */
static inline real dual_node_bin_theta_width(
        const dual_node_search_context *context, real distance)
{
    if (context->cmd->useLogHist)
        return context->bin_theta * distance;
    return context->bin_theta;
}

