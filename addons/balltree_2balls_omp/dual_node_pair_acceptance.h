/* Shared exact/approximate radial contract; included after the backend context.
 * Independent enumerations remain in tests/make_tests. */
static inline bool dual_node_pair_outside_range(
        const dual_node_search_context *context,
        const fcfc_ballnode *node1, const fcfc_ballnode *node2,
        real distance2)
{
    const real size = (real)node1->radius + (real)node2->radius;
    const real minimum = context->cmd->rminHist;
    const real maximum = context->cmd->rangeN;

    if (size <= minimum && distance2 <= rsqr(minimum - size)) return TRUE;
    if (distance2 >= rsqr(maximum + size)) return TRUE;
    return FALSE;
}

/* Return the common bin, -1 when the pair must be split, or -2 when an
 * accepted approximate pair has its centre outside the histogram range. */
static int dual_node_two_ball_bin(const dual_node_search_context *context,
                                 const fcfc_ballnode *node1,
                                 const fcfc_ballnode *node2, real distance2)
{
    const struct cmdline_data *cmd = context->cmd;
    const real size = (real)node1->radius + (real)node2->radius;
    real distance;
    real lower;
    real upper;
    int center_bin;

    if (!context->use_two_balls || !(distance2 > 0.0)
        || !(cmd->theta > 0.0))
        return -1;

    if (context->use_bin_theta) {
        real fraction;
        real coordinate;

        if (size * size > context->theta2 * distance2) return -1;
        if (cmd->useLogHist) {
            const real relative_size2 = size * size / distance2;

            if (!(cmd->rminHist > 0.0)) return -1;
            if (size * size > context->bin_theta2 * distance2) {
                if (size * size > context->half_bin_plus_slop2 * distance2)
                    return -1;
                coordinate = 0.5 * rlog(
                    distance2 / context->minimum2)
                    / context->logarithmic_bin_size;
                fraction = coordinate - rfloor(coordinate);
                if (fraction > 0.5) fraction = 1.0 - fraction;
                if (size * size
                    > rsqr(fraction * context->bin_size
                          + context->bin_theta) * distance2)
                    return -1;
                fraction = coordinate - rfloor(coordinate);
                if (size * size
                    > rsqr(fraction * context->bin_size
                          + context->bin_theta - relative_size2)
                      * distance2)
                    return -1;
                center_bin = (int)coordinate + 1;
                return center_bin >= 1 && center_bin <= cmd->sizeHistN
                    ? center_bin : -2;
            }
        } else {
            if (size > context->bin_theta) {
                if (size > 0.5 * (context->bin_size
                                  + context->bin_theta)) return -1;
                distance = rsqrt(distance2);
                coordinate = (distance - cmd->rminHist)
                           / context->bin_size;
                fraction = coordinate - rfloor(coordinate);
                if (fraction > 0.5) fraction = 1.0 - fraction;
                if (size > fraction * context->bin_size
                         + context->bin_theta) return -1;
            }
        }
        center_bin = dual_node_bin_index_squared(context, distance2);
        return center_bin < 0 ? -2 : center_bin;
    }

    distance = rsqrt(distance2);
    if (size > dual_node_bin_theta_width(context, distance)) return -1;
    lower = distance - size;
    upper = distance + size;
    if (!(lower > cmd->rminHist && upper < cmd->rangeN)) return -1;
    center_bin = dual_node_bin_index(context, distance);
    if (center_bin < 0
        || dual_node_bin_index(context, lower) != center_bin
        || dual_node_bin_index(context, upper) != center_bin)
        return -1;
    return center_bin;
}

