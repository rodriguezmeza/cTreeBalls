#ifndef CTREEBALLS_TREE_CONTRACTS_H
#define CTREEBALLS_TREE_CONTRACTS_H

/*
 * Shared numerical contracts used by the tree builder, search routines, and
 * focused regression tests. Include this header after globaldefs.h.
 */
static inline bool cballs_method_needs_balls4_scan(int search_method)
{
    bool needed = false;
#ifdef BALLS4SCANLEV
#ifdef OCTREEGGGOMP
    needed |= search_method == OCTREEGGGOMPMETHOD;
#endif
#ifdef OCTREEGGGMPI
    needed |= search_method == OCTREEGGGMPIMETHOD;
#endif
#ifdef OCTREE3PCF3DOMP
    needed |= search_method == OCTREE3PCF3DOMPMETHOD;
#endif
#ifdef OCTREE3PCF3DMPI
    needed |= search_method == OCTREE3PCF3DMPIMETHOD;
#endif
#ifdef OCTREESHEAROMP
    needed |= search_method == OCTREESHEARMETHOD;
#endif
#ifdef OCTREESHEARSPHEREOMP
    needed |= search_method == OCTREESHEARSPHEREMETHOD;
#endif
#ifdef OCTREESHEARSPHERE2BALLSOMP
    needed |= search_method == OCTREESHEARSPHERE2BALLSOMPMETHOD;
#endif
#endif
#ifdef OCTREEBALLS4OMP
    needed |= search_method == OCTREEBALLS4OMPMETHOD;
#endif
#ifdef OCTREEBALLS4MPI
    needed |= search_method == OCTREEBALLS4MPIMETHOD;
#endif
    (void)search_method;
    return needed;
}

static inline bool cballs_run_needs_balls4_scan(
        const struct cmdline_data *cmd, int search_method)
{
    bool needed = cballs_method_needs_balls4_scan(search_method);
#ifdef OCTREE2BALLSOMP
    needed |= search_method == OCTREE2BALLSMETHOD
        && cmd != NULL && cballs_opt_legacy_one_ball(cmd);
#endif
#ifdef OCTREE2BALLSMPI
    needed |= search_method == OCTREE2BALLSMPIMETHOD
        && cmd != NULL && cballs_opt_legacy_one_ball(cmd);
#endif
    if (!needed)
        return false;
#ifdef OCTREESHEARSPHEREOMP
    /* The unsmoothed spherical only-2pcf dual-tree owns a fixed node-pair
     * frontier and never consumes the body-pivot scan table. Smoothing needs
     * representative body pivots, as do exact and mixed-order runs. */
    if (search_method == OCTREESHEARSPHEREMETHOD
        && cmd != NULL && cballs_opt_only_2pcf(cmd) && cmd->theta > 0.0
        && !cballs_opt_no_one_ball(cmd) && !cballs_opt_smooth_pivot(cmd))
        return false;
#endif
#ifdef OCTREESHEARSPHERE2BALLSOMP
    if (search_method == OCTREESHEARSPHERE2BALLSOMPMETHOD
        && cmd != NULL && cballs_opt_only_2pcf(cmd)
        && !cballs_opt_smooth_pivot(cmd))
        return false;
#endif
    (void)cmd;
    return true;
}

static inline bool cballs_method_uses_compact_native_octree(int search_method)
{
    bool uses_compact_tree = false;

#ifdef OCTREE2BALLSOMP
    uses_compact_tree |= search_method == OCTREE2BALLSMETHOD;
#endif
#ifdef OCTREE2BALLSMPI
    uses_compact_tree |= search_method == OCTREE2BALLSMPIMETHOD;
#endif
    (void)search_method;
    return uses_compact_tree;
}

static inline bool cballs_run_uses_compact_native_pair(
        const struct cmdline_data *cmd, int search_method)
{
    bool uses_compact_pair =
        cballs_method_uses_compact_native_octree(search_method);

#ifdef OCTREE2BALLSOMP
    if (search_method == OCTREE2BALLSMETHOD
        && cmd != NULL && cballs_opt_legacy_one_ball(cmd))
        return false;
#endif
#ifdef OCTREE2BALLSMPI
    if (search_method == OCTREE2BALLSMPIMETHOD
        && cmd != NULL && cballs_opt_legacy_one_ball(cmd))
        return false;
#endif

    if (!cballs_opt_only_2pcf(cmd)) return uses_compact_pair;
#ifdef OCTREEBALLS4OMP
    uses_compact_pair |= search_method == OCTREEBALLS4OMPMETHOD;
#endif
#ifdef OCTREEBALLS4MPI
    uses_compact_pair |= search_method == OCTREEBALLS4MPIMETHOD;
#endif
    return uses_compact_pair;
}

static inline bool cballs_cell_accumulate_child(nodeptr parent, nodeptr child,
                                                 bool read_mask,
                                                 short *cell_mask,
                                                 compute_vector center_of_mass_sum)
{
    bool child_is_valid = !read_mask || Mask(child) != MASK_NODE_MASKED;
    compute_vector weighted_position;
#ifdef THREEPCFSHEAR
    real child_weight_sum;
#endif

    Selected(parent) |= Selected(child);
    Update(parent) |= Update(child);
    if (read_mask && cell_mask != NULL)
        *cell_mask = mask_node_combine(*cell_mask, Mask(child));

    Mass(parent) += Mass(child);
    if (child_is_valid) {
#ifdef THREEPCFSHEAR
        child_weight_sum = Type(child) == CELL
            ? ShearWeightSum(child) : Weight(child);
        ShearWeightSum(parent) += child_weight_sum;
        Gamma1(parent) += child_weight_sum*Gamma1(child);
        Gamma2(parent) += child_weight_sum*Gamma2(child);
        if (Type(child) == CELL) {
            ShearWeight2(parent) += ShearWeight2(child);
            ShearGamma2Re(parent) += ShearGamma2Re(child);
            ShearGamma2Im(parent) += ShearGamma2Im(child);
            ShearGammaAbs2(parent) += ShearGammaAbs2(child);
        } else {
            const real weighted_gamma1 = Weight(child)*Gamma1(child);
            const real weighted_gamma2 = Weight(child)*Gamma2(child);
            ShearWeight2(parent) += Weight(child)*Weight(child);
            ShearGamma2Re(parent) += weighted_gamma1*weighted_gamma1
                                         - weighted_gamma2*weighted_gamma2;
            ShearGamma2Im(parent) += 2.0*weighted_gamma1*weighted_gamma2;
            ShearGammaAbs2(parent) += weighted_gamma1*weighted_gamma1
                                      + weighted_gamma2*weighted_gamma2;
        }
#endif
        Weight(parent) += Weight(child);
#ifndef NOWKAvg
        Kappa(parent) += Weight(child)*Kappa(child);
#else
        if (Type(child) == CELL)
            Kappa(parent) += Nb(child)*Kappa(child);
        else
            Kappa(parent) += Kappa(child);
#endif

        if (Type(child) == CELL)
            Nb(parent) += Nb(child);
        else if (Type(child) == BODY || Type(child) == BODY3)
            Nb(parent) += 1;
    }

    MULVS(weighted_position, Pos(child), Mass(child));
    ADDV(center_of_mass_sum, center_of_mass_sum, weighted_position);

    return !child_is_valid
        && (Kappa(child) != 0.0 || Weight(child) != 0.0);
}

static inline void cballs_cell_finalize_averages(nodeptr parent)
{
    if (Nb(parent) <= 0)
        return;

#ifdef THREEPCFSHEAR
    if (ShearWeightSum(parent) > 0.0) {
        Gamma1(parent) /= ShearWeightSum(parent);
        Gamma2(parent) /= ShearWeightSum(parent);
    } else {
        Gamma1(parent) = 0.0;
        Gamma2(parent) = 0.0;
    }
#endif
#ifndef NOWKAvg
    if (Weight(parent) > 0.0)
        Kappa(parent) /= Weight(parent);
    else
        Kappa(parent) = 0.0;
#else
    Kappa(parent) /= Nb(parent);
    Weight(parent) /= Nb(parent);
#endif
}

static inline bool cballs_accept_body_contract(struct cmdline_data *cmd,
                                                struct global_data *gd,
                                                bodyptr p, nodeptr q,
                                                real *distance,
                                                compute_vector dr)
{
    real distance_squared;

    DOTPSUBV(distance_squared, dr, Pos(p), Pos(q));
    if (cmd->usePeriodic) {
        VWrapAll(dr);
        DOTVP(distance_squared, dr, dr);
    }
    *distance = rsqrt(distance_squared);

    return *distance < gd->Rcut;
}

#endif /* !CTREEBALLS_TREE_CONTRACTS_H */
