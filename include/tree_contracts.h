#ifndef CTREEBALLS_TREE_CONTRACTS_H
#define CTREEBALLS_TREE_CONTRACTS_H

/*
 * Shared numerical contracts used by the tree builder, search routines, and
 * focused regression tests. Include this header after globaldefs.h.
 */
static inline bool cballs_method_needs_balls4_scan(int search_method)
{
#ifdef BALLS4SCANLEV
    (void)search_method;
    return true;
#else
    bool needed = false;
#ifdef OCTREEBALLS4OMP
    needed |= search_method == OCTREEBALLS4OMPMETHOD;
#endif
#ifdef OCTREEBALLS4MPI
    needed |= search_method == OCTREEBALLS4MPIMETHOD;
#endif
    (void)search_method;
    return needed;
#endif
}

static inline bool cballs_cell_accumulate_child(nodeptr parent, nodeptr child,
                                                 bool read_mask,
                                                 short *cell_mask,
                                                 compute_vector center_of_mass_sum)
{
    bool child_is_valid = !read_mask || Mask(child) != MASK_NODE_MASKED;
    compute_vector weighted_position;

    Selected(parent) |= Selected(child);
    Update(parent) |= Update(child);
    if (read_mask && cell_mask != NULL)
        *cell_mask = mask_node_combine(*cell_mask, Mask(child));

    Mass(parent) += Mass(child);
    if (child_is_valid) {
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
