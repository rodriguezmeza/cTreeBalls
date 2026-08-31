#ifndef _options_cache_h
#define _options_cache_h

/*
 * Search options are immutable while a run is active. Cache the options used
 * by production search engines so recursive tree walks do not repeatedly scan
 * the comma-separated command string. The fallback keeps direct unit calls
 * correct when they construct cmdline_data without running startup first.
 */
#define CBALLS_CACHED_OPTION_TABLE(X) \
    X(AND_CF, and_cf, "and-CF", 0) \
    X(ASYMMETRIC, asymmetric, "asymmetric", 1) \
    X(BEHAVIOR_BALL, behavior_ball, "behavior-ball", 2) \
    X(BEHAVIOR_TREE_OMP, behavior_tree_omp, "behavior-tree-omp", 3) \
    X(BH86, bh86, "bh86", 4) \
    X(CELESTIAL, celestial, "celestial", 5) \
    X(CENTER_OF_MASS, center_of_mass, "center-of-mass", 6) \
    X(COMPUTE_HISTN, compute_histn, "compute-HistN", 7) \
    X(COMPUTE_J_NO_EQ_I, compute_j_no_eq_i, "compute-j-no-eq-i", 8) \
    X(CUTE_BOX, cute_box, "cute-box", 9) \
    X(CUTE_BOX_FMT, cute_box_fmt, "cute-box-fmt", 10) \
    X(CUTE_BOX_RMIN, cute_box_rmin, "cute-box-rmin", 11) \
    X(DEFAULT_RSMOOTH, default_rsmooth, "default-rsmooth", 12) \
    X(EDGE_CORRECTIONS, edge_corrections, "edge-corrections", 13) \
    X(EDGE_EFFECTS, edge_effects, "edge-effects", 14) \
    X(FIX_RSMOOTH, fix_rsmooth, "fix-rsmooth", 15) \
    X(GGG_CORRELATION, ggg_correlation, "GGGCorrelation", 16) \
    X(KAPPA_CONSTANT, kappa_constant, "kappa-constant", 17) \
    X(KAPPA_CONSTANT_ONE, kappa_constant_one, "kappa-constant-one", 18) \
    X(KKK_CORRELATION, kkk_correlation, "KKKCorrelation", 19) \
    X(LYA_OUTPUT_EMPTY_BINS, lya_output_empty_bins, "lya-output-empty-bins", 20) \
    X(NN_LANDY_SZALAY1, nn_landy_szalay1, "NNLandySzalay1", 21) \
    X(NN_LANDY_SZALAY2, nn_landy_szalay2, "NNLandySzalay2", 22) \
    X(NN_STANDARD, nn_standard, "NNStandard", 23) \
    X(NO_NORMALIZE_HISTZETA, no_normalize_histzeta, "no-normalize-HistZeta", 24) \
    X(NO_ONE_BALL, no_one_ball, "no-one-ball", 25) \
    X(NO_OUT_HIST, no_out_hist, "no-out-Hist", 26) \
    X(NO_TWO_BALL, no_two_ball, "no-two-ball", 27) \
    X(NO_TWO_BALLS, no_two_balls, "no-two-balls", 28) \
    X(ONLY_POS, only_pos, "only-pos", 29) \
    X(OUT_HISTZETAG, out_histzetag, "out-HistZetaG", 30) \
    X(OUT_M_HISTZETA, out_m_histzeta, "out-m-HistZeta", 31) \
    X(PATCH_WITH_ALL, patch_with_all, "patch-with-all", 32) \
    X(PIVOT_LOOP, pivot_loop, "pivot-loop", 33) \
    X(PIVOT_NUMBER, pivot_number, "pivot-number", 34) \
    X(RA_REVERSED, ra_reversed, "ra-reversed", 35) \
    X(RBIN_ARCMIN, rbin_arcmin, "rbin-arcmin", 36) \
    X(RBIN_DEGREE, rbin_degree, "rbin-degree", 37) \
    X(READ_MASK, read_mask, "read-mask", 38) \
    X(SMOOTH_PIVOT, smooth_pivot, "smooth-pivot", 39) \
    X(SW94, sw94, "sw94", 40) \
    X(WEIGHTS_NORM, weights_norm, "weights-norm", 41) \
    X(ARFKEN, arfken, "arfken", 42) \
    X(NO_ARFKEN, no_arfken, "no-arfken", 43) \
    X(GALACTIC, galactic, "galactic", 44) \
    X(ECLIPTIC, ecliptic, "ecliptic", 45) \
    X(SAME_INFILES, same_infiles, "same-infiles", 46) \
    X(FULL_SKY, full_sky, "full-sky", 47) \
    X(SMOOTH, smooth, "smooth", 48) \
    X(SET_NB_NOSEL, set_nb_nosel, "set-Nb-noSel", 49) \
    X(SMOOTH_MIN_CELL, smooth_min_cell, "smooth-min-cell", 50) \
    X(SET_DEFAULT_PARAM, set_default_param, "set-default-param", 51) \
    X(REMOVE_MEAN, remove_mean, "remove-mean", 52) \
    X(NO_CHECK_EQUAL_POSITIONS, no_check_equal_positions, \
      "no-check-two-bodies-eq-pos", 53) \
    X(ONLY_2PCF, only_2pcf, "only-2pcf", 54) \
    X(GGG_FULL_WINDOW, ggg_full_window, "ggg-full-window", 55) \
    X(GGG_PROFILE, ggg_profile, "ggg-profile", 56)

#define CBALLS_OPTION_ENUM(symbol, accessor, text, bit) \
    CBALLS_OPTF_##symbol = 1ULL << (bit),
enum cballs_cached_option_flag {
    CBALLS_CACHED_OPTION_TABLE(CBALLS_OPTION_ENUM)
};
#undef CBALLS_OPTION_ENUM

static inline bool cballs_option_cached(const struct cmdline_data *cmd,
                                        unsigned long long flag,
                                        const char *option)
{
    if (cmd == NULL || cmd->options == NULL)
        return FALSE;
    if (cmd->options_cache_valid)
        return (cmd->options_cache_flags & flag) != 0;
    return scanopt(cmd->options, (string)option);
}

static inline void cballs_refresh_option_cache(struct cmdline_data *cmd)
{
    unsigned long long flags = 0;

    if (cmd == NULL)
        return;
    if (cmd->options != NULL) {
#define CBALLS_CACHE_OPTION(symbol, accessor, text, bit) \
        if (scanopt(cmd->options, text)) flags |= CBALLS_OPTF_##symbol;
        CBALLS_CACHED_OPTION_TABLE(CBALLS_CACHE_OPTION)
#undef CBALLS_CACHE_OPTION
    }
    cmd->options_cache_flags = flags;
    cmd->options_cache_valid = TRUE;
}

#define CBALLS_OPTION_ACCESSOR(symbol, accessor, text, bit) \
    static inline bool cballs_opt_##accessor(const struct cmdline_data *cmd) \
    { \
        return cballs_option_cached(cmd, CBALLS_OPTF_##symbol, text); \
    }
CBALLS_CACHED_OPTION_TABLE(CBALLS_OPTION_ACCESSOR)
#undef CBALLS_OPTION_ACCESSOR
#undef CBALLS_CACHED_OPTION_TABLE

#endif
