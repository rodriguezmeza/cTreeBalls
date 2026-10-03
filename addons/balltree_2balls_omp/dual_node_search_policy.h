/* Backend selection, scratch ownership, task execution and reductions.
 * Numerical moment primitives live in dual_node_multipole.h. */
#ifdef DUAL_NODE_ADAPTIVE_PAIR_LEAVES
static int dual_node_adaptive_pair_leaf_capacity(
        const struct cmdline_data *cmd, const INTEGER *nbody,
        int cat1, int cat2)
{
    const real catalog_size = (real)MAX(nbody[cat1], nbody[cat2]);
    const real expected_neighbors = MAX(
        1.0, 0.25 * catalog_size * rsqr(cmd->rangeN));
    real relative_bin_width;
    real target;

    if (cmd->useLogHist && cmd->rminHist > 0.0)
        relative_bin_width = rlog(cmd->rangeN / cmd->rminHist)
                           / (real)cmd->sizeHistN;
    else
        relative_bin_width = (cmd->rangeN - cmd->rminHist)
                           / ((real)cmd->sizeHistN * cmd->rangeN);

    /* Calibrated against exact pair histograms: denser catalogs and narrower
     * bins need smaller terminal cells; a larger slop budget permits more
     * bodies per leaf.  Small catalogs rarely accept enough cell pairs to
     * repay the extra build and traversal nodes of four-body leaves. Powers
     * of two keep the tree shape deterministic. */
    target = 16.0 / rsqrt(expected_neighbors)
           * rsqrt(MAX(relative_bin_width, 1.0e-6) / 0.35)
           * rsqrt(MAX(0.125, cmd->theta));
    if (target < 6.0)
        return catalog_size >= (real)262144.0 ? 4 : 8;
    if (target < 12.0) return 8;
    return 16;
}
#endif

static int dual_node_search_log_multipole(
        struct cmdline_data *cmd, struct global_data *gd,
        bodyptr *btab, INTEGER *nbody, INTEGER ipmin, INTEGER *ipmax,
        int cat1, int cat2)
{
    const bool catalog_auto_correlation = cat1 == cat2;
    const bool auto_correlation = catalog_auto_correlation
        && !DUAL_NODE_SEPARATE_PIVOT_TREE(cmd);
    const bool only_2pcf = scanopt(cmd->options, "only-2pcf");
    const bool only_3pcf = scanopt(cmd->options, "only-3pcf");
#ifdef TWOPCF
    const bool run_2pcf = !only_3pcf;
#else
    const bool run_2pcf = FALSE;
#endif
    const bool run_3pcf = !only_2pcf;
    const size_t stride = (size_t)cmd->sizeHistN + 1;
    const int orders = cmd->mChebyshev + 1;
    const INTEGER target_tasks = run_3pcf
        ? dual_node_task_target(cmd, stride, orders)
        : dual_node_pair_frontier_target(cmd);
    int leaf_capacity =
        scanopt(cmd->options, "dual-node-singleton-leaves") ? 1 : cmd->nsmooth;
#ifdef DUAL_NODE_SCALAR_TREE_CACHE
    const bool use_tree_cache = auto_correlation
        && (run_2pcf || run_3pcf)
        && !scanopt(cmd->options, "no-balltree-tree-cache");
    bool tree_cache_hit = FALSE;
#endif
#ifdef DUAL_NODE_ADAPTIVE_PAIR_LEAVES
    if (!run_3pcf
        && !scanopt(cmd->options, "dual-node-singleton-leaves")
        && !scanopt(cmd->options, "dual-node-bucket-leaves"))
        leaf_capacity = dual_node_adaptive_pair_leaf_capacity(
            cmd, nbody, cat1, cat2);
#elif defined(DUAL_NODE_PAIR_LEAF_CAPACITY)
    if (!run_3pcf
        && !scanopt(cmd->options, "dual-node-singleton-leaves")
        && !scanopt(cmd->options, "dual-node-bucket-leaves"))
        leaf_capacity = DUAL_NODE_PAIR_LEAF_CAPACITY;
#endif
    dual_node_search_context context = {0};
    fcfc_balltreeptr tree1 = NULL;
    fcfc_balltreeptr tree2 = NULL;
    INTEGER *frontier = NULL;
    INTEGER *pair_frontier2 = NULL;
    INTEGER task_count = 0;
    INTEGER pair_frontier_count2 = 0;
    real *task_histograms = NULL;
    real *task_scratch = NULL;
    INTEGER *task_body_counts = NULL;
    INTEGER *task_cell_counts = NULL;
    size_t hist_values_per_task = 0;
    size_t scratch_values_per_level = 0;
    size_t scratch_values_per_task = 0;
    int scratch_level_count = 0;
    INTEGER body_total = 0;
    INTEGER cell_total = 0;
    INTEGER pair_test_total = 0;
    INTEGER pivot_restart_total = 0;
    INTEGER pivot_finish_total = 0;
    INTEGER frontier_failure_total = 0;
    INTEGER reuse_pairs=0, reuse_represented_pairs=0, reuse_parent_reductions=0;
    INTEGER distributed_statistics[6] = {0};
    real distributed_profile_statistics[3] = {0.0, 0.0, 0.0};
    double pivot_transport_total = 0.0;
    double scratch_clear_total = 0.0;
    double multipole_product_total = 0.0;
    dual_node_phase_timers phase_timers = {0};
#ifdef DUAL_NODE_PIVOT_PROGRESS
    dual_node_pivot_progress progress = {0};
#endif
    double phase_started = 0.0;
    int operation_status;
    int reduction_status = SUCCESS;
    int status = FAILURE;
    const double cpustart = CPUTIME;

    gd->cpu_edge_correction = 0.0;
    gd->wall_edge_correction = 0.0;

    if (only_2pcf && only_3pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf and only-3pcf are mutually exclusive",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (only_2pcf && !run_2pcf) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: only-2pcf requires TWOPCFON=1",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (ipmin != 1 || ipmax[cat1] != nbody[cat1]) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires the complete pivot catalog",
                 DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }
    if (cmd->nsmooth <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s requires nsmooth > 0", DUAL_NODE_METHOD_NAME);
        return FAILURE;
    }

    if (dual_node_reuse_prepare(&context, cmd, gd, run_3pcf) == FAILURE)
        goto cleanup;
    context.profile = scanopt(cmd->options, "dual-node-profile");
    context.timers = &phase_timers;

    verb_print(cmd->verbose, "Search: Running %s", cmd->searchMethod);
#ifdef TWOPCF
    if (run_2pcf) verb_print(cmd->verbose, " with dual-node 2PCF");
#endif
    if (run_3pcf)
        verb_print(cmd->verbose, " with LogMultipole 3PCF");
    verb_print(cmd->verbose, "\n");
#ifdef DUAL_NODE_DISTRIBUTED_ENGINE
    verb_print(cmd->verbose,
               "%s: %d ranks with deterministic cyclic frontier ownership\n",
               DUAL_NODE_METHOD_NAME, DUAL_NODE_DISTRIBUTED_SIZE());
#endif
    if (cballs_opt_no_two_balls(cmd))
        verb_print(cmd->verbose,
                   "no-two-balls: exact body pivot-neighbor accumulation\n");
    else
        verb_print(cmd->verbose,
                   "two-ball radial/angular node acceptance enabled; theta=%g\n",
                   cmd->theta);
#ifdef TWOPCF
    if (run_2pcf && scanopt(cmd->options, "dual-node-bin-theta"))
        verb_print(cmd->verbose,
                   "dual-node-compatible controlled 2PCF bin theta enabled\n");
#endif
#ifdef DUAL_NODE_NATIVE_BINARY_VIEW
    verb_print(cmd->verbose,
               "native-octree binary view with exact body pivots\n");
#else
    verb_print(cmd->verbose,
               "dual-node tree leaf capacity = %d\n", leaf_capacity);
#endif
#ifdef DUAL_NODE_SCAN_LEVEL_FRONTIER
    verb_print(cmd->verbose,
               "balanced scan-level task frontier enabled by BALLS4SCANLEV\n");
#endif
    if (context.reuse_enabled)
        verb_print(cmd->verbose, "scalar-pivot-reuse: phase_budget=%g bin_theta=%g; strict radial cutoffs\n",
                   (double)gd->scalarReusePhaseBudget, (double)gd->scalarReuseBinTheta);
    if (run_3pcf)
        verb_print(cmd->verbose, "3PCF multipoles: %s\n",
                   dual_node_normalize_3pcf(cmd)
                   ? (cballs_opt_weights_norm(cmd)
                      ? "normalized by distinct-triplet weight sum"
                      : "normalized by distinct-triplet count")
                   : "raw distinct-triplet sums");
    if (context.profile)
        verb_print(cmd->verbose, "native dual-node phase timing enabled\n");

#ifdef OPENMPCODE
    ThreadCount(cmd, gd, nbody[cat1], cat1);
#endif

    operation_status = DUAL_NODE_PREPARE_PIVOTS(
        cmd, gd, btab, nbody, ipmin, ipmax, cat1, cat2);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot preparation") == FAILURE)
        goto cleanup;
    operation_status = search_init_gd_hist_sincos(cmd, gd);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball histogram initialization") == FAILURE)
        goto cleanup;
    if (context.profile) phase_started = dual_node_timer_now();
#ifdef DUAL_NODE_SCALAR_TREE_CACHE
    if (use_tree_cache)
        operation_status = fcfc_balltree_build_scalar_role_cached(
            cmd, gd, btab[cat1], nbody[cat1], leaf_capacity,
            TRUE, &tree1, &tree_cache_hit);
    else
#endif
    operation_status = DUAL_NODE_BUILD_PIVOT_TREE(
        cmd, gd, btab[cat1], nbody[cat1], leaf_capacity, &tree1);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball first tree construction") == FAILURE)
        goto cleanup;
    DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat1, tree1->nnode);
    if (auto_correlation) {
        tree2 = tree1;
    } else {
        operation_status = DUAL_NODE_BUILD_NEIGHBOR_TREE(
            cmd, gd, btab[cat2], nbody[cat2], leaf_capacity, &tree2);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "two-ball second tree construction") == FAILURE)
            goto cleanup;
        DUAL_NODE_PUBLISH_NODE_COUNT(gd, cat2, tree2->nnode);
    }
    if (context.profile)
        phase_timers.build += dual_node_timer_now() - phase_started;
#ifdef DUAL_NODE_SCALAR_TREE_CACHE
    if (use_tree_cache)
        verb_print(cmd->verbose,
                   "%s: compact-tree cache = %s\n",
                   DUAL_NODE_METHOD_NAME,
                   tree_cache_hit ? "hit" : "miss");
#endif

    dual_node_initialize_radial_context(&context, cmd, gd);
    context.use_two_balls = !cballs_opt_no_two_balls(cmd);
    context.use_bin_theta = scanopt(cmd->options, "dual-node-bin-theta");
    context.use_three_cells = context.use_two_balls;
    context.weighted = DUAL_NODE_CONTEXT_WEIGHTED(cmd);
    dual_node_initialize_angular_tolerance(&context);

    if (context.profile) phase_started = dual_node_timer_now();
    operation_status = fcfc_balltree_frontier(
        cmd, tree1, target_tasks, &frontier, &task_count);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball pivot-frontier construction") == FAILURE)
        goto cleanup;
    if (auto_correlation) {
        pair_frontier2 = frontier;
        pair_frontier_count2 = task_count;
#ifdef TWOPCF
    } else if (run_2pcf) {
        operation_status = fcfc_balltree_frontier(
            cmd, tree2, target_tasks,
            &pair_frontier2, &pair_frontier_count2);
        if (dual_node_distributed_consensus(
                cmd, operation_status,
                "two-ball neighbor-frontier construction") == FAILURE)
            goto cleanup;
#endif
    }
    if (context.profile)
        phase_timers.frontier += dual_node_timer_now() - phase_started;

#ifdef TWOPCF
    if (run_2pcf
        && dual_node_run_pair_tasks(
            &context, tree1, tree2, auto_correlation,
            frontier, task_count,
            pair_frontier2, pair_frontier_count2) == FAILURE)
        goto cleanup;
#endif

    if (run_3pcf) {
#ifdef DUAL_NODE_PERSISTENT_PARTIAL_FRONTIER
    scratch_level_count = dual_node_balltree_depth(tree1, 0) + 1;
#elif defined(DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE)
    scratch_level_count = 1;
#else
    scratch_level_count = dual_node_balltree_depth(tree1, 0) + 1;
#endif
    if (context.reuse_enabled) scratch_level_count = 1;
    operation_status = dual_node_allocate_triple_histograms(
        cmd, task_count, stride, orders,
        &task_histograms, &hist_values_per_task,
        &task_body_counts, &task_cell_counts);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball multipole histogram allocation") == FAILURE)
        goto cleanup;
    operation_status = dual_node_allocate_multipole_scratch(
        cmd, task_count, stride, orders, scratch_level_count,
        &task_scratch, &scratch_values_per_level,
        &scratch_values_per_task);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball multipole scratch allocation") == FAILURE)
        goto cleanup;

#ifdef DUAL_NODE_PIVOT_PROGRESS
    progress.enabled = cmd->verbose >= VERBOSEMININFO
        || (cmd->verbose_log >= VERBOSEMININFO && gd->outlog != NULL);
    progress.total = tree1->npoint;
    progress.interval = MAX((INTEGER)1, cmd->stepState);
    progress.batch = MIN((INTEGER)64, progress.interval);
    context.progress = &progress;
    if (progress.enabled) {
        progress.started = dual_node_timer_now();
        dual_node_print_pivot_progress(&context);
    }
#endif
    if (context.profile) phase_started = dual_node_timer_now();
#pragma omp parallel for schedule(dynamic,1) \
    reduction(+:pair_test_total,pivot_restart_total,pivot_finish_total,frontier_failure_total,pivot_transport_total,scratch_clear_total,multipole_product_total,reuse_pairs,reuse_represented_pairs,reuse_parent_reductions)
    for (INTEGER itask = 0; itask < task_count; itask++) {
        real *hist_base = task_histograms
                        + (size_t)itask * hist_values_per_task;
        real *scratch_base = task_scratch
                           + (size_t)itask * scratch_values_per_task;
        dual_node_triple_histogram hist;
        dual_node_multipole_scratch scratch;

        if (!dual_node_distributed_task_owned(itask)) continue;

        dual_node_initialize_triple_histogram(
            &hist, hist_base, stride, orders, dual_node_window_orders(cmd));
        dual_node_initialize_multipole_scratch(
            &scratch, scratch_base, stride, orders,
            dual_node_window_orders(cmd), scratch_values_per_level);
        if (context.reuse_enabled) {
            if (dual_node_reuse_task(&context,tree1,frontier[itask],tree2,
                    auto_correlation,&scratch,&hist) == FAILURE)
                frontier_failure_total++;
        } else {
#ifdef DUAL_NODE_PERSISTENT_PARTIAL_FRONTIER
        if (context.use_two_balls) {
            const INTEGER neighbor_root = 0;
            dual_node_neighbor_frontier *levels = calloc(
                (size_t)tree1->max_depth + 1, sizeof(*levels));

            if (levels == NULL) {
                frontier_failure_total++;
            } else {
                dual_node_multipole_clear_profiled(&context, &scratch);
                if (dual_node_multipole_process_pivots_partial(
                        &context, tree1, frontier[itask], tree2,
                        auto_correlation, &scratch, &scratch,
                        scratch_base, scratch_values_per_level,
                        0, &neighbor_root, 1, levels, cmd->rangeN,
                        cmd->sizeHistN + 1, &hist) == FAILURE)
                    frontier_failure_total++;
                for (int level = 0; level <= tree1->max_depth; level++)
                    free(levels[level].nodes);
                free(levels);
            }
        } else {
            dual_node_multipole_process_body_pivots(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &hist);
        }
#elif defined(DUAL_NODE_BODY_PIVOT_LOG_MULTIPOLE)
#ifdef DUAL_NODE_PERSISTENT_NEIGHBOR_FRONTIER
        if (context.use_two_balls
            && !scanopt(cmd->options, "no-balltree-persistent-frontier")) {
            const INTEGER neighbor_root = 0;
            dual_node_neighbor_frontier *levels = calloc(
                (size_t)tree1->max_depth + 1, sizeof(*levels));

            if (levels != NULL) {
                dual_node_multipole_process_body_pivots_frontier(
                    &context, tree1, frontier[itask], tree2,
                    auto_correlation, &neighbor_root, 1,
                    levels, 0, &scratch, &hist);
                for (int level = 0; level <= tree1->max_depth; level++)
                    free(levels[level].nodes);
                free(levels);
            } else {
                dual_node_multipole_process_body_pivots(
                    &context, tree1, frontier[itask], tree2,
                    auto_correlation, &scratch, &hist);
            }
        } else {
            dual_node_multipole_process_body_pivots(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &hist);
        }
#else
        dual_node_multipole_process_body_pivots(
            &context, tree1, frontier[itask], tree2,
            auto_correlation, &scratch, &hist);
#endif
#else
        if (context.use_two_balls) {
            dual_node_multipole_clear_profiled(&context, &scratch);
            dual_node_multipole_process_pivots_partial(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &scratch,
                scratch_base, scratch_values_per_level,
                0, cmd->rangeN,
                cmd->sizeHistN + 1, &hist);
        } else {
            dual_node_multipole_process_pivots(
                &context, tree1, frontier[itask], tree2,
                auto_correlation, &scratch, &hist);
        }
#endif
        }
        reuse_pairs += scratch.reuse_pairs;
        reuse_represented_pairs += scratch.reuse_represented_pairs;
        reuse_parent_reductions += scratch.reuse_parent_reductions;
#ifdef DUAL_NODE_PIVOT_PROGRESS
        dual_node_publish_pivot_progress(&context, scratch.progress_pending);
#endif
        pair_test_total += scratch.pair_tests;
        pivot_restart_total += scratch.pivot_restarts;
        pivot_finish_total += scratch.pivot_finishes;
        pivot_transport_total += scratch.pivot_transport_seconds;
        scratch_clear_total += scratch.scratch_clear_seconds;
        multipole_product_total += scratch.multipole_product_seconds;
        task_body_counts[itask] = scratch.body_visits;
        task_cell_counts[itask] = scratch.accepted_nodes;
    }
    if (context.profile)
        phase_timers.pair_traversal += dual_node_timer_now() - phase_started;

#ifdef DUAL_NODE_PIVOT_PROGRESS
    if (progress.enabled && progress.completed == progress.total)
        verb_print_min_info(
            cmd->verbose, cmd->verbose_log, gd->outlog,
            "%s: 3PCF pivots complete; reducing histograms and finalizing outputs\n",
            DUAL_NODE_METHOD_NAME);
#endif
    operation_status = frontier_failure_total == 0 ? SUCCESS : FAILURE;
    if (operation_status == FAILURE)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: persistent neighbor-frontier allocation failed",
                 DUAL_NODE_METHOD_NAME);
    if (dual_node_distributed_consensus(
            cmd, operation_status,
            "two-ball persistent neighbor-frontier traversal") == FAILURE)
        goto cleanup;

    distributed_statistics[0] = pair_test_total;
    distributed_statistics[1] = pivot_restart_total;
    distributed_statistics[2] = pivot_finish_total;
    distributed_statistics[3] = reuse_pairs;
    distributed_statistics[4] = reuse_represented_pairs;
    distributed_statistics[5] = reuse_parent_reductions;
    distributed_profile_statistics[0] = (real)pivot_transport_total;
    distributed_profile_statistics[1] = (real)scratch_clear_total;
    distributed_profile_statistics[2] = (real)multipole_product_total;
    if (context.profile) phase_started = dual_node_timer_now();
    if (dual_node_distributed_reduce_reals(
            cmd, task_histograms,
            (size_t)task_count * hist_values_per_task) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_body_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, task_cell_counts, (size_t)task_count) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_reduce_integers(
            cmd, distributed_statistics, 6) == FAILURE)
        reduction_status = FAILURE;
    if (context.profile
        && dual_node_distributed_reduce_reals(
            cmd, distributed_profile_statistics, 3) == FAILURE)
        reduction_status = FAILURE;
    if (dual_node_distributed_consensus(
            cmd, reduction_status,
            "two-ball multipole reduction") == FAILURE)
        goto cleanup;

    if (dual_node_distributed_publish()) {
    pair_test_total = distributed_statistics[0];
    pivot_restart_total = distributed_statistics[1];
    pivot_finish_total = distributed_statistics[2];
    gd->scalarReusePairs = distributed_statistics[3];
    gd->scalarReuseRepresentedPairs = distributed_statistics[4];
    gd->scalarReuseParentReductions = distributed_statistics[5];
    if (context.reuse_enabled && context.profile)
        verb_print(cmd->verbose,"scalar-pivot-reuse: radial_pairs=%" INTEGER_FMT
                   " represented_pairs=%" INTEGER_FMT " parent_reductions=%" INTEGER_FMT "\n",
                   gd->scalarReusePairs,gd->scalarReuseRepresentedPairs,gd->scalarReuseParentReductions);
    if (context.profile) {
        phase_timers.pivot_transport_thread =
            (double)distributed_profile_statistics[0];
        phase_timers.scratch_clear_thread =
            (double)distributed_profile_statistics[1];
        phase_timers.multipole_products_thread =
            (double)distributed_profile_statistics[2];
    }
    for (INTEGER itask = 0; itask < task_count; itask++) {
        const real *base = task_histograms
                         + (size_t)itask * hist_values_per_task;
        const size_t plane = stride * stride;

        for (int order = 0; order < orders; order++) {
            const real *zcos = base
                + ((size_t)DUAL_NODE_ZETA_COS * (size_t)orders
                   + (size_t)order) * plane;
            const real *zsin = base
                + ((size_t)DUAL_NODE_ZETA_SIN * (size_t)orders
                   + (size_t)order) * plane;
            const real *zsincos = base
                + ((size_t)DUAL_NODE_ZETA_SINCOS * (size_t)orders
                   + (size_t)order) * plane;
            const real *zcossin = base
                + ((size_t)DUAL_NODE_ZETA_COSSIN * (size_t)orders
                   + (size_t)order) * plane;

            for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
                for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
                    const size_t index = (size_t)n1 * stride + (size_t)n2;
                    gd->histZetaMcos[order + 1][n1][n2] += zcos[index];
                    gd->histZetaMsin[order + 1][n1][n2] += zsin[index];
                    gd->histZetaMsincos[order + 1][n1][n2] += zsincos[index];
                    gd->histZetaMcossin[order + 1][n1][n2] += zcossin[index];
                }
        }
        body_total += task_body_counts[itask];
        cell_total += task_cell_counts[itask];
    }

    }
    operation_status = dual_node_distributed_publish()
        ? dual_node_publish_edge(cmd, gd, task_histograms, task_count,
                                hist_values_per_task, stride, orders)
        : SUCCESS;
    if (dual_node_distributed_consensus(
            cmd, operation_status, "two-ball edge correction") == FAILURE)
        goto cleanup;
    if (dual_node_distributed_publish()) {
    if (dual_node_normalize_3pcf(cmd)) {
        for (int n1 = 1; n1 <= cmd->sizeHistN; n1++)
            for (int n2 = 1; n2 <= cmd->sizeHistN; n2++) {
            const size_t index = (size_t)n1 * stride + (size_t)n2;
            real denominator = 0.0;

            for (INTEGER itask = 0; itask < task_count; itask++) {
                const real *base = task_histograms
                    + (size_t)itask * hist_values_per_task;
                denominator += base[
                    DUAL_NODE_ZETA_COMPONENTS
                    * (size_t)orders * stride * stride + index];
            }
            for (int order = 1; order <= orders; order++) {
                gd->histZetaMcos[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcos[order][n1][n2], denominator);
                gd->histZetaMsin[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsin[order][n1][n2], denominator);
                gd->histZetaMsincos[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMsincos[order][n1][n2], denominator);
                gd->histZetaMcossin[order][n1][n2] = cballs_normalize_or_zero(
                    gd->histZetaMcossin[order][n1][n2], denominator);
            }
        }
    }
    }
    }
    if (context.profile && run_3pcf)
        phase_timers.reduction += dual_node_timer_now() - phase_started;

    gd->cpusearch = CPUTIME - cpustart;
#ifdef TWOPCF
    if (run_2pcf)
        verb_print(cmd->verbose,
                   "%s: nbbcalc = %" INTEGER_FMT
                   ", nbccalc = %" INTEGER_FMT
                   ", frontier tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, gd->nbbcalc, gd->nbccalc, task_count);
#endif
    if (run_3pcf) {
        verb_print(cmd->verbose,
                   "%s: neighbor body visits = %" INTEGER_FMT
                   ", accepted cell orientations = %" INTEGER_FMT
                   ", pivot tasks = %" INTEGER_FMT "\n",
                   cmd->searchMethod, body_total, cell_total, task_count);
        verb_print(cmd->verbose,
                   "%s: pair tests = %" INTEGER_FMT
                   ", restarted pivot scans = %" INTEGER_FMT
                   ", finished pivots = %" INTEGER_FMT "\n",
                   cmd->searchMethod, pair_test_total,
                   pivot_restart_total, pivot_finish_total);
    }
    if (context.profile && dual_node_distributed_publish())
        verb_print(TRUE,
                   "%s phase-timers: build_wall=%.9g, frontier_wall=%.9g, "
                   "pair_traversal_wall=%.9g, pivot_transport_thread=%.9g, "
                   "scratch_clear_thread=%.9g, multipole_products_thread=%.9g, "
                   "reduction_wall=%.9g\n",
                   DUAL_NODE_METHOD_NAME, phase_timers.build,
                   phase_timers.frontier, phase_timers.pair_traversal,
                   phase_timers.pivot_transport_thread,
                   phase_timers.scratch_clear_thread,
                   phase_timers.multipole_products_thread,
                   phase_timers.reduction);
    verb_print(cmd->verbose, "Going out: CPU time = %lf\n", gd->cpusearch);
    status = SUCCESS;

cleanup:
    free(task_scratch);
    free(task_cell_counts);
    free(task_body_counts);
    free(task_histograms);
    if (!auto_correlation) free(pair_frontier2);
    free(frontier);
    if (tree2 != tree1) DUAL_NODE_RELEASE_TREE(tree2);
    DUAL_NODE_RELEASE_TREE(tree1);
    return status;
}
