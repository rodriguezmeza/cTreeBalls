/* Engine lookup and help share the generated capability table. */
#include "globaldefs.h"
#include "engine_registry_generated.h"

global int cballs_search_method_id(const char *method)
{
    size_t i;
    if (method == NULL) return -1;
    if (strnull(method)) return SEARCHNULL;
    for (i = 0; i < sizeof(cballs_engine_registry)/sizeof(cballs_engine_registry[0]); ++i)
        if (strcmp(method, cballs_engine_registry[i].name) == 0)
            return cballs_engine_registry[i].id;
    for (i = 0; cballs_engine_aliases[i].name != NULL; ++i)
        if (strcmp(method, cballs_engine_aliases[i].name) == 0)
            return cballs_engine_aliases[i].id;
    return -1;
}

global int cballs_print_search_methods(struct cmdline_data* cmd,
                               struct global_data* gd)
{
    const size_t method_count =
        sizeof(cballs_engine_registry) / sizeof(cballs_engine_registry[0]);
    size_t available_count = 0;
    size_t i;
    int method_id;

    (void)gd;

    for (i = 0; i < method_count; ++i) {
        method_id = cballs_engine_registry[i].id;
        if (method_id >= 0) ++available_count;
    }

    verb_print_zero(cmd->verbose,
                    "\nSearching methods registered in this executable (%zu):\n",
                    available_count);
    verb_print_zero(cmd->verbose,
                    "Select one with search=<name> or searchMethod=<name>.\n");
    verb_print_zero(cmd->verbose,
                    "Common controls: rangeN, rminHist, sizeHistN, useLogHist, "
                    "theta, numberThreads, and options.\n");
    verb_print_zero(cmd->verbose,
                    "Only methods enabled by this executable's build profile are shown.\n");

    for (i = 0; i < method_count; ++i) {
        method_id = cballs_engine_registry[i].id;
        if (method_id < 0) continue;
        verb_print_zero(cmd->verbose, "\n- %s (id=%d)\n",
                        cballs_engine_registry[i].name, method_id);
        verb_print_zero(cmd->verbose, "  geometry: %s\n",
                        cballs_engine_registry[i].geometry);
        verb_print_zero(cmd->verbose, "  computes: %s\n",
                        cballs_engine_registry[i].correlations);
        verb_print_zero(cmd->verbose, "  use: %s\n",
                        cballs_engine_registry[i].usage);
        if (strcmp(cballs_engine_registry[i].name, "octree-2balls-omp") == 0
            || strcmp(cballs_engine_registry[i].name, "octree-2balls-mpi") == 0) {
#ifdef SMOOTHPIVOT
            verb_print_zero(cmd->verbose,
                "  smooth-pivot: default-on only with options=legacy-one-ball; add no-smooth-pivot to disable it\n");
#else
            verb_print_zero(cmd->verbose,
                "  smooth-pivot: legacy-one-ball supports it, but SMOOTHPIVOTON=1 is not compiled\n");
#endif
        } else if (cballs_method_supports_smooth_pivot(cballs_engine_registry[i].name)) {
#ifdef SMOOTHPIVOT
            verb_print_zero(cmd->verbose,
                "  smooth-pivot: default-on; add options=no-smooth-pivot to disable\n");
#else
            verb_print_zero(cmd->verbose,
                "  smooth-pivot: supported but not compiled; set SMOOTHPIVOTON=1 and rebuild\n");
#endif
        } else {
            verb_print_zero(cmd->verbose,
                "  smooth-pivot: unsupported; SMOOTHPIVOTON does not change this engine\n");
        }
        if (strcmp(cballs_engine_registry[i].name,
                   "lya-1d-tree-same-los-2pcf-omp") == 0)
            verb_print_zero(cmd->verbose,
                "  input: one x y z delta weight forest_id catalog (lya-ascii), or "
                "cyballs.set_forest_catalog(positions, delta, weights, forest_ids). "
                "Only pairs within the same forest are accepted. Each occupied "
                "forest/bin is normalized first, then forests are averaged equally. "
                "DEFDIMENSION=3 and usePeriodic=false are required.\n");
        else if (strncmp(cballs_engine_registry[i].name, "lya-", 4) == 0)
            verb_print_zero(cmd->verbose,
                "  input: one x y z delta weight forest_id catalog (lya-ascii), or "
                "cyballs.set_forest_catalog(positions, delta, weights, forest_ids). "
                "DEFDIMENSION=3, observer-centered comoving coordinates and "
                "usePeriodic=false are required even for radial searches. "
                "Pairs exclude the same quasar; triplets require three distinct "
                "quasars. Histograms are weight-normalized; empty bins are zero. "
                "Use tests/python/lya_corr_all_engines.py for DESI/eBOSS FITS, NPZ or ASCII, "
                "one-time loading and MPI broadcasting. No smooth-pivot.\n");
        if (strncmp(cballs_engine_registry[i].name, "octree-3pcf-3d-", 15) == 0)
            verb_print_zero(cmd->verbose,
                "  modes: only-2pcf-3d, only-3pcf-3d, or "
                "compute-2pcf-3d,compute-3pcf-3d. survey-estimator-3d "
                "uses data/random catalogs and window correction. "
                "exclude-same-los excludes pivot LOS matches only. "
                "exclude-all-same-los (alias lya-distinct-forests) requires "
                "three distinct LOS/forest IDs in every triplet.\n");
        if (strstr(cballs_engine_registry[i].name, "2balls") != NULL
            && strstr(cballs_engine_registry[i].name, "shear") == NULL)
            verb_print_zero(cmd->verbose,
                "  edge correction: add edge-corrections,no-normalize-HistZeta "
                "for complex scalar 3PCF window deconvolution. Window modes "
                "extend through 2*mChebyshev. weights-norm weights both signal "
                "and window; unsupported corrected bins are NaN; window diagnostics distinguish empty, singular, and valid bins. 2PCF is unchanged.\n");
    }

    verb_print_zero(cmd->verbose,
                    "\nScalar angular engines preserve the observer frame in 3D and use "
                    "tangent-plane angles with Euclidean chord bins. Coincident/radial/antipodal "
                    "legs are excluded from angular multipoles, not from ordinary pair counts.\n"
                    "tests/python/kappa_corr_all_engines.py reuses one in-memory catalog across active native "
                    "engines and writes ordinary and "
                    "flattened radial-bin 3PCF plots. tests/python/shear_corr_all_engines.py provides "
                    "active full-sky spin-2 comparisons. All three drivers record native MainLoop "
                    "and complete Python-call timings separately. tests/python/lya_corr_all_engines.py provides "
                    "the corresponding forest workflow and rejects spin-2 engines explicitly.\n"
                    "Use options=make-info to inspect the build profile and "
                    "options=print-options for the full option list.\n\n");
    return SUCCESS;
}
