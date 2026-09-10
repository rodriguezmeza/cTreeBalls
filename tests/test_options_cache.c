#include "globaldefs.h"

static int require_true(bool condition, const char *message)
{
    if (condition)
        return SUCCESS;
    fprintf(stderr, "FAIL: %s\n", message);
    return FAILURE;
}

int main(void)
{
    struct cmdline_data first;
    struct cmdline_data second;

    memset(&first, 0, sizeof(first));
    memset(&second, 0, sizeof(second));

    first.searchMethod = "octree-ggg-omp";
    second.searchMethod = "octree-2balls-omp";
    first.options = "read-mask,no-one-ball,smooth-pivot,full-sky";
    second.options = "behavior-ball,edge-corrections";

    if (require_true(cballs_opt_read_mask(&first),
                     "uncached fallback did not find read-mask") == FAILURE)
        return EXIT_FAILURE;
    if (require_true(!cballs_opt_behavior_ball(&first),
                     "uncached fallback reported an absent option") == FAILURE)
        return EXIT_FAILURE;

    cballs_refresh_option_cache(&first);
    cballs_refresh_option_cache(&second);
    if (require_true(cballs_opt_read_mask(&first)
                     && cballs_opt_no_one_ball(&first)
                     && cballs_opt_smooth_pivot_requested(&first)
                     && cballs_opt_full_sky(&first)
                     && !cballs_opt_smooth(&first)
                     && !cballs_opt_no_check_equal_positions(&first),
                     "first cache lost a present option") == FAILURE)
        return EXIT_FAILURE;
#ifdef SMOOTHPIVOT
    if (require_true(cballs_opt_smooth_pivot(&first),
                     "compiled supported method did not enable smoothing") == FAILURE)
        return EXIT_FAILURE;
#else
    if (require_true(!cballs_opt_smooth_pivot(&first),
                     "uncompiled smoothing became effective") == FAILURE)
        return EXIT_FAILURE;
#endif
    if (require_true(!cballs_opt_read_mask(&second)
                     && cballs_opt_behavior_ball(&second)
                     && cballs_opt_edge_corrections(&second)
                     && !cballs_opt_smooth_pivot(&second),
                     "independent caches contaminated each other") == FAILURE)
        return EXIT_FAILURE;

    first.options = "behavior-ball,no-smooth-pivot";
    cballs_refresh_option_cache(&first);
    if (require_true(!cballs_opt_read_mask(&first)
                     && !cballs_opt_no_one_ball(&first)
                     && cballs_opt_behavior_ball(&first)
                     && cballs_opt_no_smooth_pivot(&first)
                     && !cballs_opt_smooth_pivot(&first),
                     "refresh retained stale option bits") == FAILURE)
        return EXIT_FAILURE;

    first.options = "";
    cballs_refresh_option_cache(&first);
#ifdef SMOOTHPIVOT
    if (require_true(cballs_opt_smooth_pivot(&first),
                     "supported method did not default to smoothing") == FAILURE)
        return EXIT_FAILURE;
#else
    if (require_true(!cballs_opt_smooth_pivot(&first),
                     "compiled-off method defaulted to smoothing") == FAILURE)
        return EXIT_FAILURE;
#endif

    second.options = "smooth,no-check-two-bodies-eq-pos,ggg-full-window,ggg-profile,legacy-one-ball";
    cballs_refresh_option_cache(&second);
    if (require_true(cballs_opt_smooth(&second)
                     && cballs_opt_no_check_equal_positions(&second)
                     && cballs_opt_ggg_full_window(&second)
                     && cballs_opt_ggg_profile(&second)
                     && cballs_opt_legacy_one_ball(&second)
                     && !cballs_opt_smooth_pivot_requested(&second)
                     && !cballs_opt_smooth_pivot(&second),
                     "high cache bits or exact token matching failed") == FAILURE)
        return EXIT_FAILURE;

    puts("PASS: search option cache fallback, refresh, and ownership");
    return EXIT_SUCCESS;
}
