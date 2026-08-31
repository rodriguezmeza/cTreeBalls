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
                     && cballs_opt_smooth_pivot(&first)
                     && cballs_opt_full_sky(&first)
                     && !cballs_opt_smooth(&first)
                     && !cballs_opt_no_check_equal_positions(&first),
                     "first cache lost a present option") == FAILURE)
        return EXIT_FAILURE;
    if (require_true(!cballs_opt_read_mask(&second)
                     && cballs_opt_behavior_ball(&second)
                     && cballs_opt_edge_corrections(&second),
                     "independent caches contaminated each other") == FAILURE)
        return EXIT_FAILURE;

    first.options = "behavior-ball";
    cballs_refresh_option_cache(&first);
    if (require_true(!cballs_opt_read_mask(&first)
                     && !cballs_opt_no_one_ball(&first)
                     && cballs_opt_behavior_ball(&first),
                     "refresh retained stale option bits") == FAILURE)
        return EXIT_FAILURE;

    second.options = "smooth,no-check-two-bodies-eq-pos,ggg-full-window,ggg-profile";
    cballs_refresh_option_cache(&second);
    if (require_true(cballs_opt_smooth(&second)
                     && cballs_opt_no_check_equal_positions(&second)
                     && cballs_opt_ggg_full_window(&second)
                     && cballs_opt_ggg_profile(&second)
                     && !cballs_opt_smooth_pivot(&second),
                     "high cache bits or exact token matching failed") == FAILURE)
        return EXIT_FAILURE;

    puts("PASS: search option cache fallback, refresh, and ownership");
    return EXIT_SUCCESS;
}
