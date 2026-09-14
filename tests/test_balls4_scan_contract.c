#include <assert.h>

#include "globaldefs.h"
#include "tree_contracts.h"

bool scanopt(string options, string key)
{
    (void)options;
    (void)key;
    assert(!"option-cache contract test unexpectedly reparsed options");
    return false;
}

int main(void)
{
#ifdef BALLS4SCANLEV
#ifdef OCTREE3PCF3DOMP
    assert(cballs_method_needs_balls4_scan(OCTREE3PCF3DOMPMETHOD));
#endif
#ifdef OCTREE3PCF3DMPI
    assert(cballs_method_needs_balls4_scan(OCTREE3PCF3DMPIMETHOD));
#endif
#ifdef OCTREESHEARSPHERE2BALLSOMP
    {
        struct cmdline_data cmd = {0};

        assert(cballs_method_needs_balls4_scan(
            OCTREESHEARSPHERE2BALLSOMPMETHOD));
        cmd.searchMethod = "octree-shear-sphere-2balls-omp";
        cmd.options = "";
        cmd.options_cache_valid = TRUE;
        cmd.options_cache_flags = CBALLS_OPTF_ONLY_2PCF
                                | CBALLS_OPTF_NO_SMOOTH_PIVOT;
        assert(!cballs_run_needs_balls4_scan(
            &cmd, OCTREESHEARSPHERE2BALLSOMPMETHOD));
        cmd.options_cache_flags = CBALLS_OPTF_ONLY_2PCF;
#ifdef SMOOTHPIVOT
        assert(cballs_run_needs_balls4_scan(
            &cmd, OCTREESHEARSPHERE2BALLSOMPMETHOD));
#endif
    }
#endif
#else
    assert(!cballs_method_needs_balls4_scan(66));
#endif
#ifdef KDTREE2BALLSOMP
    assert(!cballs_method_needs_balls4_scan(KDTREE2BALLSOMPMETHOD));
#endif
#ifdef KDTREE2BALLSMPI
    assert(!cballs_method_needs_balls4_scan(KDTREE2BALLSMPIMETHOD));
#endif
    return 0;
}
