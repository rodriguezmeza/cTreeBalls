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
#ifdef OCTREEGGGOMP
    assert(cballs_method_needs_balls4_scan(OCTREEGGGOMPMETHOD));
#endif
#ifdef OCTREEGGGMPI
    assert(cballs_method_needs_balls4_scan(OCTREEGGGMPIMETHOD));
#endif
#ifdef OCTREE3PCF3DOMP
    assert(cballs_method_needs_balls4_scan(OCTREE3PCF3DOMPMETHOD));
#endif
#ifdef OCTREE3PCF3DMPI
    assert(cballs_method_needs_balls4_scan(OCTREE3PCF3DMPIMETHOD));
#endif
#ifdef OCTREESHEAROMP
    assert(cballs_method_needs_balls4_scan(OCTREESHEARMETHOD));
#endif
#ifdef OCTREESHEARSPHEREOMP
    {
        struct cmdline_data cmd = {0};

        assert(cballs_method_needs_balls4_scan(OCTREESHEARSPHEREMETHOD));
        cmd.searchMethod = "octree-shear-sphere-omp";
        cmd.theta = 1.0;
        cmd.options = "";
        cmd.options_cache_valid = TRUE;
        cmd.options_cache_flags = CBALLS_OPTF_ONLY_2PCF;
#ifdef SMOOTHPIVOT
        assert(cballs_run_needs_balls4_scan(
            &cmd, OCTREESHEARSPHEREMETHOD));
        cmd.options_cache_flags = CBALLS_OPTF_ONLY_2PCF
                                | CBALLS_OPTF_NO_SMOOTH_PIVOT;
#endif
        assert(!cballs_run_needs_balls4_scan(
            &cmd, OCTREESHEARSPHEREMETHOD));
        cmd.options_cache_flags = CBALLS_OPTF_ONLY_2PCF
                                | CBALLS_OPTF_NO_ONE_BALL;
        assert(cballs_run_needs_balls4_scan(
            &cmd, OCTREESHEARSPHEREMETHOD));
    }
#endif
#else
    assert(!cballs_method_needs_balls4_scan(66));
#endif
#ifdef KDTREEOMP
    assert(!cballs_method_needs_balls4_scan(KDTREEOMPMETHOD));
#endif
#ifdef KDTREEMPI
    assert(!cballs_method_needs_balls4_scan(KDTREEMPIMETHOD));
#endif
#ifdef KDTREE2BALLSOMP
    assert(!cballs_method_needs_balls4_scan(KDTREE2BALLSOMPMETHOD));
#endif
#ifdef KDTREE2BALLSMPI
    assert(!cballs_method_needs_balls4_scan(KDTREE2BALLSMPIMETHOD));
#endif
#ifdef OCTREEBALLS4OMP
    assert(cballs_method_needs_balls4_scan(OCTREEBALLS4OMPMETHOD));
#endif
#ifdef OCTREEBALLS4MPI
    assert(cballs_method_needs_balls4_scan(OCTREEBALLS4MPIMETHOD));
#endif

    return 0;
}
