#include <assert.h>

#include "globaldefs.h"
#include "tree_contracts.h"

int main(void)
{
#ifdef BALLS4SCANLEV
    assert(cballs_method_needs_balls4_scan(66));
#else
    assert(!cballs_method_needs_balls4_scan(66));
#endif
#ifdef OCTREEBALLS4OMP
    assert(cballs_method_needs_balls4_scan(OCTREEBALLS4OMPMETHOD));
#endif
#ifdef OCTREEBALLS4MPI
    assert(cballs_method_needs_balls4_scan(OCTREEBALLS4MPIMETHOD));
#endif

    return 0;
}
