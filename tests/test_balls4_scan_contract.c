#include <assert.h>

#include "globaldefs.h"
#include "tree_contracts.h"

int main(void)
{
#ifdef BALLS4SCANLEV
    assert(cballs_method_needs_balls4_scan(66));
#elif defined(OCTREEKKKBALLS4OMP)
    assert(cballs_method_needs_balls4_scan(OCTREEKKKBALLS4OMPMETHOD));
    assert(!cballs_method_needs_balls4_scan(66));
#else
    assert(!cballs_method_needs_balls4_scan(66));
#endif

    return 0;
}
