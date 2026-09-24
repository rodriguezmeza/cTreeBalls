/* Direct ownership and legacy-adapter isolation without estimator mocks. */
#define global
#include "globaldefs.h"
#include <assert.h>

int main(void)
{
    cballs_runtime_state *original = cballs_runtime_current();
    cballs_runtime_state *a = cballs_runtime_create(), *b = cballs_runtime_create();
    real io_a = 3, io_b = 7;
    body cat_a = {0}, cat_b = {0};
    assert(a && b && a != b);
    assert(cballs_runtime_activate(a) == SUCCESS);
    inout_xval = &io_a;
    bodytable[0] = &cat_a;
    a->mpi[0].rank = 3;
    a->tree_workspace.inode = 27;
#ifdef USEGSL
    r_gsl = gsl_rng_alloc(gsl_rng_mt19937);
    assert(r_gsl);
    gsl_rng_set(r_gsl, 17);
    double first = gsl_rng_uniform(r_gsl), second = gsl_rng_uniform(r_gsl);
    gsl_rng_set(r_gsl, 17);
    assert(gsl_rng_uniform(r_gsl) == first);
#endif
    assert(cballs_runtime_activate(b) == SUCCESS);
    assert(inout_xval == NULL && bodytable[0] == NULL && b->mpi[0].rank == 0 && b->tree_workspace.inode == 0);
    inout_xval = &io_b;
    bodytable[0] = &cat_b;
#ifdef USEGSL
    assert(r_gsl == NULL);
    r_gsl = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r_gsl, 91);
#endif
    assert(cballs_runtime_activate(a) == SUCCESS);
    assert(inout_xval == &io_a && bodytable[0] == &cat_a && a->mpi[0].rank == 3 && a->tree_workspace.inode == 27);
#ifdef USEGSL
    assert(gsl_rng_uniform(r_gsl) == second);
    gsl_rng_free(r_gsl); r_gsl = NULL;
#endif
    assert(cballs_runtime_activate(b) == SUCCESS);
    assert(inout_xval == &io_b && bodytable[0] == &cat_b);
#ifdef USEGSL
    gsl_rng_free(r_gsl); r_gsl = NULL;
#endif
    cballs_runtime_destroy(a); /* Destroying an inactive context preserves b. */
    assert(cballs_runtime_current() == b && inout_xval == &io_b);
    cballs_runtime_destroy(b);
    assert(cballs_runtime_current() == original);
    assert(inout_xval == NULL && bodytable[0] == NULL);
    puts("PASS: context-owned RNG/I-O/MPI views and interleaved catalog adapter");
    return 0;
}
