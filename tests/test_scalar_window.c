/* Analytic windows and allocation failures, independent of any tree traversal. */
#define global
#include "globaldefs.h"
#include <unistd.h>

#define DUAL_NODE_ZETA_COMPONENTS 4
#define DUAL_NODE_METHOD_NAME "scalar-window-test"
static int dual_node_window_orders(struct cmdline_data *cmd)
{ return 2*cmd->mChebyshev+1; }
#include "dual_node_edge_correction.h"

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #c); return 1; \
} } while (0)

static int analytic_solver(void)
{
    double ratio;
    double complex identity[9] = {1,0,0,0,1,0,0,0,1}, zero[3] = {0};
    CHECK(dual_node_edge_solve(identity, zero, 3, &ratio) == CBALLS_WINDOW_VALID);
    CHECK(ratio == 1 && zero[0] == 0 && zero[1] == 0 && zero[2] == 0);
    double complex singular[9] = {1,1,1,1,1,1,1,1,1};
    CHECK(dual_node_edge_solve(singular, zero, 3, &ratio) == CBALLS_WINDOW_SINGULAR);
    CHECK(ratio == 0);
    double complex narrow[9] = {1,0,0,0,1,0,0,0,1e-16};
    CHECK(dual_node_edge_solve(narrow, zero, 3, &ratio) == CBALLS_WINDOW_SINGULAR);
    CHECK(ratio == 0); /* The original pivot threshold must not be loosened. */
    double complex nonfinite[1] = {NAN};
    CHECK(dual_node_edge_solve(nonfinite, zero, 1, &ratio) == CBALLS_WINDOW_NONFINITE);
    CHECK(isnan(ratio));
    double complex matrix[9] = {2,I,0,-I,3,1,0,1,4};
    const double complex expected[3] = {1+2*I, -2+I, .5-I};
    double complex rhs[3] = {0};
    for (int i=0; i<3; i++) for (int j=0; j<3; j++) rhs[i] += matrix[3*i+j]*expected[j];
    CHECK(dual_node_edge_solve(matrix, rhs, 3, &ratio) == CBALLS_WINDOW_VALID);
    CHECK(ratio > 0 && ratio <= 1);
    for (int i=0; i<3; i++) CHECK(cabs(rhs[i]-expected[i]) < 1e-14);
    puts("PASS: valid zero, complex solution, singular, near-singular, nonfinite solver");
    return 0;
}

static real ***histogram(void)
{
    real ***p = calloc(3, sizeof(*p));
    for (int m=0; m<3; m++) {
        p[m] = calloc(3, sizeof(*p[m]));
        for (int i=0; i<3; i++) p[m][i] = calloc(3, sizeof(*p[m][i]));
    }
    return p;
}

static void free_histogram(real ***p)
{
    for (int m=0; m<3; m++) {
        for (int i=0; i<3; i++) free(p[m][i]);
        free(p[m]);
    }
    free(p);
}

static int publication(void)
{
    struct cmdline_data cmd = {0};
    struct global_data gd = {0};
    cmd.sizeHistN = 2; cmd.mChebyshev = 1;
    cmd.options = "edge-corrections,no-out-Hist";
    for (int fail=0; fail<3; fail++) {
        cballs_test_fail_allocation_after(fail);
        CHECK(cballs_scalar_window_begin(&cmd, &gd) == FAILURE);
        CHECK(!gd.scalar_window_status && !gd.scalar_window_w0
              && !gd.scalar_window_pivot_ratio && !gd.scalar_window_ready);
        cballs_test_reset_allocation_failure();
    }
    real ***histograms[6];
    for (int i=0; i<6; i++) histograms[i] = histogram();
    gd.histZetaMcos = histograms[0]; gd.histZetaMsin = histograms[1];
    gd.histZetaMsincos = histograms[2]; gd.histZetaMcossin = histograms[3];
    gd.histZetaM_EE = histograms[4]; gd.histZetaM_EE_Im = histograms[5];
    real tasks[135] = {0};
    real *window = tasks + 81;
    window[4] = 4; /* Pure W0: identity system, supported zero field. */
    for (int m=0; m<3; m++) window[9*m+5] = window[9*m+8] = 1;
    for (int repeat=0; repeat<3; repeat++) {
        CHECK(dual_node_publish_edge(&cmd, &gd, tasks, 1, 135, 3, 2) == SUCCESS);
        CHECK(gd.scalar_window_ready && gd.scalar_window_bins == 2);
        CHECK(gd.scalar_window_status[0] == CBALLS_WINDOW_VALID);
        CHECK(gd.scalar_window_status[1] == CBALLS_WINDOW_SINGULAR);
        CHECK(gd.scalar_window_status[2] == CBALLS_WINDOW_EMPTY);
        CHECK(gd.scalar_window_status[3] == CBALLS_WINDOW_SINGULAR);
        CHECK(gd.scalar_window_w0[0] == 4 && gd.scalar_window_w0[2] == 0);
        for (int m=1; m<=2; m++) {
            CHECK(gd.histZetaM_EE[m][1][1] == 0 && gd.histZetaM_EE_Im[m][1][1] == 0);
            CHECK(isnan(gd.histZetaM_EE[m][1][2]) && isnan(gd.histZetaM_EE_Im[m][1][2]));
            CHECK(isnan(gd.histZetaM_EE[m][2][1]) && isnan(gd.histZetaM_EE_Im[m][2][1]));
        }
    }
    cballs_scalar_window_free(&gd);
    cballs_scalar_window_free(&gd);
    CHECK(!gd.scalar_window_ready && gd.scalar_window_bins == 0);
    for (int i=0; i<6; i++) free_histogram(histograms[i]);
    puts("PASS: three allocation failures, published validity/NaN/zero, repeated cleanup");
    return 0;
}

int main(void)
{ return analytic_solver() || publication() ? EXIT_FAILURE : EXIT_SUCCESS; }
