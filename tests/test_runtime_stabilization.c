/* Small analytic and fault-injection cases; no large catalog allocation. */
/* Define the runtime globals here, so the archive need not pull in main.o. */
#define global
#include "globaldefs.h"
#include "kdtree.h"
#include <errno.h>
#include <limits.h>
#include <unistd.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #c); return 1; \
} } while (0)

static int normalization(void)
{
    struct cmdline_data cmd = {0};
    struct global_data gd = {0};
    real nn[4], cf[4];
    const INTEGER counts[] = {46340, 46341, 65536,
#ifdef LONGINT
                               (INTEGER)INT_MAX + 17,
#endif
                               8};
    /* Independent, explicit shell endpoints: linear, log, zero-cutoff log.
       The zero-cutoff first bin reflects C truncation in existing producers. */
    const double lower[3][3] = {{1, 2, 3}, {1, 2, 4}, {.0001, .01, .1}};
    const double upper[3][3] = {{2, 3, 4}, {2, 4, 8}, {.01, .1, 1}};
    const double xi[3] = {-.25, .5, 0};
    cmd.options = "compute-HistN,and-CF";
    cmd.sizeHistN = 3;
    gd.histNN = nn;
    gd.histCF = cf;
    for (int k = 0; k < NDIM; k++) gd.Box[k] = 10;
    for (int mode = 0; mode < 3; mode++) {
        cmd.useLogHist = mode != 0;
        cmd.rminHist = mode == 2 ? 0 : 1;
        cmd.rangeN = mode == 0 ? 4 : mode == 1 ? 8 : 1;
        cmd.logHistBinsPD = 1;
        gd.deltaR = mode == 1 ? log10(2.) : 1;
        for (size_t j = 0; j < sizeof(counts)/sizeof(counts[0]); j++) {
            double ordered[3];
            for (int n = 0; n < 3; n++) {
                double shell = NDIM == 3
                    ? 4*PI/3*(pow(upper[mode][n], 3)-pow(lower[mode][n], 3))
                    : PI*(pow(upper[mode][n], 2)-pow(lower[mode][n], 2));
                ordered[n] = (double)counts[j]*(double)counts[j]
                    * shell/pow(10., NDIM)*(1+xi[n]);
                nn[n+1] = ordered[n];
                cf[n+1] = 123;
            }
            CHECK(search_compute_HistN(&cmd, &gd, counts[j]) == SUCCESS);
            for (int n = 1; n <= 3; n++) {
                CHECK(fabs(cf[n]-xi[n-1]) < 1e-12);
                CHECK(fabs(nn[n]/ordered[n-1]-.5) < 1e-12);
            }
        }
    }
    for (int n = 1; n <= 3; n++) nn[n] = 0;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == SUCCESS);
    CHECK(cf[1] == -1 && cf[2] == -1 && cf[3] == -1);
    CHECK(search_compute_HistN(&cmd, &gd, 0) == FAILURE);
    CHECK(strstr(cmd.error_message, "body count") != NULL);
    gd.Box[0] = 0;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == FAILURE);
    gd.Box[0] = 10;
    nn[1] = NAN;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == FAILURE);
    nn[1] = -1;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == FAILURE);
    nn[1] = 0;
    gd.deltaR = 0;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == FAILURE);
    cmd.options = "compute-HistN";
    nn[1] = 12; cf[1] = 123;
    CHECK(search_compute_HistN(&cmd, &gd, 8) == SUCCESS);
    CHECK(nn[1] == 6 && cf[1] == 123);
    puts("PASS: analytic normalization, integer boundaries, empty/invalid bins");
    return 0;
}

static int constructors(void)
{
    struct cmdline_data cmd = {0};
    struct global_data gd = {0};
    body catalog[8] = {0};
    cmd.options = "";
    cmd.theta = 1;
    for (int i = 0; i < 8; i++) {
        for (int k = 0; k < NDIM; k++) Pos(&catalog[i])[k] = (i >> k)&1;
        Mass(&catalog[i]) = Weight(&catalog[i]) = Kappa(&catalog[i]) = 1;
    }
    CHECK(init_kdtree(&cmd, &gd, NULL, 8) == NULL);
    CHECK(init_kdtree(&cmd, &gd, catalog, 0) == NULL);
    CHECK(build_kdtree(&cmd, &gd, NULL, 2) == FAILURE);
    finish_kdtree(NULL);
    for (int fail = 0; fail < 4; fail++) {
        ballxptr kd;
        cballs_test_fail_allocation_after(fail);
        kd = init_kdtree(&cmd, &gd, catalog, 8);
        if (fail < 2) CHECK(kd == NULL);
        else {
            CHECK(kd != NULL);
            CHECK(build_kdtree(&cmd, &gd, kd, 2) == FAILURE);
        }
        CHECK(strstr(cmd.error_message, "allocation failed") != NULL);
        finish_kdtree(kd);
        cballs_test_reset_allocation_failure();
    }
    for (int repeat = 0; repeat < 64; repeat++) {
        ballxptr kd = init_kdtree(&cmd, &gd, catalog, 8);
        CHECK(kd != NULL);
        CHECK(build_kdtree(&cmd, &gd, kd, 0) == FAILURE);
        CHECK(build_kdtree(&cmd, &gd, kd, 2) == SUCCESS);
        CHECK(kd->body_order != NULL && kd->ntab != NULL);
        finish_kdtree(kd);
    }
    puts("PASS: KD context/pointer/node/order failures and repeated destruction");
    return 0;
}

static int directories(void)
{
    char root[] = "/tmp/cballs-mkdir-XXXXXX", path[512], nested[512];
    FILE *stream;
    CHECK(mkdir_p(NULL, 0700) == -1 && errno == EINVAL);
    CHECK(mkdir_p("", 0700) == -1 && errno == EINVAL);
    CHECK(mkdir_p("/", 0700) == 0);
    CHECK(mkdir_p(".", 0700) == 0);
    CHECK(mkdtemp(root) != NULL);
    snprintf(path, sizeof(path), "%s/file", root);
    stream = fopen(path, "w"); CHECK(stream != NULL); CHECK(fclose(stream) == 0);
    CHECK(mkdir_p(path, 0700) == -1 && errno == ENOTDIR);
    snprintf(nested, sizeof(nested), "%s/file/child", root);
    CHECK(mkdir_p(nested, 0700) == -1 && errno == ENOTDIR);
    CHECK(unlink(path) == 0);
    snprintf(path, sizeof(path), "%s/literal $(printf bad) 'quote'", root);
    CHECK(mkdir_p(path, 0700) == 0);
    CHECK(mkdir_p(path, 0700) == 0);
    CHECK(rmdir(path) == 0);
    CHECK(rmdir(root) == 0);
    puts("PASS: literal mkdir, existing directories, invalid input/file collisions");
    return 0;
}

int main(void)
{
    if (normalization() || constructors() || directories()) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
