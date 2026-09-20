#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "globaldefs.h"
#ifdef KDTREEOMP
#include "kdtree.h"
#endif
#ifdef BALLTREEOMP
#include "fcfc_balltree.h"
#endif

/* Keep the standalone main.o out of libcballs.a when this focused test links. */
real *inout_xval = NULL;
real *inout_yval = NULL;
real *inout_zval = NULL;
real *inout_uval = NULL;
real *inout_vval = NULL;

_Static_assert(sizeof(real) == sizeof(double),
               "cTreeBalls computation precision must remain double");
_Static_assert(sizeof(((node *)0)->mass) == sizeof(real),
               "body mass must use computation precision");
_Static_assert(sizeof(((node *)0)->kappa) == sizeof(real),
               "body fields must use computation precision");
_Static_assert(sizeof(((node *)0)->weight) == sizeof(real),
               "body weights must use computation precision");
_Static_assert(sizeof(((node *)0)->pos[0]) == sizeof(cballs_storage_real),
               "positions must use geometry storage precision");
_Static_assert(sizeof(((node *)0)->radius) == sizeof(cballs_storage_real),
               "node radii must use geometry storage precision");
_Static_assert(sizeof(((cell *)0)->size) == sizeof(cballs_storage_real),
               "cell sizes must use geometry storage precision");

#ifdef SINGLEP
_Static_assert(sizeof(cballs_storage_real) == sizeof(float),
               "SINGLEP must select float geometry storage");
#ifdef KDTREEOMP
_Static_assert(sizeof(((kd_leaf_point *)0)->pos[0]) == sizeof(float),
               "packed KD coordinates must use SINGLEP storage");
_Static_assert(sizeof(((kd_leaf_point *)0)->kappa) == sizeof(real),
               "packed KD fields must retain computation precision");
#endif
#ifdef BALLTREEOMP
_Static_assert(sizeof(((fcfc_ballpoint *)0)->pos[0]) == sizeof(float),
               "packed ball coordinates must use SINGLEP storage");
_Static_assert(sizeof(((fcfc_ballpoint *)0)->kappa) == sizeof(real),
               "packed ball fields must retain computation precision");
#endif
#else
_Static_assert(sizeof(cballs_storage_real) == sizeof(double),
               "the default profile must preserve double geometry storage");
#endif

static void test_search_bounds(void)
{
    const real value = 0.10000000000000001;
    const real upper = (real)cballs_store_upper_bound(value);
    const real search = (real)cballs_store_search_bound(value);

    assert(upper >= value);
    assert(search >= value);
#ifdef SINGLEP
    assert(search > value);
#else
    assert(upper == value);
    assert(search == value);
#endif
}

static void test_binary_vector_format(void)
{
    const real source[NDIM] = {
        0.123456789012345,
#if NDIM > 1
        -2.987654321098765,
#endif
#if NDIM > 2
        17.00000011920929,
#endif
    };
    real disk[NDIM];
    vector stored;
    vector restored;
    char errmsg[256];
    FILE *file;
    int k;

    DO_COORD(k) {
        stored[k] = (cballs_storage_real)source[k];
        restored[k] = 0;
        disk[k] = 0;
    }

    file = tmpfile();
    assert(file != NULL);
    assert(out_vector_bin_checked(file, stored, "test_mixed_precision",
                                  "tmpfile", errmsg, sizeof(errmsg)) == SUCCESS);
    assert(fflush(file) == 0);
    assert(ftell(file) == (long)(NDIM * sizeof(real)));

    rewind(file);
    assert(fread(disk, sizeof(real), NDIM, file) == NDIM);
    DO_COORD(k)
        assert(disk[k] == (real)stored[k]);

    rewind(file);
    in_vector_bin(file, restored);
    DO_COORD(k)
        assert(restored[k] == stored[k]);
    assert(fclose(file) == 0);

    file = tmpfile();
    assert(file != NULL);
    assert(fwrite(source, sizeof(real), NDIM, file) == NDIM);
    rewind(file);
    in_vector_bin(file, restored);
    DO_COORD(k)
        assert(restored[k] == (cballs_storage_real)source[k]);
    assert(fclose(file) == 0);
}

int main(void)
{
    test_search_bounds();
    test_binary_vector_format();
    printf("PASS: %s; binary vectors remain %zu-byte doubles\n",
           CballsStoragePrecision, sizeof(real));
    return 0;
}
