#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "chealpix.h"

#define CHECK(condition, message)                     \
    do {                                              \
        if (!(condition)) {                           \
            fprintf(stderr, "FAIL: %s\n", message); \
            goto cleanup;                             \
        }                                             \
    } while (0)

int main(void)
{
    const long nside = 2;
    const long npix = nside2npix(nside);
    char coordsys[32] = {0};
    char filename[160] = {0};
    char ordering[32] = {0};
    float *nested = NULL;
    float *map = NULL;
    float *original_map;
    long file_nside = 0;
    long file_npix = 0;
    long ipnest;
    long ipring;
    int status;
    int result = EXIT_FAILURE;

    nested = (float *)malloc((size_t)npix * sizeof(*nested));
    CHECK(nested != NULL, "could not allocate test map");

    for (ipnest = 0; ipnest < npix; ipnest++) {
        nest2ring(nside, ipnest, &ipring);
        CHECK(ipring >= 0 && ipring < npix, "nest2ring rejected valid pixel");
        nested[ipnest] = (float)ipring + 0.25f;
    }

    snprintf(filename, sizeof(filename),
             "/tmp/cballs-healpix-ordering-%ld.fits", (long)getpid());
    unlink(filename);
    status = write_healpix_map_status(nested, nside, filename, 1, "G");
    CHECK(status == 0, "could not write NESTED FITS map");

    status = get_fits_size_status(filename, &file_nside, ordering,
                                  &file_npix);
    CHECK(status == 0, "could not read FITS map dimensions");
    CHECK(file_nside == nside && file_npix == npix,
          "FITS map dimensions changed");

    status = read_healpix_map_status(filename, &file_nside, coordsys,
                                     ordering, &map);
    CHECK(status == 0 && map != NULL, "could not read NESTED FITS map");
    original_map = map;
    status = healpix_map_to_ring_status(&map, file_nside, ordering);
    CHECK(status == 0, "NESTED-to-RING conversion failed");
    CHECK(map != original_map, "NESTED conversion did not replace the map");

    for (ipring = 0; ipring < npix; ipring++)
        CHECK(fabsf(map[ipring] - ((float)ipring + 0.25f)) < 1.0e-6f,
              "converted map has a value at the wrong RING pixel");

    original_map = map;
    status = healpix_map_to_ring_status(&map, file_nside, "  ring  ");
    CHECK(status == 0 && map == original_map,
          "RING input should remain unchanged");

    status = healpix_map_to_ring_status(&map, file_nside, "NUNIQ");
    CHECK(status != 0 && map == original_map,
          "unknown ordering should fail without taking ownership");
    status = healpix_map_to_ring_status(&map, 3, "NESTED");
    CHECK(status != 0 && map == original_map,
          "invalid NESTED NSIDE should fail without taking ownership");

    result = EXIT_SUCCESS;

cleanup:
    if (filename[0] != '\0')
        unlink(filename);
    free(map);
    free(nested);
    if (result == EXIT_SUCCESS)
        puts("PASS: CFITSIO NESTED maps are normalized to RING");
    return result;
}
