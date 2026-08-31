#ifndef _lya_forest_defs_h
#define _lya_forest_defs_h

#define LyaForestId(x) (((bodyptr)(x))->lya_forest_id)
#define LyaDistance(x) (((bodyptr)(x))->lya_distance)
#define LyaLOS(x)      (((bodyptr)(x))->lya_los)

#define INLYAASCII 18

/* Shared family classification for validation and MPI dispatch. */
static inline int lya_forest_method_kind(const char *name)
{
    if (name == NULL) return -1;
    if (strcmp(name, "lya-2pcf-omp") == 0 || strcmp(name, "lya-2pcf-mpi") == 0) return 0;
    if (strcmp(name, "lya-3pcf-omp") == 0 || strcmp(name, "lya-3pcf-mpi") == 0) return 1;
    if (strcmp(name, "lya-2pcf-3pcf-omp") == 0 || strcmp(name, "lya-2pcf-3pcf-mpi") == 0) return 2;
    if (strcmp(name, "lya-1d-2pcf-omp") == 0 || strcmp(name, "lya-1d-2pcf-mpi") == 0) return 3;
    if (strcmp(name, "lya-1d-3pcf-omp") == 0 || strcmp(name, "lya-1d-3pcf-mpi") == 0) return 4;
    if (strcmp(name, "lya-1d-2pcf-3pcf-omp") == 0 || strcmp(name, "lya-1d-2pcf-3pcf-mpi") == 0) return 5;
    if (strcmp(name, "lya-1d-tree-2pcf-omp") == 0 || strcmp(name, "lya-1d-tree-2pcf-mpi") == 0) return 6;
    return -1;
}

static inline int lya_forest_is_mpi_method(const char *name)
{
    if (lya_forest_method_kind(name) < 0) return 0;
    return strcmp(strrchr(name, '-'), "-mpi") == 0;
}


#endif
