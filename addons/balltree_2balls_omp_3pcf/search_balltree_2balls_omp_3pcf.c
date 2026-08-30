/*
 * Compile the shared TreeCorr numerical kernels with the task-frontier 3PCF
 * scheduler and a distinct public entry point.
 */
#define TREECORR_TASK_FRONTIER_ENGINE 1
#define TREECORR_LOG_MULTIPOLE_ENGINE 1
#define TREECORR_METHOD_NAME "balltree-2balls-omp_3pcf"
#define searchcalc_balltree_2balls_omp searchcalc_balltree_2balls_omp_3pcf
#ifdef TWOPCF
#undef TWOPCF
#endif
#include "../balltree_2balls_omp/search_balltree_2balls_omp.c"
