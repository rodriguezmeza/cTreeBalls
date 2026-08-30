#ifndef _startrun_lya_forest_omp_11_h
#define _startrun_lya_forest_omp_11_h

if (strcmp(method_str, "lya-2pcf-omp") == 0) *method_int = 169;
if (strcmp(method_str, "lya-3pcf-omp") == 0) *method_int = 170;
if (strcmp(method_str, "lya-2pcf-3pcf-omp") == 0) *method_int = 171;
if (strcmp(method_str, "lya-1d-2pcf-omp") == 0)
    *method_int = LYAFOREST1D2PCFMETHOD;
if (strcmp(method_str, "lya-1d-3pcf-omp") == 0)
    *method_int = LYAFOREST1D3PCFMETHOD;
if (strcmp(method_str, "lya-1d-2pcf-3pcf-omp") == 0)
    *method_int = LYAFOREST1D2PCF3PCFMETHOD;
if (strcmp(method_str, "lya-1d-tree-2pcf-omp") == 0)
    *method_int = LYAFOREST1DTREE2PCFMETHOD;

#endif
