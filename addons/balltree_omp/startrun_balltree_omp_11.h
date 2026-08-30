#ifndef _startrun_balltree_omp_11_h
#define _startrun_balltree_omp_11_h

if (strcmp(method_str, "balltree-omp") == 0)
    *method_int = 167;
#ifdef BALLTREEMPI
if (strcmp(method_str, "balltree-mpi") == 0)
    *method_int = 168;
#endif

#endif /* !_startrun_balltree_omp_11_h */
