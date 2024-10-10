// Use:
//#include "cballs_octree_kkk_omp_00.h"

#ifndef _cballs_octree_kkk_omp_00_h
#define _cballs_octree_kkk_omp_00_h

#ifdef NMultipoles
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd);
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
                               struct  global_data* gd);
#ifdef NONORMHIST
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
                                           struct  global_data* gd);
// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_normalized(struct  cmdline_data* cmd,
                                            struct  global_data* gd);
#endif
#endif // ! NMultipoles

#endif	// ! _cballs_octree_kkk_omp_00_h
