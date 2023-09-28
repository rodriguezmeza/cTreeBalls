// Use:
//NOLSST:
//#include "protodefs_patch.h"

#ifndef _protodefs_patch_h
#define _protodefs_patch_h

//B Tree utilities
global int search_init(gdhistptr hist);
global int search_free(gdhistptr hist);
global int computeBodyProperties(bodyptr p, int nbody, gdhistptr hist);

global int search_init_3pcfbf(gdhistptr_3pcfbf hist);
global int search_free_3pcfbf(gdhistptr_3pcfbf hist);
global int computeBodyProperties_3pcfbf(bodyptr p, int nbody, gdhistptr_3pcfbf hist);
//E

#endif	// ! _protodefs_patch_h
