// Use:
//NOLSST:
//#include "tpcf_patch.h"

#ifndef _tpcf_patch_h
#define _tpcf_patch_h

        case DIRECTSIMPLE:          // search=direct-simple
            verb_print(cmd.verbose, "\n\tevalHist: direct method simple\n\n");
            searchcalc_direct_simple(bodytab, cmd.nbody, 1, cmd.nbody); break;

case DIRECTSIMPLESINCOS:
    verb_print(cmd.verbose, "\n\tevalHist: direct method simple (base sincos)\n\n");
    evalHistograms_direct_simple_sincos(); break;

#endif	// ! _tpcf_02_h
