// Use:
//#include "cballsio_include_10.h"

#ifndef _cballsio_include_10_h
#define _cballsio_include_10_h

//#ifdef OCTREEKKKOMP
//#include "cballsio_octree_kkk_omp_10.h"
//#endif

/*
// problems with OCTREEKKKOMPON = 0
free_dvector(gd->histN2pcf,1,cmd->sizeHistN);
// 2pcf
//B kappa Avg Rmin
free_dvector(gd->histNNSubN2pcftotal,1,cmd->sizeHistN);
//E
free_dvector(gd->histNNSubN2pcf,1,cmd->sizeHistN);
//
*/

/*
 Add your addon item here
 */

#ifdef OCTREESHEAROMP
free(gd->histShearGammaIm);
gd->histShearGammaIm = NULL;
free(gd->histShearGammaRe);
gd->histShearGammaRe = NULL;
free(gd->histShearDenominatorIm);
gd->histShearDenominatorIm = NULL;
free(gd->histShearDenominatorRe);
gd->histShearDenominatorRe = NULL;
free(gd->histShearGammaMultipoleIm);
gd->histShearGammaMultipoleIm = NULL;
free(gd->histShearGammaMultipoleRe);
gd->histShearGammaMultipoleRe = NULL;
free(gd->histShearGammaNumeratorIm);
gd->histShearGammaNumeratorIm = NULL;
free(gd->histShearGammaNumeratorRe);
gd->histShearGammaNumeratorRe = NULL;
free(gd->histShearXiWeight);
gd->histShearXiWeight = NULL;
free(gd->histShearXiMinusIm);
gd->histShearXiMinusIm = NULL;
free(gd->histShearXiMinusRe);
gd->histShearXiMinusRe = NULL;
free(gd->histShearXiPlusIm);
gd->histShearXiPlusIm = NULL;
free(gd->histShearXiPlusRe);
gd->histShearXiPlusRe = NULL;
gd->shearMultipoleMax = 0;
gd->shearAngularBins = 0;
#endif

#endif	// ! _cballsio_include_10_h
