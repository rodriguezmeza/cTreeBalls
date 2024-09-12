// Use:
//#include "cballsio_kdtree_omp_10.h"

#ifndef _cballsio_kdtree_omp_10_h
#define _cballsio_kdtree_omp_10_h

#ifdef NMultipoles
if (cmd->computeTPCF) {
    free_dmatrix3D(gd->NhistZetaGmIm,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gd->NhistZetaGmRe,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    free_dmatrix3D(gd->NhistZetaMcossin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gd->NhistZetaMsincos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gd->NhistZetaMsin,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gd->NhistZetaMcos,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(gd->NhistZetaM,
                   1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix(gd->NhistXisin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(gd->NhistXicos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    free_dmatrix(gd->NhistXi,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
}

 free_dvector(gd->NhistXi2pcf,1,cmd->sizeHistN);

free_dvector(gd->NhistNNN,1,cmd->sizeHistN);
// 2pcf
//B kappa Avg Rmin
free_dvector(gd->NhistNNSubXi2pcftotal,1,cmd->sizeHistN);
//E
free_dvector(gd->NhistNNSubXi2pcf,1,cmd->sizeHistN);
//
free_dvector(gd->NhistNNSub,1,cmd->sizeHistN);
free_dvector(gd->NhistCF,1,cmd->sizeHistN);
free_dvector(gd->NhistNN,1,cmd->sizeHistN);
#endif

#endif	// ! _cballsio_kdtree_omp_10_h
