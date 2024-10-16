// Use:
//#include "startrun_octree_kkk_omp_10.h"

#ifndef _startrun_octree_kkk_omp_10_h
#define _startrun_octree_kkk_omp_10_h

#ifdef NMultipoles
//B Already allocated in startrun_kdtree_omp_10.h
#ifndef KDTREEOMP
    gd->NhistNNSub = dvector(1,cmd->sizeHistN);
    bytes_tot_local += cmd->sizeHistN*sizeof(real);

    gd->NhistXi = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    gd->NhistXicos = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    gd->NhistXisin = dmatrix(1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    bytes_tot_local += 3*(cmd->mChebyshev+1)*cmd->sizeHistN*sizeof(real);
    gd->NhistZetaM = dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,
                              1,cmd->sizeHistN);
    bytes_tot_local +=
            (cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);

    gd->NhistZetaMcos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gd->NhistZetaMsin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gd->NhistZetaMsincos =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    gd->NhistZetaMcossin =
            dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            4*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
    gd->NhistZetaGmRe =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    gd->NhistZetaGmIm =
                dmatrix3D(1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
    bytes_tot_local +=
            2*(cmd->mChebyshev+1)*cmd->sizeHistN*cmd->sizeHistN*sizeof(real);
#endif
//E
#endif // ! NMultipoles

#endif	// ! _startrun_octree_kkk_omp_10_h
