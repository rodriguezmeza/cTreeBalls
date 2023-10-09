// Use:
//NOLSST:
//#include "treeutils_patch.h"

#ifndef _treeutils_patch_h
#define _treeutils_patch_h

global int search_init(gdhistptr hist)
{
    int n;

#ifdef TPCF
    hist->Chebs = dvector(1,cmd.mchebyshev+1);
    hist->xiOUTVP = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
    hist->histZetaMtmp = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
#endif
    hist->histXi2pcfsub = dvector(1,cmd.sizeHistN);

    for (n = 1; n <= cmd.sizeHistN; n++) {
        hist->histXi2pcfsub[n] = 0.0;
    }

    return _SUCCESS_;
}

global int search_free(gdhistptr hist)
{
//    verb_print_debug(1, "\nAqui voy (00)\n");
    free_dvector(hist->histXi2pcfsub,1,cmd.sizeHistN);
#ifdef TPCF
    free_dmatrix(hist->histZetaMtmp,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(hist->xiOUTVP,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dvector(hist->Chebs,1,cmd.mchebyshev+1);
#endif
//    verb_print_debug(1, "\nAqui voy (000)\n");

    return _SUCCESS_;
}

global int search_init_3pcfbf(gdhistptr_3pcfbf hist)
{
    int n, m, l;

#ifdef TPCF
    hist->Chebs = dvector(1,cmd.mchebyshev+1);
//    hist->xiOUTVP = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
//    hist->histZetaMtmp = dmatrix(1,cmd.sizeHistN,1,cmd.sizeHistN);
#endif
    hist->histXi2pcfsub = dvector(1,cmd.sizeHistN);

#ifdef TPCF
    hist->histNNN = dvector(1,cmd.sizeHistN);
    hist->histNNNSub = dmatrix3D(1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    hist->histXi3pcf = dmatrix3D(1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    for (n = 1; n <= cmd.sizeHistN; n++) {
        hist->histNNN[n] = 0.0;
        for (m = 1; m <= cmd.sizeHistN; m++)
            for (l = 1; l <= cmd.sizeHistTheta; l++)
                hist->histXi3pcf[n][m][l] = 0.0;
    }
#endif

    for (n = 1; n <= cmd.sizeHistN; n++) {
        hist->histXi2pcfsub[n] = 0.0;
    }

    return _SUCCESS_;
}

global int search_free_3pcfbf(gdhistptr_3pcfbf hist)
{
#ifdef TPCF
    free_dmatrix3D(hist->histXi3pcf,1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix3D(hist->histNNNSub,1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dvector(hist->histNNN,1,cmd.sizeHistN);
#endif
    free_dvector(hist->histXi2pcfsub,1,cmd.sizeHistN);
#ifdef TPCF
//    free_dmatrix(hist->histZetaMtmp,1,cmd.sizeHistN,1,cmd.sizeHistN);
//    free_dmatrix(hist->xiOUTVP,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dvector(hist->Chebs,1,cmd.mchebyshev+1);
#endif

    return _SUCCESS_;
}


global int computeBodyProperties(bodyptr p, int nbody, gdhistptr hist)
{
    int n;
    int m;
    real xi, xi_2p;

    xi = Kappa(p)/nbody;
    xi_2p = Kappa(p);
#ifdef TPCF
    for (m=1; m<=cmd.mchebyshev+1; m++)
        for (n=1; n<=cmd.sizeHistN; n++)
            gd.histXi[m][n] /= MAX(gd.histNSub[n],1.0);
    for (m=1; m<=cmd.mchebyshev+1; m++){
        OUTVP_ext(hist->xiOUTVP, gd.histXi[m], gd.histXi[m],cmd.sizeHistN);
        CLRM_ext(hist->histZetaMtmp,cmd.sizeHistN);
        MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd.sizeHistN);
        ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],hist->histZetaMtmp,cmd.sizeHistN);
    }
#endif

    for (n=1; n<=cmd.sizeHistN; n++) {
        gd.histXi2pcf[n] += xi_2p*hist->histXi2pcfsub[n];
    }

    return _SUCCESS_;
}

global int computeBodyProperties_3pcfbf(bodyptr p, int nbody, gdhistptr_3pcfbf hist)
{
    int n;
    int m;
    real xi, xi_2p;

    xi = Kappa(p)/nbody;
    xi_2p = Kappa(p);
#ifdef TPCF
/*
    for (m=1; m<=cmd.mchebyshev+1; m++)
        for (n=1; n<=cmd.sizeHistN; n++)
            gd.histXi[m][n] /= MAX(gd.histNSub[n],1.0);
    for (m=1; m<=cmd.mchebyshev+1; m++){
        OUTVP_ext(hist->xiOUTVP, gd.histXi[m], gd.histXi[m],cmd.sizeHistN);
        CLRM_ext(hist->histZetaMtmp,cmd.sizeHistN);
        MULMS_ext(hist->histZetaMtmp,hist->xiOUTVP,xi,cmd.sizeHistN);
        ADDM_ext(gd.histZetaM[m],gd.histZetaM[m],hist->histZetaMtmp,cmd.sizeHistN);
    }
*/
#endif

    for (n=1; n<=cmd.sizeHistN; n++) {
        gd.histXi2pcf[n] += xi_2p*hist->histXi2pcfsub[n];
    }

    return _SUCCESS_;
}

#endif	// ! _treeutils_patch_h
