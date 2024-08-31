// Use:
//#include "cballs_pxd_05.h"

#ifndef _cballs_pxd_05_h
#define _cballs_pxd_05_h

int get_theta(struct  cmdline_data* cmd, real *theta)
{
    *theta = cmd->theta;
    return SUCCESS;
}

int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN)
{
    *sizeHistN = cmd->sizeHistN;
    return SUCCESS;
}

int get_rBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        gd->rBins[n] = rBin;
    }

    return SUCCESS;
}

// get matrix ZetaM for each m multipole
#define _COS_         1
#define _SIN_         2
#define _SINCOS_      3
#define _COSSIN_      4
int get_HistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd,
                                int m, int type, ErrorMsg errmsg)
{
    int n1, n2;
    int n;

    class_test((m <= 0 || m > cmd->mChebyshev + 1),
               errmsg,
               "\nget_HistZetaM_sincos: not allowed value of m = %d\n",
               m);

    switch(type) {
        case _COS_:
            n = 1;
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->histZetaMFlatten[n++] = gd->histZetaMcos[m][n1][n2];
                }
            }
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                        "\nget_HistZetaM_sincos: n = %d vs sqr(sizeHistN) = %d\n",
                           n,cmd->sizeHistN);
            break;
        case _SIN_:
            n = 1;
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->histZetaMFlatten[n++] = gd->histZetaMsin[m][n1][n2];
                }
            }
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                        "\nget_HistZetaM_sincos: n = %d vs sqr(sizeHistN) = %d\n",
                           n,cmd->sizeHistN);
            break;
        case _SINCOS_:
            n = 1;
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->histZetaMFlatten[n++] = gd->histZetaMsincos[m][n1][n2];
                }
            }
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                        "\nget_HistZetaM_sincos: n = %d vs sqr(sizeHistN) = %d\n",
                           n,cmd->sizeHistN);
            break;
        case _COSSIN_:
            n = 1;
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->histZetaMFlatten[n++] = gd->histZetaMcossin[m][n1][n2];
                }
            }
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                        "\nget_HistZetaM_sincos: n = %d vs sqr(sizeHistN) = %d\n",
                           n,cmd->sizeHistN);
            break;
        default:
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                           "\nget_HistZetaM_sincos: default type...\n");
            n = 1;
            for (n1=1; n1<=cmd->sizeHistN; n1++) {
                for (n2=1; n2<=cmd->sizeHistN; n2++) {
                    gd->histZetaMFlatten[n++] = gd->histZetaMcos[m][n1][n2];
                }
            }
            if (cmd->verbose>=3)
                verb_print(cmd->verbose,
                        "\nget_HistZetaM_sincos: n = %d vs sqr(sizeHistN) = %d\n",
                           n,cmd->sizeHistN);
            break;
    }

    return SUCCESS;
}
#undef _COS_
#undef _SIN_
#undef _SINCOS_
#undef _COSSIN_

//E PXD functions


#endif	// ! _cballs_pxd_05_h
