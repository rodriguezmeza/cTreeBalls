// Use:
//#include "cballs_pxd_05.h"

// included in addons/addons_include/source/cballs_include_05.h

#ifndef _cballs_pxd_05_h
#define _cballs_pxd_05_h

//B parameters section

//B flags
int get_tree_allocated(struct global_data* gd, short *value)
{
    *value = gd->tree_allocated;
    return SUCCESS;
}

int get_allocated_2(struct global_data* gd, short *value)
{
    *value = gd->gd_allocated_2;
    return SUCCESS;
}

int get_bodytable_allocated(struct global_data* gd, short *value)
{
    *value = gd->bodytable_allocated;
    return SUCCESS;
}

int get_histograms_allocated(struct global_data* gd, short *value)
{
    *value = gd->histograms_allocated;
    return SUCCESS;
}

int get_gd_allocated(struct global_data* gd, short *value)
{
    *value = gd->gd_allocated;
    return SUCCESS;
}

int get_cmd_allocated(struct global_data* gd, short *value)
{
    *value = gd->cmd_allocated;
    return SUCCESS;
}
//E

int get_nthreads(struct  cmdline_data* cmd, int *value)
{
    *value = cmd->numthreads;
    return SUCCESS;
}

//B from codex
// This (the getters) approach is safer than accessing self.cmd.scanLevel directly from Cython,
//  because the C helper is compiled with the same macros as the C library
//  and can guard unsupported builds cleanly.
int get_scanLevel(struct cmdline_data* cmd, int *value)
{
#if defined(OCTREESMOOTHING) || defined(BALLS)
    *value = cmd->scanLevel;
    return SUCCESS;
#else
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "scanLevel is only available when OCTREESMOOTHING or BALLS is enabled");
    return FAILURE;
#endif
}

int get_scanLevelRoot(struct cmdline_data* cmd, int *value)
{
#if defined(OCTREESMOOTHING) || defined(BALLS)
    *value = cmd->scanLevelRoot;
    return SUCCESS;
#else
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "scanLevelRoot is only available when OCTREESMOOTHING or BALLS is enabled");
    return FAILURE;
#endif
}

int get_scanLevelMin(struct cmdline_data* cmd, int *value)
{
#if defined(OCTREESMOOTHING) || defined(BALLS)
    *value = cmd->scanLevelMin;
    return SUCCESS;
#else
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "scanLevelMin is only available when OCTREESMOOTHING or BALLS is enabled");
    return FAILURE;
#endif
}
//E


//B version 1.0.1
int get_nmultipoles(struct  cmdline_data* cmd, int *value)
{
    *value = cmd->mChebyshev;
    return SUCCESS;
}
//E

int get_theta(struct  cmdline_data* cmd, real *theta)
{
    *theta = cmd->theta;
    return SUCCESS;
}

int get_rsmooth(struct  global_data* gd, real *value)
{
    *value = gd->rsmooth[0];
    return SUCCESS;
}

int get_cputime(struct  global_data* gd, real *cputime)
{
    *cputime = gd->cpusearch;
    return SUCCESS;
}

int get_sizeHistN(struct  cmdline_data* cmd, int *sizeHistN)
{
    *sizeHistN = cmd->sizeHistN;
    return SUCCESS;
}

// use same version as is in cmdline_defs.h and/or in addons/class_lib/common.h
//  see also setup.py
int get_version(struct  cmdline_data* cmd, char *param)
{
    sprintf(param,"%s","1.0.1");
    return SUCCESS;
}


int get_rootDir(struct cmdline_data* cmd, char *value)
{
    if (cmd->rootDir == NULL) {
        sprintf(value, "%s", "");
        return SUCCESS;
    }

    snprintf(value, MAXLENGTHOFFILES, "%s", cmd->rootDir);
    return SUCCESS;
}
//E parameters section

//B global parameter section

//B
int get_nbody(struct cmdline_data* cmd, struct global_data* gd, INTEGER *value)
{
    (void)gd;

    if (cmd == NULL || value == NULL)
        return FAILURE;

    *value = cmd->nbody;
    return SUCCESS;
}
//E

int get_computeTPCF(struct  cmdline_data* cmd, struct  global_data* gd, short *value)
{
    *value = gd->computeTPCF;

    return SUCCESS;
}
//E global parameter section


//B added by cBalls
local int require_live_hist_vector(struct cmdline_data* cmd,
                                   struct global_data* gd,
                                   string routineName,
                                   realptr src,
                                   string srcName)
{
    if (cmd->sizeHistN <= 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: invalid sizeHistN=%d", routineName, cmd->sizeHistN);
        return FAILURE;
    }

    if (gd->histograms_allocated != TRUE) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: histogram arrays are not allocated", routineName);
        return FAILURE;
    }

    if (gd->vecPXD == NULL || src == NULL) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "%s: missing PXD vector or source histogram %s", routineName, srcName);
        return FAILURE;
    }

    return SUCCESS;
}
//E

//B histograms section
int get_rBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;

    if (require_live_hist_vector(cmd, gd, "get_rBins", gd->rBins, "rBins") == FAILURE)
        return FAILURE;

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

//B added by cBalls
int get_HistNN(struct cmdline_data* cmd, struct global_data* gd)
{
    int n;

    if (require_live_hist_vector(cmd, gd, "get_HistNN", gd->histNN, "histNN") == FAILURE)
        return FAILURE;

    for (n=1; n<=cmd->sizeHistN; n++)
        gd->vecPXD[n] = gd->histNN[n];

    return SUCCESS;
}
//E


int get_HistCF(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    if (require_live_hist_vector(cmd, gd, "get_HistCF", gd->histCF, "histCF") == FAILURE)
        return FAILURE;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histCF[n];
    }

    return SUCCESS;
}

int get_HistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    if (require_live_hist_vector(cmd, gd, "get_HistXi2pcf", gd->histXi2pcf, "histXi2pcf") == FAILURE)
        return FAILURE;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histXi2pcf[n];
    }

    return SUCCESS;
}

//B cross
int get_HistXi2pcf12(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    if (require_live_hist_vector(cmd, gd, "get_HistXi2pcf12", gd->histXi2pcf12, "histXi2pcf12") == FAILURE)
        return FAILURE;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histXi2pcf12[n];
    }

    return SUCCESS;
}

int get_HistXi2pcf13(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n;

    if (require_live_hist_vector(cmd, gd, "get_HistXi2pcf13", gd->histXi2pcf13, "histXi2pcf13") == FAILURE)
        return FAILURE;

    for (n=1; n<=cmd->sizeHistN; n++) {
        gd->vecPXD[n] = gd->histXi2pcf13[n];
    }

    return SUCCESS;
}
//E

// get matrix ZetaM for each m multipole
#define _COS_         1
#define _SIN_         2
#define _SINCOS_      3
#define _COSSIN_      4

int get_HistZetaMsincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int m, int type, ErrorMsg errmsg)
{
    string routineName = "get_HistZetaMsincos";
    int n1, n2;
    real ***histZetaMptr = NULL;
    string histName = "histZetaMcos";

    if (gd->computeTPCF==TRUE) {
        class_test((m <= 0 || m > cmd->mChebyshev + 1),
                   errmsg,"\n%s: not allowed value of m = %d\n", routineName, m);
        class_test((gd->matPXD == NULL),
                   errmsg,"\n%s: matPXD is not allocated; call cyballs Run(level=[\"MainLoop\"]) before PXD getters\n", routineName);

        switch(type) {
            case _COS_:
                histZetaMptr = gd->histZetaMcos;
                histName = "histZetaMcos";
                break;
            case _SIN_:
                histZetaMptr = gd->histZetaMsin;
                histName = "histZetaMsin";
                break;
            case _SINCOS_:
                histZetaMptr = gd->histZetaMsincos;
                histName = "histZetaMsincos";
                break;
            case _COSSIN_:
                histZetaMptr = gd->histZetaMcossin;
                histName = "histZetaMcossin";
                break;
            default:
                class_test(TRUE,
                           errmsg,"\n%s: not allowed value of type = %d (use 1=cos, 2=sin, 3=sincos, 4=cossin)\n", routineName, type);
        }

        class_test((histZetaMptr == NULL),
                   errmsg,"\n%s: %s is not allocated\n", routineName, histName);
        class_test((histZetaMptr[m] == NULL),
                   errmsg,"\n%s: %s[%d] is not allocated\n", routineName, histName, m);

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            class_test((histZetaMptr[m][n1] == NULL),
                       errmsg,"\n%s: %s[%d][%d] is not allocated\n", routineName, histName, m, n1);
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                gd->matPXD[n1-1][n2-1] = histZetaMptr[m][n1][n2];
            }
        }
    }

    return SUCCESS;
}

#undef _COS_
#undef _SIN_
#undef _SINCOS_
#undef _COSSIN_

int get_HistZetaM_EE(struct  cmdline_data* cmd,
                     struct  global_data* gd,
                     int m, ErrorMsg errmsg)
{
    string routineName = "get_HistZetaM_EE";
    int n1, n2;

    if (gd->computeTPCF==TRUE) {
        class_test((m <= 0 || m > cmd->mChebyshev + 1),
                   errmsg,"\n%s: not allowed value of m = %d\n", routineName, m);
        class_test((gd->matPXD == NULL),
                   errmsg,"\n%s: matPXD is not allocated; call cyballs Run(level=[\"MainLoop\"]) before PXD getters\n", routineName);
        class_test((gd->histZetaM_EE == NULL),
                   errmsg,"\n%s: histZetaM_EE is not allocated\n", routineName);
        class_test((gd->histZetaM_EE[m] == NULL),
                   errmsg,"\n%s: histZetaM_EE[%d] is not allocated\n", routineName, m);

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            class_test((gd->histZetaM_EE[m][n1] == NULL),
                       errmsg,"\n%s: histZetaM_EE[%d][%d] is not allocated\n", routineName, m, n1);
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                gd->matPXD[n1-1][n2-1] = gd->histZetaM_EE[m][n1][n2];
            }
        }
    }

    return SUCCESS;
}

int get_HistZetaM_EE_Im(struct cmdline_data* cmd,
                        struct global_data* gd,
                        int m, ErrorMsg errmsg)
{
    string routineName = "get_HistZetaM_EE_Im";
    int n1, n2;

    if (gd->computeTPCF==TRUE) {
        class_test((m <= 0 || m > cmd->mChebyshev + 1),
                   errmsg,"\n%s: not allowed value of m = %d\n", routineName, m);
        class_test((gd->matPXD == NULL),
                   errmsg,"\n%s: matPXD is not allocated; call cyballs Run(level=[\"MainLoop\"]) before PXD getters\n", routineName);
        class_test((gd->histZetaM_EE_Im == NULL),
                   errmsg,"\n%s: histZetaM_EE_Im is not allocated\n", routineName);
        class_test((gd->histZetaM_EE_Im[m] == NULL),
                   errmsg,"\n%s: histZetaM_EE_Im[%d] is not allocated\n", routineName, m);

        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            class_test((gd->histZetaM_EE_Im[m][n1] == NULL),
                       errmsg,"\n%s: histZetaM_EE_Im[%d][%d] is not allocated\n",
                       routineName, m, n1);
            for (n2=1; n2<=cmd->sizeHistN; n2++)
                gd->matPXD[n1-1][n2-1] = gd->histZetaM_EE_Im[m][n1][n2];
        }
    }
    return SUCCESS;
}

//E histograms section

//E PXD functions


#endif	// ! _cballs_pxd_05_h
