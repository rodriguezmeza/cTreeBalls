// Use:
//#include "cballs_octree_kkk_omp_05.h"

#ifndef _cballs_octree_kkk_omp_05_h
#define _cballs_octree_kkk_omp_05_h

local int PrintHistN2pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    sprintf(namebuf, "%s%s%s", gd->fpfnamehistN2pcfFileName,
            cmd->suffixOutFiles, EXTFILES);
    verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
    outstr = stropen(namebuf, "w!");


    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            if (cmd->rminHist==0) {
                rbinlog = ((real)(n-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
            } else {
                rbinlog = rlog10(cmd->rminHist) + ((real)(n)-0.5)*gd->deltaR;
            }
            rBin=rpow(10.0,rbinlog);
        } else {
            rBin = cmd->rminHist + ((real)n-0.5)*gd->deltaR;
        }
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->histN2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}

#endif	// ! _cballs_octree_kkk_omp_05_h
