// Use:
//#include "startrun_save_restore_05b.h"

#ifndef _startrun_save_restore_05b_h
#define _startrun_save_restore_05b_h

} else { // if !strnull(cmd.restorefile)

    // We must check the order of memory allocation and deallocation
    
    restorestate(cmd, gd, cmd->restorefile);
    startrun_ParamStat(cmd, gd);

    gd->stopflag = 0;                           // is 1 after stop-state
    
    setFilesDirs(cmd, gd);
    setFilesDirs_log(cmd, gd);
    strcpy(gd->mode,"a");                       // set at restorestate
    if(!(gd->outlog=fopen(gd->logfilePath, gd->mode)))
        error("\nstart_Common: error opening file '%s' \n",gd->logfilePath);
    
    verb_log_print(cmd->verbose_log, gd->outlog,
                   "\n\nAdded after restart from restart file\n");
    verb_log_print(cmd->verbose_log, gd->outlog, "\nnbody=%d\n",cmd->nbody);
    for (ifile=0; ifile<gd->ninfiles; ifile++) {
        verb_log_print(cmd->verbose_log, gd->outlog, "\nnbody=%d\n",gd->nbodyTable[ifile]);
    }
    
    class_call_cballs(random_init(cmd, gd, cmd->seed), errmsg, errmsg);

    if (!strnull(cmd->restorefile) && !strnull(cmd->infile)) {
        fprintf(stdout,
            "\nCheckParameters: Warning! : %s\n\n",
            "You are using options restorefile and infile at the same time");
    }

    class_call_cballs(CheckParameters(cmd, gd), errmsg, errmsg);
    
    PrintState(cmd, gd, cmd->restorefile);

    //B with this routine histograms arrays are allocated...
    //      it is located in restorestate
    //  class_call_cballs(startrun_memoryAllocation(cmd, gd), errmsg, errmsg);
    //      contains allocations for:
    //  rBins,
    //  histZetaMFlatten,
    //  histNN,
    //  histCF,
    //  histNNSub,
    //  histNNSubXi2pcf,
    //  histNNSubXi2pcftotal,
    //  histNNN,
    //  histXi2pcf,
    //  if computeTPCF:
    //      histXi,
    //      histXicos,
    //      histXisin,
    //      histZetaM,
    //      histZetaMcos,
    //      histZetaMsin,
    //      histZetaMsincos,
    //      histZetaMcossin,
    //      ifdef USEGSL
    //          histXi_gsl,
    //          histZetaMatrix,
    //          histZetaMatrix[m].histZetaM,
    //      end ifdef
    //      histZetaGmRe,
    //      histZetaGmIm,
    //  end if
    //  if computeShearCF:
    //      histXitt,
    //      histXixx,
    //      histXitx,
    //  end if
    //#ifdef ADDONS
    //#include "startrun_include_10.h"              // should be sync with
    //  "cballsio_include_10.h"
    //#endif
    //E
    
}

#endif	// ! _startrun_save_restore_05b_h
