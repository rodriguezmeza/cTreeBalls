/*==============================================================================
 MODULE: cballs.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose:
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

local int PrintHistNN(struct cmdline_data* cmd, struct  global_data* gd);
local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistXi3pcf(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM(struct  cmdline_data* cmd, struct  global_data* gd);
local int PrintHistZetaM_exp(struct  cmdline_data* cmd,
                             struct  global_data* gd);
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd);
// Seems obsolete. Delete. Only appears in balls_method addon
local int PrintHistZeta_theta2_fix(struct  cmdline_data* cmd,
                                   struct  global_data* gd);
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                               struct  global_data* gd);
local int PrintEvalHist(struct  cmdline_data* cmd, struct  global_data* gd);

local int PrintHistZetaG(struct  cmdline_data* cmd,
                                       struct  global_data* gd);
local int PrintHistZetaGm_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd);
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                      struct  global_data* gd);

#ifdef ADDONS
#include "cballs_include_01.h"
#endif


int MainLoop(struct  cmdline_data* cmd, struct  global_data* gd)
{
    bodyptr p,q;
    real kavg;
    INTEGER in;

    int ifile = 0;

#ifdef ADDONS
#include "cballs_include_02.h"
#endif

    gd->flagSmoothCellMin = FALSE;
    gd->flagSmooth = FALSE;
    gd->flagSetNbNoSel = FALSE;

    if (scanopt(cmd->options, "smooth-min-cell")) {
        verb_print(cmd->verbose,
                   "\n\tMainLoop: smooth cell min: try making tree...\n\n");
        DO_BODY(p,bodytable[gd->iCatalogs[0]],
                bodytable[gd->iCatalogs[0]]+gd->nbodyTable[gd->iCatalogs[0]])
            Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[gd->iCatalogs[0]],
                 gd->nbodyTable[gd->iCatalogs[0]], gd->iCatalogs[0]);
//B
        free(bodytable[gd->iCatalogs[0]]);
        gd->bytes_tot -= gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body);
        gd->nbodyTable[gd->iCatalogs[0]] = gd->nnodescanlevTable[gd->iCatalogs[0]];
        bodytable[ifile] =
                    (bodyptr) allocate(gd->nnodescanlevTable[gd->iCatalogs[0]]
                    * sizeof(body));
        gd->bytes_tot += gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body);
        verb_print(cmd->verbose,
"Allocated %g MByte for final (smoothCellMin) particle (%ld) storage (%ld).\n",
                   gd->nnodescanlevTable[gd->iCatalogs[0]]*sizeof(body)*INMB,
                   gd->nnodescanlevTable[gd->iCatalogs[0]],
                   gd->nnodescanlevTable[gd->iCatalogs[0]]);
        kavg = 0;
        q = nodetable;
        bodytable[gd->iCatalogs[0]] = nodetable;
        in = 0;
        DO_BODY(p,bodytable[gd->iCatalogs[0]],
            bodytable[gd->iCatalogs[0]]+gd->nnodescanlevTable[gd->iCatalogs[0]]) {
            kavg += Kappa(p);
            in++;
        }
        verb_print(cmd->verbose,"%ld %d %ld %lg %ld %lg\n",
                   in,Type(p-1),Id(p-1),Weight(p-1),Nb(p-1),Kappa(p-1));
//E
        verb_print(cmd->verbose,
                   "smoothCellMin: %ld particles in nodetable\n", in);
        verb_print(cmd->verbose,
                   "smoothCellMin: Average of kappa (%ld particles) = %le\n",
                   gd->nnodescanlevTable[gd->iCatalogs[0]],
                   kavg/((real)gd->nnodescanlevTable[gd->iCatalogs[0]]) );
        gd->flagSmoothCellMin = TRUE;
    } // ! smooth-min-cell

    if (scanopt(cmd->options, "smooth")
        && !scanopt(cmd->options, "set-Nb-noSel")) {
        verb_print(cmd->verbose, "\n\tMainLoop: smooth: try making tree...\n\n");
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
            Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], cmd->nbody, 0);
//B
        free(bodytable[ifile]);
        gd->bytes_tot -= cmd->nbody*sizeof(body);
        cmd->nbody = gd->nbodysm;
            bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
        gd->bytes_tot += cmd->nbody*sizeof(body);
        verb_print(cmd->verbose,
            "Allocated %g MByte for final (smooth) particle (%ld) storage.\n",
            cmd->nbody*sizeof(body)*INMB, cmd->nbody);
        kavg = 0;
        q = bodytabsm;
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            SETV(Pos(p),Pos(q));
// BODY3
            Type(p) = Type(q);
#ifdef BODY3ON
            Nbb(p) = Nbb(q);
#endif
            Weight(p) = Weight(q);
            Kappa(p) = Kappa(q);
            Id(p) = p-bodytable[ifile]+1;
            Update(p) = TRUE;
            q++;
            kavg += Kappa(p);
        }
        free(bodytabsm);
        gd->bytes_tot -= gd->nbodysm*sizeof(body);
//E
        verb_print(cmd->verbose,
                   "smooth: Average of kappa (%ld particles) = %le\n",
                   cmd->nbody, kavg/((real)cmd->nbody) );
        gd->flagSmooth = TRUE;
    } // ! smooth && set-Nb-noSel

    if ( scanopt(cmd->options, "smooth")
        && scanopt(cmd->options, "set-Nb-noSel") ) {
        verb_print(cmd->verbose,
               "\n\tMainLoop: smooth & set-Nb-noSel: try making tree...\n\n");
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
            Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], cmd->nbody, 0);
//B
        free(bodytable[ifile]);
        gd->bytes_tot -= cmd->nbody*sizeof(body);
        cmd->nbody = gd->nbodySel;
        bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
        gd->bytes_tot += cmd->nbody*sizeof(body);
        verb_print(cmd->verbose,
        "Allocated %g MByte for final (smooth-noSel) particle (%ld) storage.\n",
        cmd->nbody*sizeof(body)*INMB, cmd->nbody);
        kavg = 0;
        q = bodytabSel;
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
            SETV(Pos(p),Pos(q));
// BODY3
            Type(p) = Type(q);
#ifdef BODY3ON
            Nbb(p) = Nbb(q);
#endif
            Weight(p) = Weight(q);
            Kappa(p) = Kappa(q);
            Id(p) = p-bodytable[ifile]+1;
            Update(p) = TRUE;
            q++;
            kavg += Kappa(p);
        }
        free(bodytabSel);
        gd->bytes_tot -= gd->nbodySel*sizeof(body);
//E
        verb_print(cmd->verbose,
            "smooth-set-Nb-noSel: Average of kappa (%ld particles) = %le\n",
            cmd->nbody, kavg/((real)cmd->nbody) );
        gd->flagSmooth = TRUE;
        gd->flagSetNbNoSel = TRUE;
    } // ! smooth && set-Nb-noSel

    if ( scanopt(cmd->options, "make-tree") ) {
    verb_print(cmd->verbose, "\n\tMainLoop: make-tree: try making tree...\n\n");
        DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
            Update(p) = TRUE;
        MakeTree(cmd, gd, bodytable[ifile], cmd->nbody, 0);
    }

    if (scanopt(cmd->options, "stop")) {
        if (!strnull(cmd->outfile))
            Output(cmd, gd, bodytable, gd->nbodyTable, ifile);
        verb_print(cmd->verbose, "\n\tMainLoop: stopping...\n\n");
        exit(1);
    }

    EvalHist(cmd, gd);

    if (!gd->stopflag)
        PrintEvalHist(cmd, gd);

#ifdef DEBUG
    if (!strnull(cmd->outfile))
        Output(cmd, gd);
#endif

//B Post-processing:
    double cpustart;
    char buf[200];
    if (scanopt(cmd->options, "post-processing")) {
        cpustart = CPUTIME;
        sprintf(buf,"%s",cmd->script);
        verb_print(cmd->verbose, 
                   "\npost-processing: executing '%s'...", cmd->script);
        system(buf);
        verb_print(cmd->verbose, "done.\n");
        gd->cputotal += CPUTIME - cpustart;
        verb_print(cmd->verbose, "cpu time expended in this script %g\n\n",
                   CPUTIME - cpustart);
    }
//E

    return SUCCESS;
}

int EvalHist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    bodyptr p;
    int ifile;

    switch(gd->searchMethod_int) {
        case TREEOMPMETHODSINCOS:                   // search=tree-omp-sincos
            verb_print(cmd->verbose,
                    "\n\tEvalHist: with normal tree method (sincos-omp)\n\n");
            for (ifile=0; ifile<gd->ninfiles; ifile++) {
                DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
                Update(p) = TRUE;
                MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            }
            class_call_cballs(searchcalc_normal_sincos(cmd, gd, bodytable,
                            gd->nbodyTable, 1, gd->nbodyTable, gd->iCatalogs[0],
                            gd->iCatalogs[1]), errmsg, errmsg);
            break;

#ifdef ADDONS
#include "cballs_include_03.h"
#endif

        case SEARCHNULL:
            verb_print(cmd->verbose, "\n\tEvalHist: null search method.\n");
            verb_print(cmd->verbose,
                    "\n\tevalHist: with normal tree method (sincos-omp)\n\n");
            for (ifile=0; ifile<gd->ninfiles; ifile++) {
                DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
                Update(p) = TRUE;
                MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            }
            searchcalc_normal_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
                        gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
        default:
            verb_print(cmd->verbose, "\n\tEvalHist: dafault search method.\n");
            verb_print(cmd->verbose,
                    "\n\tevalHist: with normal tree method (sincos-omp)\n\n");
            for (ifile=0; ifile<gd->ninfiles; ifile++) {
                DO_BODY(p,bodytable[ifile],bodytable[ifile]+gd->nbodyTable[ifile])
                Update(p) = TRUE;
                MakeTree(cmd, gd, bodytable[ifile], gd->nbodyTable[ifile], ifile);
            }
            searchcalc_normal_sincos(cmd, gd, bodytable, gd->nbodyTable, 1,
                        gd->nbodyTable, gd->iCatalogs[0], gd->iCatalogs[1]);
            break;
    }

    return SUCCESS;
}

local int PrintEvalHist(struct  cmdline_data* cmd, struct  global_data* gd)
{
    double cpustart = CPUTIME;
    if (!scanopt(cmd->options, "no-out-Hist")) {
    switch(gd->searchMethod_int) {
        case TREEOMPMETHODSINCOS:
            verb_print(cmd->verbose,
            "\n\tPrintEvalHist: printing normal tree method (omp-sincos)\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            if (cmd->computeTPCF) {
                PrintHistZetaM_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-m-HistZeta"))
                    PrintHistZetaMm_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
            }
            
            break;
        case SEARCHNULL:
            verb_print(cmd->verbose,
                       "\n\tPrintEvalHist: printing null search method.\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            if (cmd->computeTPCF) {
                PrintHistZetaM_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-m-HistZeta"))
                    PrintHistZetaMm_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
            }
            break;
        default:
            verb_print(cmd->verbose,
                       "\n\tPrintEvalHist: printing dafault search method.\n\n");
            if (scanopt(cmd->options, "compute-HistN")) PrintHistNN(cmd, gd);
            PrintHistrBins(cmd, gd);
            PrintHistXi2pcf(cmd, gd);
            if (cmd->computeTPCF) {
                PrintHistZetaM_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-m-HistZeta"))
                    PrintHistZetaMm_sincos(cmd, gd);
                if (scanopt(cmd->options, "out-HistZetaG")) {
                    PrintHistZetaMZetaGm_sincos(cmd, gd);
                }
            }
            break;

#ifdef ADDONS
#include "cballs_include_04.h"
#endif

        }
    } // ! scanoptions no-out-Hist

    gd->cputotalinout += CPUTIME - cpustart;

    return SUCCESS;
}

local int PrintHistNN(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistNNFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistNNFileName);

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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->histNN[n]);
    }
    fclose(outstr);

    if (scanopt(cmd->options, "and-CF"))
        PrintHistCF(cmd, gd);

    return SUCCESS;
}

local int PrintHistCF(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistCFFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistCFFileName);

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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->histCF[n]);
    }
    fclose(outstr);

    return SUCCESS;
}

local int PrintHistrBins(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistrBinsFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistrBinsFileName);

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
        fprintf(outstr,"%16.8e\n",rBin);
    }
    fclose(outstr);

    return SUCCESS;
}

local int PrintHistXi2pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistXi2pcfFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistXi2pcfFileName);

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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->histXi2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}

//B Update or delete. Seems the same as above routine
local int PrintHistXi3pcf(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

    outstr = stropen(gd->fpfnamehistXi2pcfFileName, "w!");

    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistXi2pcfFileName);

    for (n=1; n<=cmd->sizeHistN; n++) {
        if (cmd->useLogHist) {
            rbinlog = ((real)(n)-0.5)*gd->deltaR;
            rBin=rpow(10,rbinlog);
        } else {
            rBin = ((int)n-0.5)*gd->deltaR;
        }
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->histXi2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}
//E


local int PrintHistZetaM(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin1, rBin2;
    real rbinlog1, rbinlog2;

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s_%d%s", gd->fpfnamehistZetaMFileName, m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                if (scanopt(cmd->options, "xr1r2")) {
                    if (cmd->useLogHist) {
                        if (cmd->rminHist==0) {
                            rbinlog1 = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
                            rbinlog2 = ((real)(n2-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
                        } else {
                            rbinlog1 = rlog10(cmd->rminHist)
                            + ((real)(n1)-0.5)*gd->deltaR;
                            rbinlog2 = rlog10(cmd->rminHist)
                            + ((real)(n2)-0.5)*gd->deltaR;
                        }
                        rBin1=rpow(10.0,rbinlog1);
                        rBin2=rpow(10.0,rbinlog2);
                    } else {
                        rBin1 = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
                        rBin2 = cmd->rminHist + ((real)n2-0.5)*gd->deltaR;
                    }
                    fprintf(outstr,"%16.8e ",
                            rBin1*rBin2*gd->histZetaM[m][n1][n2]);
                } else {
                    fprintf(outstr,"%16.8e ",gd->histZetaM[m][n1][n2]);
                }
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

#ifdef USEGSL
local int PrintHistZetaM_exp(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int n1, n2, m;
    stream outstr;

    outstr = stropen(gd->fpfnamehistZetaMFileName, "w!");
    verb_print_q(2, cmd->verbose,
               "Printing : to a file %s ...\n",gd->fpfnamehistZetaMFileName);

    gsl_complex Atmp;

    for (m=0; m<=cmd->mChebyshev; m++) {
        fprintf(outstr,"\n\nm=%d\n",m);
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                Atmp = gsl_matrix_complex_get(histZetaMatrix[m].histZetaM,n1,n1);
                fprintf(outstr,"%16.8e %16.8e ",GSL_REAL(Atmp),GSL_IMAG(Atmp));
            }
            fprintf(outstr,"\n");
        }
    }
    fclose(outstr);

    return SUCCESS;
}
#endif // ! USEGSL

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaMcos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaMsin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sincos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaMsincos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cossin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaMcossin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// PXD function
int saveHistZetaM_sincos(struct  cmdline_data* cmd,
                         struct  global_data* gd,
                         int mChebyshev, int sizeHistN,
                         real ***histZetaMcos)
{
    int n1, n2, m;

    verb_print_debug(1, "\nAqui voy (0)\n");

//    for (m=1; m<=mChebyshev+1; m++) {
    for (m=1; m<=mChebyshev; m++) {
        printf("%d\n",m);
        for (n1=1; n1<=sizeHistN; n1++) {
            printf("%d\n",n1);
            for (n2=1; n2<=sizeHistN; n2++) {
                printf("%d %g\n",n2, gd->histZetaMcos[m][n1][n2]);
//                histZetaMcos[m-1][n1-1][n2-1] = gd->histZetaMcos[m][n1][n2];
            }
        }
    }

    verb_print_debug(1, "\nAqui voy (1)\n");

    return SUCCESS;
}

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

// Seems obsolete. Delete. Only appears in balls_method addon
local int PrintHistZeta_theta2_fix(struct  cmdline_data* cmd,
                                   struct  global_data* gd)
{
    real rBin, rbinlog;
    int n1, m;
    stream outstr;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    int Nbins;
    char namebuf[256];
    Nbins = cmd->sizeHistN;
#ifdef USEGSL
    gsl_complex Atmp;                               // Use in cballs_04.h
#endif

#ifdef ADDONS
#include "cballs_include_05.h"
#endif

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s_%d%s", gd->fpfnamemhistZetaMFileName, m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta = gd->histZetaM[m][n1][n1];
            Zeta2 = gd->histZetaM[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaM[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaM[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaM[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}
//E

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos(struct  cmdline_data* cmd,
                                struct  global_data* gd)
{
    real rBin, rbinlog;
    int n1, m;
    stream outstr;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    int Nbins;
    char namebuf[256];

    Nbins = cmd->sizeHistN;

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta = gd->histZetaMcos[m][n1][n1];
            Zeta2 = gd->histZetaMcos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                    + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
                Zeta = gd->histZetaMsin[m][n1][n1];
                Zeta2 = gd->histZetaMsin[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gd->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gd->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gd->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                        + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta = gd->histZetaMsincos[m][n1][n1];
            Zeta2 = gd->histZetaMsincos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cossin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            if (cmd->useLogHist) {
                if (cmd->rminHist==0) {
                    rbinlog = ((real)(n1-cmd->sizeHistN))/cmd->logHistBinsPD
                            + rlog10(cmd->rangeN);
                } else {
                    rbinlog = rlog10(cmd->rminHist) + ((real)(n1)-0.5)*gd->deltaR;
                }
                rBin=rpow(10.0,rbinlog);
            } else {
                rBin = cmd->rminHist + ((real)n1-0.5)*gd->deltaR;
            }
            Zeta = gd->histZetaMcossin[m][n1][n1];
            Zeta2 = gd->histZetaMcossin[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, full correlation function at each phi bins
local int PrintHistZetaG(struct  cmdline_data* cmd,
                                       struct  global_data* gd)
{
    int n1, n2, l;
    stream outstr;
    char namebuf[256];

    fclose(outstr);
    for (l=1; l<=cmd->sizeHistTheta; l++) {
        sprintf(namebuf, "%s_%d%s", gd->fpfnamehistZetaGFileName,
                l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histXi3pcf[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM as obtained from full 3pcf (bf) ZetaG, for each m multipole
local int PrintHistZetaGm_sincos(struct  cmdline_data* cmd,
                                 struct  global_data* gd)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Re", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaGmRe[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Im", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaGmIm[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}


// Saves matrix ZetaG, real and imaginary parts, obtained from ZetaM multipoles
//  also saves full 3pcf ZetaG matrix for each phi bins obtained from inverse FFT
local int PrintHistZetaMZetaGm_sincos(struct  cmdline_data* cmd,
                                     struct  global_data* gd)
{
    int n1, n2, m, l;
    stream outstr;
    char namebuf[256];

    int NP = 2*(cmd->mChebyshev+1);
    double ***histZetaG;
    double ***histZetaG_Im;
    histZetaG = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    histZetaG_Im = dmatrix3D(1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

#ifdef USEGSL
    double *data;
    gsl_fft_real_wavetable * real;
    gsl_fft_real_workspace * work;
    gsl_fft_halfcomplex_wavetable * hc;

    data=dvector(0,NP-1);
    work = gsl_fft_real_workspace_alloc (NP);
    real = gsl_fft_real_wavetable_alloc (NP);
    hc = gsl_fft_halfcomplex_wavetable_alloc (NP);
#else
    double *data;
    data=dvector(1,NP);
#endif

    //B Sum cos^2 + sin^2 and sincos - sincos
    // mchebyshev + 1 < sizeHistTheta/2
    // and mchebyshev + 1 must be a power of 2 also
    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                gd->histZetaGmRe[m][n1][n2] =
                        gd->histZetaMcos[m][n1][n2]
                        + gd->histZetaMsin[m][n1][n2];
                gd->histZetaGmIm[m][n1][n2] =
                        gd->histZetaMsincos[m][n1][n2]
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                        - gd->histZetaMcossin[m][n1][n2];
            }
#ifdef USEGSL
            for (m=0; m<cmd->mChebyshev+1; m++) {
                data[2*m] = gd->histZetaGmRe[m+1][n1][n2];
                data[2*m+1] = gd->histZetaGmIm[m+1][n1][n2];
            }

            gsl_fft_complex_radix2_inverse (data, 1, NP/2);
            for (l=0; l<NP; l++) {                 // l denote angular separation
                histZetaG[l+1][n1][n2] = data[l];
            }
#else
            int isign = -1;                         // sign in imaginary unit
            for (m=1; m<=cmd->mChebyshev+1; m++) {
                data[2*m-1] = gd->histZetaGmRe[m][n1][n2];
                data[2*m] = gd->histZetaGmIm[m][n1][n2];
            }
            dfour1(data,NP/2,isign);                // Inverse Fourier transform
                                                    // data has Re and Im parts
            for (l=1; l<=NP; l++) {                 // l denote angular separation
                histZetaG[l][n1][n2] = (2.0/(double)NP)*data[l];
            }
#endif
        }
    }
    //E

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Re", m, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaGmRe[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGmFileName,
                "_Im", m, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->histZetaGmIm[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    for (l=1; l<=cmd->mChebyshev+1; l++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                "_fftinv_Re",l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[2*l-1][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (l=1; l<=cmd->mChebyshev+1; l++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaGFileName,
                "_fftinv_Im",l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[2*l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    //B Sum cos^2 + sin^2 and sincos - sincos
    // mchebyshev + 1 < sizeHistTheta/2
    // and mchebyshev + 1 must be a power of 2 also
    double deltaTheta = TWOPI/((double)NP);
    for (n1=1; n1<=cmd->sizeHistN; n1++) {
        for (n2=1; n2<=cmd->sizeHistN; n2++) {
            for (l=1; l<=NP; l++) {              // l denote angular sep.
                histZetaG[l][n1][n2] = gd->histZetaMcos[1][n1][n2]
                                        + gd->histZetaMsin[1][n1][n2];
                histZetaG_Im[l][n1][n2] = gd->histZetaMsincos[1][n1][n2]
                // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                                            - gd->histZetaMcossin[1][n1][n2];
                for (m=2; m<=cmd->mChebyshev+1; m++) {
                    histZetaG[l][n1][n2] += 2.0*(gd->histZetaMcos[m][n1][n2]
                                            + gd->histZetaMsin[m][n1][n2])
                                            *rcos(((double)(m*l))*deltaTheta);
                    histZetaG_Im[l][n1][n2] +=
                                    2.0*(gd->histZetaMsincos[m][n1][n2]
                    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
                                    - gd->histZetaMcossin[m][n1][n2])
                                    *rcos(((double)(m*l))*deltaTheta);
                }
            }
        }
    }
    //E
    for (l=1; l<=NP; l++) {
        sprintf(namebuf, "%s_%s_%d%s",
                gd->fpfnamehistZetaGFileName, "Re", l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (l=1; l<=NP; l++) {
        sprintf(namebuf, "%s_%s_%d%s",
                gd->fpfnamehistZetaGFileName, "Im", l, EXTFILES);
        verb_print_q(2, cmd->verbose,
                    "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",histZetaG_Im[l][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

#ifdef USEGSL
    gsl_fft_halfcomplex_wavetable_free (hc);
    gsl_fft_real_workspace_free (work);
    gsl_fft_real_wavetable_free (real);
    free_dvector(data,0,NP-1);
#else
    free_dvector(data,1,NP);
#endif

    free_dmatrix3D(histZetaG_Im,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dmatrix3D(histZetaG,1,NP,1,cmd->sizeHistN,1,cmd->sizeHistN);

    return SUCCESS;
}

#undef MHISTZETA

