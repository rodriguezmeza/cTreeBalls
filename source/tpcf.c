/*==============================================================================
 MODULE: tpcf.c				[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose:
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/

#include "globaldefs.h"

local int printHistN(void);
local int printHistCF(void);
local int printHistXi2pcf(void);
local int printHistXi3pcf(void);
#ifdef TPCF
local int printHistZetaM(void);
local int printHistZetaM_exp(void);
local int printHistZetaM_sincos(void);
local int printHistZeta(void);
local int printHistZeta_sincos(void);
#endif
local int evalHist(void);
local int printEvalHist(void);


int MainLoop(void)
{
    bodyptr p,q;
    real kavg;


    gd.flagSmooth = FALSE;
    gd.flagSetNbNoSel = FALSE;

    if (scanopt(cmd.options, "smooth") && !scanopt(cmd.options, "set-Nb-noSel")) {
        verb_print(cmd.verbose, "\n\tMainLoop: smooth: try making tree...\n\n");
        DO_BODY(p,bodytab,bodytab+cmd.nbody)
            Update(p) = TRUE;
        maketree(bodytab, cmd.nbody);
//B
        free(bodytab);
        gd.bytes_tot -= cmd.nbody*sizeof(body);
        cmd.nbody = gd.nbodysm;
        bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
        gd.bytes_tot += cmd.nbody*sizeof(body);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for final (smooth) particle (%ld) storage.\n",
                   cmd.nbody*sizeof(body)*INMB, cmd.nbody);
        kavg = 0;
        q = bodytabsm;
        DO_BODY(p, bodytab, bodytab+cmd.nbody) {
            SETV(Pos(p),Pos(q));
// BODY3
//            Type(q) = BODY;
            Type(p) = Type(q);
            Nbb(p) = Nbb(q);
//            verb_print_debug(1, "\nAqui voy (0):: %ld %d\n", Type(p), Nbb(p));
            Weight(p) = Weight(q);
            Kappa(p) = Kappa(q);
            Id(p) = p-bodytab+1;
            Update(p) = TRUE;
            q++;
            kavg += Kappa(p);
        }
        free(bodytabsm);
        gd.bytes_tot -= gd.nbodysm*sizeof(body);
//E
        verb_print(cmd.verbose, "smooth: Average of kappa (%ld particles) = %le\n",
                   cmd.nbody, kavg/((real)cmd.nbody) );
        gd.flagSmooth = TRUE;
    }

    if ( scanopt(cmd.options, "smooth") && scanopt(cmd.options, "set-Nb-noSel") ) {
    verb_print(cmd.verbose, "\n\tMainLoop: smooth & set-Nb-noSel: try making tree...\n\n");
        DO_BODY(p,bodytab,bodytab+cmd.nbody)
            Update(p) = TRUE;
        maketree(bodytab, cmd.nbody);
//B
        free(bodytab);
        gd.bytes_tot -= cmd.nbody*sizeof(body);
        cmd.nbody = gd.nbodySel;
        bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
        gd.bytes_tot += cmd.nbody*sizeof(body);
        verb_print(cmd.verbose,
                   "Allocated %g MByte for final (smooth-noSel) particle (%ld) storage.\n",
                   cmd.nbody*sizeof(body)*INMB, cmd.nbody);
        kavg = 0;
        q = bodytabSel;
        DO_BODY(p, bodytab, bodytab+cmd.nbody) {
            SETV(Pos(p),Pos(q));
// BODY3
//            Type(q) = BODY;
            Type(p) = Type(q);
            Nbb(p) = Nbb(q);
            Weight(p) = Weight(q);
            Kappa(p) = Kappa(q);
            Id(p) = p-bodytab+1;
            Update(p) = TRUE;
            q++;
            kavg += Kappa(p);
        }
        free(bodytabSel);
        gd.bytes_tot -= gd.nbodySel*sizeof(body);
//E
        verb_print(cmd.verbose,
                   "smooth-set-Nb-noSel: Average of kappa (%ld particles) = %le\n",
                   cmd.nbody, kavg/((real)cmd.nbody) );
        gd.flagSmooth = TRUE;
        gd.flagSetNbNoSel = TRUE;
    }

    if (scanopt(cmd.options, "stop")) {
        if (!strnull(cmd.outfile))
            output();
        verb_print(cmd.verbose, "\n\tMainLoop: stopping...\n\n");
        exit(1);
    }

// To see the hists::
#ifdef DEBUG
//    DO_BODY(p,bodytab,bodytab+cmd.nbody)
//        Update(p) = TRUE;
//    maketree(bodytab, cmd.nbody);
#endif
//    if (!strnull(cmd.outfile))
//        output();

    evalHist();

    if (!gd.stopflag)
        printEvalHist();

#ifdef DEBUG
    if (!strnull(cmd.outfile))
        output();
#endif

    return _SUCCESS_;
}

local int evalHist(void)
{
    bodyptr p;

    switch(gd.searchMethod_int) {
        case DIRECT3PCFOMP:             // search=direct-3pcf-omp
            verb_print(cmd.verbose, "\n\tevalHist: direct method simple (3pcf-omp)\n\n");
            search_direct_3pcf_omp(bodytab, cmd.nbody, 1, cmd.nbody); break;
        case TREEOMPMETHOD:         // search=tree-omp
            verb_print(cmd.verbose, "\n\tevalHist: with normal tree method (omp)\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
//            verb_print_debug(1, "\nAqui voy (0)\n");
            maketree(bodytab, cmd.nbody);
//            if ( scanopt(cmd.options, "smooth") )
//                free(bodytabsm);
//            if ( scanopt(cmd.options, "set-Nb-noSel") )
//                free(bodytabSel);
//            verb_print_debug(1, "\nAqui voy (1)\n");
            searchcalc_normal_omp(bodytab, cmd.nbody, 1, cmd.nbody);
            break;
        case TREE3PCFBFOMPMETHOD:         // search=tree-3pcf-direct
            verb_print(cmd.verbose, "\n\tevalHist: with normal tree method (tree-3pcf-direct-omp)\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
            maketree(bodytab, cmd.nbody);
        searchcalc_normal_3pcf_direct_omp(bodytab, cmd.nbody, 1, cmd.nbody);
            break;

#ifdef BALLS
        case BALLSOMPMETHOD:
            verb_print(cmd.verbose, "\n\tevalHist: with balls tree-omp method\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
            maketree(bodytab, cmd.nbody);
            searchcalc_balls_omp(bodytab, cmd.nbody, 1, cmd.nbody);
            break;
#endif

        case TREEOMPMETHODSINCOS:   // search=tree-omp-sincos
            verb_print(cmd.verbose, "\n\tevalHist: with normal tree method (sincos-omp)\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
            maketree(bodytab, cmd.nbody);
            searchcalc_normal_omp_sincos(bodytab, cmd.nbody, 1, cmd.nbody);
            break;

        case SEARCHNULL:
            verb_print(cmd.verbose, "\n\tevalHist: null search method.\n");
            verb_print(cmd.verbose, "\tevalHist: with normal tree method (omp)\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
            maketree(bodytab, cmd.nbody);
            searchcalc_normal_omp(bodytab, cmd.nbody, 1, cmd.nbody);
            break;
        default:
            verb_print(cmd.verbose, "\n\tevalHist: dafault search method.\n");
            verb_print(cmd.verbose, "\tevalHist: with normal tree method (omp)\n\n");
            DO_BODY(p,bodytab,bodytab+cmd.nbody)
                Update(p) = TRUE;
            maketree(bodytab, cmd.nbody);
            searchcalc_normal_omp(bodytab, cmd.nbody, 1, cmd.nbody);
            break;



    }

    return _SUCCESS_;
}

local int printEvalHist(void)
{
#ifdef MPICODE
    if (ThisTask==0) {
#endif
        double cpustart = CPUTIME;
        if (!scanopt(cmd.options, "no-out-Hist")) {
    switch(gd.searchMethod_int) {
        case DIRECT3PCFOMP:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing direct method simple\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
            printHistXi3pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif
            break;
        case TREEOMPMETHOD:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing normal tree method (omp)\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif
            break;
        case TREE3PCFBFOMPMETHOD:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing normal tree method (tree-3pcf-direct-omp)\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
        #ifdef TPCF
            printHistZetaM();
            printHistZeta();
        #endif
            break;

#ifdef BALLS
        case BALLSOMPMETHOD:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing  balls tree-omp method\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif
            break;
#endif

        case TREEOMPMETHODSINCOS:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing normal tree method (omp-sincos)\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
        #ifdef TPCF
            printHistZetaM_sincos();
            printHistZeta_sincos();
        #endif
            break;

        case SEARCHNULL:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing null search method.\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif
        default:
            verb_print(cmd.verbose, "\n\tprintEvalHist: printing dafault search method.\n\n");
            if (scanopt(cmd.options, "compute-HistN")) printHistN();
            printHistXi2pcf();
#ifdef TPCF
            printHistZetaM();
            printHistZeta();
#endif



            }
        }

        gd.cputotalinout += CPUTIME - cpustart;

#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local int printHistN(void)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    outstr = stropen(gd.fpfnamehistNFileName, "w!");

    verb_print(cmd.verbose, "Printing : to a file %s ...\n",gd.fpfnamehistNFileName);

    for (n=1; n<=cmd.sizeHistN; n++) {
#ifdef LOGHIST
        if (cmd.rminHist==0) {
            rbinlog = ((real)(n-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
        } else {
            rbinlog = rlog10(cmd.rminHist) + ((real)(n))*gd.deltaR;
        }
        rBin=rpow(10.0,rbinlog);
#else
        rBin = cmd.rminHist + ((real)n-0.5)*gd.deltaR;
#endif
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd.histN[n]);
    }
    fclose(outstr);
#ifdef MPICODE
    }
#endif

    if (scanopt(cmd.options, "and-CF"))
        printHistCF();

    return _SUCCESS_;
}

local int printHistCF(void)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    outstr = stropen(gd.fpfnamehistCFFileName, "w!");

    verb_print(cmd.verbose, "Printing : to a file %s ...\n",gd.fpfnamehistCFFileName);

    for (n=1; n<=cmd.sizeHistN; n++) {
#ifdef LOGHIST
        if (cmd.rminHist==0) {
            rbinlog = ((real)(n-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
        } else {
            rbinlog = rlog10(cmd.rminHist) + ((real)(n))*gd.deltaR;
        }
        rBin=rpow(10.0,rbinlog);
#else
        rBin = cmd.rminHist + ((real)n-0.5)*gd.deltaR;
#endif
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd.histCF[n]);
    }
    fclose(outstr);
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local int printHistXi2pcf(void)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    outstr = stropen(gd.fpfnamehistXi2pcfFileName, "w!");

    verb_print(cmd.verbose, "Printing : to a file %s ...\n",gd.fpfnamehistXi2pcfFileName);

    for (n=1; n<=cmd.sizeHistN; n++) {
#ifdef LOGHIST
        if (cmd.rminHist==0) {
            rbinlog = ((real)(n-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
        } else {
            rbinlog = rlog10(cmd.rminHist) + ((real)(n))*gd.deltaR;
        }
        rBin=rpow(10.0,rbinlog);
#else
        rBin = cmd.rminHist + ((real)n-0.5)*gd.deltaR;
#endif
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd.histXi2pcf[n]);
    }
    fclose(outstr);
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local int printHistXi3pcf(void)
{
    real rBin, rbinlog;
    int n;
    stream outstr;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    outstr = stropen(gd.fpfnamehistXi2pcfFileName, "w!");

    verb_print(cmd.verbose, "Printing : to a file %s ...\n",gd.fpfnamehistXi2pcfFileName);

    for (n=1; n<=cmd.sizeHistN; n++) {
#ifdef LOGHIST
        rbinlog = ((real)(n-0.5))*gd.deltaR;
        rBin=rpow(10,rbinlog);
#else
        rBin = ((int)n-0.5)*gd.deltaR;
#endif
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd.histXi2pcf[n]);
    }
    fclose(outstr);
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}



#ifdef TPCF

local int printHistZetaM(void)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];
    real rBin1, rBin2;
    real rbinlog1, rbinlog2;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s_%d%s", gd.fpfnamehistZetaMFileName, m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            for (n2=1; n2<=cmd.sizeHistN; n2++) {
                if (scanopt(cmd.options, "xr1r2")) {
#ifdef LOGHIST
//                    rbinlog1 = ((real)(n1-0.5))*gd.deltaR;
//                    rBin1=rpow(10,rbinlog1);
//                    rbinlog2 = ((real)(n2-0.5))*gd.deltaR;
//                    rBin2=rpow(10,rbinlog2);
                    if (cmd.rminHist==0) {
                        rbinlog1 = ((real)(n1-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
                        rbinlog2 = ((real)(n2-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
                    } else {
                        rbinlog1 = rlog10(cmd.rminHist) + ((real)(n1))*gd.deltaR;
                        rbinlog2 = rlog10(cmd.rminHist) + ((real)(n2))*gd.deltaR;
                    }
                    rBin1=rpow(10.0,rbinlog1);
                    rBin2=rpow(10.0,rbinlog2);
#else
//                    rBin1 = ((int)n1-0.5)*gd.deltaR;
//                    rBin2 = ((int)n2-0.5)*gd.deltaR;
                    rBin1 = cmd.rminHist + ((real)n1-0.5)*gd.deltaR;
                    rBin2 = cmd.rminHist + ((real)n2-0.5)*gd.deltaR;
#endif
//                if (scanopt(cmd.options, "xr1r2")) {
                    fprintf(outstr,"%16.8e ",rBin1*rBin2*gd.histZetaM[m][n1][n2]);
                } else {
                    fprintf(outstr,"%16.8e ",gd.histZetaM[m][n1][n2]);
                }
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local int printHistZetaM_exp(void)
{
    int n1, n2, m;
    stream outstr;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    outstr = stropen(gd.fpfnamehistZetaMFileName, "w!");
    verb_print(cmd.verbose, "Printing : to a file %s ...\n",gd.fpfnamehistZetaMFileName);

    gsl_complex Atmp;

    for (m=0; m<=cmd.mchebyshev; m++) {
        fprintf(outstr,"\n\nm=%d\n",m);
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            for (n2=1; n2<=cmd.sizeHistN; n2++) {
                Atmp = gsl_matrix_complex_get(histZetaMatrix[m].histZetaM,n1,n1);
//                fprintf(outstr,"%16.8e ",Atmp);
                fprintf(outstr,"%16.8e %16.8e ",GSL_REAL(Atmp),GSL_IMAG(Atmp));
            }
            fprintf(outstr,"\n");
        }
    }
        fclose(outstr);
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local int printHistZetaM_sincos(void)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd.fpfnamehistZetaMFileName, "_cos", m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            for (n2=1; n2<=cmd.sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd.histZetaMcos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd.fpfnamehistZetaMFileName, "_sin", m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            for (n2=1; n2<=cmd.sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd.histZetaMsin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd.fpfnamehistZetaMFileName, "_sincos", m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            for (n2=1; n2<=cmd.sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd.histZetaMsincos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

#define MHISTZETA \
"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n"

local int printHistZeta(void)
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

    m = cmd.mToPlot;
    Nbins = cmd.sizeHistN;
    
    gsl_complex Atmp;

#ifdef MPICODE
    if (ThisTask==0) {
#endif


    for (m = 1; m <= cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s_%d%s", gd.fpfnamemhistZetaFileName, m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
#ifdef LOGHIST
//            rbinlog = ((real)(n1-0.5))*gd.deltaR;
//            rBin=rpow(10,rbinlog);
            if (cmd.rminHist==0) {
                rbinlog = ((real)(n1-cmd.sizeHistN))/NLOGBINPD + rlog10(cmd.rangeN);
            } else {
                rbinlog = rlog10(cmd.rminHist) + ((real)(n1))*gd.deltaR;
            }
            rBin=rpow(10.0,rbinlog);
#else
//          rBin = ((int)n1-0.5)*gd.deltaR;
            rBin = cmd.rminHist + ((real)n1-0.5)*gd.deltaR;
#endif
            Zeta = gd.histZetaM[m][n1][n1];
            Zeta2 = gd.histZetaM[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd.histZetaM[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd.histZetaM[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd.histZetaM[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

// Correct this to have logscale
local int printHistZeta_sincos(void)
{
    real rBin;
    int n1, m;
    stream outstr;
    real Zeta;
    real Zeta2;
    real Zeta3;
    real Zeta4;
    real Zeta5;
    int Nbins;
    char namebuf[256];

    Nbins = cmd.sizeHistN;

#ifdef MPICODE
    if (ThisTask==0) {
#endif
    for (m = 1; m <= cmd.mchebyshev+1; m++) {
        sprintf(namebuf, "%s_%d%s", gd.fpfnamemhistZetaFileName, m, EXTFILES);
        verb_print(cmd.verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd.sizeHistN; n1++) {
            rBin = ((int)n1-0.5)*cmd.rangeN/cmd.sizeHistN;
            Zeta = gd.histZetaMcos[m][n1][n1];
            Zeta2 = gd.histZetaMcos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd.histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd.histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd.histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

#undef MHISTZETA
#endif // ! TPCF
