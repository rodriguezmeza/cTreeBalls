// Use:
//#include "cballs_include_05_NMultipoles.h"

#ifndef _cballs_include_05_NMultipoles_h
#define _cballs_include_05_NMultipoles_h

#ifdef NMultipoles
// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_N(struct  cmdline_data* cmd,
                                struct  global_data* gd)
{
    int n1, n2, m;
    stream outstr;
    char namebuf[256];

    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->NhistZetaMcos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sin_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->NhistZetaMsin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_sincos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->NhistZetaMsincos[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }
    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m=1; m<=cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamehistZetaMFileName,
                "_cossin_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        for (n1=1; n1<=cmd->sizeHistN; n1++) {
            for (n2=1; n2<=cmd->sizeHistN; n2++) {
                fprintf(outstr,"%16.8e ",gd->NhistZetaMcossin[m][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_N(struct  cmdline_data* cmd,
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
                "_cos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta = gd->NhistZetaMcos[m][n1][n1];
            Zeta2 = gd->NhistZetaMcos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->NhistZetaMcos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->NhistZetaMcos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->NhistZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sin_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
                Zeta = gd->NhistZetaMsin[m][n1][n1];
                Zeta2 = gd->NhistZetaMsin[m][n1][(int)(Nbins/4.0)];
                Zeta3 = gd->NhistZetaMsin[m][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gd->NhistZetaMsin[m][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gd->NhistZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos_N", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta = gd->NhistZetaMsincos[m][n1][n1];
            Zeta2 = gd->NhistZetaMsincos[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->NhistZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->NhistZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->NhistZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cossin_", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta = gd->NhistZetaMcossin[m][n1][n1];
            Zeta2 = gd->NhistZetaMcossin[m][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->NhistZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->NhistZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->NhistZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}
//#endif


//B Additional definitions... delete ASAP

//#ifdef NMultipoles
local int PrintHistXi2pcf_N(struct  cmdline_data* cmd, struct  global_data* gd)
{
    real rBin, rbinlog;
    int n;
    stream outstr;
    char namebuf[256];

    sprintf(namebuf, "%s%s%s%s", gd->fpfnamehistXi2pcfFileName,
            "_N", cmd->suffixOutFiles, EXTFILES);
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
        fprintf(outstr,"%16.8e %16.8e\n",rBin,gd->NhistXi2pcf[n]);
    }
    fclose(outstr);

    return SUCCESS;
}


#ifdef NONORMHIST

// Saves matrix ZetaM for each m multipole
local int PrintHistZetaM_sincos_normalized(struct  cmdline_data* cmd,
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
                fprintf(outstr,"%16.8e ",
                        gd->histZetaMcos[m][n1][n2]/gd->NhistZetaMcos[1][n1][n2]);
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
                fprintf(outstr,"%16.8e ",
                        gd->histZetaMsin[m][n1][n2]/gd->NhistZetaMcos[1][n1][n2]);
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
                fprintf(outstr,"%16.8e ",
                        gd->histZetaMsincos[m][n1][n2]/gd->NhistZetaMcos[1][n1][n2]);
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
                fprintf(outstr,"%16.8e ",
                        gd->histZetaMcossin[m][n1][n2]/gd->NhistZetaMcos[1][n1][n2]);
            }
            fprintf(outstr,"\n");
        }
        fclose(outstr);
    }

    return SUCCESS;
}

// Saves matrix ZetaM for each m multipole at a set of theta2 angles
local int PrintHistZetaMm_sincos_normalized(struct  cmdline_data* cmd,
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

    real Norm;

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_cos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta  = gd->histZetaMcos[m][n1][n1]
                  /gd->NhistZetaMcos[1][n1][n1];
            Zeta2 = gd->histZetaMcos[m][n1][(int)(Nbins/4.0)]
                  /gd->NhistZetaMcos[1][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMcos[m][n1][(int)(2.0*Nbins/4.0)]
                  /gd->NhistZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMcos[m][n1][(int)(3.0*Nbins/4.0)]
                  /gd->NhistZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMcos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                  /gd->NhistZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }
        
    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sin", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
                Zeta = gd->histZetaMsin[m][n1][n1]
                     /gd->NhistZetaMcos[1][n1][n1];
                Zeta2 = gd->histZetaMsin[m][n1][(int)(Nbins/4.0)]
                      /gd->NhistZetaMcos[1][n1][(int)(Nbins/4.0)];
                Zeta3 = gd->histZetaMsin[m][n1][(int)(2.0*Nbins/4.0)]
                      /gd->NhistZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
                Zeta4 = gd->histZetaMsin[m][n1][(int)(3.0*Nbins/4.0)]
                      /gd->NhistZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
                Zeta5 = gd->histZetaMsin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                      /gd->NhistZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
                fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    for (m = 1; m <= cmd->mChebyshev+1; m++) {
        sprintf(namebuf, "%s%s_%d%s", gd->fpfnamemhistZetaMFileName,
                "_sincos", m, EXTFILES);
        verb_print_q(2, cmd->verbose, "Printing : to a file %s ...\n",namebuf);
        outstr = stropen(namebuf, "w!");
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta = gd->histZetaMsincos[m][n1][n1]
                    /gd->NhistZetaMcos[1][n1][n1];
            Zeta2 = gd->histZetaMsincos[m][n1][(int)(Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMsincos[m][n1][(int)(2.0*Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMsincos[m][n1][(int)(3.0*Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMsincos[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
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
        fprintf(outstr,MHISTZETAHEADER);
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
            Zeta = gd->histZetaMcossin[m][n1][n1]
                    /gd->NhistZetaMcos[1][n1][n1];
            Zeta2 = gd->histZetaMcossin[m][n1][(int)(Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(Nbins/4.0)];
            Zeta3 = gd->histZetaMcossin[m][n1][(int)(2.0*Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(2.0*Nbins/4.0)];
            Zeta4 = gd->histZetaMcossin[m][n1][(int)(3.0*Nbins/4.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(3.0*Nbins/4.0)];
            Zeta5 = gd->histZetaMcossin[m][n1][(int)(4.0*Nbins/4.0 - 1.0)]
                     /gd->NhistZetaMcos[1][n1][(int)(4.0*Nbins/4.0 - 1.0)];
            fprintf(outstr,MHISTZETA,rBin,Zeta,Zeta2,Zeta3,Zeta4,Zeta5);
        }
        fclose(outstr);
    }

    return SUCCESS;
}


#endif // ! NONORMHIST

#endif // ! NMultipoles


/*
 Add your addon item here
 */

#endif	// ! _cballs_include_05_NMultipoles_h
