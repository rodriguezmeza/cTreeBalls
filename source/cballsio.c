/*==============================================================================
 MODULE: cballsio.c		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: Routines to drive input and output data
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

local int inputdata_ascii(struct cmdline_data*, struct  global_data*,
                          string filename, int);
local int inputdata_bin(struct cmdline_data*, struct  global_data*,
                        string filename, int);
local int inputdata_takahasi(struct cmdline_data*, struct  global_data*,
                             string filename, int);

local int outputdata(struct cmdline_data*, struct  global_data*,
                     bodyptr, INTEGER nbody);
local int outputdata_ascii(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
local int outputdata_bin(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
local int EndRun_FreeMemory(struct cmdline_data*, struct  global_data*);

#ifdef ADDONS
#include "cballsio_include_00.h"
#endif

local int outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;

int InputData(struct cmdline_data* cmd, 
              struct  global_data* gd, string filename, int ifile)
{
    double cpustart = CPUTIME;

    switch(gd->infilefmt_int) {
        case INCOLUMNS:
            verb_print(cmd->verbose,"\n\tInput in columns (ascii) format...\n");
            inputdata_ascii(cmd, gd, filename, ifile); break;
        case INNULL:
            verb_print(cmd->verbose,"\n\t(Null) Input in columns (ascii) format...\n");
            inputdata_ascii(cmd, gd, filename, ifile); break;
        case INCOLUMNSBIN:
            verb_print(cmd->verbose,"\n\tInput in binary format...\n");
            inputdata_bin(cmd, gd, filename, ifile); break;
        case INTAKAHASI:
            verb_print(cmd->verbose,"\n\tInput in takahasi format...\n");
            class_call_cballs(inputdata_takahasi(cmd, gd, filename, ifile),
                       errmsg, errmsg);
            break;

#ifdef ADDONS
#include "cballsio_include_01.h"
#endif

        default:
            verb_print(cmd->verbose,
                       "\n\tInput: Unknown input format (%s)...",cmd->infilefmt);
            verb_print(cmd->verbose,
                       "\n\tInput in default columns (ascii) format...\n");
            class_call_cballs(inputdata_ascii(cmd, gd, filename, ifile),
                                              errmsg, errmsg);
            break;
    }

    gd->cputotalinout += CPUTIME - cpustart;
    verb_print(cmd->verbose, "\n\tinputdata :: reading time = %f\n",
               CPUTIME - cpustart);

    return SUCCESS;
}

local int inputdata_ascii(struct cmdline_data* cmd, struct  global_data* gd,
                           string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real weight=1;

    gd->input_comment = "Column form input file";

    instr = stropen(filename, "r");

    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd->nbody);
    if (cmd->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd->nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);

    gd->nbodyTable[ifile] = cmd->nbody;

// Check the center of the box!!!
#if NDIM == 3
    real Lx, Ly, Lz;
#ifdef SINGLEP
    in_real_double(instr, &Lx);
    in_real_double(instr, &Ly);
    in_real_double(instr, &Lz);
#else
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    in_real(instr, &Lz);
#endif
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
    gd->Box[2] = Lz;
#else
    real Lx, Ly;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
#endif

    verb_print(cmd->verbose, "\tInput: nbody and ndim: %d %d...\n", cmd->nbody, ndim);

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_vector(instr, Pos(p));
        in_real(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
    }

    fclose(instr);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
        kavg += Kappa(p);
    }
    kavg /= ((real)cmd->nbody);
    real kstd;
    real sum=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        sum += rsqr(Kappa(p) - kavg);
    }
    kstd = rsqrt( sum/((real)cmd->nbody - 1.0) );
    verb_print(cmd->verbose,
               "inputdata_ascii: average and std dev of kappa ");
    verb_print(cmd->verbose,
               "(%ld particles) = %le %le\n", cmd->nbody, kavg, kstd);

//B Locate particles with same position
    if (scanopt(cmd->options, "check-eq-pos")) {
    bodyptr q;
    real dist2;
    vector distv;
    bool flag=0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody-1)
        DO_BODY(q, p+1, bodytable[ifile]+cmd->nbody)
            if (p != q) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                if (dist2 == 0.0) {
                    flag=1;
                }
            }
    if (flag) error("inputdata_ascii: at least two bodies have same position\n");
    }
//E

    return SUCCESS;
}


local int inputdata_bin(struct cmdline_data* cmd, struct  global_data* gd,
                         string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    real weight=1;

    gd->input_comment = "Binary input file";

    instr = stropen(filename, "r");
    in_int_bin_long(instr, &cmd->nbody);
    verb_print(cmd->verbose, "\tInput: nbody %d\n", cmd->nbody);
    if (cmd->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd->nbody);
    in_int_bin(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    verb_print(cmd->verbose, "\tInput: nbody and ndim: %d %d...\n", cmd->nbody, ndim);

#ifdef SINGLEP
    in_real_bin_double(instr, &gd->Box[0]);
    in_real_bin_double(instr, &gd->Box[1]);
#if NDIM == 3
    in_real_bin_double(instr, &gd->Box[2]);
#endif
#else
    in_real_bin(instr, &gd->Box[0]);
    in_real_bin(instr, &gd->Box[1]);
#if NDIM == 3
    in_real_bin(instr, &gd->Box[2]);
#endif
#endif

#if NDIM == 3
    verb_print(cmd->verbose, "\tInput: Box: %g %g %g\n", 
               gd->Box[0], gd->Box[1], gd->Box[2]);
#else
    verb_print(cmd->verbose, "\tInput: Box: %g %g %g\n", gd->Box[0], gd->Box[1]);
#endif

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        in_vector_bin(instr, Pos(p));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_real_bin(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
    }
    fclose(instr);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
    }

    return SUCCESS;
}


//B BEGIN:: Reading Takahasi simulations
//From Takahasi web page. Adapted to our needs

#include<math.h>
#include<stdio.h>
#include<stdlib.h>


void pix2ang(long pix, int nside, double *theta, double *phi);

local int Takahasi_region_selection(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat, int ifile);
local int Takahasi_region_selection_3d_all(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax,
            real *zmin, real *zmax, int ifile);
local int Takahasi_region_selection_3d(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax,
            real *zmin, real *zmax, int ifile);

#if NDIM == 2
local int Takahasi_region_selection_2d(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax, int ifile);
#endif

local int inputdata_takahasi(struct cmdline_data* cmd, struct  global_data* gd,
                             string filename, int ifile)
{
    FILE *fp;
    long i,j,npix,dummy;
    long jj[6]={536870908,1073741818,1610612728,2147483638,2684354547,3221225457};
    int negi,nside;
//    double theta,phi;
//    char file[200];

    gd->input_comment = "Takahasi input file";

//E Begin reading Takahasi file
    fp = stropen(filename, "rb");

    fread(&negi, sizeof(int), 1, fp);
    fread(&nside, sizeof(int), 1, fp);
    fread(&npix, sizeof(long), 1, fp);
    fread(&dummy, sizeof(long), 1, fp);

    float(*conv);     // convergence
    conv=(float *)malloc(sizeof(float)*npix);
    float(*shear1);   // shear 1
    shear1=(float *)malloc(sizeof(float)*npix);
    float(*shear2);   // shear 2
    shear2=(float *)malloc(sizeof(float)*npix);
    float(*rotat);    // rotation
    rotat=(float *)malloc(sizeof(float)*npix);

    verb_print(cmd->verbose, "\nAllocated %g MByte for pixel storage.\n",
               sizeof(float)*npix*4/(1024.0*1024.0));

    for(j=0;j<npix;j++){
      fread(&conv[j], sizeof(float), 1, fp);
      for(i=0;i<6;i++) if(j==jj[i]) fread(&dummy, sizeof(long), 1, fp);
    }

    fread(&dummy, sizeof(long), 1, fp);
    for(j=0;j<npix;j++){
      fread(&shear1[j], sizeof(float), 1, fp);
      for(i=0;i<6;i++) if(j==jj[i]) fread(&dummy, sizeof(long), 1, fp);
    }

    fread(&dummy, sizeof(long), 1, fp);
    for(j=0;j<npix;j++){
      fread(&shear2[j], sizeof(float), 1, fp);
      for(i=0;i<6;i++) if(j==jj[i]) fread(&dummy, sizeof(long), 1, fp);
    }

    fread(&dummy, sizeof(long), 1, fp);
    for(j=0;j<npix;j++){
      fread(&rotat[j], sizeof(float), 1, fp);
      for(i=0;i<6;i++) if(j==jj[i]) fread(&dummy, sizeof(long), 1, fp);
    }

    fclose(fp);             // Close Takahasi file.

    verb_print(cmd->verbose, "\n\tinputdata_takahasi: total read points = %ld\n",npix);

//E End reading Takahasi file

    Takahasi_region_selection(cmd, gd,
                              nside, npix, conv, shear1, shear2, rotat, ifile);

    free(conv);
    free(shear1);
    free(shear2);
    free(rotat);

    return SUCCESS;
}

local int Takahasi_region_selection(struct cmdline_data* cmd, 
                                    struct  global_data* gd,
                                    int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat, int ifile)
{
    double theta,phi;
    long i;

//B Computing min and max of theta and phi:
    real theta_min, theta_max;
    real phi_min, phi_max;
    pix2ang(0,nside,&theta_min,&phi_min);
    theta_max = theta_min;
    phi_max = phi_min;

    for(i=1;i<npix;i++){   // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
        theta_min = MIN(theta_min,theta);
        theta_max = MAX(theta_max,theta);
        phi_min = MIN(phi_min,phi);
        phi_max = MAX(phi_max,phi);
    }
    verb_print(cmd->verbose, "\n\tinputdata_takahasi: min and max of theta = %f %f\n",theta_min, theta_max);
    verb_print(cmd->verbose, "\tinputdata_takahasi: min and max of phi = %f %f\n",phi_min, phi_max);
//E

//B Selection of a region: the (center) lower edge is random... or not
    real rphi, rtheta;
//B Change selection to random or fix, given by thetaL, phiL, thetaR, phiR
    if (scanopt(cmd->options, "random-point")) {
        rphi    = 2.0 * PI * xrandom(0.0, 1.0);
        rtheta    = racos(1.0 - 2.0 * xrandom(0.0, 1.0));
        verb_print(cmd->verbose,
                   "\n\tinputdata_takahasi: random theta and phi = %f %f\n",rtheta, rphi);
    } else {
// The radius of the region is:
//        rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL)))
// ... and the center:
        rphi = 0.5*(cmd->phiR + cmd->phiL);
        rtheta = 0.5*(cmd->thetaR + cmd->thetaL);
        verb_print(cmd->verbose,
                   "\n\tinputdata_takahasi: fix theta and phi = %f %f\n",rtheta, rphi);
    }
//E

    if (cmd->lengthBox>PI) {
        fprintf(stdout,"\ninputdata_takahasi: Warning! %s\n%s\n\n",
                "lengthBox is greater than one of the angular ranges...",
                "using length = PI");
        cmd->lengthBox = PI;
    }

//B Set rotation dtheta and dphi and theta_rot and phi_rot
    real dtheta_rot, dphi_rot;
    real rtheta_rot, rphi_rot;

    if (scanopt(cmd->options, "rotation")) {
        dtheta_rot = rtheta - PI/2,
        dphi_rot = rphi - PI;
        rtheta_rot = rtheta - dtheta_rot;
        rphi_rot = rphi - dphi_rot;
        verb_print(cmd->verbose,
                   "\tinputdata_takahasi: rotated theta and phi = %f %f\n",
                   rtheta_rot, rphi_rot);
    } else {
        dtheta_rot = 0.0,
        dphi_rot = 0.0;
        rtheta_rot = rtheta;
        rphi_rot = rphi;
        verb_print(cmd->verbose,
                "\tinputdata_takahasi: theta and phi (no-rotation) = %f %f\n",
                rtheta_rot, rphi_rot);
    }
//E

    real thetaL, thetaR;
    real phiL, phiR;
// Here we chose for the box, left and right values
// Default to select-region
// Fix-center is given by thetaL, phiL, thetaR, phiR.
// But it doesn´t use L and R. It use instead the size of the box, lBox
    if ( scanopt(cmd->options, "rotation")
        || scanopt(cmd->options, "fix-center") ) {
// Center is at rotated chosen angles
        thetaL = rtheta_rot - 0.5*cmd->lengthBox;
        thetaR = rtheta_rot + 0.5*cmd->lengthBox;
        phiL = rphi_rot - 0.5*cmd->lengthBox;
        phiR = rphi_rot + 0.5*cmd->lengthBox;
        verb_print(cmd->verbose,
"\tinputdata_takahasi: theta and phi of the center of the selected region = %lf %lf\n",
                   0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
        verb_print(cmd->verbose,
                "\tinputdata_takahasi: radius of the selected region = %lf\n",
                rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL))) );
    } else {
// Fix-center is given by thetaL, phiL, thetaR, phiR.
// But it does use L and R. It doesn´t use the size of the box, lBox
            thetaL = cmd->thetaL;
            thetaR = cmd->thetaR;
            phiL = cmd->phiL;
            phiR = cmd->phiR;
            dtheta_rot = 0.0;
            dphi_rot = 0.0;
            verb_print(cmd->verbose,
"\tinputdata_takahasi: theta and phi of the center of the selected region = %f %f\n",
                       0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
            verb_print(cmd->verbose,
                "\tinputdata_takahasi: radius of the selected region = %lf\n",
                rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL))) );
    }

    verb_print(cmd->verbose,
               "\tinputdata_takahasi: left and right theta = %f %f\n", 
               thetaL, thetaR);
    verb_print(cmd->verbose,
               "\tinputdata_takahasi: left and right phi = %f %f\n", 
               phiL, phiR);
    verb_print(cmd->verbose,
               "\tinputdata_takahasi: theta and phi d_rotation = %f %f\n",
               dtheta_rot, dphi_rot);
//E

#if NDIM == 3
    real xmin, ymin, zmin;
    real xmax, ymax, zmax;

    if (scanopt(cmd->options, "all")) {
        Takahasi_region_selection_3d_all(cmd, gd,
                                         nside, npix, conv, shear1, shear2, rotat,
            dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
            &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, ifile);
    } else {
        Takahasi_region_selection_3d(cmd, gd,
                                     nside, npix, conv, shear1, shear2, rotat,
                dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, ifile);
    }
#else   // ! TREEDIM
    real xmin, ymin;
    real xmax, ymax;

    Takahasi_region_selection_2d(cmd, gd,
                nside, npix, conv, shear1, shear2, rotat,
                dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                &xmin, &xmax, &ymin, &ymax, ifile);
#endif

#if NDIM == 3
    gd->Box[0] = xmax-xmin;
    gd->Box[1] = ymax-ymin; gd->Box[2] = zmax-zmin;
#else
    gd->Box[0] = xmax-xmin;
    gd->Box[1] = ymax-ymin;
#endif

    return SUCCESS;
}

#if NDIM == 3
local int Takahasi_region_selection_3d_all(struct cmdline_data* cmd, 
                                           struct  global_data* gd,
                                           int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
    real *xmin, real *xmax, real *ymin, real *ymax, real *zmin, real *zmax, int ifile)
{
    long i;
    bodyptr p;
    real weight = 1.0;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;

    cmd->nbody = npix;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for all particle (%ld) storage.\n",
               cmd->nbody*sizeof(body)/(1024.0*1024.0),cmd->nbody);

    *xmin=0., *ymin=0., *zmin=0.;
    *xmax=0., *ymax=0., *zmax=0.;

    for(i=0;i<npix;i++){                                // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
//        printf("%ld %f %f %f %f \n", 
//                i, conv[i], shear1[i], shear2[i], rotat[i]);
//        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytable[ifile]+i;
        iselect++;
        if (scanopt(cmd->options, "rotation")) {
            theta_rot = theta - dtheta_rot;
            phi_rot = phi - dphi_rot;
        } else {
            theta_rot = theta;
            phi_rot = phi;
        }

        spherical_to_cartesians(cmd, gd, theta_rot, phi_rot, Pos(p));

        if (!scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = conv[i];
        else
            Kappa(p) = 2.0;
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+iselect;

        *xmin = Pos(p)[0];
        *ymin = Pos(p)[1];
        *zmin = Pos(p)[2];
        *xmax = Pos(p)[0];
        *ymax = Pos(p)[1];
        *zmax = Pos(p)[2];

        Update(p) = TRUE;

    } // ! end for

    real kavg = 0;
    for(i=0;i<npix;i++){
        p = bodytable[ifile] +i;
        *xmin = MIN(*xmin,Pos(p)[0]);
        *ymin = MIN(*ymin,Pos(p)[1]);
        *zmin = MIN(*zmin,Pos(p)[2]);
        *xmax = MAX(*xmax,Pos(p)[0]);
        *ymax = MAX(*ymax,Pos(p)[1]);
        *zmax = MAX(*zmax,Pos(p)[2]);
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose, 
               "\n\tinputdata_takahasi: min and max of x = %f %f\n",*xmin, *xmax);
    verb_print(cmd->verbose, 
               "\tinputdata_takahasi: min and max of y = %f %f\n",*ymin, *ymax);
    verb_print(cmd->verbose, 
               "\tinputdata_takahasi: min and max of z = %f %f\n",*zmin, *zmax);

    verb_print(cmd->verbose,
        "\n\tinputdata_takahasi: selected all read points and nbody: %ld %ld\n",
        iselect, cmd->nbody);

    verb_print(cmd->verbose, 
               "inputdata_takahasi: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

local int Takahasi_region_selection_3d(struct cmdline_data* cmd, 
                                       struct  global_data* gd,
                                       int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
    real *xmin, real *xmax, real *ymin, real *ymax, real *zmin, real *zmax, int ifile)
{
    long i;
    bodyptr p;
    real weight = 1.0;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;
    
    bodyptr bodytabtmp;
    cmd->nbody = npix;
    bodytabtmp = (bodyptr) allocate(cmd->nbody * sizeof(body));
    verb_print(cmd->verbose, "\nAllocated %g MByte for particle (%ld) storage.\n",
               cmd->nbody*sizeof(body)/(1024.0*1024.0),cmd->nbody);

    *xmin=0., *ymin=0., *zmin=0.;
    *xmax=0., *ymax=0., *zmax=0.;

    for(i=0;i<npix;i++){                                // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
//        printf("%ld %f %f %f %f \n", 
//                i, conv[i], shear1[i], shear2[i], rotat[i]);
//        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;

        if (scanopt(cmd->options, "all")) {
            iselect++;
            if (scanopt(cmd->options, "rotation")) {
                theta_rot = theta - dtheta_rot;
                phi_rot = phi - dphi_rot;
            } else {
                theta_rot = theta;
                phi_rot = phi;
            }

            spherical_to_cartesians(cmd, gd, theta_rot, phi_rot, Pos(p));

            if (!scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = conv[i];
            else
                Kappa(p) = 2.0;
            Type(p) = BODY;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;

            *xmin = Pos(p)[0];
            *ymin = Pos(p)[1];
            *zmin = Pos(p)[2];
            *xmax = Pos(p)[0];
            *ymax = Pos(p)[1];
            *zmax = Pos(p)[2];

            Update(p) = TRUE;

        } else {
            if (thetaL < theta - dtheta_rot && theta - dtheta_rot < thetaR) {
                if (phiL < phi - dphi_rot && phi - dphi_rot < phiR) {
                    iselect++;
                    if (scanopt(cmd->options, "rotation")) {
                        theta_rot = theta - dtheta_rot;
                        phi_rot = phi - dphi_rot;
                    } else {
                        theta_rot = theta;
                        phi_rot = phi;
                    }

                    spherical_to_cartesians(cmd, gd, theta_rot, phi_rot, Pos(p));

                    if (!scanopt(cmd->options, "kappa-fix"))
                        Kappa(p) = conv[i];
                    else
                        Kappa(p) = 2.0;
                    Type(p) = BODY;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;

                    *xmin = Pos(p)[0];
                    *ymin = Pos(p)[1];
                    *zmin = Pos(p)[2];
                    *xmax = Pos(p)[0];
                    *ymax = Pos(p)[1];
                    *zmax = Pos(p)[2];

                    Update(p) = TRUE;
                }
            }
        }
    } // ! end for

    bodyptr q;
    if (!scanopt(cmd->options, "all"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
//    bodytab = (bodyptr) allocate(cmd->nbody * sizeof(body));
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    real kavg = 0;
    INTEGER ij=0;
    for(i=0;i<npix;i++){
        q = bodytabtmp+i;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+i;
            *xmin = MIN(*xmin,Pos(p)[0]);
            *ymin = MIN(*ymin,Pos(p)[1]);
            *zmin = MIN(*zmin,Pos(p)[2]);
            *xmax = MAX(*xmax,Pos(p)[0]);
            *ymax = MAX(*ymax,Pos(p)[1]);
            *zmax = MAX(*zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
        }
    }
    verb_print(cmd->verbose, "\n\tinputdata_takahasi: min and max of x = %f %f\n",*xmin, *xmax);
    verb_print(cmd->verbose, "\tinputdata_takahasi: min and max of y = %f %f\n",*ymin, *ymax);
    verb_print(cmd->verbose, "\tinputdata_takahasi: min and max of z = %f %f\n",*zmin, *zmax);
    free(bodytabtmp);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
                   "\n\tinputdata_takahasi: selected read points and nbody: %ld %ld\n",
                   iselect, cmd->nbody);
    else
        verb_print(cmd->verbose, 
                   "\n\tinputdata_takahasi: selected read points = %ld\n",iselect);

    verb_print(cmd->verbose, 
               "inputdata_takahasi: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}

#else

local int Takahasi_region_selection_2d(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax, int ifile)
{
    long i;
    bodyptr p;
    real weight = 1;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;

    bodyptr bodytabtmp;
    cmd->nbody = npix;
    bodytabtmp = (bodyptr) allocate(cmd->nbody * sizeof(body));
    verb_print(cmd->verbose, "\nAllocated %g MByte for particle storage.\n",
               cmd->nbody*sizeof(body)/(1024.0*1024.0));

    *xmin=0., *ymin=0.;
    *xmax=0., *ymax=0.;

    real ra, dec;

    for(i=0;i<npix;i++){   // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
    //        printf("%ld %f %f %f %f \n", i, conv[i], shear1[i], shear2[i], rotat[i]);
    //        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;

        if (scanopt(cmd->options, "all")) {
            iselect++;
            if (scanopt(cmd->options, "rotation")) {
                theta_rot = theta - dtheta_rot;
                phi_rot = phi - dphi_rot;
            } else {
                theta_rot = theta;
                phi_rot = phi;
            }

//            spherical_to_cartesians(theta_rot, phi_rot, Pos(p));
// Better use theta, phi as x, y

            ra = phi_rot;
            dec = theta_rot;
            Pos(p)[0] = ra;
            Pos(p)[1] = dec;

            Kappa(p) = conv[i];
            Type(p) = BODY;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;

            *xmin = Pos(p)[0];
            *ymin = Pos(p)[1];
            *xmax = Pos(p)[0];
            *ymax = Pos(p)[1];

            Update(p) = TRUE;

        } else {

            if (thetaL < theta - dtheta_rot && theta - dtheta_rot < thetaR) {
                if (phiL < phi - dphi_rot && phi - dphi_rot < phiR) {
                    iselect++;
                    if (scanopt(cmd->options, "rotation")) {
                        theta_rot = theta - dtheta_rot;
                        phi_rot = phi - dphi_rot;
                    } else {
                        theta_rot = theta;
                        phi_rot = phi;
                    }

// Better use theta, phi as x, y

                    ra = phi_rot;
                    dec = theta_rot;
                    Pos(p)[0] = ra;
                    Pos(p)[1] = dec;

                    Kappa(p) = conv[i];
                    Type(p) = BODY;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;

                    *xmin = Pos(p)[0];
                    *ymin = Pos(p)[1];
                    *xmax = Pos(p)[0];
                    *ymax = Pos(p)[1];

                    Update(p) = TRUE;
                }
            }
        } // ! all
    } // ! end loop i

    bodyptr q;
    if (!scanopt(cmd->options, "all"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for all particle (%ld) storage.\n",
               cmd->nbody*sizeof(body)/(1024.0*1024.0),cmd->nbody);


    INTEGER ij=0;
    for(i=0;i<npix;i++){
        q = bodytabtmp+i;
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+i;
            *xmin = MIN(*xmin,Pos(p)[0]);
            *ymin = MIN(*ymin,Pos(p)[1]);
            *xmax = MAX(*xmax,Pos(p)[0]);
            *ymax = MAX(*ymax,Pos(p)[1]);
            ij++;
        }
    }
    verb_print(cmd->verbose, 
               "\n\tinputdata_takahasi: min and max of x = %f %f\n",
               *xmin, *xmax);
    verb_print(cmd->verbose, 
               "\tinputdata_takahasi: min and max of y = %f %f\n",
               *ymin, *ymax);

    free(bodytabtmp);
    
    verb_print(cmd->verbose, 
               "\n\tinputdata_takahasi: selected read points = %ld\n",iselect);

    return SUCCESS;
}

#endif

// convert a pixel index (pix) to an angular position (theta,phi)[rad]
//  in spherical coordinates
void pix2ang(long pix, int nside, double *theta, double *phi)
{
  long npix=12*nside*(long)nside, ncap=2*(long)nside*(nside-1);
  long i,j,pp,s;
  double ph,z;

  if(pix<ncap){ // North polar cap
    ph=0.5*(pix+1.);
    i=(long)(sqrt(ph-sqrt((double)((long)ph))))+1;
    j=pix+1-2*i*(i-1);

    z=1.-i*i/(3.*nside*nside);
    *phi=0.5*M_PI/i*(j-0.5);
  }
  else if(pix<(npix-ncap)){ // Equatorial belt
    pp=pix-ncap;
    i=(long)(0.25*pp/nside)+nside;
    j=pp%(4*nside)+1;
    s=(i+nside)%2+1;

    z=4./3.-2.*i/(3.*nside);
    *phi=0.5*M_PI/nside*(j-0.5*s);
  }
  else{ // South polar cap
    ph=0.5*(npix-pix);
    i=(long)(sqrt(ph-sqrt((double)((long)ph))))+1;
    j=4*i+1-(npix-pix-2*i*(i-1));

    z=-1.0+(i*i)/(3.*nside*nside);
    *phi=0.5*M_PI/i*(j-0.5);
    }

  *theta=acos(z);
}

//E End:: Reading Takahasi simulations


int StartOutput(struct cmdline_data *cmd)
{
    outfilefmt_string_to_int(cmd->outfilefmt, &outfilefmt_int);

    if (! strnull(cmd->options))
        verb_print(cmd->verbose, "\n\toptions: %s\n", cmd->options);

    return SUCCESS;
}

int Output(struct cmdline_data* cmd, struct  global_data* gd,
           bodyptr *btable, INTEGER *nbody, int ifile)
{
    double cpustart = CPUTIME;
    if (! strnull(cmd->outfile)) {
        outputdata(cmd, gd, btable[ifile], nbody[ifile]);
    }
    gd->cputotalinout += CPUTIME - cpustart;

    return SUCCESS;
}

local int outputdata(struct cmdline_data* cmd, struct  global_data* gd,
                     bodyptr btable, INTEGER nbody)
{
    switch(outfilefmt_int) {
        case OUTCOLUMNS:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii(cmd, gd, btable, nbody); break;
        case OUTCOLUMNSBIN:
            verb_print(cmd->verbose, "\n\tbinary format output\n");
            outputdata_bin(cmd, gd, btable, nbody); break;
        case OUTNULL:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii(cmd, gd, btable, nbody); break;

#ifdef ADDONS
#include "cballsio_include_03.h"
#endif

        default:
            verb_print(cmd->verbose, 
                    "\n\toutput: Unknown output format...\n\tprinting in default format (columns-ascii)...\n");
                outputdata_ascii(cmd, gd, btable, nbody); break;
    }

    return SUCCESS;
}

local int outputdata_ascii(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, gd->fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
    fprintf(outstr,"#   nbody NDIM\n# %ld %d ",cmd->nbody,NDIM);
#if NDIM == 3
    fprintf(outstr,"%lf %lf %lf\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    fprintf(outstr,"%lf %lf\n",gd->Box[0],gd->Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+cmd->nbody) {
        out_vector_mar(outstr, Pos(p));
        out_real_mar(outstr, Kappa(p));
//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
        out_bool_mar(outstr, HIT(p));
//        out_bool_mar(outstr, Selected(p));
//        out_bool_mar(outstr, Update(p));
#endif
//E
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    verb_print(cmd->verbose, "\tdata output to file %s\n", namebuf);

    return SUCCESS;
}

local int outputdata_bin(struct cmdline_data* cmd, struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody)
{
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, gd->fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
    out_int_bin_long(outstr, nbody);
    out_int_bin(outstr, NDIM);
    out_real_bin(outstr, gd->Box[0]);
    out_real_bin(outstr, gd->Box[1]);
#if NDIM == 3
    out_real_bin(outstr, gd->Box[2]);
#endif
    DO_BODY(p, bodytab, bodytab+nbody)
        out_vector_bin(outstr, Pos(p));
    DO_BODY(p, bodytab, bodytab+nbody)
        out_real_bin(outstr, Kappa(p));
    fclose(outstr);
    verb_print(cmd->verbose, "\tdata output to file %s\n", namebuf);

    return SUCCESS;
}


global int infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(infmt_str,"columns-ascii") == 0)         *infmt_int = INCOLUMNS;
    if (strnull(infmt_str))                             *infmt_int = INNULL;
    if (strcmp(infmt_str,"binary") == 0)                *infmt_int = INCOLUMNSBIN;
    if (strcmp(infmt_str,"takahasi") == 0)              *infmt_int = INTAKAHASI;

#ifdef ADDONS
#include "cballsio_include_08.h"
#endif

    return SUCCESS;
}

local int outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"columns-ascii") == 0)    *outfmt_int = OUTCOLUMNS;
    if (strnull(outfmt_str))                        *outfmt_int = OUTNULL;
    if (strcmp(outfmt_str,"binary") == 0)           *outfmt_int = OUTCOLUMNSBIN;

#ifdef ADDONS
#include "cballsio_include_09.h"
#endif

    return SUCCESS;
}

// I/O directories:
global void setFilesDirs_log(struct cmdline_data* cmd, struct  global_data* gd)
{
    char buf[200];

    sprintf(gd->tmpDir,"%s/%s",cmd->rootDir,"tmp");

    double cpustart = CPUTIME;
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd->tmpDir,gd->tmpDir);
    system(buf);
    gd->cputotalinout += CPUTIME - cpustart;

    sprintf(gd->logfilePath,"%s/cballs%s.log",gd->tmpDir,cmd->suffixOutFiles);
}

global void setFilesDirs(struct cmdline_data* cmd, struct  global_data* gd)
{
    char buf[200];

    double cpustart = CPUTIME;

    int ndefault = 0;
    int ipos;
    char *dp1, *dp2;
    int lenDir = strlen(cmd->rootDir);
    for (int i; i< lenDir; i++) {
        if(cmd->rootDir[i] == '/') {
            ipos = i+1;
            ndefault++;
        }
    }
    if (ndefault>1)
        error("setFilesDirs: more '/' than 1 in 'rootDir=%s'. Use only 1 or none\n",
              cmd->rootDir);

    if (ndefault == 0) {
        sprintf(gd->outputDir,cmd->rootDir);
        sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd->outputDir,gd->outputDir);
        verb_print_q(3, cmd->verbose_log,"\nsystem: %s\n",buf);
        system(buf);
    } else {
        dp1 = (char*) malloc((ipos-1)*sizeof(char));
        strncpy(dp1, cmd->rootDir, ipos-1);
        dp2 = (char*) malloc((lenDir-ipos)*sizeof(char));
        strncpy(dp2, cmd->rootDir + ipos, lenDir-ipos);
        verb_print_q(2,cmd->verbose,
                    "setFilesDirs: '/' counts %d pos %d and %s %s\n",
                    ndefault, ipos, dp1, dp2);

        sprintf(gd->outputDir,cmd->rootDir);
        sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",dp1,dp1);
        verb_print_q(3,cmd->verbose_log,"\nsystem: %s\n",buf);
        system(buf);
        sprintf(buf,"if [ ! -d %s ]; then mkdir %s/%s; fi",gd->outputDir,dp1,dp2);
        verb_print_q(3,cmd->verbose_log,"\nsystem: %s\n",buf);
        system(buf);
    }
    gd->cputotalinout += CPUTIME - cpustart;

    sprintf(gd->fpfnameOutputFileName,"%s/%s%s%s",
            cmd->rootDir,cmd->outfile,cmd->suffixOutFiles,EXTFILES);
    sprintf(gd->fpfnamehistNNFileName,"%s/%s%s%s",
            cmd->rootDir,cmd->histNNFileName,cmd->suffixOutFiles,EXTFILES);
// Will be gd->histNNFileName. Set it in gd structure
    sprintf(gd->fpfnamehistCFFileName,"%s/%s%s%s",
            cmd->rootDir,"histCF",cmd->suffixOutFiles,EXTFILES);
    sprintf(gd->fpfnamehistrBinsFileName,"%s/%s%s%s",
            cmd->rootDir,"rbins",cmd->suffixOutFiles,EXTFILES);
    sprintf(gd->fpfnamehistXi2pcfFileName,"%s/%s%s%s",
        cmd->rootDir,cmd->histXi2pcfFileName, cmd->suffixOutFiles,EXTFILES);
    sprintf(gd->fpfnamehistZetaGFileName,"%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"G",cmd->suffixOutFiles);
    sprintf(gd->fpfnamehistZetaGmFileName,"%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"G",cmd->suffixOutFiles);
    sprintf(gd->fpfnamehistZetaMFileName,"%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"M",cmd->suffixOutFiles);
    sprintf(gd->fpfnamemhistZetaMFileName,"%s/%s%s%s%s",
            cmd->rootDir,"m",cmd->histZetaFileName,"M",cmd->suffixOutFiles);
    sprintf(gd->fpfnameCPUFileName,"%s/cputime%s%s",
            cmd->rootDir,cmd->suffixOutFiles,EXTFILES);
}


int EndRun(struct cmdline_data* cmd, struct  global_data* gd)
{
    stream outstr;

	fclose(gd->outlog);

    if (cmd->verbose >= 2) {
        printf("\nrSize \t\t= %lf\n", gd->rSizeTable[0]);
        printf("nbbcalc \t= %ld\n", gd->nbbcalc);
        printf("nbccalc \t= %ld\n", gd->nbccalc);
// BALLS
        printf("ncccalc \t= %ld\n", gd->ncccalc);
//
        printf("tdepth \t\t= %d\n", gd->tdepthTable[0]);
        printf("ncell\t\t= %ld\n", gd->ncellTable[0]);
        verb_print_q(3,cmd->verbose,"sameposcount \t= %ld\n",gd->sameposcount);
#ifdef OPENMPCODE
        printf("cpusearch \t= %lf (be aware of the number of threads)\n", gd->cpusearch*cmd->numthreads);     // in minutes
#else
        printf("cpusearch \t= %lf\n", gd->cpusearch);     // in minutes
#endif
        printf("cputotalinout \t= %lf\n", gd->cputotalinout);
    }

    if (cmd->verbose > 2) {
#ifdef OPENMPCODE
        real cpuTotal = gd->cpusearch*cmd->numthreads;
#else
        real cpuTotal = gd->cpusearch;
#endif
        real m1, m2, m3, m4;
        m1 =1.0;
        m2 =1.0e-5;
        m3 =1.0e-9;
        m4 =1.0e-12;
        outstr = stropen(gd->fpfnameCPUFileName, "a");
        fprintf(outstr,"%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                (double) cmd->nbody, cpuTotal,
                            m1*(double) cmd->nbody,
                            m2*((double) cmd->nbody)*rlog10((double) cmd->nbody),
                            m3*rsqr((double) cmd->nbody),
                            m4*rpow(((double) cmd->nbody), 3.0));
        fclose(outstr);
    }

    if (cmd->verbose > 0) {
        printf("\nFinal CPU time : %lf %s\n", CPUTIME - gd->cpuinit,PRNUNITOFTIMEUSED);// in minutes
        printf("Final real time: %ld", (rcpu_time()-gd->cpurealinit)); // in minutes
        printf(" %s\n\n", PRNUNITOFTIMEUSED);                          // Only work this way
    }

    EndRun_FreeMemory(cmd, gd);

    return SUCCESS;
}

//
// We must check the order of memory allocation and dealocation!!!
//
local int EndRun_FreeMemory(struct cmdline_data* cmd, struct  global_data* gd)
{
    int m;

    if (cmd->computeShearCF) {
        free_dvector(gd->histXitx,1,cmd->sizeHistN);
        free_dvector(gd->histXixx,1,cmd->sizeHistN);
        free_dvector(gd->histXitt,1,cmd->sizeHistN);
    }

    if (cmd->computeTPCF) {
        free_dmatrix3D(gd->histZetaGmIm,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaGmRe,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
#ifdef USEGSL
        for (m=1; m<=cmd->mChebyshev+1; m++)
            gsl_matrix_complex_free(histZetaMatrix[m].histZetaM);
        free(histZetaMatrix);
        gsl_matrix_complex_free(gd->histXi_gsl);
#endif
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        free_dmatrix3D(gd->histZetaMcossin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMsincos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaM,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,cmd->sizeHistN);
        free_dmatrix(gd->histXisin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(gd->histXicos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(gd->histXi,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
    }

     free_dvector(gd->histXi2pcf,1,cmd->sizeHistN);
     
#ifdef ADDONS
#include "cballsio_include_10.h"
#endif
     
     free_dvector(gd->histNNN,1,cmd->sizeHistN);
     // 2pcf
     //B kappa Avg Rmin
     free_dvector(gd->histNNSubXi2pcftotal,1,cmd->sizeHistN);
     //E
     free_dvector(gd->histNNSubXi2pcf,1,cmd->sizeHistN);
     //
     free_dvector(gd->histNNSub,1,cmd->sizeHistN);
     free_dvector(gd->histCF,1,cmd->sizeHistN);
     free_dvector(gd->histNN,1,cmd->sizeHistN);

//B Set gsl uniform random :: If not needed globally this line have to go to testdata
#ifdef USEGSL
    gsl_rng_free (gd->r);
#endif

    return SUCCESS;
}

#ifdef ADDONS
#include "cballsio_include_11a.h"
#endif


#ifdef ADDONS
#include "cballsio_include_11b.h"
#endif

