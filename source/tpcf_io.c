/*==============================================================================
 MODULE: tpcf_io.c		[tpcf]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	april 2023
 Purpose: Routines to drive input and output data
 Language: C
 Use:
 Major revisions:
 ==============================================================================*/


#include "globaldefs.h"

local void inputdata_ascii(void);
#if NDIM == 3
local void inputdata_ascii_2d_to_3d(void);
#endif
local void inputdata_bin(void);
local int inputdata_takahasi(void);
local void outputdata(void);
local void outputdata_ascii(void);
local void outputdata_bin(void);
//local void outputdata_smooth(void);
//local void outputdata_bin_smooth(void);
//local void outputdata_ascii_smooth(void);


local void outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;

int inputdata(void)
{
    double cpustart = CPUTIME;

//    int ifile;

//    for (ifile=0; ifile<gd.nfiles; ifile++) {

#ifdef MPICODE
    if (ThisTask==0) {                      // Only task 0 can get points catalog
#endif
    switch(gd.infilefmt_int) {
        case INCOLUMNS:
            printf("\n\tInput in columns (ascii) format...\n");
            inputdata_ascii(); break;
#if NDIM == 3
        case INCOLUMNS2DTO3D:
            printf("\n\tInput in columns (ascii-2d-to-3d) format...\n");
            inputdata_ascii_2d_to_3d(); break;
#endif
        case INNULL:
            printf("\n\t(Null) Input in columns (ascii) format...\n");
            inputdata_ascii(); break;
        case INCOLUMNSBIN:
            printf("\n\tInput in binary format...\n");
            inputdata_bin(); break;
        case INTAKAHASI:
            printf("\n\tInput in takahasi format...\n");
            inputdata_takahasi(); break;


        default:
            printf("\n\tInput: Unknown input format...");
            printf("\n\tInput in default columns (ascii) format...\n");
            inputdata_ascii(); break;
    }
#ifdef MPICODE
    }
#endif


//    } // end of ifile

    gd.cputotalinout += CPUTIME - cpustart;
    verb_print(cmd.verbose, "\n\tinputdata :: reading time = %f\n",CPUTIME - cpustart);

    return _SUCCESS_;
}

local void inputdata_ascii(void)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real weight=1;

    gd.model_comment = "Column form input file";

//    instr = stropen(cmd.infile, "r");
    instr = stropen(gd.infilenames[0], "r");

    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd.nbody);
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);

// Check the center of the box!!!
#if NDIM == 3
    real Lx, Ly, Lz;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    in_real(instr, &Lz);
    gd.Box[0] = Lx;
    gd.Box[1] = Ly;
    gd.Box[2] = Lz;
#else
    real Lx, Ly;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    gd.Box[0] = Lx;
    gd.Box[1] = Ly;
#endif

    verb_print(cmd.verbose, "\tInput: nbody and ndim: %d %d...\n", cmd.nbody, ndim);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        in_vector(instr, Pos(p));
        in_real(instr, &Kappa(p));
        if (scanopt(cmd.options, "kappa-constant"))
            Kappa(p) = 2.0;
    }

    fclose(instr);

    real kavg;
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytab+1;
        kavg += Kappa(p);
    }
    verb_print(cmd.verbose, "inputdata_ascii: average of kappa (%ld particles) = %le\n",
               cmd.nbody, kavg/((real)cmd.nbody) );

//B Locate particles with same position
    if (scanopt(cmd.options, "check-eq-pos")) {
    bodyptr q;
    real dist2;
    vector distv;
//    INTEGER ip, iq;
    bool flag=0;
    DO_BODY(p, bodytab, bodytab+cmd.nbody-1)
        DO_BODY(q, p+1, bodytab+cmd.nbody)
            if (p != q) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                if (dist2 == 0.0) {
                    flag=1;
                }
            }
    if (flag) error("inputdata_ascii: at least two bodies have same position\n");
    }
//E
}

#if NDIM == 3
local void inputdata_ascii_2d_to_3d(void)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real weight=1;

    gd.model_comment = "Column form input file (2d-to-3d)";

    instr = stropen(cmd.infile, "r");
    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd.nbody);
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int(instr, &ndim);
    if (ndim != 2)
        error("inputdata: ndim = %d; expected 2\n", ndim);

// Check the center of the box!!!
    real Lx, Ly;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    gd.Box[0] = Lx;
    gd.Box[1] = Ly;
// Added this line to set lbox in z direction. Check if it es necessary
    gd.Box[2] = Ly;

    verb_print(cmd.verbose, "\tInput: nbody and ndim: %d %d...\n", cmd.nbody, ndim);
    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

//    real tmp;
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        in_real(instr, &Pos(p)[0]);
        in_real(instr, &Pos(p)[1]);
        in_real(instr, &Kappa(p));
        if (scanopt(cmd.options, "kappa-constant"))
            Kappa(p) = 2.0;
    }

    fclose(instr);

//B Find MIN and MAX
    real theta, phi;
    real theta_min, theta_max;
    real phi_min, phi_max;
    p = bodytab;
    theta_max = theta_min = Pos(p)[0];
    phi_max = phi_min  = Pos(p)[1];

//    for(i=1;i<npix;i++){   // Healpix ring scheme
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        theta = Pos(p)[0];
        phi = Pos(p)[1];
        theta_min = MIN(theta_min,theta);
        theta_max = MAX(theta_max,theta);
        phi_min = MIN(phi_min,phi);
        phi_max = MAX(phi_max,phi);
    }
    verb_print(cmd.verbose, "\n\tinputdata_AA: min and max of theta = %f %f\n",
               theta_min, theta_max);
    verb_print(cmd.verbose, "\tinputdata_AA: min and max of phi = %f %f\n",
               phi_min, phi_max);
//E

    real ra, dec;       // phi, theta from pix2ang :: column 2, column 1, respectively
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        ra = Pos(p)[0];
        dec = Pos(p)[1];
//        Pos(p)[2] = rsin(dec)*rcos(ra);
//        Pos(p)[0] = rsin(dec)*rsin(ra);
//        Pos(p)[1] = rcos(dec);
//
        Pos(p)[0] = rsin(dec)*rcos(ra);
        Pos(p)[1] = rsin(dec)*rsin(ra);
        Pos(p)[2] = rcos(dec);
    }

    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytab+1;
    }
}
#endif

local void inputdata_bin(void)
{
    stream instr;
    int ndim;
    bodyptr p;
    real weight=1;

    gd.model_comment = "Binary input file";

    instr = stropen(cmd.infile, "r");
    in_int_bin_long(instr, &cmd.nbody);
    verb_print(cmd.verbose, "\tInput: nbody %d\n", cmd.nbody);
    if (cmd.nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd.nbody);
    in_int_bin(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    verb_print(cmd.verbose, "\tInput: nbody and ndim: %d %d...\n", cmd.nbody, ndim);

    in_real_bin(instr, &gd.Box[0]);
    in_real_bin(instr, &gd.Box[1]);
#if NDIM == 3
    in_real_bin(instr, &gd.Box[2]);
#endif

#if NDIM == 3
    verb_print(cmd.verbose, "\tInput: Box: %g %g %g\n", gd.Box[0], gd.Box[1], gd.Box[2]);
#else
    verb_print(cmd.verbose, "\tInput: Box: %g %g %g\n", gd.Box[0], gd.Box[1]);
#endif

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    DO_BODY(p, bodytab, bodytab+cmd.nbody)
        in_vector_bin(instr, Pos(p));
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        in_real_bin(instr, &Kappa(p));
        if (scanopt(cmd.options, "kappa-constant"))
            Kappa(p) = 2.0;
    }
    fclose(instr);

    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
        Type(p) = BODY;
        Weight(p) = weight;
        Id(p) = p-bodytab+1;
    }
}





//B BEGIN:: Reading Takahasi simulations
//From Takahasi web page. Adapted to our needs

#include<math.h>
#include<stdio.h>
#include<stdlib.h>


void pix2ang(long pix, int nside, double *theta, double *phi);

local int Takahasi_region_selection(int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat);
local int Takahasi_region_selection_3d(int nside, int npix,
    float *conv, float *shear1, float *shear2, float *rotat,
    real dtheta_rot, real thetaL, real thetaR,
    real dphi_rot, real phiL, real phiR,
real *xmin, real *xmax, real *ymin, real *ymax, real *zmin, real *zmax);
#if NDIM == 2
local int Takahasi_region_selection_2d(int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
                                       real *xmin, real *xmax, real *ymin, real *ymax);
#endif

local int inputdata_takahasi(void)
{
    FILE *fp;
    long i,j,npix,dummy;
    long jj[6]={536870908,1073741818,1610612728,2147483638,2684354547,3221225457};
    int negi,nside;
//    double theta,phi;
//    char file[200];

    gd.model_comment = "Takahasi input file";

//E Begin reading Takahasi file
    fp = stropen(cmd.infile, "rb");

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

    verb_print(cmd.verbose, "\nAllocated %g MByte for pixel storage.\n",
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

    verb_print(cmd.verbose, "\n\tinputdata_takahasi: total read points = %ld\n",npix);

//E End reading Takahasi file

    Takahasi_region_selection(nside, npix, conv, shear1, shear2, rotat);

    free(conv);
    free(shear1);
    free(shear2);
    free(rotat);

    return 0;
}

local int Takahasi_region_selection(int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat)
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
    verb_print(cmd.verbose, "\n\tinputdata_takahasi: min and max of theta = %f %f\n",theta_min, theta_max);
    verb_print(cmd.verbose, "\tinputdata_takahasi: min and max of phi = %f %f\n",phi_min, phi_max);
//E

//B Selection of a region: the (center) lower edge is random... or not
    real rphi, rtheta;
//B Change selection to random or fix, given by thetaL, phiL, thetaR, phiR
    if (scanopt(cmd.options, "random-point")) {
        rphi    = 2.0 * PI * xrandom(0.0, 1.0);
        rtheta    = racos(1.0 - 2.0 * xrandom(0.0, 1.0));
        verb_print(cmd.verbose,
                   "\n\tinputdata_takahasi: random theta and phi = %f %f\n",rtheta, rphi);
    } else {
// The radius of the region is:
//        rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL)))
// ... and the center:
        rphi = 0.5*(cmd.phiR + cmd.phiL);
        rtheta = 0.5*(cmd.thetaR + cmd.thetaL);
        verb_print(cmd.verbose,
                   "\n\tinputdata_takahasi: fix theta and phi = %f %f\n",rtheta, rphi);
    }
//E

    if (cmd.lengthBox>PI) {
        fprintf(stdout,"\ninputdata_takahasi: Warning! %s\n%s\n\n",
                "lengthBox is greater than one of the angular ranges...",
                "using length = PI");
        cmd.lengthBox = PI;
    }

//B Set rotation dtheta and dphi and theta_rot and phi_rot
    real dtheta_rot, dphi_rot;
    real rtheta_rot, rphi_rot;

    if (scanopt(cmd.options, "rotation")) {
        dtheta_rot = rtheta - PI/2,
        dphi_rot = rphi - PI;
        rtheta_rot = rtheta - dtheta_rot;
        rphi_rot = rphi - dphi_rot;
        verb_print(cmd.verbose,
                   "\tinputdata_takahasi: rotated theta and phi = %f %f\n",
                   rtheta_rot, rphi_rot);
    } else {
        dtheta_rot = 0.0,
        dphi_rot = 0.0;
        rtheta_rot = rtheta;
        rphi_rot = rphi;
        verb_print(cmd.verbose,
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
//    if (scanopt(cmd.options, "rotation")) {
    if ( scanopt(cmd.options, "rotation") || scanopt(cmd.options, "fix-center") ) {
// Center is at rotated chosen angles
        thetaL = rtheta_rot - 0.5*cmd.lengthBox;
        thetaR = rtheta_rot + 0.5*cmd.lengthBox;
        phiL = rphi_rot - 0.5*cmd.lengthBox;
        phiR = rphi_rot + 0.5*cmd.lengthBox;
        verb_print(cmd.verbose,
"\tinputdata_takahasi: theta and phi of the center of the selected region = %lf %lf\n",
                   0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
        verb_print(cmd.verbose,
"\tinputdata_takahasi: radius of the selected region = %lf\n",
                   rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL))) );
    } else {
// Fix-center is given by thetaL, phiL, thetaR, phiR.
// But it does use L and R. It doesn´t use the size of the box, lBox
            thetaL = cmd.thetaL;
            thetaR = cmd.thetaR;
            phiL = cmd.phiL;
            phiR = cmd.phiR;
            dtheta_rot = 0.0;
            dphi_rot = 0.0;
            verb_print(cmd.verbose,
"\tinputdata_takahasi: theta and phi of the center of the selected region = %f %f\n",
                       0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
            verb_print(cmd.verbose,
"\tinputdata_takahasi: radius of the selected region = %lf\n",
                       rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL))) );
    }

    verb_print(cmd.verbose,
               "\tinputdata_takahasi: left and right theta = %f %f\n", thetaL, thetaR);
    verb_print(cmd.verbose,
               "\tinputdata_takahasi: left and right phi = %f %f\n", phiL, phiR);
    verb_print(cmd.verbose,
               "\tinputdata_takahasi: theta and phi d_rotation = %f %f\n",
               dtheta_rot, dphi_rot);
//E

#if NDIM == 3

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;

    Takahasi_region_selection_3d(nside, npix, conv, shear1, shear2, rotat,
            dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
            &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);

#else   // ! TREEDIM

    real xmin, ymin;
    real xmax, ymax;

    Takahasi_region_selection_2d(nside, npix, conv, shear1, shear2, rotat,
            dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                                 &xmin, &xmax, &ymin, &ymax);

#endif

#if NDIM == 3
    gd.Box[0] = xmax-xmin;
    gd.Box[1] = ymax-ymin; gd.Box[2] = zmax-zmin;
#else
    gd.Box[0] = xmax-xmin;
    gd.Box[1] = ymax-ymin;
#endif

    return _SUCCESS_;
}

#if NDIM == 3
local int Takahasi_region_selection_3d(int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
    real *xmin, real *xmax, real *ymin, real *ymax, real *zmin, real *zmax)
{
    long i;
    bodyptr p;
    real weight = 1.0;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;
    
    bodyptr bodytabtmp;
    cmd.nbody = npix;
    bodytabtmp = (bodyptr) allocate(cmd.nbody * sizeof(body));
    verb_print(cmd.verbose, "\nAllocated %g MByte for particle (%ld) storage.\n",
               cmd.nbody*sizeof(body)/(1024.0*1024.0),cmd.nbody);

    *xmin=0., *ymin=0., *zmin=0.;
    *xmax=0., *ymax=0., *zmax=0.;

    for(i=0;i<npix;i++){   // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
//        printf("%ld %f %f %f %f \n", i, conv[i], shear1[i], shear2[i], rotat[i]);
//        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;

        if (scanopt(cmd.options, "all")) {
            iselect++;
            if (scanopt(cmd.options, "rotation")) {
                theta_rot = theta - dtheta_rot;
                phi_rot = phi - dphi_rot;
            } else {
                theta_rot = theta;
                phi_rot = phi;
            }

            spherical_to_cartesians(theta_rot, phi_rot, Pos(p));

            if (!scanopt(cmd.options, "kappa-constant"))
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
                    if (scanopt(cmd.options, "rotation")) {
                        theta_rot = theta - dtheta_rot;
                        phi_rot = phi - dphi_rot;
                    } else {
                        theta_rot = theta;
                        phi_rot = phi;
                    }

                    spherical_to_cartesians(theta_rot, phi_rot, Pos(p));

                    if (!scanopt(cmd.options, "kappa-fix"))
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
    if (!scanopt(cmd.options, "all"))
        cmd.nbody = iselect;

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    real kavg = 0;
    INTEGER ij=0;
    for(i=0;i<npix;i++){
        q = bodytabtmp+i;
        if(Update(q)) {
            p = bodytab+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Weight(p) = weight;
            Id(p) = p-bodytab+i;
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
    verb_print(cmd.verbose, "\n\tinputdata_takahasi: min and max of x = %f %f\n",*xmin, *xmax);
    verb_print(cmd.verbose, "\tinputdata_takahasi: min and max of y = %f %f\n",*ymin, *ymax);
    verb_print(cmd.verbose, "\tinputdata_takahasi: min and max of z = %f %f\n",*zmin, *zmax);
    free(bodytabtmp);

    if (scanopt(cmd.options, "all"))
        verb_print(cmd.verbose,
                   "\n\tinputdata_takahasi: selected read points and nbody: %ld %ld\n",
                   iselect, cmd.nbody);
    else
        verb_print(cmd.verbose, "\n\tinputdata_takahasi: selected read points = %ld\n",iselect);

    verb_print(cmd.verbose, "inputdata_takahasi: average of kappa (%ld particles) = %le\n",
               cmd.nbody, kavg/((real)cmd.nbody) );

    return _SUCCESS_;
}

#else

local int Takahasi_region_selection_2d(int nside, int npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax)
{
    long i;
    bodyptr p;
    real weight = 1;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;

    bodyptr bodytabtmp;
    cmd.nbody = npix;
    bodytabtmp = (bodyptr) allocate(cmd.nbody * sizeof(body));
    verb_print(cmd.verbose, "\nAllocated %g MByte for particle storage.\n",
               cmd.nbody*sizeof(body)/(1024.0*1024.0));

    *xmin=0., *ymin=0.;
    *xmax=0., *ymax=0.;

    real ra, dec;

    for(i=0;i<npix;i++){   // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
    //        printf("%ld %f %f %f %f \n", i, conv[i], shear1[i], shear2[i], rotat[i]);
    //        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;

        if (scanopt(cmd.options, "all")) {
            iselect++;
            if (scanopt(cmd.options, "rotation")) {
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
//            *zmin = Pos(p)[2];
            *xmax = Pos(p)[0];
            *ymax = Pos(p)[1];
//            *zmax = Pos(p)[2];

            Update(p) = TRUE;

        } else {

        if (thetaL < theta - dtheta_rot && theta - dtheta_rot < thetaR) {
            if (phiL < phi - dphi_rot && phi - dphi_rot < phiR) {
                iselect++;
                if (scanopt(cmd.options, "rotation")) {
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
//                Id(p) = p-bodytab+iselect;
                Id(p) = p-bodytabtmp+iselect;

                *xmin = Pos(p)[0];
                *ymin = Pos(p)[1];
                *xmax = Pos(p)[0];
                *ymax = Pos(p)[1];

                Update(p) = TRUE;
            }
        }
    }
    }

    bodyptr q;
    if (!scanopt(cmd.options, "all"))
        cmd.nbody = iselect;
//    cmd.nbody = iselect;

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));
    INTEGER ij=0;
    for(i=0;i<npix;i++){
        q = bodytabtmp+i;
        if(Update(q)) {
            p = bodytab+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Weight(p) = weight;
            Id(p) = p-bodytab+i;
            *xmin = MIN(*xmin,Pos(p)[0]);
            *ymin = MIN(*ymin,Pos(p)[1]);
            *xmax = MAX(*xmax,Pos(p)[0]);
            *ymax = MAX(*ymax,Pos(p)[1]);
            ij++;
        }
    }
    verb_print(cmd.verbose, "\n\tinputdata_takahasi: min and max of x = %f %f\n",*xmin, *xmax);
    verb_print(cmd.verbose, "\tinputdata_takahasi: min and max of y = %f %f\n",*ymin, *ymax);
    free(bodytabtmp);
    
    verb_print(cmd.verbose, "\n\tinputdata_takahasi: selected read points = %ld\n",iselect);

    return _SUCCESS_;
}

#endif


void pix2ang(long pix, int nside, double *theta, double *phi)  // convert a pixel index (pix) to an angular position (theta,phi)[rad] in spherical coordinates
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


int StartOutput(void)
{
    outfilefmt_string_to_int(cmd.outfilefmt, &outfilefmt_int);

    if (! strnull(cmd.options))
        verb_print(cmd.verbose, "\n\toptions: %s\n", cmd.options);

    return _SUCCESS_;
}

int output(void)
{
    double cpustart = CPUTIME;
#ifdef MPICODE
    if (ThisTask==0) {
#endif
    if (! strnull(cmd.outfile)) {
        outputdata();

//        if (scanopt(cmd.options, "smooth"))
//            outputdata_smooth();


    }
    gd.cputotalinout += CPUTIME - cpustart;
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local void outputdata(void)
{
    switch(outfilefmt_int) {
        case OUTCOLUMNS:
            verb_print(cmd.verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii(); break;
        case OUTCOLUMNSBIN:
            verb_print(cmd.verbose, "\n\tcolumns-bin format output\n");
            outputdata_bin(); break;
        case OUTNULL:
            verb_print(cmd.verbose, "\n\tcolumns-bin format output\n");
            outputdata_ascii(); break;
        default:
            verb_print(cmd.verbose, "\n\toutput: Unknown output format...\n\tprinting in default format (columns-ascii)...\n");
                outputdata_ascii(); break;
    }
}

/*
local void outputdata_smooth(void)
{
    switch(outfilefmt_int) {
        case OUTCOLUMNS:
            verb_print(cmd.verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii_smooth(); break;
        case OUTCOLUMNSBIN:
            verb_print(cmd.verbose, "\n\tcolumns-bin format output\n");
            outputdata_bin_smooth(); break;
        case OUTNULL:
            verb_print(cmd.verbose, "\n\tcolumns-bin format output\n");
            outputdata_ascii_smooth(); break;
        default:
            verb_print(cmd.verbose, "\n\toutput: Unknown output format...\n\tprinting in default format (columns-ascii)...\n");
            outputdata_ascii_smooth(); break;
    }
}
*/

local void outputdata_ascii(void)
{
    char namebuf[256];
//    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, gd.fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
    fprintf(outstr,"#   nbody NDIM\n# %ld %d ",cmd.nbody,NDIM);
#if NDIM == 3
    fprintf(outstr,"%lf %lf %lf\n",gd.Box[0],gd.Box[1],gd.Box[2]);
#else
    fprintf(outstr,"%lf %lf\n",gd.Box[0],gd.Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+cmd.nbody) {
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
    verb_print(cmd.verbose, "\tdata output to file %s\n", namebuf);
}

/*
local void outputdata_ascii_smooth(void)
{
    char namebuf[256];
//    struct stat buf;
    stream outstr;
    bodyptr p;

//B Locate particles with same position
    int MAXNUMPOINTS=10;
    int *checklist;
    checklist = ivector(1,MAXNUMPOINTS);

    bodyptr q;
    real dist2;
    vector distv;
//    INTEGER ip, iq;
    bool flag=0;
    int ifound=0;
    int i;
    if (scanopt(cmd.options, "check-eq-pos")) {
    DO_BODY(p, bodytabsm, bodytabsm+gd.nbodysm-1)
        DO_BODY(q, p+1, bodytabsm+gd.nbodysm) {
            DOTPSUBV(dist2, distv, Pos(p), Pos(q));
            if (dist2 == 0.0) {
                ifound++;
                checklist[ifound] = Id(q);
                    flag=1;
                }
        }
    verb_print(cmd.verbose,
               "outputdata_ascii_smooth: found %d pair of bodies with same position\n",
               ifound);
    }
//E
    sprintf(namebuf,"%s/%s%s_smooth%s",cmd.rootDir,cmd.outfile,cmd.suffixOutFiles,EXTFILES);
    outstr = stropen(namebuf, "w!");
    fprintf(outstr,"#   nbody NDIM\n# %ld %d ",gd.nbodysm-ifound,NDIM);
#if NDIM == 3
    fprintf(outstr,"%lf %lf %lf\n",gd.Box[0],gd.Box[1],gd.Box[2]);
#else
    fprintf(outstr,"%lf %lf\n",gd.Box[0],gd.Box[1]);
#endif

    DO_BODY(p, bodytabsm, bodytabsm+gd.nbodysm) {
        if (scanopt(cmd.options, "check-eq-pos")) {
        flag=1;
        for (i=1; i<=ifound; i++) {
            if (Id(p)==checklist[i]) {
                flag=0;
                continue;
            }
        }
        if (!flag) continue;
        }
        out_vector_mar(outstr, Pos(p));
        out_real_mar(outstr, Kappa(p));
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    verb_print(cmd.verbose, "\toutputdata_ascii_smooth: data output to file %s\n", namebuf);
}
*/

local void outputdata_bin(void)
{
    char namebuf[256];
//    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf, gd.fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
    out_int_bin_long(outstr, cmd.nbody);
    out_int_bin(outstr, NDIM);
    out_real_bin(outstr, gd.Box[0]);
    out_real_bin(outstr, gd.Box[1]);
#if NDIM == 3
    out_real_bin(outstr, gd.Box[2]);
#endif
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_vector_bin(outstr, Pos(p));
    DO_BODY(p, bodytab, bodytab+cmd.nbody)
        out_real_bin(outstr, Kappa(p));
    fclose(outstr);
    verb_print(cmd.verbose, "\tdata output to file %s\n", namebuf);
}

/*
local void outputdata_bin_smooth(void)
{
    char namebuf[256];
//    struct stat buf;
    stream outstr;
    bodyptr p;

    sprintf(namebuf,"%s/%s%s_smooth%s",cmd.rootDir,cmd.outfile,cmd.suffixOutFiles,EXTFILES);
    outstr = stropen(namebuf, "w!");
    out_int_bin_long(outstr, gd.nbodysm);
    out_int_bin(outstr, NDIM);
    out_real_bin(outstr, gd.Box[0]);
    out_real_bin(outstr, gd.Box[1]);
#if NDIM == 3
    out_real_bin(outstr, gd.Box[2]);
#endif
    DO_BODY(p, bodytabsm, bodytabsm+gd.nbodysm)
        out_vector_bin(outstr, Pos(p));
    DO_BODY(p, bodytabsm, bodytabsm+gd.nbodysm)
        out_real_bin(outstr, Kappa(p));
    fclose(outstr);
    verb_print(cmd.verbose, "\tdata output to file %s\n", namebuf);
}
*/

local void outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"columns-ascii") == 0)    *outfmt_int = OUTCOLUMNS;
    if (strnull(outfmt_str))                        *outfmt_int = OUTNULL;
    if (strcmp(outfmt_str,"binary") == 0)       *outfmt_int = OUTCOLUMNSBIN;
}

// I/O directories:
global void setFilesDirs_log(void)
{
    char buf[200];

#ifdef MPICODE
        if(ThisTask==0) {
#endif
    sprintf(gd.tmpDir,"%s/%s",cmd.rootDir,"tmp");

    double cpustart = CPUTIME;
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);
    gd.cputotalinout += CPUTIME - cpustart;

    sprintf(gd.logfilePath,"%s/cballs%s.log",gd.tmpDir,cmd.suffixOutFiles);
#ifdef MPICODE
        }
#endif
}

global void setFilesDirs(void)
{
    char buf[200];

#ifdef MPICODE
        if(ThisTask==0) {
#endif
    sprintf(gd.outputDir,cmd.rootDir);
    double cpustart = CPUTIME;
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.outputDir,gd.outputDir);
//    fprintf(gd.outlog,"\nsystem: %s\n",buf);
#ifdef DEBUG
    if (cmd.verbose_log>=3)
    verb_log_print(cmd.verbose_log,gd.outlog,"\nsystem: %s\n",buf);
#endif
    system(buf);
    gd.cputotalinout += CPUTIME - cpustart;

    sprintf(gd.fpfnameOutputFileName,"%s/%s%s%s",cmd.rootDir,cmd.outfile,cmd.suffixOutFiles,EXTFILES);
    sprintf(gd.fpfnamehistNFileName,"%s/%s%s%s",cmd.rootDir,cmd.histNFileName,cmd.suffixOutFiles,EXTFILES);
// Will be gd.histNFileName. Set it in gd structure
    sprintf(gd.fpfnamehistCFFileName,"%s/%s%s%s",cmd.rootDir,"histCF",cmd.suffixOutFiles,EXTFILES);
    sprintf(gd.fpfnamehistrBinsFileName,"%s/%s%s%s",cmd.rootDir,"rbins",cmd.suffixOutFiles,EXTFILES);
    sprintf(gd.fpfnamehistXi2pcfFileName,"%s/%s%s%s",cmd.rootDir,cmd.histXi2pcfFileName,cmd.suffixOutFiles,EXTFILES);
    sprintf(gd.fpfnamehistZetaMFileName,"%s/%s%s",cmd.rootDir,cmd.histZetaMFileName,cmd.suffixOutFiles);
    sprintf(gd.fpfnamemhistZetaFileName,"%s/%s%s",cmd.rootDir,cmd.mhistZetaFileName,cmd.suffixOutFiles);
    sprintf(gd.fpfnameCPUFileName,"%s/cputime%s%s",cmd.rootDir,cmd.suffixOutFiles,EXTFILES);
#ifdef MPICODE
        }
#endif
}

//
// We must check the order of memory allocation and dealocation!!!
//
int EndRun(void)
{
//    char   buf[200];
    stream outstr;

// Here freeing memory lines should go

    int m;

//
// Freeing allocated must consider the allocated memory by others, Task != 0. See StartRun.
//
//#ifdef TPCF
//    for (m=1; m<=cmd.mchebyshev+1; m++)
//        gsl_matrix_complex_free(histZetaMatrix[m].histZetaM);
/*
    free(histZetaMatrix);
    gsl_matrix_complex_free(gd.histXi_gsl);
    free_dmatrix3D(gd.histZetaMsincos,1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix3D(gd.histZetaMsin,1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix3D(gd.histZetaMcos,1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix3D(gd.histZetaM,1,cmd.mchebyshev+1,1,cmd.sizeHistN,1,cmd.sizeHistN);
    free_dmatrix(gd.histXisin,1,cmd.mchebyshev+1,1,cmd.sizeHistN);
    free_dmatrix(gd.histXicos,1,cmd.mchebyshev+1,1,cmd.sizeHistN);
    free_dmatrix(gd.histXi,1,cmd.mchebyshev+1,1,cmd.sizeHistN);
#endif
    free_dmatrix3D(gd.histXi3pcf,1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    free_dvector(gd.histXi2pcf,1,cmd.sizeHistN);
    free_dmatrix3D(gd.histNNNSub,1,cmd.sizeHistN,1,cmd.sizeHistN,1,cmd.sizeHistTheta);
    free_dvector(gd.histNNN,1,cmd.sizeHistN);
// 2pcf
    free_dvector(gd.histNSubXi2pcf,1,cmd.sizeHistN);
//
    free_dvector(gd.histNSub,1,cmd.sizeHistN);
    free_dvector(gd.histCF,1,cmd.sizeHistN);
    free_dvector(gd.histN,1,cmd.sizeHistN);
*/
//B Set gsl uniform random :: If not needed globally this line have to go to testdata
//    gsl_rng_free (gd.r);

#ifdef MPICODE
        if(ThisTask==0) {
#endif
	fclose(gd.outlog);
#ifdef MPICODE
        }
#endif

#ifdef MPICODE
        if(ThisTask==0) {
#endif
    if (cmd.verbose >= 2) {
        printf("\nrSize \t\t= %lf\n", gd.rSize);
        printf("nbbcalc \t= %ld\n", gd.nbbcalc);
        printf("nbccalc \t= %ld\n", gd.nbccalc);
// BALLS
        printf("ncccalc \t= %ld\n", gd.ncccalc);
//
        printf("tdepth \t\t= %d\n", gd.tdepth);
        printf("ncell \t\t= %ld\n", gd.ncell);
        printf("cpusearch \t= %lf (be aware of the number of threads)\n", gd.cpusearch*cmd.numthreads);     // in minutes
        printf("cputotalinout \t= %lf\n", gd.cputotalinout);
    }

    if (cmd.verbose > 2) {
        real cpuTotal = gd.cpusearch*cmd.numthreads;
        real m1, m2, m3, m4;
        m1 =1.0;
        m2 =1.0e-5;
        m3 =1.0e-9;
        m4 =1.0e-12;
        outstr = stropen(gd.fpfnameCPUFileName, "a");
        fprintf(outstr,"%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                (double) cmd.nbody, cpuTotal,
                            m1*(double) cmd.nbody,
                            m2*((double) cmd.nbody)*rlog10((double) cmd.nbody),
                            m3*rsqr((double) cmd.nbody),
                            m4*rpow(((double) cmd.nbody), 3.0));
        fclose(outstr);
    }

    if (cmd.verbose > 0) {
        printf("\nFinal CPU time : %lf %s\n", CPUTIME - gd.cpuinit,PRNUNITOFTIMEUSED);// in minutes
//        printf("Final real time: %lf", (rcpu_time()-gd.cpurealinit));    // in minutes
        printf("Final real time: %lf", (rcpu_time()-gd.cpurealinit));    // in minutes
        printf(" %s\n\n", PRNUNITOFTIMEUSED);                               // Only work this way
    }
#ifdef MPICODE
        }
#endif

    return _SUCCESS_;
}

