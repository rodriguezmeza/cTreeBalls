/* ==============================================================================
!	MODULE: testdata.c				[cTreeBalls]					    		!
!	Written by: M.A. Rodriguez-Meza 											!
!	Date: april 2023				            								!
!	Purpose: Setting some test models											!
!	Language: C																	!
!	Use: 'testdata();'															!
!	Major revisions:															!
!==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
local int testdata_sc_random(struct  cmdline_data* cmd, struct  global_data* gd);
local int testdata_sc(struct  cmdline_data* cmd, struct  global_data* gd);
#if defined(THREEDIM)
local int testdata_unit_sphere_random(struct  cmdline_data* cmd, struct  global_data* gd);
#endif

local int Compute_nbody(struct  cmdline_data* cmd, struct  global_data* gd);
local int Compute_Box_Units(struct  cmdline_data* cmd, struct  global_data* gd);

local vector N;                                     // Needed in testdata_sc
local real mass;
local real weight;

#define SIMPLECUBICRANDOM   0
#define SIMPLECUBIC         1
#define UNKNOWN             2
#define NULLMODEL           3
#define UNITSPHERE          4

local int model_string_to_int(struct  cmdline_data* cmd, struct  global_data* gd,
                              string model_str,int *model_int);

#ifdef NR3
#include "xran.h"
#endif

int TestData(struct  cmdline_data* cmd, struct  global_data* gd)
{
    int model_int;

    model_string_to_int(cmd, gd, cmd->testmodel, &model_int);
    switch (model_int){
        case SIMPLECUBIC:
            verb_print(cmd->verbose, "\nUsing simple cubic test model\n");
            verb_print(cmd->verbose, " with nbody = %ld", cmd->nbody);
            testdata_sc(cmd, gd);
            break;
        case SIMPLECUBICRANDOM:
            verb_print(cmd->verbose, "\nUsing simple cubic random test model");
            verb_print(cmd->verbose, " with nbody = %ld\n", cmd->nbody);
            testdata_sc_random(cmd, gd);
            break;
#if defined(THREEDIM)
        case UNITSPHERE:
            verb_print(cmd->verbose, "\nUsing unit sphere random test model\n");
            verb_print(cmd->verbose, " with nbody = %ld\n", cmd->nbody);
#ifdef NR3
            testdata_unit_sphere_random_nr3(cmd, gd);
#else
            testdata_unit_sphere_random(cmd, gd);
#endif
            break;
#endif
        case NULLMODEL:
            verb_print(cmd->verbose, 
                       "\nNull test model. Using default test model");
            verb_print(cmd->verbose, " with nbody = %ld", cmd->nbody);
//            testdata_sc_random(cmd, gd);
            if (cmd->usePeriodic==TRUE) {
                verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                testdata_sc_random(cmd, gd);
            } else {
                if (cmd->computeTPCF==TRUE) {
                    verb_print(cmd->verbose, "\nUsing unit sphere random test model\n");
                    testdata_unit_sphere_random(cmd, gd);
                } else {
                    verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                    testdata_sc_random(cmd, gd);
                }
            }
            break;
        case UNKNOWN:
            verb_print(cmd->verbose,
                       "\nUnknown test model. Using default test model");
            verb_print(cmd->verbose, " with nbody = %ld", cmd->nbody);
//            testdata_sc_random(cmd, gd);
            verb_print(cmd->verbose, "\nDefault test model type %s", cmd->testmodel);
            if (cmd->usePeriodic==TRUE) {
                verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                testdata_sc_random(cmd, gd);
            } else {
                if (cmd->computeTPCF==TRUE) {
                    verb_print(cmd->verbose, "\nUsing unit sphere random test model");
                    verb_print(cmd->verbose, " with nbody = %ld\n", cmd->nbody);
                    testdata_unit_sphere_random(cmd, gd);
                } else {
                    verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                    testdata_sc_random(cmd, gd);
                }
            }
            break;
        default:
            verb_print(cmd->verbose, "\nDefault test model type %s", cmd->testmodel);
            verb_print(cmd->verbose, " with nbody = %ld", cmd->nbody);
            if (cmd->usePeriodic==TRUE) {
                verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                testdata_sc_random(cmd, gd);
            } else {
                if (cmd->computeTPCF==TRUE) {
                    verb_print(cmd->verbose, "\nUsing unit sphere random test model\n");
                    testdata_unit_sphere_random(cmd, gd);
                } else {
                    verb_print(cmd->verbose, "\nUsing simple cubic random test model\n");
                    testdata_sc_random(cmd, gd);
                }
            }
    }

    return SUCCESS;
}

local int model_string_to_int(struct  cmdline_data* cmd,
                              struct  global_data* gd,
                              string model_str, int *model_int)
{
    *model_int = -1;
    if (strcmp(model_str,"simple-cubic") == 0)          *model_int=SIMPLECUBIC;
    if (strcmp(model_str,"simple-cubic-random") == 0)
        *model_int=SIMPLECUBICRANDOM;
#if defined(THREEDIM)
    if (strcmp(model_str,"unit-sphere-random") == 0)
        *model_int=UNITSPHERE;
#endif
    if (strnull(model_str)) {
        *model_int = NULLMODEL;
        gd->model_comment =
               "null test model ... running deafult (simple-cubic-random)";
        verb_log_print(cmd->verbose_log, gd->outlog,
        "\n\tmodel_string_to_int: default test model (simple-cubic-random)...\n");
    }
    if (*model_int == -1) {
        *model_int=UNKNOWN;
        gd->model_comment =
        "Unknown test model ... running deafult model";
//        "Unknown test model ... running deafult (simple-cubic-random)";
        verb_log_print(cmd->verbose_log, gd->outlog,
            "\n\nmodel_string_to_int: Unknown test model... %s",cmd->testmodel);
//        verb_log_print(cmd->verbose_log, gd->outlog,
//            "\n\trunning default test model (simple-cubic-random)...\n");
        verb_log_print(cmd->verbose_log, gd->outlog,
            "\n\trunning default test model...\n");
    }

    return SUCCESS;
}

#undef SIMPLECUBICRANDOM
#undef SIMPLECUBIC
#undef UNKNOWN
#undef NULLMODEL
#undef UNITSPHERE

// If we change seed or the list of randoms, there will be changes in
//  the configurations of points... and results. That is way line 123 below.
local int testdata_sc_random(struct  cmdline_data* cmd,
                             struct  global_data* gd)
{
    bodyptr p;
    real Ekin;
    real tmass;
	int k;

    gd->model_comment = "Random Cubic Box Model";

    Compute_nbody(cmd,gd);
    Compute_Box_Units(cmd,gd);

    int ifile=0;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    mass = 1.0;
    weight = 1.0;
    tmass=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
		Id(p) = p-bodytable[ifile]+1;
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        DO_COORD(k) {
            Pos(p)[k] = xrandom(0.0, gd->Box[k]);   // Box lower edge is (0,0,0)
//            Pos(p)[k] = xrandom(-0.5*gd->Box[k], 0.5*gd->Box[k]);
                                                    // Box center is (0,0,0)
        }
		DO_COORD(k)
            Ekin = grandom(0.0,1.0);                // Dummy line to reproduce
                                                    //  results
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        else
            Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd->Box[0]/20.0))
                            * rsin(2.0*PI*Pos(p)[1]/(gd->Box[1]/20.0));

		tmass += Mass(p);
    }
    
    verb_log_print(cmd->verbose_log, gd->outlog, "\nCreated bodies = %d",cmd->nbody);

    return SUCCESS;
}

local int testdata_sc(struct  cmdline_data* cmd, struct  global_data* gd)
{
    bodyptr p;
    real tmass;
	int k, i,j,l;
	vector delta, c;
	real f, os;

    gd->model_comment = "Simple Cubic Box Model";

    int ifile=0;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    mass = 1.0;
    weight = 1.0;

    Compute_nbody(cmd, gd);
    Compute_Box_Units(cmd, gd);

	os =0.05*gd->Box[0];						    // offset
	f = 1.0 - 2.0*os/gd->Box[0];				    // scaling delta
                                                    //  to 1 - 2 os/Lk
	DO_COORD(k)
		delta[k]=f*gd->Box[k]/((real) (N[k]-1));

	if (NDIM == 3)
		fprintf(gd->outlog,"\ndelta: %g %g %g\n",
				delta[0],delta[1],delta[2]);
	else if (NDIM == 2)
		fprintf(gd->outlog,"\ndelta: %g %g\n",
				delta[0],delta[1]);
	else error("\n\nWrong NDIM!\n\n");

	p = bodytable[ifile];
	if (NDIM == 3)
		for (i=1; i<=N[0]; i++) {
			c[0] = ((real) (i-1))*delta[0] + os;
			for (j=1; j<=N[1]; j++) {
				c[1] = ((real) (j-1))*delta[1] + os;
				for (l=1; l<=N[2]; l++) {
					c[2] = ((real) (l-1))*delta[2] + os;
					SETV(Pos(p),c);
					++p;
				}
			}
		}
	else if (NDIM == 2)
		for (i=1; i<=N[0]; i++) {
			c[0] = ((real) (i-1))*delta[0] + os;
			for (j=1; j<=N[1]; j++) {
				c[1] = ((real) (j-1))*delta[1] + os;
				SETV(Pos(p),c);
				++p;
			}
		}
	else error("\n\nWrong NDIM!\n\n");

    verb_log_print(cmd->verbose_log, gd->outlog, "\nCreated bodies = %d",cmd->nbody);

    tmass=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
		Id(p) = p-bodytable[ifile]+1;
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        else
            Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd->Box[0]/20.0))
                            * rsin(2.0*PI*Pos(p)[1]/(gd->Box[1]/20.0));
        tmass += Mass(p);
    }

    return SUCCESS;
}

#if defined(THREEDIM)
local int testdata_unit_sphere_random(struct  cmdline_data* cmd, 
                                      struct  global_data* gd)
{
// Set of standard parameters to define histogram bins for a unit-sphere
//  rangeN =    0.0633205   # 8 arcmin aprox
//  rminHist =  0.00213811  # 200 arcmin aprox
//  sizeHistN = 20
//
    bodyptr p;
    real tmass;
    real phi, theta;
    real ra, dec;

    gd->model_comment = "Random unit sphere Model";

    gd->Box[0] = 2.0;
    gd->Box[1] = 2.0;
    gd->Box[2] = 2.0;

    verb_print(cmd->verbose,
               "using sizeHistN, rangeN and rminHist for a unit-sphere:");
    verb_print(cmd->verbose,
               "%d %g %g\n",
               cmd->sizeHistN, cmd->rangeN, cmd->rminHist);
    //E

    int ifile=0;
    int nbody = cmd->nbody;

    mass = 1.0;
    weight = 1.0;
    tmass=0.0;
    
    bodyptr bodytabtmp;
    bodytabtmp = (bodyptr) allocate(nbody * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for temporal particle (%ld) storage.\n",
               nbody*sizeof(body)*INMB, nbody);

    real xmin, ymin, zmin;
    real xmax, ymax, zmax;
    xmin=0., ymin=0., zmin=0.;
    xmax=0., ymax=0., zmax=0.;

    INTEGER i;
    INTEGER iselect = 0;

    DO_BODY(p, bodytabtmp, bodytabtmp+nbody) {
        Update(p) = FALSE;
        phi    = 2.0 * PI_D * xrandom(0.0, 1.0);
        theta    = racos(1.0 - 2.0 * xrandom(0.0, 1.0));
        if (scanopt(cmd->options, "patch")) {
            //B
            if (cmd->thetaL < theta && theta < cmd->thetaR) {
                if (cmd->phiL < phi && phi < cmd->phiR) {
                    iselect++;
//                    spherical_to_cartesians(cmd, gd, theta, phi, Pos(p));
                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));
                    if (!scanopt(cmd->options, "kappa-constant")) {
                        // a takahasi simulation with nside1024 gives:
                        //  inputdata_ascii: average and std dev of kappa
                        //  (12582912 particles) = 2.936728e-08 6.528639e-03
                        Kappa(p) = grandom(2.936728e-08, 6.528639e-03);
//                        Kappa(p) = grandom(1.0, 0.25);
                    } else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }
                    Type(p) = BODY;
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;
                    Update(p) = TRUE;
                    xmin = Pos(p)[0];
                    ymin = Pos(p)[1];
                    zmin = Pos(p)[2];
                    xmax = Pos(p)[0];
                    ymax = Pos(p)[1];
                    zmax = Pos(p)[2];
                } // ! in phi range
            } // in ! in theta range
            //E
        } else { // ! all
            //B
            iselect++;
//            spherical_to_cartesians(cmd, gd, theta, phi, Pos(p));
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
            if (!scanopt(cmd->options, "kappa-constant")) {
                // a takahasi simulation with nside1024 gives:
                //  inputdata_ascii: average and std dev of kappa
                //  (12582912 particles) = 2.936728e-08 6.528639e-03
                Kappa(p) = grandom(2.936728e-08, 6.528639e-03);
                //                        Kappa(p) = grandom(1.0, 0.25);
            } else {
                Kappa(p) = 2.0;
                if (scanopt(cmd->options, "kappa-constant-one"))
                    Kappa(p) = 1.0;
            }
            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;
            Update(p) = TRUE;
            xmin = Pos(p)[0];
            ymin = Pos(p)[1];
            zmin = Pos(p)[2];
            xmax = Pos(p)[0];
            ymax = Pos(p)[1];
            zmax = Pos(p)[2];
            //E
        } // ! all
    } // ! end loop p

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    real kavg = 0;
    INTEGER ij=0;
    DO_BODY(q, bodytabtmp, bodytabtmp+nbody) {
        if(Update(q)) {
            p = bodytable[ifile]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytable[ifile]+1;
            xmin = MIN(xmin,Pos(p)[0]);
            ymin = MIN(ymin,Pos(p)[1]);
            zmin = MIN(zmin,Pos(p)[2]);
            xmax = MAX(xmax,Pos(p)[0]);
            ymax = MAX(ymax,Pos(p)[1]);
            zmax = MAX(zmax,Pos(p)[2]);
            ij++;
            kavg += Kappa(p);
            tmass += Mass(p);
        }
    }

    verb_print(cmd->verbose,
               "\n\ttestdata_unit_sphere_random: min and max of x = %f %f\n",
               xmin, xmax);
    verb_print(cmd->verbose,
               "\ttestdata_unit_sphere_random: min and max of y = %f %f\n",
               ymin, ymax);
    verb_print(cmd->verbose,
               "\ttestdata_unit_sphere_random: min and max of z = %f %f\n",
               zmin, zmax);

    free(bodytabtmp);
    verb_print(cmd->verbose,
            "\nFreed %g MByte for temporal particle (%ld) storage.\n",
               nbody*sizeof(body)*INMB,nbody);

    if (scanopt(cmd->options, "all"))
        verb_print(cmd->verbose,
            "\n\ttestdata_unit_sphere_random: selected read points and nbody: %ld %ld\n",
            iselect, cmd->nbody);
    else
        verb_print(cmd->verbose,
            "\n\ttestdata_unit_sphere_random: selected read points = %ld\n",iselect);

    verb_print(cmd->verbose,
            "testdata_unit_sphere_random: average of kappa (%ld particles) = %le\n",
            cmd->nbody, kavg/((real)cmd->nbody) );


    verb_log_print(cmd->verbose_log, gd->outlog, "\nCreated bodies = %d",cmd->nbody);

    return SUCCESS;
}
#endif

local int Compute_nbody(struct  cmdline_data* cmd, struct  global_data* gd)
{
    fprintf(gd->outlog,"\n\nInitial nbody=%ld",cmd->nbody);
    if (NDIM == 3) {
        N[0] = (int) rpow(cmd->nbody,1./3.);        // nbody = N^3
        N[1] = N[0]; N[2] = N[0];
        cmd->nbody = N[0]*N[1]*N[2];
    } else if (NDIM == 2) {
        N[0] = (int) rsqrt(cmd->nbody);             // nbody = N^2
        N[1] = N[0];
        cmd->nbody = N[0]*N[1];
    } else error("\n\nWrong NDIM!\n\n");

    return SUCCESS;
}

local int Compute_Box_Units(struct  cmdline_data* cmd, struct  global_data* gd)
{
    mass = 1.0;                                     // Units system
    if (NDIM == 3) {                                // Size of box side x
        gd->Box[0] = cmd->lengthBox;
        gd->Box[1] = gd->Box[0]; gd->Box[2] = gd->Box[0];
    } else if (NDIM == 2) {                         // Size of box side x
        gd->Box[0] = cmd->lengthBox;
        gd->Box[1] = gd->Box[0];
    } else error("\n\nWrong NDIM!\n\n");

    return SUCCESS;
}


// Seams that this routine is not used... delete it!!!
//  or move it to cballsutils...
//
// This function makes sure that all bodies coordinates (Pos) are
// mapped onto the interval [0, Box].
#ifdef PERIODIC
global int doBoxWrapping(struct  cmdline_data* cmd, struct  global_data* gd)
{
    bodyptr p;
    int j;

    int ifile=0;

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        DO_COORD(j) {
            while(Pos(p)[j] < 0)
                Pos(p)[j] += gd->Box[j];
            
            while(Pos(p)[j] > gd->Box[j])
                Pos(p)[j] -= gd->Box[j];
        }
    }

    return SUCCESS;
}
#endif


//B socket:
#ifdef ADDONS
#include "testdata_include.h"
#endif
//E

//E
