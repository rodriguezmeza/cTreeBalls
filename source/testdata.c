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
            verb_print(cmd->verbose, "\nUsing simple cubic test model");
            testdata_sc(cmd, gd);
            break;
        case SIMPLECUBICRANDOM:
            verb_print(cmd->verbose, "\nUsing simple cubic random test model");
            testdata_sc_random(cmd, gd);
            break;
#if defined(THREEDIM)
        case UNITSPHERE:
            verb_print(cmd->verbose, "\nUsing unit sphere random test model");
#ifdef NR3
            testdata_unit_sphere_random_nr3(cmd, gd);
#else
            testdata_unit_sphere_random(cmd, gd);
#endif
            break;
#endif
        case NULLMODEL:
            verb_print(cmd->verbose, 
                       "\nNull test model. Using simple cubic random test model");
            testdata_sc_random(cmd, gd);
            break;
        case UNKNOWN:
            verb_print(cmd->verbose, "\nUnknown test model. Using simple cubic random test model");
            testdata_sc_random(cmd, gd);
            break;
        default:
            verb_print(cmd->verbose, "\nUnknown test model type %s", cmd->testmodel);
            verb_print(cmd->verbose, "\nUsing simple cubic random test model");
            testdata_sc_random(cmd, gd);
    }

    return SUCCESS;
}

local int model_string_to_int(struct  cmdline_data* cmd, struct  global_data* gd,
                               string model_str,int *model_int)
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
        *model_int=SIMPLECUBICRANDOM;
        gd->model_comment =
           "Unknown test model ... running deafult (simple-cubic-random)";
        verb_log_print(cmd->verbose_log, gd->outlog,
            "\n\nmodel_string_to_int: Unknown test model... %s",cmd->testmodel);
        verb_log_print(cmd->verbose_log, gd->outlog,
            "\n\trunning default test model (simple-cubic-random)...\n");
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
    bodyptr p;
    real tmass;
    real phi, theta;
    real ra, dec;

    gd->model_comment = "Random unit sphere Model";

    gd->Box[0] = 2.0;
    gd->Box[1] = 2.0;
    gd->Box[2] = 2.0;

    int ifile=0;
    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));

    mass = 1.0;
    weight = 1.0;
    tmass=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        phi    = 2.0 * PI_D * xrandom(0.0, 1.0);
        theta    = racos(1.0 - 2.0 * xrandom(0.0, 1.0));
        if (scanopt(cmd->options, "no-arfken")) {
            ra = theta;
            dec = PIO2 - phi;
            Pos(p)[0] = rcos(dec)*rcos(ra);
            Pos(p)[1] = rcos(dec)*rsin(ra);
            Pos(p)[2] = rsin(dec);
        } else {
            ra = phi;
            dec = theta;
            // Standard transformation. See Arfken
            Pos(p)[0] = rsin(dec)*rcos(ra);
            Pos(p)[1] = rsin(dec)*rsin(ra);
            Pos(p)[2] = rcos(dec);
        }
        Id(p) = p-bodytable[ifile]+1;
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;                         // use kapp-constant
        else {
            // a takahasi simulation with nside1024 gives:
            //  inputdata_ascii: average and std dev of kappa
            //  (12582912 particles) = 2.936728e-08 6.528639e-03
            Kappa(p) = grandom(2.936728e-08, 6.528639e-03);
//            Kappa(p) = grandom(1.0, 0.25);
        }
        tmass += Mass(p);
    }
    
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


#ifdef ADDONS
#include "testdata_include.h"
#endif

//E
