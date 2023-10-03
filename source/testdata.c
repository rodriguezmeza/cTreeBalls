/* ==============================================================================
!	MODULE: testdata.c				[cTreeBalls]					    		!
!	Written by: M.A. Rodriguez-Meza 											!
!	Date: april 2023				            								!
!	Purpose: Setting some test models											!
!	Language: C																	!
!	Use: 'testdata();'															!
!	Major revisions:															!
!==============================================================================*/

#include "globaldefs.h"
local void testdata_sc_random(void);
local void testdata_sc(void);

local void Compute_nbody(void);
local void Compute_Box_Units(void);
//local void Compute_Parameters(void);
//void Compute_Parameters(void);

local vector N;                             // Needed in testdata_sc
local real weight;

#define SIMPLECUBICRANDOM   0
#define SIMPLECUBIC         1
#define UNKNOWN             2
#define NULLMODEL           3

local void model_string_to_int(string model_str,int *model_int);

int testdata(void)
{
    int model_int;
    model_string_to_int(cmd.testmodel, &model_int);
#ifdef MPICODE
    if (ThisTask==0) {                      // Only task 0 can get points catalog
#endif
    switch (model_int){
        case SIMPLECUBIC:
            verb_print(cmd.verbose, "\nUsing simple cubic test model");
            testdata_sc();
            break;
        case SIMPLECUBICRANDOM:
            verb_print(cmd.verbose, "\nUsing simple cubic random test model");
            testdata_sc_random();
            break;
        case NULLMODEL:
            verb_print(cmd.verbose, "\nNull test model. Using simple cubic random test model");
            testdata_sc_random();
            break;
        case UNKNOWN:
            verb_print(cmd.verbose, "\nUnknown test model. Using simple cubic random test model");
            testdata_sc_random();
            break;
        default:
            verb_print(cmd.verbose, "\nUnknown test model type %s", cmd.testmodel);
            verb_print(cmd.verbose, "\nUsing simple cubic random test model");
            testdata_sc_random();
    }
#ifdef MPICODE
    }
#endif

    return _SUCCESS_;
}

local void model_string_to_int(string model_str,int *model_int)
{
    *model_int = -1;
    if (strcmp(model_str,"simple-cubic") == 0)          *model_int=SIMPLECUBIC;
    if (strcmp(model_str,"simple-cubic-random") == 0)   *model_int=SIMPLECUBICRANDOM;
    if (strnull(model_str)) {
        *model_int = NULLMODEL;
        gd.model_comment =
               "null test model ... running deafult (simple-cubic-random)";
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\n\tmodel_string_to_int: default test model (simple-cubic-random)...\n");
    }
    if (*model_int == -1) {
        *model_int=SIMPLECUBICRANDOM;
        gd.model_comment =
           "Unknown test model ... running deafult (simple-cubic-random)";
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\n\nmodel_string_to_int: Unknown test model... %s",cmd.testmodel);
        verb_log_print(cmd.verbose_log, gd.outlog,
                       "\n\trunning default test model (simple-cubic-random)...\n");
    }
}

#undef SIMPLECUBICRANDOM
#undef SIMPLECUBIC
#undef UNKNOWN
#undef NULLMODEL

// If we change seed or the list of randoms, there will be changes in
//  the configurations of points... and results. That is way line 119 below.
local void testdata_sc_random(void)
{
    bodyptr p;
    real Ekin;
    real tweight;
	int k;

    gd.model_comment = "Random Cubic Box Model";

    Compute_nbody();
    Compute_Box_Units();

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    tweight=0.0;
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
		Id(p) = p-bodytab+1;
        Type(p) = BODY;
        Weight(p) = weight;
        DO_COORD(k) {
            Pos(p)[k]    = xrandom(0.0, gd.Box[k]);    // Box lower edge is (0,0,0)
//            Pos(p)[k]    = xrandom(-0.5*gd.Box[k], 0.5*gd.Box[k]);    // Box center is (0,0,0)
        }
		DO_COORD(k)
            Ekin = grandom(0.0,1.0);            // Dummy line to reproduce results
        if (scanopt(cmd.options, "kappa-constant"))
            Kappa(p) = 2.0;
        else
            Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd.Box[0]/20.0))
                            * rsin(2.0*PI*Pos(p)[1]/(gd.Box[1]/20.0));

		tweight += Weight(p);
    }

    verb_log_print(cmd.verbose_log, gd.outlog, "\nCreated bodies = %d",cmd.nbody);
}

local void testdata_sc(void)
{
    bodyptr p;
    real tweight;
//    int n;
	int k, i,j,l;
	vector delta, c;
	real f, os;

    gd.model_comment = "Simple Cubic Box Model";

    bodytab = (bodyptr) allocate(cmd.nbody * sizeof(body));

    Compute_nbody();
    Compute_Box_Units();

	os =0.05*gd.Box[0];						// offset
	f = 1.0 - 2.0*os/gd.Box[0];				// scaling delta to 1 - 2 os/Lk

	DO_COORD(k)
		delta[k]=f*gd.Box[k]/((real) (N[k]-1));

	if (NDIM == 3)
		fprintf(gd.outlog,"\ndelta: %g %g %g\n",
				delta[0],delta[1],delta[2]);
	else if (NDIM == 2)
		fprintf(gd.outlog,"\ndelta: %g %g\n",
				delta[0],delta[1]);
	else error("\n\nWrong NDIM!\n\n");

	p = bodytab;
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

    verb_log_print(cmd.verbose_log, gd.outlog, "\nCreated bodies = %d",cmd.nbody);

    tweight=0.0;
	DO_BODY(p, bodytab, bodytab+cmd.nbody) {
		Id(p) = p-bodytab+1;
        Type(p) = BODY;
        Weight(p) = weight;
//        Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/N[0])
//                        * rsin(2.0*PI*Pos(p)[1]/N[1]);
        if (scanopt(cmd.options, "kappa-constant"))
            Kappa(p) = 2.0;
        else
            Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd.Box[0]/20.0))
                            * rsin(2.0*PI*Pos(p)[1]/(gd.Box[1]/20.0));
        tweight += Weight(p);
    }
}

local void Compute_nbody(void)
{
#ifdef MPICODE
        if(ThisTask==0) {
#endif
    fprintf(gd.outlog,"\n\nInitial nbody=%ld",cmd.nbody);
#ifdef MPICODE
        }
#endif
    if (NDIM == 3) {
        N[0] = (int) rpow(cmd.nbody,1./3.);        // nbody = N^3
        N[1] = N[0]; N[2] = N[0];
        cmd.nbody = N[0]*N[1]*N[2];
    } else if (NDIM == 2) {
        N[0] = (int) rsqrt(cmd.nbody);            // nbody = N^2
        N[1] = N[0];
        cmd.nbody = N[0]*N[1];
    } else error("\n\nWrong NDIM!\n\n");
}

local void Compute_Box_Units(void)
{
    weight = 1.0;        // Units system
    if (NDIM == 3) {                        // Size of box side x
        gd.Box[0] = cmd.lengthBox;
        gd.Box[1] = gd.Box[0]; gd.Box[2] = gd.Box[0];
    } else if (NDIM == 2) {                    // Size of box side x
        gd.Box[0] = cmd.lengthBox;
        gd.Box[1] = gd.Box[0];
    } else error("\n\nWrong NDIM!\n\n");
}


// This function makes sure that all bodies coordinates (Pos) are
// mapped onto the interval [0, Box].
#ifdef PERIODIC
global void doBoxWrapping(void)
{
    bodyptr p;
    int j;

    DO_BODY(p, bodytab, bodytab+cmd.nbody)
        DO_COORD(j) {
            while(Pos(p)[j] < 0)
                Pos(p)[j] += gd.Box[j];

            while(Pos(p)[j] > gd.Box[j])
                Pos(p)[j] -= gd.Box[j];
      }
}
#endif
