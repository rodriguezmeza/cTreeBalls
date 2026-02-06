/*==============================================================================
 HEADER: global_data.h		[cTreeBalls]
 Written by: Mario A. Rodriguez-Meza
 Starting date: april 2023
 Purpose: Definitions of global variables and parameters
 Language: C
 Use: '#include "global_data.h"
 Major revisions:
 ==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#ifndef _global_data_h
#define _global_data_h

#include "vectdefs.h"
#include "datastruc_defs.h"
#include "common_defs.h"

struct global_data{

    //B PXD functions
    // need to be public, not use "#ifdef PXD" directive
    real *rBins;
    real *vecPXD;
    real **matPXD;
    real *histZetaMFlatten;
    //E PXD functions

    real cpuinit;
    long cpurealinit;                               // get time of the day

	string headline0;
	string headline1;
	string headline2;
	string headline3;

    FILE *outlog;

	char mode[2];

//B Tree
    INTEGER ncell;
    int tdepth;
    INTEGER actmax;
    INTEGER sameposcount;
// BALLS
    INTEGER ncccalc;
    FILE *outnodelev;
    FILE *outbodylev;
    INTEGER nsmoothcount;
//
    INTEGER nbccalc;
    INTEGER nbbcalc;
    real rSize;                                     // Maximum r of the box,needed?

    real rSizeTable[MAXITEMS];                      // Maximum r of the box
    INTEGER ncellTable[MAXITEMS];
    INTEGER nbodyTable[MAXITEMS];
    int tdepthTable[MAXITEMS];
    INTEGER nnodescanlevTable[MAXITEMS];
    INTEGER nnodescanlev_rootTable[MAXITEMS];
//E

//B BALLS4
    INTEGER nnodescanlevTableB4[MAXITEMS];
//E

    real cputree;

// Tree:
     real Rcut;                                     // Cutoff radius
     real RcutSq;
    real cpusearch;
//

    int infilefmt_int;

#ifdef USEGSL
    gsl_rng * r;            // It is used r_gsl in globaldefs.h. Check!!!
#endif

#ifdef SINGLEP
    double Box[NDIM];
#else
    vector Box;
#endif

// -----------------------------------
    //B Histogram arrays PXD versions
    realptr histNNPXD;
    //E Histogram arrays PXD versions

//B Histogram arrays
    realptr histNN;
    realptr histCF;
    realptr histNNSub;
// 2pcf
    realptr histNNSubXi2pcf;
#ifdef SMOOTHPIVOT
    realptr histNNSubXi2pcftotal;
#endif
    real *histXi2pcf;
//B TPCF
    real ***histZetaMcos;
    real ***histZetaMsin;
    real ***histZetaMsincos;
// Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    real ***histZetaMcossin;
//E
    real ***histZetaM;
//
    real *histNNN;
    real ***histNNNSub;
    real *histXi2pcf_omp;
    real ***histXi3pcf;
//
    real **histXicos;
    real **histXisin;
#ifdef USEGSL
    gsl_matrix_complex *histXi_gsl;
#endif

//B To save total 3pcf
    //B look for where these two arrays are allocated...
    real **histZetaGcos;
    real **histZetaGsin;
    //E
    real ***histZetaGmRe;
    real ***histZetaGmIm;
//E

//E Histograms arrays
// -----------------------------------

    int searchMethod_int;

    real deltaR;
    real deltaRmin;
    real deltaRmax;
    real *deltaRV;
    real *ddeltaRV;

    real deltaPhi;

//B File pointers:
    char logfilePath[MAXLENGTHOFFILES];
    char tmpDir[MAXLENGTHOFFILES];

    char fpfnameOutputFileName[MAXLENGTHOFFILES];
    char fpfnamehistNNFileName[MAXLENGTHOFFILES];
    char fpfnamehistCFFileName[MAXLENGTHOFFILES];
    char fpfnamehistrBinsFileName[MAXLENGTHOFFILES];
    char fpfnamehistXi2pcfFileName[MAXLENGTHOFFILES];
//B TPCF
    char fpfnamehistZetaGFileName[MAXLENGTHOFFILES];
    char fpfnamehistZetaGmFileName[MAXLENGTHOFFILES];
    char fpfnamehistZetaMFileName[MAXLENGTHOFFILES];
    char fpfnamemhistZetaMFileName[MAXLENGTHOFFILES];
    char fpfnameCPUFileName[MAXLENGTHOFFILES];
//E
//E
    string model_comment;
    string input_comment;
    string output_comment;

//B save-restore
    int stopflag;
    INTEGER ip;
//E
    real cputotalinout;
    real cputotal;

    INTEGER bytes_tot;
    INTEGER bytes_tot_cells;

//B To see the bodies belonging to a cell:
    INTEGER nbodySel;
//E
    bool bh86, sw94;

    real i_deltaR;

    // 2pcf
    realptr histNNSubN2pcf;
#ifdef SMOOTHPIVOT
    realptr histNNSubN2pcftotal;
#endif
    real *histN2pcf;
#ifdef ADDONS
    char fpfnamehistN2pcfFileName[MAXLENGTHOFFILES];
#endif
    //

    char fnameData_kd[128];
    char fnameOut_kd[128];
    int input_format_kd;
    int use_tree_kd;
    int max_tree_order_kd;
    int max_tree_nparts_kd;
    int use_pm_kd;
    INTEGER n_objects_kd;
    float l_box_kd;
    float l_box_half_kd;

    int ninfiles;
    char infilenames[MAXITEMS][MAXLENGTHOFFILES];
    char infilefmtname[MAXITEMS][MAXLENGTHOFFMTFILES];
    int iCatalogs[MAXITEMS];

    int nsmooth[MAXITEMS];          // deprecated, will be deleted

    //B to control memory allocation/deallocation
    bool cmd_allocated;
    bool random_allocated;
    bool gd_allocated;
    bool gd_allocated_2;
    bool histograms_allocated;
    bool tree_allocated;
    bool bodytable_allocated;
    //E

    INTEGER nnode;
    INTEGER rnnode;

#define MAXLEVEL  32
    real Rcell[MAXLEVEL];           // used only in treeload's scanLevel routine
#undef MAXLEVEL

#ifdef BALLS4SCANLEV
    bool flagBalls4Scanlevel;
#endif

    real rsmooth[MAXITEMS];
    bool rsmoothFlag;

    int irsmooth;

    //B to control memory allocation/deallocation
    bool flagPrint;
    bool rootDirFlag;
    bool searchMethodFlag;
    bool infilefmtFlag;
    bool infileFlag;
    bool outfilefmtFlag;
    bool outfileFlag;
    bool rsmoothFlagFree;
    bool rootDirFlagFree;
    bool iCatalogsFlag;
    bool histNNFileNameFlag;
    bool histXi2pcfFileNameFlag;
    bool histZetaFileNameFlag;
    bool suffixOutFilesFlag;
    bool testmodelFlag;
    bool preScriptFlag;
    bool posScriptFlag;
    bool optionsFlag;
    //E

//B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    INTEGER pivotCount;
#endif
//E

    //B correction 2025-04-06
    INTEGER pivotNumber;
    //E

    INTEGER stepState;

    short coordTag;

//B socket:
#ifdef ADDONS
#include "globaldefs_include_03.h"
#endif
//E

};

#endif // ! _global_data_h

