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

#ifndef _global_data_h
#define _global_data_h

#include "vectdefs.h"
#include "datastruc_defs.h"
#include "common_defs.h"

struct global_data{

    //B PXD functions
    real *rBins;
    real *histZetaMFlatten;
    //E

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
    real rSize;                                     // Maximum r of the box

    real rSizeTable[MAXITEMS];                      // Maximum r of the box
    INTEGER ncellTable[MAXITEMS];
    INTEGER nbodyTable[MAXITEMS];
    int tdepthTable[MAXITEMS];
    INTEGER nnodescanlevTable[MAXITEMS];
    INTEGER nnodescanlev_rootTable[MAXITEMS];
//E

    real cputree;

// Tree:
     real Rcut;                                     // Cutoff radius
     real RcutSq;
    real cpusearch;
//
#ifdef CELLMETHOD
// Cell search
     vectorI cells;
    INTEGER *cellList;
//
#endif

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
//B Histogram arrays
    realptr histNN;
    realptr histCF;
    realptr histNNSub;
// 2pcf
    realptr histNNSubXi2pcf;
//B kappa Avg Rmin
    realptr histNNSubXi2pcftotal;
//E
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
// TPCF
    real **histXi;
    real **histXicos;
    real **histXisin;
#ifdef USEGSL
    gsl_matrix_complex *histXi_gsl;
#endif

//B To save total 3pcf
    real **histZetaGcos;
    real **histZetaGsin;
    real ***histZetaGmRe;
    real ***histZetaGmIm;
//E
    
//B To save total 3pcf shear
    real *histXitt;
    real *histXixx;
    real *histXitx;
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
    char outputDir[MAXLENGTHOFFILES];
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
    INTEGER nbodysm;
    INTEGER nbodybf;

    bool bh86, sw94;

    real i_deltaR;

#ifdef ADDONS
#include "global_data_include.h"
#endif

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
    char *infilenames[MAXITEMS];
    char *infilefmtname[MAXITEMS];
    int iCatalogs[MAXITEMS];

    int nsmooth[MAXITEMS];
    INTEGER nnode;
    INTEGER rnnode;

    int scanLevelMin[MAXITEMS];

    char nodesfilePath[MAXLENGTHOFFILES];
    int nnodescanlev;
// Root nodes:
    int nnodescanlev_root;
#define MAXLEVEL  32
    real Rcell[MAXLEVEL];           // used only in treeload's scanLevel routine
#undef MAXLEVEL
//
    char bodiesfilePath[MAXLENGTHOFFILES];

    bool flagSmoothCellMin;
    bool flagSmooth;
    bool flagSetNbNoSel;
//B BUCKET
    real rminCell[2];               // used only in treeload's scanLevel routine
                                    //  and in search_balls_omp;
                                    //      search_balls_kk_omp;
                                    //      search_octree_kk_balls4_omp;
                                    //      search_octree_kkk_balls4_omp;
//E

    real rsmooth[MAXITEMS];
    bool rsmoothFlag;

    int irsmooth;

#ifdef ADDONS
#include "globaldefs_include_03.h"
#endif

};

#endif // ! _global_data_h

