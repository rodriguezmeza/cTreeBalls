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

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"

local int inputdata_ascii(struct cmdline_data*, struct  global_data*,
                          string filename, int);
local int inputdata_ascii_all(struct cmdline_data*, struct  global_data*,
                          string filename, int);
local int inputdata_bin(struct cmdline_data*, struct  global_data*,
                        string filename, int);
local int inputdata_bin_all(struct cmdline_data*, struct  global_data*,
                        string filename, int);
local int inputdata_takahashi(struct cmdline_data*, struct  global_data*,
                             string filename, int);

local int outputdata(struct cmdline_data*, struct  global_data*,
                     bodyptr, INTEGER nbody);
local int outputdata_ascii(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
local int outputdata_ascii_all(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
local int outputdata_bin(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);
local int outputdata_bin_all(struct cmdline_data*, struct  global_data*,
                         bodyptr, INTEGER);

//B socket:
#ifdef ADDONS
#include "cballsio_include_00.h"
#endif
//E

local int outfilefmt_string_to_int(string,int *);
local int outfilefmt_int;


/*
 InputData routine:

 To be called by StartRun_Common in startrun.c:
    InputData(cmd, gd, gd->infilenames[ifile], ifile);

 This routine is in charge of reading catalog of data
    to be analyzed

 Arguments:
    * `cmd`:        Input: structure cmdline_data pointer
    * `gd`:         Input: structure global_data pointer
    * `filename`:   Input: catalog of data filename
    * `ifile`:      Input: catalog file tag
 Return (the error status):
    int SUCCESS or FAILURE
 */
int InputData(struct cmdline_data* cmd,
              struct  global_data* gd, string filename, int ifile)
{
    string routineName = "InputData";
    double cpustart = CPUTIME;

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\n%s: reading data catalog...\n", routineName);
    switch(gd->infilefmt_int) {
        case INCOLUMNS:
            verb_print_normal_info(cmd->verbose,
                                cmd->verbose_log, gd->outlog,
                                "\n\tInput in columns (ascii) format...\n");
            inputdata_ascii(cmd, gd, filename, ifile); break;
        case INCOLUMNSALL:
            verb_print_normal_info(cmd->verbose,
                            cmd->verbose_log, gd->outlog,
                            "\tInput in columns (ascii) all format...\n");
            inputdata_ascii_all(cmd, gd, filename, ifile); break;
        case INNULL:
            verb_print_normal_info(cmd->verbose,
                        cmd->verbose_log, gd->outlog,
                        "\n\t(Null) Input in columns (ascii) format...\n");
            inputdata_ascii(cmd, gd, filename, ifile); break;
        case INCOLUMNSBIN:
            verb_print_normal_info(cmd->verbose,
                                   cmd->verbose_log, gd->outlog,
                                   "\n\tInput in binary format...\n");
            inputdata_bin(cmd, gd, filename, ifile); break;
        case INCOLUMNSBINALL:
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "\n\tInput in binary-all format...\n");
            inputdata_bin_all(cmd, gd, filename, ifile); break;
        case INTAKAHASHI:
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "\n\tInput in takahashi format...\n");
            class_call_cballs(inputdata_takahashi(cmd, gd, filename, ifile),
                       errmsg, errmsg);
            break;

//B socket:
#ifdef ADDONS
#include "cballsio_include_01.h"
#endif
//E

        default:
            verb_print(cmd->verbose,
                       "\n\tInput: Unknown input format (%s)...",cmd->infilefmt);
            if (scanopt(cmd->infilefmt, "fits")) {
                verb_print(cmd->verbose,
                "\n\tInput: set CFITSIOON = 1 in ");
                verb_print(cmd->verbose,
                "addons/Makefile_addons_settings file...");
                verb_print(cmd->verbose,
                "\n\t\t and compile again ($ make clean; make).");
                error("\n\tgoing out...\n");
            }
            verb_print(cmd->verbose,
                       "\n\tInput in default columns (ascii) format...\n");
            class_call_cballs(inputdata_ascii(cmd, gd, filename, ifile),
                                              errmsg, errmsg);
            break;
    }
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\tdone reading.\n");

    gd->cputotalinout += CPUTIME - cpustart;

    // DEBUG WARNING!!
    //B There is a Segmentation fault: 11 if run as:
    //  cballs in=./scripts/Abraham/kappa_nres12_zs9NS256r000.txt \
    //  options=header-info
    // (and works with 'options=0')
#ifdef DEBUG
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
               "\tinputdata :: reading time = %f\n",CPUTIME-cpustart);
#else
    // but if comment above line and use this... works (?)
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                "\tinputdata :: reading time = %f %s\n",
                CPUTIME - cpustart, PRNUNITOFTIMEUSED);
#endif
    // seems it is associated to the design of
    //      void verb_print(int verbose, string fmt, ...)
    // need to check this carefully
    //E

    return SUCCESS;
}

//B gives treeload expandbox: rSize = 4.000000
//#define EPSILON 1.0E-7
//E
#define EPSILON 1.0E-8

/*
 InputData_all_in_one routine:

 To be called by StartRun_Common in startrun.c

 This routine is in charge of reading catalogs of data
    to be analyzed and then combine all of them in one catalog

 Arguments:
    * `cmd`:        Input: structure cmdline_data pointer
    * `gd`:         Input: structure global_data pointer
 Return (the error status):
    int SUCCESS or FAILURE
 */
global int InputData_all_in_one(struct cmdline_data* cmd,
                               struct  global_data* gd)
{
    string routineName = "InputData_all_in_one";
    bodyptr p,q;
    INTEGER i, l, ij;
    int j;
    int k;
    bodyptr bodytabtmp;

    cmd->nbody = 0;
    for (j=0; j<gd->ninfiles; j++)
        cmd->nbody += gd->nbodyTable[j];
    bodytabtmp = (bodyptr) allocate(cmd->nbody * sizeof(body));

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
    "\n\t%s: Allocated %g MByte for tmp bodytable storage (total bodies=%ld).\n",
    routineName, cmd->nbody*sizeof(body)*INMB, cmd->nbody);

    INTEGER iselect = 0;
    l=0;
    ij=0;
    for (j=0; j<gd->ninfiles; j++) {
        i = 0;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\t%s: processing file %d... with ",
                               routineName, j);
        DO_BODY(q, bodytable[j], bodytable[j]+gd->nbodyTable[j]) {
            p = bodytabtmp + ij;
                DO_COORD(k){
                    Pos(p)[k] = Pos(q)[k];
                    Pos(p)[k] +=
                        EPSILON*grandom(0.0, 1.0);  // use 0.01*gd->rSizeTable[j]
                }
            Kappa(p) = Kappa(q);
            if (scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = 2.0;
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
            
            Mass(p) = Mass(q);
            
#ifdef THREEPCFSHEAR
            Gamma1(p) = Gamma1(q);
            Gamma2(p) = Gamma2(q);
#endif
            Weight(p) = Weight(q);
            Type(p) = BODY;
            Id(p) = p-bodytabtmp+1;
            Mask(p) = Mask(q);
            if (Mask(p) == 0) {
                iselect++;
            }
            i++;
            ij++;
        }
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "%ld bodies\n", i);
        l += i;
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: masked pixels = %ld\n",
                                  routineName, iselect);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: unmasked pixels = %ld\n",
                                  routineName, cmd->nbody-iselect);

    if (l!=cmd->nbody || ij!=cmd->nbody)
        error("\n%s: nbody (%ld) not equal to read bodies (%ld, %ld)\n\n",
              routineName, cmd->nbody, i, ij);

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,"\n");
    for (j=0; j<gd->ninfiles; j++) {
        free(bodytable[j]);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "\tfreed %g %s (catalog %d with %ld bodies).\n",
                            gd->nbodyTable[j]*sizeof(body)*INMB,
                            "MByte for particle storage", j, gd->nbodyTable[j]);
        gd->bytes_tot -= gd->nbodyTable[j]*sizeof(body);
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tallocating %ld bodies (%ld)... ",
                           l, cmd->nbody-iselect);
    gd->nbodyTable[0] = cmd->nbody-iselect;
    bodytable[0] = (bodyptr) allocate(gd->nbodyTable[0] * sizeof(body));
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "done allocating.\n");
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
"\n\t%s: Allocated %g MByte for final bodytable storage (total bodies=%ld).\n",
        routineName, gd->nbodyTable[0]*sizeof(body)*INMB, gd->nbodyTable[0]);
    gd->bytes_tot += gd->nbodyTable[0]*sizeof(body);

    real kavg = 0;
    ij=0;
    for(i=0;i<cmd->nbody;i++){
        q = bodytabtmp+i;
        if (Mask(q) == 1) {
            p = bodytable[0]+ij;
            Pos(p)[0] = Pos(q)[0];
            Pos(p)[1] = Pos(q)[1];
            Pos(p)[2] = Pos(q)[2];
            Kappa(p) = Kappa(q);
            
#ifdef THREEPCFSHEAR
            Gamma1(p) = Gamma1(q);
            Gamma2(p) = Gamma2(q);
#endif
            
            Type(p) = Type(q);
            Mass(p) = Mass(q);
            Weight(p) = Weight(q);
            Id(p) = p-bodytable[0]+i;
            kavg += Kappa(p);
            Update(p) = Update(q);
            Mask(p) = Mask(q);
            ij++;
        }
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t%s: final unmasked pixels = %ld\n",
                        routineName, ij);

    if (gd->nbodyTable[0]>0) {
        kavg /= ((real)gd->nbodyTable[0]);
        real kstd;
        real sum=0.0;
        DO_BODY(p, bodytable[0], bodytable[0]+gd->nbodyTable[0]) {
            sum += rsqr(Kappa(p) - kavg);
        }
        kstd = rsqrt( sum/((real)gd->nbodyTable[0] - 1.0) );
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\t%s: average and std dev of kappa ",
                               routineName);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "(%ld particles) = %le %le\n",
                               gd->nbodyTable[0], kavg, kstd);
    } else {
        error("%s: no unmasked bodies (nbody=%ld) were given... exiting...\n",
              routineName, gd->nbodyTable[0]);
    }

    free(bodytabtmp);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                    "\n\tfreed %g MByte for tmp storage (%ld bodies).\n\n",
                    cmd->nbody * sizeof(body)*INMB, cmd->nbody);

    for (i=0; i<gd->ninfiles; i++) {
        (gd->iCatalogs[i]) = 0;
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                            "\t%s: iCatalogs final values: %d\n",
                            routineName, gd->iCatalogs[i]);
    }

    return SUCCESS;
}
#undef EPSILON

local int inputdata_ascii(struct cmdline_data* cmd, struct  global_data* gd,
                           string filename, int ifile)
{
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real mass=1;
    real weight=1;

    gd->input_comment = "Column form input file";

    instr = stropen(filename, "r");

    if (scanopt(cmd->options, "header-info")){
        InputData_check_file(filename);
        fgets(firstline,200,instr);
        verb_print(cmd->verbose, "\n\tinputdata_ascii: header of %s\n", filename);
        verb_print(cmd->verbose, "\t1st line: %s", firstline);
        fgets(firstline,200,instr);
        verb_print(cmd->verbose, "\t2nd line: %s\n", firstline);
        rewind(instr);
        if (scanopt(cmd->options, "stop")) {
            fclose(instr);
            exit(1);
        }
    }

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

    verb_print(cmd->verbose,
               "\tinputdata_ascii: nbody and ndim: %d %d...\n",
               cmd->nbody, ndim);
    verb_print(cmd->verbose,
               "\tinputdata_ascii: lbox dimensions: ");
    int k;
    DO_COORD(k)
    verb_print(cmd->verbose,
               "%g ", gd->Box[k]);
    verb_print(cmd->verbose,"\n\n");

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_vector(instr, Pos(p));
        in_real(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        if (scanopt(cmd->options, "kappa-constant-one"))
            Kappa(p) = 1.0;

#ifdef THREEPCFSHEAR
        //B 3pcf shear
        Gamma1(p) = 1.0;
        Gamma2(p) = 1.0;
        //E
#endif

    }

    fclose(instr);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Mask(p) = TRUE;                             // initialize body's Mask
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

local int inputdata_ascii_all(struct cmdline_data* cmd, struct  global_data* gd,
                           string filename, int ifile)
{
    string routineName = "inputdata_ascii_all";
    stream instr;
    int ndim;
    bodyptr p;
    char gato[1], firstline[20];
    real mass=1;
    real weight=1;

    gd->input_comment = "Column form input file all";

    instr = stropen(filename, "r");

    if (scanopt(cmd->options, "header-info")){
        InputData_check_file(filename);
        fgets(firstline,200,instr);
        verb_print(cmd->verbose, "\n\t%s: header of %s\n", routineName,filename);
        verb_print(cmd->verbose, "\t1st line: %s", firstline);
        fgets(firstline,200,instr);
        verb_print(cmd->verbose, "\t2nd line: %s\n", firstline);
        rewind(instr);
        if (scanopt(cmd->options, "stop")) {
            fclose(instr);
            exit(1);
        }
    }

    fgets(firstline,200,instr);
    fscanf(instr,"%1s",gato);
    in_int_long(instr, &cmd->nbody);
    if (cmd->nbody < 1)
        error("%s: nbody = %d is absurd\n", routineName, cmd->nbody);
    in_int(instr, &ndim);
    if (ndim != NDIM)
        error("%s: ndim = %d; expected %d\n", routineName, ndim, NDIM);

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
#else // ! NDIM
    real Lx, Ly;
    in_real(instr, &Lx);
    in_real(instr, &Ly);
    gd->Box[0] = Lx;
    gd->Box[1] = Ly;
#endif

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t%s: nbody and ndim: %d %d...\n",
                        routineName, cmd->nbody, ndim);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "\t%s: lbox dimensions: ", routineName);
    int k;
    DO_COORD(k)
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                        "%g ", gd->Box[k]);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog, "\n");

    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    INTEGER iselect = 0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_vector(instr, Pos(p));
        in_real(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        if (scanopt(cmd->options, "kappa-constant-one"))
            Kappa(p) = 1.0;
        in_real(instr, &Weight(p));
        in_short(instr, &Mask(p));
        if (Mask(p) == 0) {
            iselect++;
        }
#ifdef THREEPCFSHEAR
        Gamma1(p) = 1.0;
        Gamma2(p) = 1.0;
#endif
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: masked pixels = %ld\n",
                           routineName, iselect);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: unmasked pixels = %ld\n",
                           routineName, gd->nbodyTable[ifile]-iselect);

    fclose(instr);

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
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
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\t%s: average and std dev of kappa ", routineName);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
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
                    if (dist2 == 0.0)
                        flag=1;
                }
        if (flag)
            error("%s: at least two bodies have same position\n", routineName);
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
    real mass=1;
    real weight=1;

    gd->input_comment = "Binary input file";

    instr = stropen(filename, "r");

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose, "\n\tinputdata_bin: header of %s\n", 
                   filename);
        in_int_bin_long(instr, &cmd->nbody);
        verb_print(cmd->verbose, "\t1st line: %d\n", cmd->nbody);
        in_int_bin(instr, &ndim);
        verb_print(cmd->verbose, "\t2nd line: %d\n", ndim);
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
        verb_print(cmd->verbose, "\tinputdata_bin: Box: %g %g %g\n\n",
                   gd->Box[0], gd->Box[1], gd->Box[2]);
#else
        verb_print(cmd->verbose, "\tinputdata_bin: Box: %g %g %g\n\n",
                   gd->Box[0], gd->Box[1]);
#endif
        rewind(instr);
        if (scanopt(cmd->options, "stop")) {
            fclose(instr);
            exit(1);
        }
    }

    in_int_bin_long(instr, &cmd->nbody);
    verb_print(cmd->verbose, "\tInput: nbody %d\n", cmd->nbody);
    if (cmd->nbody < 1)
        error("inputdata: nbody = %d is absurd\n", cmd->nbody);
    in_int_bin(instr, &ndim);
    if (ndim != NDIM)
        error("inputdata: ndim = %d; expected %d\n", ndim, NDIM);
    verb_print(cmd->verbose, "\tInput: nbody and ndim: %d %d...\n", 
               cmd->nbody, ndim);

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
        if (scanopt(cmd->options, "kappa-constant-one"))
            Kappa(p) = 1.0;
    }
    fclose(instr);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Mask(p) = TRUE;                             // initialize body's Mask
        Id(p) = p-bodytable[ifile]+1;
    }

    return SUCCESS;
}

local int inputdata_bin_all(struct cmdline_data* cmd, struct  global_data* gd,
                         string filename, int ifile)
{
    string routineName = "inputdata_bin_all";
    stream instr;
    int ndim;
    bodyptr p;
    real mass=1;
    real weight=1;

    gd->input_comment = "Binary-all input file";

    instr = stropen(filename, "r");

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose, "\n\t%s: header of %s\n",
                   routineName, filename);
        in_int_bin_long(instr, &cmd->nbody);
        verb_print(cmd->verbose, "\t1st line: %d\n", cmd->nbody);
        in_int_bin(instr, &ndim);
        verb_print(cmd->verbose, "\t2nd line: %d\n", ndim);
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
        verb_print(cmd->verbose, "\t%s: Box: %g %g %g\n\n",
                   routineName, gd->Box[0], gd->Box[1], gd->Box[2]);
#else
        verb_print(cmd->verbose, "\t%s: Box: %g %g %g\n\n",
                   routineName, gd->Box[0], gd->Box[1]);
#endif
        rewind(instr);
        if (scanopt(cmd->options, "stop")) {
            fclose(instr);
            exit(1);
        }
    }

    in_int_bin_long(instr, &cmd->nbody);
    verb_print(cmd->verbose, "\t%s: nbody %d\n", routineName, cmd->nbody);
    if (cmd->nbody < 1)
        error("%s: nbody = %d is absurd\n", routineName, cmd->nbody);
    in_int_bin(instr, &ndim);
    if (ndim != NDIM)
        error("%s: ndim = %d; expected %d\n", routineName, ndim, NDIM);
    verb_print(cmd->verbose, "\t%s: nbody and ndim: %d %d...\n",
               routineName, cmd->nbody, ndim);

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
    verb_print(cmd->verbose, "\t%s: Box: %g %g %g\n",
               routineName, gd->Box[0], gd->Box[1], gd->Box[2]);
#else
    verb_print(cmd->verbose, "\t%s: Box: %g %g %g\n",
               routineName, gd->Box[0], gd->Box[1]);
#endif

    gd->nbodyTable[ifile] = cmd->nbody;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bytes_tot += cmd->nbody*sizeof(body);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        in_vector_bin(instr, Pos(p));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        in_real_bin(instr, &Kappa(p));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        if (scanopt(cmd->options, "kappa-constant-one"))
            Kappa(p) = 1.0;
    }
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        in_real_bin(instr, &Weight(p));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody)
        in_short_bin(instr, &Mask(p));
    fclose(instr);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Type(p) = BODY;
        Mass(p) = mass;
        Id(p) = p-bodytable[ifile]+1;
    }

    return SUCCESS;
}


//B BEGIN:: Reading Takahasi simulations
//From Takahashi web page. Adapted to our needs

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

local int inputdata_takahashi(struct cmdline_data* cmd, struct  global_data* gd,
                             string filename, int ifile)
{
    string routinename = "inputdata_takahashi";
    FILE *fp;
    long i,j,npix,dummy;
    long jj[6]={536870908,1073741818,1610612728,2147483638,2684354547,3221225457};
    int negi,nside;
//    double theta,phi;
//    char file[200];

    gd->input_comment = "Takahasi input file";

//E Begin reading Takahashi file
    fp = stropen(filename, "rb");

    if (scanopt(cmd->options, "header-info")){
        verb_print(cmd->verbose, "\n\t%s: header of %s... ",
                   routinename, filename);
        verb_print(cmd->verbose, "not available yet... sorry!\n\n");
        if (scanopt(cmd->options, "stop")) {
            fclose(fp);
            exit(1);
        }
    }

    
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

    verb_print(cmd->verbose,
               "\nAllocated %g MByte for pixel storage.\n",
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

    fclose(fp);                                     // Close Takahashi file.

    verb_print(cmd->verbose,
               "\n\t%s: total read points = %ld\n",routinename, npix);

//E End reading Takahashi file

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
    string routinename = "Takahasi_region_selection";
    double theta,phi;
    long i;

//B Computing min and max of theta and phi:
    real theta_min, theta_max;
    real phi_min, phi_max;
    pix2ang(0,nside,&theta_min,&phi_min);
    theta_max = theta_min;
    phi_max = phi_min;

    for(i=1;i<npix;i++){                            // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
        theta_min = MIN(theta_min,theta);
        theta_max = MAX(theta_max,theta);
        phi_min = MIN(phi_min,phi);
        phi_max = MAX(phi_max,phi);
    }
    verb_print(cmd->verbose,
               "\t%s: min and max of theta = %f %f\n",
               routinename, theta_min, theta_max);
    verb_print(cmd->verbose,
               "\t%s: min and max of phi = %f %f\n",
               routinename, phi_min, phi_max);
//E

//B Selection of a region: the (center) lower edge is random... or not
    real rphi, rtheta;
//B Change selection to random or fix, given by thetaL, phiL, thetaR, phiR
    if (scanopt(cmd->options, "random-point")) {
        rphi    = 2.0 * PI * xrandom(0.0, 1.0);
        rtheta    = racos(1.0 - 2.0 * xrandom(0.0, 1.0));
        verb_print(cmd->verbose,
                   "\t%s: random theta and phi = %f %f\n",
                   routinename, rtheta, rphi);
    } else {
// The radius of the region is:
//        rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL)))
// ... and the center:
        rphi = 0.5*(cmd->phiR + cmd->phiL);
        rtheta = 0.5*(cmd->thetaR + cmd->thetaL);
        verb_print(cmd->verbose,
                   "\t%s: fix theta and phi = %f %f\n",
                   routinename, rtheta, rphi);
    }
//E

    if (cmd->lengthBox>PI) {
        fprintf(stdout,"\n%s: Warning! %s\n%s\n\n",
                "lengthBox is greater than one of the angular ranges...",
                "using length = PI",
                routinename);
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
                   "\t%s: rotated theta and phi = %f %f\n",
                   routinename, rtheta_rot, rphi_rot);
    } else {
        dtheta_rot = 0.0,
        dphi_rot = 0.0;
        rtheta_rot = rtheta;
        rphi_rot = rphi;
        verb_print(cmd->verbose,
                "\t%s: theta and phi (no-rotation) = %f %f\n",
                   routinename, rtheta_rot, rphi_rot);
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
        "\t%s: theta and phi of the center of the selected region = %lf %lf\n",
                   routinename, 0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
        verb_print(cmd->verbose,
                "\t%s: radius of the selected region = %lf\n",
                   routinename,
                   rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL)))
                   );
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
"\tinputdata_takahashi: theta and phi of the center of the selected region = %f %f\n",
                       0.5*(thetaR + thetaL), 0.5*(phiR + phiL));
            verb_print(cmd->verbose,
                    "\t%s: radius of the selected region = %lf\n",
                    routinename,
                    rsqrt(rsqr(0.5*(thetaR - thetaL)) + rsqr(0.5*(phiR - phiL)))
                       );
    }

    verb_print(cmd->verbose,
               "\t%s: left and right theta = %f %f\n",
               routinename, thetaL, thetaR);
    verb_print(cmd->verbose,
               "\t%s: left and right phi = %f %f\n",
               routinename, phiL, phiR);
    verb_print(cmd->verbose,
               "\t%s: theta and phi d_rotation = %f %f\n",
               routinename, dtheta_rot, dphi_rot);
//E

#if NDIM == 3
    real xmin, ymin, zmin;
    real xmax, ymax, zmax;

    if (scanopt(cmd->options, "patch")) {
        Takahasi_region_selection_3d(cmd, gd,
                                     nside, npix, conv, shear1, shear2, rotat,
                dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, ifile);
    } else {
        Takahasi_region_selection_3d_all(cmd, gd,
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
    string routinename = "Takahasi_region_selection_3d_all";
    long i;
    bodyptr p;
    real mass = 1.0;
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

    for(i=0;i<npix;i++){                            // Healpix ring scheme
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

        coordinate_transformation(cmd, gd, theta, phi, Pos(p));

        if (!scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = conv[i];
        else {
            Kappa(p) = 2.0;
            if (scanopt(cmd->options, "kappa-constant-one"))
                Kappa(p) = 1.0;
        }

#ifdef THREEPCFSHEAR
        //B 3pcf shear
        Gamma1(p) = shear1[i];
        Gamma2(p) = shear2[i];
        //E
#endif

        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+iselect;

        *xmin = Pos(p)[0];
        *ymin = Pos(p)[1];
        *zmin = Pos(p)[2];
        *xmax = Pos(p)[0];
        *ymax = Pos(p)[1];
        *zmax = Pos(p)[2];

        Update(p) = TRUE;
        Mask(p) = TRUE;                             // initialize body's Mask

        //B correction 2025-05-03 :: look for edge-effects
        // activate a flag for this catalog, that it is using patch-with-all
        //  and use it in EvalHist routine...
#if defined(NMultipoles) && defined(NONORMHIST)
        if (scanopt(cmd->options, "patch-with-all")) {
            UpdatePivot(p) = TRUE;
            if (thetaL < theta - dtheta_rot && theta - dtheta_rot < thetaR) {
                if (phiL < phi - dphi_rot && phi - dphi_rot < phiR) {
                    UpdatePivot(p) = TRUE;
                    gd->pivotCount += 1;
                } else {
                    UpdatePivot(p) = FALSE;
                }
            } else {
                UpdatePivot(p) = FALSE;
            }
        }
#endif
        //E
    } // ! end for

    //B correction 2025-05-03 :: look for edge-effects
#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
            "\nsearchcalc_tc_kkk_omp: total number of pixels to be pivots: %ld\n",
            gd->pivotCount);
    }
#endif
    //E
    
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
               "\n\t%s: min and max of x = %f %f\n",
               routinename, *xmin, *xmax);
    verb_print(cmd->verbose,
               "\t%s: min and max of y = %f %f\n",
               routinename, *ymin, *ymax);
    verb_print(cmd->verbose,
               "\t%s: min and max of z = %f %f\n",
               routinename, *zmin, *zmax);

    verb_print(cmd->verbose,
        "\n\t%s: selected all read points and nbody: %ld %ld\n",
               routinename, iselect, cmd->nbody);

    verb_print(cmd->verbose, 
               "\t%s: average of kappa (%ld particles) = %le\n",
               routinename, cmd->nbody, kavg/((real)cmd->nbody) );

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
    real mass = 1.0;
    real weight = 1.0;

    real theta, phi;
    real theta_rot, phi_rot;
    INTEGER iselect = 0;
    
    bodyptr bodytabtmp;
    cmd->nbody = npix;
    bodytabtmp = (bodyptr) allocate(cmd->nbody * sizeof(body));
    verb_print(cmd->verbose,
               "\nAllocated %g MByte for particle (%ld) storage.\n",
               cmd->nbody*sizeof(body)/(1024.0*1024.0),cmd->nbody);

    *xmin=0., *ymin=0., *zmin=0.;
    *xmax=0., *ymax=0., *zmax=0.;

    for(i=0;i<npix;i++){                            // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
//        printf("%ld %f %f %f %f \n", 
//                i, conv[i], shear1[i], shear2[i], rotat[i]);
//        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;
        Mask(p) = TRUE;                             // initialize body's Mask

        if (scanopt(cmd->options, "patch")) {
            //B
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

                    coordinate_transformation(cmd, gd, theta, phi, Pos(p));

                    if (!scanopt(cmd->options, "kappa-constant"))
                        Kappa(p) = conv[i];
                    else {
                        Kappa(p) = 2.0;
                        if (scanopt(cmd->options, "kappa-constant-one"))
                            Kappa(p) = 1.0;
                    }

#ifdef THREEPCFSHEAR
                    //B 3pcf shear
                    Gamma1(p) = shear1[i];
                    Gamma2(p) = shear2[i];
                    //E
#endif

                    Type(p) = BODY;
                    Mass(p) = mass;
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
            //E
        } else { // ! all
            //B
            iselect++;
            if (scanopt(cmd->options, "rotation")) {
                theta_rot = theta - dtheta_rot;
                phi_rot = phi - dphi_rot;
            } else {
                theta_rot = theta;
                phi_rot = phi;
            }

            coordinate_transformation(cmd, gd, theta, phi, Pos(p));

            if (!scanopt(cmd->options, "kappa-constant"))
                Kappa(p) = conv[i];
            else {
                Kappa(p) = 2.0;
                if (scanopt(cmd->options, "kappa-constant-one"))
                    Kappa(p) = 1.0;
            }

#ifdef THREEPCFSHEAR
            //B 3pcf shear
            Gamma1(p) = shear1[i];
            Gamma2(p) = shear2[i];
            //E
#endif

            Type(p) = BODY;
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;

            *xmin = Pos(p)[0];
            *ymin = Pos(p)[1];
            *zmin = Pos(p)[2];
            *xmax = Pos(p)[0];
            *ymax = Pos(p)[1];
            *zmax = Pos(p)[2];

            Update(p) = TRUE;
            //E
        } // ! all
    } // ! end for

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
        cmd->nbody = iselect;

    gd->nbodyTable[ifile] = cmd->nbody;
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

#ifdef THREEPCFSHEAR
            //B 3pcf shear
            Gamma1(p) = Gamma1(q);
            Gamma2(p) = Gamma2(q);
            //E
#endif

            Type(p) = Type(q);
            Mass(p) = mass;
            Weight(p) = weight;
            Mask(p) = Mask(q);
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
    verb_print(cmd->verbose, "\n\tinputdata_takahashi: min and max of x = %f %f\n",*xmin, *xmax);
    verb_print(cmd->verbose, "\tinputdata_takahashi: min and max of y = %f %f\n",*ymin, *ymax);
    verb_print(cmd->verbose, "\tinputdata_takahashi: min and max of z = %f %f\n",*zmin, *zmax);
    free(bodytabtmp);

    if (scanopt(cmd->options, "patch"))
        verb_print(cmd->verbose,
                   "\n\tinputdata_takahashi: selected read points = %ld\n",iselect);
    else
        verb_print(cmd->verbose,
                   "\n\tinputdata_takahashi: selected read points and nbody: %ld %ld\n",
                   iselect, cmd->nbody);

    verb_print(cmd->verbose, 
               "inputdata_takahashi: average of kappa (%ld particles) = %le\n",
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
    real mass = 1;
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

    for(i=0;i<npix;i++){                            // Healpix ring scheme
        pix2ang(i,nside,&theta,&phi);
    //        printf("%ld %f %f %f %f \n", i, conv[i], shear1[i], shear2[i], rotat[i]);
    //        printf("%ld %f %f %f \n", i, theta, phi, conv[i]);
        p = bodytabtmp+i;
        Update(p) = FALSE;
        Mask(p) = TRUE;                             // initialize body's Mask

        if (scanopt(cmd->options, "patch")) {
            //B
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
                    Mass(p) = mass;
                    Weight(p) = weight;
                    Id(p) = p-bodytabtmp+iselect;

                    *xmin = Pos(p)[0];
                    *ymin = Pos(p)[1];
                    *xmax = Pos(p)[0];
                    *ymax = Pos(p)[1];

                    Update(p) = TRUE;
                }
            }
            //E
        } else { // ! all
            //B
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
            Mass(p) = mass;
            Weight(p) = weight;
            Id(p) = p-bodytabtmp+iselect;

            *xmin = Pos(p)[0];
            *ymin = Pos(p)[1];
            *xmax = Pos(p)[0];
            *ymax = Pos(p)[1];

            Update(p) = TRUE;
            //E
        } // ! all
    } // ! end loop i

    bodyptr q;
    if (scanopt(cmd->options, "patch"))
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
            Mass(p) = mass;
            Weight(p) = weight;
            Mask(p) = Mask(q);
            Id(p) = p-bodytable[ifile]+i;
            *xmin = MIN(*xmin,Pos(p)[0]);
            *ymin = MIN(*ymin,Pos(p)[1]);
            *xmax = MAX(*xmax,Pos(p)[0]);
            *ymax = MAX(*ymax,Pos(p)[1]);
            ij++;
        }
    }
    verb_print(cmd->verbose, 
               "\n\tinputdata_takahashi: min and max of x = %f %f\n",
               *xmin, *xmax);
    verb_print(cmd->verbose, 
               "\tinputdata_takahashi: min and max of y = %f %f\n",
               *ymin, *ymax);

    free(bodytabtmp);
    
    verb_print(cmd->verbose, 
               "\n\tinputdata_takahashi: selected read points = %ld\n",iselect);

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

//E End:: Reading Takahashi simulations


int StartOutput(struct cmdline_data *cmd, struct  global_data* gd)
{
    //B clear some char arrays
    gd->logfilePath[0] = '\0';
    gd->fpfnameOutputFileName[0] = '\0';
    gd->fnameData_kd[0] = '\0';
    gd->fnameOut_kd[0] = '\0';
    //E

    outfilefmt_string_to_int(cmd->outfilefmt, &outfilefmt_int);

    if (cmd->verbose>=VERBOSEMININFO)
        if (! strnull(cmd->options))
            verb_print(cmd->verbose, "\n\toptions: %s\n", cmd->options);

    return SUCCESS;
}

/*
 OutputData routine:

 To be called by MainLoop in cballs.c:
    OutputData(cmd, gd, bodytable, gd->nbodyTable, ifile);

 This routine is in charge of saving a catalog of data

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
    * `btable`: Input: a body pointer structure array
    * `nbody`: Input: number of points in table array
    * `ifile`: Input: catalog file tag
 Return (the error status):
    int SUCCESS or FAILURE
 */
int OutputData(struct cmdline_data* cmd, struct  global_data* gd,
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
        case OUTCOLUMNSALL:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii_all(cmd, gd, btable, nbody); break;
        case OUTCOLUMNSBIN:
            verb_print(cmd->verbose, "\n\tbinary format output\n");
            outputdata_bin(cmd, gd, btable, nbody); break;
        case OUTCOLUMNSBINALL:
            verb_print(cmd->verbose, "\n\tbinary-all format output\n");
            outputdata_bin_all(cmd, gd, btable, nbody); break;
        case OUTNULL:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            outputdata_ascii(cmd, gd, btable, nbody); break;

//B socket:
#ifdef ADDONS
#include "cballsio_include_03.h"
#endif
//E

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
#if NDIM == 3
    fprintf(outstr,"# nbody NDIM Lx Ly Lz\n# %ld %d ",nbody,NDIM);
    fprintf(outstr,"%lf %lf %lf\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    fprintf(outstr,"# nbody NDIM Lx Ly\n# %ld %d ",nbody,NDIM);
    fprintf(outstr,"%lf %lf\n",gd->Box[0],gd->Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        out_vector_mar(outstr, Pos(p));
        out_real_mar(outstr, Kappa(p));
//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
        out_bool_mar(outstr, HIT(p));
#endif
//E
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    verb_print(cmd->verbose, "\tdata output to file %s\n", namebuf);

    return SUCCESS;
}

local int outputdata_ascii_all(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    string routineName = "outputdata_ascii_all";
    char namebuf[256];
    stream outstr;
    bodyptr p;

    sprintf(namebuf, gd->fpfnameOutputFileName);
    outstr = stropen(namebuf, "w!");
#if NDIM == 3
    fprintf(outstr,"# nbody NDIM Lx Ly Lz\n# %ld %d ",nbody,NDIM);
    fprintf(outstr,"%lf %lf %lf\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    fprintf(outstr,"# nbody NDIM Lx Ly\n# %ld %d ",nbody,NDIM);
    fprintf(outstr,"%lf %lf\n",gd->Box[0],gd->Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        out_vector_mar(outstr, Pos(p));
        out_real_mar(outstr, Kappa(p));
        out_real_mar(outstr, Weight(p));
        out_short_mar(outstr, Mask(p));
//B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
        out_bool_mar(outstr, HIT(p));
#endif
//E
        fprintf(outstr,"\n");
    }
    fclose(outstr);
    verb_print(cmd->verbose, "\t%s: data output to file %s\n",
               routineName, namebuf);

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

local int outputdata_bin_all(struct cmdline_data* cmd, struct  global_data* gd,
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
    DO_BODY(p, bodytab, bodytab+nbody)
        out_real_bin(outstr, Weight(p));
    DO_BODY(p, bodytab, bodytab+nbody)
        out_short_bin(outstr, Mask(p));
    //B BALLS :: DIAGNOSTICS (DEBUG)
#ifdef DEBUG
    DO_BODY(p, bodytab, bodytab+nbody)
        out_bool_mar(outstr, HIT(p));
#endif
    //E
    fclose(outstr);
    verb_print(cmd->verbose, "\tdata output to file %s\n", namebuf);

    return SUCCESS;
}


global int infilefmt_string_to_int(string infmt_str,int *infmt_int)
{
    *infmt_int=-1;
    if (strcmp(infmt_str,"columns-ascii") == 0)     *infmt_int = INCOLUMNS;
    if (strcmp(infmt_str,"columns-ascii-all") == 0) *infmt_int = INCOLUMNSALL;
    if (strnull(infmt_str))                         *infmt_int = INNULL;
    if (strcmp(infmt_str,"binary") == 0)            *infmt_int = INCOLUMNSBIN;
    if (strcmp(infmt_str,"binary-all") == 0)        *infmt_int = INCOLUMNSBINALL;
    if (strcmp(infmt_str,"takahashi") == 0)          *infmt_int = INTAKAHASHI;

//B socket:
#ifdef ADDONS
#include "cballsio_include_08.h"
#endif
//E

    return SUCCESS;
}

local int outfilefmt_string_to_int(string outfmt_str,int *outfmt_int)
{
    *outfmt_int=-1;
    if (strcmp(outfmt_str,"columns-ascii") == 0)     *outfmt_int = OUTCOLUMNS;
    if (strcmp(outfmt_str,"columns-ascii-all") == 0) *outfmt_int = OUTCOLUMNSALL;
    if (strnull(outfmt_str))                         *outfmt_int = OUTNULL;
    if (strcmp(outfmt_str,"binary") == 0)           *outfmt_int = OUTCOLUMNSBIN;
    if (strcmp(outfmt_str,"binary-all") == 0)      *outfmt_int = OUTCOLUMNSBINALL;

//B socket:
#ifdef ADDONS
#include "cballsio_include_09.h"
#endif
//E

    return SUCCESS;
}

// I/O directories:
global void setFilesDirs_log(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "setFilesDirs_log";
    char buf[BUFFERSIZE];

    debug_tracking_s("001", routineName);
    if (cmd->verbose_log>0) {           // gd->logfilePath is defined
        debug_tracking_s("002", cmd->rootDir);
        sprintf(gd->tmpDir,"%s/%s",cmd->rootDir,"tmp");
        double cpustart = CPUTIME;
        debug_tracking_s("003", gd->tmpDir);
        sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",
                gd->tmpDir,gd->tmpDir);
        system(buf);
        debug_tracking("004");
        gd->cputotalinout += CPUTIME - cpustart;
        sprintf(gd->logfilePath,"%s/cballs%s.log",
                gd->tmpDir,cmd->suffixOutFiles);
    }
    debug_tracking("005... final");
}

global void setFilesDirs(struct cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "setFilesDirs";
    char buf[BUFFERSIZE];

    char outputDir[MAXLENGTHOFFILES];

    double cpustart = CPUTIME;

    int ndefault = 0;
    int *ipos;
    char *dp1, *dp2;
    int lenDir = strlen(cmd->rootDir);
    int i;

    debug_tracking_s("001", routineName);
    debug_tracking_s("002: init", cmd->rootDir);

    if (gd->rootDirFlag==TRUE) {
        
        int nslashs = MAXNSLASHS;
        ipos = (int*) malloc((nslashs)*sizeof(int));
        dp1 = (char*) malloc((lenDir)*sizeof(char));

        for (i=0; i< lenDir; i++) {
            if(cmd->rootDir[i] == '/') {
                ipos[ndefault] = i+1;
                ndefault++;
            }
        }
        if (ndefault>nslashs)
            error("%s: more '/' than %d in 'rootDir=%s'. Use only %d or none\n",
                  routineName, nslashs, cmd->rootDir, nslashs);
        
        if (ndefault == 0) {
            sprintf(outputDir,cmd->rootDir);
            sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",
                    outputDir,outputDir);
            if (cmd->verbose >= 3)
                verb_print_q(3, cmd->verbose_log,"\nsystem: %s\n",buf);
            system(buf);
        } else {
            for (i=0; i<ndefault; i++) {
                debug_tracking("003");
                strncpy(dp1, cmd->rootDir, ipos[i]-1);
                sprintf(buf,"if [ ! -d %s ]; then mkdir -p %s; fi",dp1,dp1);
                verb_print_q(3,cmd->verbose_log,"\nsystem: %d: %s\n",i,buf);
                system(buf);
                debug_tracking("004");
            }
            strncpy(dp1, cmd->rootDir, lenDir);
            sprintf(buf,"if [ ! -d %s ]; then mkdir -p %s; fi",dp1,dp1);
            verb_print_q(3,cmd->verbose_log,"\nsystem: %d: %s\n",i,buf);
            system(buf);
        }
        gd->cputotalinout += CPUTIME - cpustart;
        
        debug_tracking_s("005", cmd->rootDir);

        sprintf(gd->fpfnameOutputFileName,"%s/%s%s%s",
                cmd->rootDir,cmd->outfile,cmd->suffixOutFiles,EXTFILES);
        sprintf(gd->fpfnamehistNNFileName,"%s/%s%s%s",
                cmd->rootDir,cmd->histNNFileName,cmd->suffixOutFiles,EXTFILES);
        sprintf(gd->fpfnamehistCFFileName,"%s/%s%s%s",
                cmd->rootDir,"histCF",cmd->suffixOutFiles,EXTFILES);
        sprintf(gd->fpfnamehistrBinsFileName,"%s/%s%s%s",
                cmd->rootDir,"rbins",cmd->suffixOutFiles,EXTFILES);
        sprintf(gd->fpfnamehistXi2pcfFileName,"%s/%s",
                cmd->rootDir,cmd->histXi2pcfFileName);
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
        
        free(ipos);
        

        //B socket:
#ifdef ADDONS
#include "cballsio_include_09b.h"
#endif
        //E
        free(dp1);
    } // ! rootDirFlag
    debug_tracking_s("006: final", cmd->rootDir);
}


/*
 EndRun routine:

 To be called in main:
    EndRun(&cmd, &gd);

 This routine is in charge of closing log file, printing a summary
    of the run and freeing the allocated memory

 Arguments:
    * `cmd`: Input: structure cmdline_data pointer
    * `gd`: Input: structure global_data pointer
 Return (the error status):
    int SUCCESS or FAILURE
 */
int EndRun(struct cmdline_data* cmd, struct  global_data* gd)
{
    stream outstr;

    if (cmd->verbose_log>0)
        fclose(gd->outlog);

    if (cmd->verbose >= VERBOSENORMALINFO) {
        //B only catalog 0 is considered... modify to include others
        printf("\nrSize \t\t= %lf\n", gd->rSizeTable[0]);
        printf("nbbcalc \t= %ld\n", gd->nbbcalc);
        printf("nbccalc \t= %ld\n", gd->nbccalc);
        printf("ncccalc \t= %ld\n", gd->ncccalc);
        printf("tdepth \t\t= %d\n", gd->tdepthTable[0]);
        // Consider other tree cells, like kdtree
        printf("ncell\t\t= %ld\n", gd->ncellTable[0]);
        //E
        verb_print_q(3,cmd->verbose,"sameposcount \t= %ld\n",gd->sameposcount);
#ifdef OPENMPCODE
        printf("cpusearch \t= %lf %s\n",
               gd->cpusearch, PRNUNITOFTIMEUSED);
#else
        printf("cpusearch \t= %lf %s\n",
               gd->cpusearch, PRNUNITOFTIMEUSED);
#endif
        printf("cputotalinout \t= %lf %s\n",
               gd->cputotalinout, PRNUNITOFTIMEUSED);
    }

    if (cmd->verbose > VERBOSENOINFO) {
        real cpuTotal = CPUTIME - gd->cpuinit;
        printf("\nFinal CPU time : %lf %s\n",
               cpuTotal, PRNUNITOFTIMEUSED);
        if (scanopt(cmd->options, "measure-cputime")) {
                outstr = stropen(gd->fpfnameCPUFileName, "a");
            fprintf(outstr,"%g %g %g\n",
                    (double) cmd->nbody, cpuTotal, gd->cpusearch);
                fclose(outstr);
        }
        printf("Final real time: %ld",
               (rcpu_time()-gd->cpurealinit));
        printf(" %s\n\n", PRNUNITOFTIMEUSED);       // Only work this way
    }

    EndRun_FreeMemory(cmd, gd);

    return SUCCESS;
}

//
// We must check the order of memory allocation and deallocation!!!
//
global int EndRun_FreeMemory(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
//    printf("\n Flags: %d, %d, %d, %d, %d, %d\n", gd->tree_allocated,
//           gd->gd_allocated_2,
//           gd->bodytable_allocated,
//           gd->histograms_allocated,
//           gd->gd_allocated,
//           gd->cmd_allocated);
    
    if (gd->tree_allocated == TRUE)
        EndRun_FreeMemory_tree(cmd, gd);

    if (gd->gd_allocated_2 == TRUE)
        EndRun_FreeMemory_gd_2(cmd, gd);

    if (gd->bodytable_allocated == TRUE)
        EndRun_FreeMemory_bodytable(cmd, gd);

    if (gd->histograms_allocated == TRUE)
        EndRun_FreeMemory_histograms(cmd, gd);

    if (gd->gd_allocated == TRUE)
        EndRun_FreeMemory_gd(cmd, gd);
    if (gd->cmd_allocated == TRUE)
        EndRun_FreeMemory_cmd(cmd, gd);

    return SUCCESS;
}

global int EndRun_FreeMemory_tree(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    int ifile;

#ifdef BALLS4SCANLEV
    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        free(nodetablescanlevB4[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++)
            free(nodetablescanlevB4[ifile]);
    }
#endif

    if (!scanopt(cmd->searchMethod, "kdtree-omp")
        && !scanopt(cmd->searchMethod, "kdtree-box-omp") ) {
        freeTree(cmd, gd);
    }

    gd->tree_allocated = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_bodytable(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    int ifile;

    if (scanopt(cmd->options, "read-mask")) {
        ifile=0;
        free(bodytable[ifile]);
    } else {
        for (ifile=0; ifile<gd->ninfiles; ifile++)
            free(bodytable[ifile]);
    }

    gd->bodytable_allocated = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_histograms(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    free_dvector(gd->histN2pcf,1,cmd->sizeHistN);
    // 2pcf
#ifdef SMOOTHPIVOT
    free_dvector(gd->histNNSubN2pcftotal,1,cmd->sizeHistN);
#endif
    free_dvector(gd->histNNSubN2pcf,1,cmd->sizeHistN);
    //E


//B socket:
#ifdef ADDONS
#include "cballsio_include_10.h"                    // this is empty and can
                                                    //be remove these 3 lines
#endif
//E

#ifdef TPCF
        free_dmatrix3D(gd->histZetaGmIm,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaGmRe,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
        free_dmatrix3D(gd->histZetaMcossin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMsincos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMsin,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaMcos,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix3D(gd->histZetaM,
                       1,cmd->mChebyshev+1,1,cmd->sizeHistN,1,
                       cmd->sizeHistN);
        free_dmatrix(gd->histXisin,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
        free_dmatrix(gd->histXicos,1,cmd->mChebyshev+1,1,cmd->sizeHistN);
#endif

    free_dvector(gd->histXi2pcf,1,cmd->sizeHistN);
     
    free_dvector(gd->histNNN,1,cmd->sizeHistN);
    // 2pcf
#ifdef SMOOTHPIVOT
    free_dvector(gd->histNNSubXi2pcftotal,1,cmd->sizeHistN);
#endif
    free_dvector(gd->histNNSubXi2pcf,1,cmd->sizeHistN);
    //
    free_dvector(gd->histNNSub,1,cmd->sizeHistN);
    free_dvector(gd->histCF,1,cmd->sizeHistN);
    free_dvector(gd->histNN,1,cmd->sizeHistN);

    //B Histogram arrays PXD versions
#ifdef PXD
    free_dvector(gd->histZetaMFlatten,1,cmd->sizeHistN*cmd->sizeHistN);
    free_dvector(gd->rBins,1,cmd->sizeHistN);
    free_dmatrix(gd->matPXD,1,cmd->sizeHistN,1,cmd->sizeHistN);
    free_dvector(gd->vecPXD,1,cmd->sizeHistN);
#endif
    //E Histogram arrays PXD versions

    gd->histograms_allocated = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_gd(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "EndRun_FreeMemory_gd";
    int ifile;

    debug_tracking_s("001", routineName);
    //B Set gsl uniform random :: If not needed globally
    //      this line have to go to testdata
    #ifdef USEGSL
        gsl_rng_free (gd->r);           // allocated by random_init
                                        //  routine in startrun.c
    #endif
    //E

    gd->gd_allocated = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_gd_2(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    if (cmd->useLogHist) {
        if (cmd->rminHist==0) {
            free_dvector(gd->deltaRV,1,cmd->sizeHistN);
        } else {
            free_dvector(gd->deltaRV,1,cmd->sizeHistN);
            free_dvector(gd->ddeltaRV,1,cmd->sizeHistN);
        }
    }

    gd->gd_allocated_2 = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_cmd(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "EndRun_FreeMemory_cmd";
    //B 2026-01: added to follow up: freeing several cmd strings variables:
    // first one must be version, check!!
    debug_tracking_s("001", routineName);
#ifdef SAVERESTORE
    if (cmd->restorefile!=NULL)
    free(cmd->restorefile);
    if (cmd->statefile!=NULL)
    free(cmd->statefile);           // it uses strdup() function, check!!
#endif
#ifdef IOLIB
    if (gd->columnsFlag==TRUE)
        free(cmd->columns);
#endif
    debug_tracking("002");
    if (cmd->options!=NULL)
    free(cmd->options);
    debug_tracking("003");
//    free(cmd->posScript);
    debug_tracking("004");
//    free(cmd->preScript);
//    free(cmd->testmodel);
//    free(cmd->suffixOutFiles);
    debug_tracking("005");
// they all have assigned a string: "histNN", ...
//    free(cmd->histZetaFileName);
//    free(cmd->histXi2pcfFileName);
//    free(cmd->histNNFileName);
//    free(cmd->outfilefmt);
//    free(cmd->outfile);
    if (gd->rootDirFlagFree==TRUE)
        free(cmd->rootDir);
    if (gd->iCatalogsFlag==TRUE)
        free(cmd->iCatalogs);
    if (gd->infilefmtFlag==TRUE)
        free(cmd->infilefmt);
    if (gd->infileFlag==TRUE)
        free(cmd->infile);
    if (gd->rsmoothFlagFree==TRUE)
        free(cmd->rsmooth);
    if (gd->searchMethodFlag==TRUE)
        free(cmd->searchMethod);
        debug_tracking("006");
    // last one must be paramfile, check!!

    gd->cmd_allocated = FALSE;

    return SUCCESS;
}


//B socket:
#ifdef ADDONS
#include "cballsio_include_11a.h"
#endif
//E

//B socket:
#ifdef ADDONS
#include "cballsio_include_11b.h"
#endif
//E
