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
#include "input_contracts.h"
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include <stdint.h>

#ifdef CLASSLIB
#define cBALLS_FAIL(cmd, ...)                                           \
    do {                                                                \
        snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__);    \
        return FAILURE;                                                 \
    } while (0)
#else
#define cBALLS_FAIL(cmd, ...) error(__VA_ARGS__)
#endif


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
local int InputData_local(struct cmdline_data* cmd,
                          struct global_data* gd, string filename, int ifile);

typedef struct {
    struct cmdline_data *cmd;
    struct global_data *gd;
    string filename;
    int ifile;
} inputdata_context;

local int inputdata_guarded_local(void *argument)
{
    inputdata_context *context = argument;
    return InputData_local(context->cmd, context->gd,
                           context->filename, context->ifile);
}

int InputData(struct cmdline_data* cmd,
              struct global_data* gd, string filename, int ifile)
{
    inputdata_context context = {cmd, gd, filename, ifile};
    /* Catch local allocation failures before the collective, not outside it.
     * Every rank must enter this checkpoint, including successful readers. */
    int status = cballs_allocation_guard(inputdata_guarded_local, &context,
                                          cmd->error_message, _ERRORMSGSIZE_);
    if (status == FAILURE && strstr(cmd->error_message, filename) == NULL) {
        char detail[_ERRORMSGSIZE_];
        snprintf(detail, sizeof(detail), "%s", cmd->error_message);
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "catalog input '%s': %.1024s", filename, detail);
    }
#ifdef CBALLS_MPI_ENABLED
    status = cballs_mpi_consensus(cmd, status, "MPI catalog input");
#endif
    return status;
}

local int InputData_local(struct cmdline_data* cmd,
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
            class_call_cballs(inputdata_ascii(cmd, gd, filename, ifile), errmsg, errmsg);
            break;

            
        case INCOLUMNSALL:
            verb_print_normal_info(cmd->verbose,
                            cmd->verbose_log, gd->outlog,
                            "\tInput in columns (ascii) all format...\n");
            class_call_cballs(inputdata_ascii_all(cmd, gd, filename, ifile), errmsg, errmsg);
            break;
        case INNULL:
            verb_print_normal_info(cmd->verbose,
                        cmd->verbose_log, gd->outlog,
                        "\n\t(Null) Input in columns (ascii) format...\n");
            class_call_cballs(inputdata_ascii(cmd, gd, filename, ifile), errmsg, errmsg);
            break;
        case INCOLUMNSBIN:
            verb_print_normal_info(cmd->verbose,
                                   cmd->verbose_log, gd->outlog,
                                   "\n\tInput in binary format...\n");
            class_call_cballs(inputdata_bin(cmd, gd, filename, ifile), errmsg, errmsg);
            break;
        case INCOLUMNSBINALL:
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                   "\n\tInput in binary-all format...\n");
            class_call_cballs(inputdata_bin_all(cmd, gd, filename, ifile), errmsg, errmsg);
            break;
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
                cBALLS_FAIL(cmd, "\n\tgoing out...\n");
            }
            verb_print(cmd->verbose,
                       "\n\tInput in default columns (ascii) format...\n");
            class_call_cballs(inputdata_ascii(cmd, gd, filename, ifile),
                                              errmsg, errmsg);
            break;
    }
    /* A final common boundary also checks values produced by coordinate
     * conversion. Header-only and mask-only readers may not publish a body. */
    const int separate_mask = scanopt(cmd->options, "read-mask") && ifile == 1;
    if (!gd->inputHeaderFlag && !gd->stopflag && !separate_mask
        && cballs_input_validate_bodies(cmd, filename, bodytable[ifile],
                                         gd->nbodyTable[ifile]) == FAILURE)
        return FAILURE;
    if (!gd->inputHeaderFlag)
        for (int k = 0; k < NDIM; ++k)
            if (cballs_input_finite(cmd, filename, 0, "box dimension",
                                    gd->Box[k]) == FAILURE)
                return FAILURE;
    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,
            "\tdone reading.\n");

    gd->cputotalinout += CPUTIME - cpustart;

    // DEBUG WARNING!!
    //B There is a Segmentation fault: 11 if run as:
    //  cballs in=./scripts/Abraham/kappa_nres12_zs9NS256r000.txt
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
#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
            LyaForestId(p) = LyaForestId(q);
            LyaDistance(p) = LyaDistance(q);
            SETV(LyaLOS(p), LyaLOS(q));
#endif
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
        cBALLS_FAIL(cmd,
                    "\n%s: nbody (%ld) not equal to read bodies (%ld, %ld)\n\n",
                    routineName, cmd->nbody, i, ij);

    verb_print_min_info(cmd->verbose, cmd->verbose_log, gd->outlog,"\n");
    
    //B
    for (j = 0; j < gd->ninfiles && j < MAXITEMS; j++) {
        if (bodytable[j] != NULL) {
            verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                                "\tfreed %g %s (catalog %d with %ld bodies).\n",
                                gd->nbodyTable[j]*sizeof(body)*INMB,
                                "MByte for particle storage", j, gd->nbodyTable[j]);
            gd->bytes_tot -= gd->nbodyTable[j] * sizeof(body);
            free(bodytable[j]);
            bodytable[j] = NULL;
            gd->nbodyTable[j] = 0;
        }
    }

    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "\n\tallocating %ld bodies (%ld)... ",
                           l, cmd->nbody-iselect);
    gd->nbodyTable[0] = cmd->nbody-iselect;
    bodytable[0] = (bodyptr) allocate(gd->nbodyTable[0] * sizeof(body));
    gd->bodytable_allocated = TRUE;
    //E
    
    
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
#if defined(LYAFORESTOMP) || defined(LYAFORESTMPI)
            LyaForestId(p) = LyaForestId(q);
            LyaDistance(p) = LyaDistance(q);
            SETV(LyaLOS(p), LyaLOS(q));
#endif
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
        cBALLS_FAIL(cmd,
                "%s: no unmasked bodies (nbody=%ld) were given... exiting...\n",
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

/* Native ASCII catalogs are also read inside Python.  Never call the legacy
 * in_* helpers here: their conversion failures terminate the host process. */
local int ascii_input_error(struct cmdline_data *cmd, const char *filename,
                            INTEGER row, const char *field, const char *reason)
{
    if (row == 0)
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "inputdata_ascii: %s: header %s: %s", filename, field, reason);
    else
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "inputdata_ascii: %s: data row %" INTEGER_FMT ", %s: %s",
                 filename, row, field, reason);
    return FAILURE;
}

local int ascii_token(struct cmdline_data *cmd, stream instr,
                      const char *filename, INTEGER row, const char *field,
                      char token[256])
{
    int next;
    if (fscanf(instr, "%255s", token) != 1)
        return ascii_input_error(cmd, filename, row, field,
                                 ferror(instr) ? "read error" : "unexpected end of file");
    /* A width-limited scanf must not accept a valid prefix of a longer token. */
    next = fgetc(instr);
    if (next != EOF && !isspace((unsigned char)next))
        return ascii_input_error(cmd, filename, row, field, "token exceeds 255 characters");
    if (ferror(instr))
        return ascii_input_error(cmd, filename, row, field, "read error");
    return SUCCESS;
}

local int ascii_integer(struct cmdline_data *cmd, stream instr,
                        const char *filename, INTEGER row, const char *field,
                        long minimum, long maximum, long *value)
{
    char token[256], *end;
    if (ascii_token(cmd, instr, filename, row, field, token) == FAILURE)
        return FAILURE;
    errno = 0;
    *value = strtol(token, &end, 10);
    if (errno == ERANGE || end == token || *end != '\0'
        || *value < minimum || *value > maximum)
        return ascii_input_error(cmd, filename, row, field, "invalid or out-of-range integer");
    return SUCCESS;
}

local int ascii_real(struct cmdline_data *cmd, stream instr,
                     const char *filename, INTEGER row, const char *field,
                     real *value)
{
    char token[256], *end;
    double parsed;
    if (ascii_token(cmd, instr, filename, row, field, token) == FAILURE)
        return FAILURE;
    errno = 0;
    parsed = strtod(token, &end);
    if (end == token || *end != '\0')
        return ascii_input_error(cmd, filename, row, field, "invalid real number");
    if (!isfinite(parsed))
        return ascii_input_error(cmd, filename, row, field, "non-finite real number");
    if (errno == ERANGE)
        return ascii_input_error(cmd, filename, row, field, "out-of-range real number");
    *value = (real)parsed;
    if (!isfinite(*value))
        return ascii_input_error(cmd, filename, row, field, "non-finite stored value");
    return SUCCESS;
}

/* Consume a whole header/comment line, including comments longer than 199
 * characters.  Header inspection uses the same checked path as normal input. */
local int ascii_header_line(struct cmdline_data *cmd, stream instr,
                            const char *filename, const char *field, int show)
{
    char line[200];
    do {
        if (fgets(line, sizeof(line), instr) == NULL)
            return ascii_input_error(cmd, filename, 0, field,
                                     ferror(instr) ? "read error" : "unexpected end of file");
        if (show) verb_print(cmd->verbose, "%s", line);
    } while (strchr(line, '\n') == NULL && !feof(instr));
    return SUCCESS;
}

local int inputdata_ascii_columns(struct cmdline_data *cmd,
                                   struct global_data *gd, string filename,
                                   int ifile, int layout)
{
    const int all_columns = layout == 1;
    const int positions_only = layout == 2;
    const int angles = layout == 3;
    const int file_ndim = angles ? 2 : NDIM;
    stream instr = NULL;
    bodyptr catalog = NULL, p;
    INTEGER nbody, iselect = 0;
    long count, ndim, mask;
    real box[NDIM], kavg = 0.0, sum = 0.0;
    int marker;
    int k;
#ifdef LONGINT
    const long integer_max = LONG_MAX;
#else
    const long integer_max = INT_MAX;
#endif
    const char *routineName = all_columns ? "inputdata_ascii_all" : "inputdata_ascii";

    gd->input_comment = all_columns ? "Column form input file all" : "Column form input file";
    OPEN_OUTPUT_OR_FAIL(instr, filename, "r");

    if (scanopt(cmd->options, "header-info")) {
        verb_print(cmd->verbose, "\n\t%s: header of %s\n", routineName, filename);
        if (ascii_header_line(cmd, instr, filename, "comment", TRUE) == FAILURE
            || ascii_header_line(cmd, instr, filename, "dimensions", TRUE) == FAILURE)
            goto fail;
        rewind(instr);
        if (scanopt(cmd->options, "stop")) {
            fclose(instr);
            gd->inputHeaderFlag = TRUE;
            gd->stopflag = TRUE;
            return SUCCESS;
        }
    }

    if (ascii_header_line(cmd, instr, filename, "comment", FALSE) == FAILURE)
        goto fail;
    /* The legacy format permits both "#123" and "# 123". Consume only the
     * marker character; the checked integer reader validates the full count. */
    do {
        marker = fgetc(instr);
    } while (marker != EOF && isspace((unsigned char)marker));
    if (marker != '#') {
        ascii_input_error(cmd, filename, 0, "marker",
                          marker == EOF ? (ferror(instr) ? "read error"
                                                       : "unexpected end of file")
                                        : "expected #");
        goto fail;
    }
    if (ascii_integer(cmd, instr, filename, 0, "nbody", 1, integer_max, &count) == FAILURE
        || ascii_integer(cmd, instr, filename, 0, "ndim", file_ndim, file_ndim, &ndim) == FAILURE)
        goto fail;
    /* bytes_tot uses INTEGER; check before multiplication, allocation or cast. */
    if ((uintmax_t)count > (uintmax_t)integer_max / sizeof(body)
        || (uintmax_t)count > (uintmax_t)SIZE_MAX / sizeof(body)
        || gd->bytes_tot > integer_max - count * (long)sizeof(body)) {
        ascii_input_error(cmd, filename, 0, "nbody", "catalog byte size is not representable");
        goto fail;
    }
    nbody = (INTEGER)count;
    for (k = 0; k < file_ndim; ++k) {
        if (ascii_real(cmd, instr, filename, 0, "box dimension", &box[k]) == FAILURE)
            goto fail;
    }
    if (angles) box[NDIM - 1] = box[1];
    if (cballs_calloc_checked((void **)&catalog, (size_t)nbody, sizeof(body),
                              "native ASCII catalog", cmd->error_message,
                              _ERRORMSGSIZE_) == FAILURE)
        goto fail;

    DO_BODY(p, catalog, catalog + nbody) {
        INTEGER row = (INTEGER)(p - catalog) + 1;
        for (k = 0; k < file_ndim; ++k) {
            real coordinate;
            if (ascii_real(cmd, instr, filename, row, "position", &coordinate) == FAILURE)
                goto fail;
            Pos(p)[k] = coordinate;
            if (!isfinite(Pos(p)[k])) {
                ascii_input_error(cmd, filename, row, "position", "non-finite stored value");
                goto fail;
            }
        }
        Kappa(p) = 2.0;
        if (!positions_only
            && ascii_real(cmd, instr, filename, row, "kappa", &Kappa(p)) == FAILURE)
            goto fail;
#if NDIM == 3
        if (angles) {
            real theta = Pos(p)[0], phi = Pos(p)[1];
            coordinate_transformation(cmd, gd, theta, phi, Pos(p));
        }
#endif
        /* Validate the input even when a constant-field override is requested. */
        if (!positions_only) {
            if (scanopt(cmd->options, "kappa-constant")) Kappa(p) = 2.0;
            if (scanopt(cmd->options, "kappa-constant-one")
                && (!angles || scanopt(cmd->options, "kappa-constant")))
                Kappa(p) = 1.0;
        }
        Weight(p) = 1.0;
        Mask(p) = TRUE;
        if (all_columns) {
            if (ascii_real(cmd, instr, filename, row, "weight", &Weight(p)) == FAILURE
                || ascii_integer(cmd, instr, filename, row, "mask", SHRT_MIN, SHRT_MAX, &mask) == FAILURE)
                goto fail;
            Mask(p) = (short)mask;
            if (Mask(p) == 0) iselect++;
        }
        Type(p) = BODY;
        Mass(p) = 1.0;
        Id(p) = row;
#ifdef THREEPCFSHEAR
        Gamma1(p) = 1.0;
        Gamma2(p) = 1.0;
#endif
        kavg += Kappa(p);
    }
    if (cballs_input_validate_bodies(cmd, filename, catalog, nbody) == FAILURE)
        goto fail;
    if (fclose(instr) != 0) {
        instr = NULL;
        ascii_input_error(cmd, filename, 0, "stream", "close error");
        goto fail;
    }
    instr = NULL;

    if (scanopt(cmd->options, "check-eq-pos")) {
        bodyptr q;
        real dist2;
        vector distv;
        DO_BODY(p, catalog, catalog + nbody - 1)
            DO_BODY(q, p + 1, catalog + nbody) {
                DOTPSUBV(dist2, distv, Pos(p), Pos(q));
                if (dist2 == 0.0) {
                    ascii_input_error(cmd, filename, (INTEGER)(q - catalog) + 1,
                                      "position", "at least two bodies have same position");
                    goto fail;
                }
            }
    }

    /* Publish only a completely read and validated catalog.  Earlier catalogs
     * remain owned by the caller's normal cleanup if this one fails. */
    bodytable[ifile] = catalog;
    gd->bodytable_allocated = TRUE;
    gd->nbodyTable[ifile] = cmd->nbody = nbody;
    gd->bytes_tot += (INTEGER)((size_t)nbody * sizeof(body));
    DO_COORD(k) gd->Box[k] = box[k];

    verb_print(cmd->verbose, "\t%s: nbody and ndim: %" INTEGER_FMT " %ld...\n",
               routineName, nbody, ndim);
    verb_print(cmd->verbose, "\t%s: lbox dimensions: ", routineName);
    DO_COORD(k) verb_print(cmd->verbose, "%g ", box[k]);
    verb_print(cmd->verbose, "\n");
    if (all_columns) {
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\t%s: masked pixels = %" INTEGER_FMT "\n", routineName, iselect);
        verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                               "\t%s: unmasked pixels = %" INTEGER_FMT "\n", routineName, nbody - iselect);
    }
    kavg /= (real)nbody;
    DO_BODY(p, catalog, catalog + nbody) sum += rsqr(Kappa(p) - kavg);
    verb_print_normal_info(cmd->verbose, cmd->verbose_log, gd->outlog,
                           "%s: average and std dev of kappa (%" INTEGER_FMT " particles) = %le %le\n",
                           routineName, nbody, kavg, nbody > 1 ? rsqrt(sum / ((real)nbody - 1.0)) : 0.0);
    return SUCCESS;

fail:
    if (instr != NULL) fclose(instr);
    free(catalog);
    return FAILURE;
}

local int inputdata_ascii(struct cmdline_data *cmd, struct global_data *gd,
                          string filename, int ifile)
{
    return inputdata_ascii_columns(cmd, gd, filename, ifile, FALSE);
}

local int inputdata_ascii_all(struct cmdline_data *cmd, struct global_data *gd,
                              string filename, int ifile)
{
    return inputdata_ascii_columns(cmd, gd, filename, ifile, TRUE);
}

/* Keep the native binary layout: INTEGER count, int dimension, real box,
 * position arrays, scalar arrays, and optional real weights / short masks. */
local int binary_input_read(struct cmdline_data *cmd, stream instr,
                             const char *filename, INTEGER row,
                             const char *field, void *value,
                             size_t item_size, size_t count)
{
    if (fread(value, item_size, count, instr) == count) return SUCCESS;
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "inputdata_binary: %s: %s row %" INTEGER_FMT ", %s: %s",
             filename, row ? "data" : "header", row, field,
             ferror(instr) ? "read error" : "unexpected end of file");
    return FAILURE;
}

local int inputdata_binary_columns(struct cmdline_data *cmd,
                                    struct global_data *gd, string filename,
                                    int ifile, int all_columns)
{
    stream instr = NULL;
    bodyptr catalog = NULL, p;
    INTEGER nbody;
    int ndim, k;
    real box[NDIM], disk_position[NDIM];
#ifdef LONGINT
    const uintmax_t integer_max = LONG_MAX;
#else
    const uintmax_t integer_max = INT_MAX;
#endif
    gd->input_comment = all_columns ? "Binary-all input file" : "Binary input file";
    if (stropen_checked(filename, "rb", &instr, cmd->error_message,
                         _ERRORMSGSIZE_) == FAILURE)
        return FAILURE;
#define BINARY_READ(row, field, value, size, count) \
    do { if (binary_input_read(cmd, instr, filename, row, field, value, size, count) \
              == FAILURE) goto fail; } while (0)
    BINARY_READ(0, "nbody", &nbody, sizeof(nbody), 1);
    BINARY_READ(0, "ndim", &ndim, sizeof(ndim), 1);
    if (nbody < 1 || (uintmax_t)nbody > integer_max / sizeof(body)
        || (uintmax_t)nbody > SIZE_MAX / sizeof(body)
        || gd->bytes_tot < 0
        || (uintmax_t)gd->bytes_tot > integer_max - (uintmax_t)nbody * sizeof(body)
        || ndim != NDIM) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "inputdata_binary: %s: invalid header dimensions or catalog byte size",
                 filename);
        goto fail;
    }
    BINARY_READ(0, "box dimension", box, sizeof(real), NDIM);
    DO_COORD(k)
        if (cballs_input_finite(cmd, filename, 0, "box dimension", box[k]) == FAILURE)
            goto fail;
    if (scanopt(cmd->options, "header-info")) {
        verb_print(cmd->verbose, "\nBinary header of %s: %" INTEGER_FMT " %d\n",
                   filename, nbody, ndim);
        DO_COORD(k) verb_print(cmd->verbose, " %g", box[k]);
        verb_print(cmd->verbose, "\n");
        if (scanopt(cmd->options, "stop")) {
            if (fclose(instr) != 0) { instr = NULL; goto close_fail; }
            gd->inputHeaderFlag = TRUE;
            gd->stopflag = TRUE;
            return SUCCESS;
        }
    }
    if (cballs_calloc_checked((void **)&catalog, (size_t)nbody, sizeof(body),
                              "native binary catalog", cmd->error_message,
                              _ERRORMSGSIZE_) == FAILURE)
        goto fail;
    DO_BODY(p, catalog, catalog + nbody) {
        INTEGER row = p - catalog + 1;
        BINARY_READ(row, "position", disk_position, sizeof(real), NDIM);
        DO_COORD(k) {
            if (cballs_input_finite(cmd, filename, (size_t)row,
                                    "position", disk_position[k]) == FAILURE)
                goto fail;
            Pos(p)[k] = (cballs_storage_real)disk_position[k];
        }
        Type(p) = BODY; Mass(p) = 1.0; Weight(p) = 1.0;
        Mask(p) = TRUE; Id(p) = row;
#ifdef THREEPCFSHEAR
        Gamma1(p) = Gamma2(p) = 1.0;
#endif
    }
    DO_BODY(p, catalog, catalog + nbody) {
        INTEGER row = p - catalog + 1;
        BINARY_READ(row, "kappa", &Kappa(p), sizeof(real), 1);
        if (cballs_input_finite(cmd, filename, (size_t)row, "kappa", Kappa(p)) == FAILURE)
            goto fail;
        if (scanopt(cmd->options, "kappa-constant")) Kappa(p) = 2.0;
        if (scanopt(cmd->options, "kappa-constant-one")) Kappa(p) = 1.0;
    }
    if (all_columns) {
        DO_BODY(p, catalog, catalog + nbody)
            BINARY_READ(p - catalog + 1, "weight", &Weight(p), sizeof(real), 1);
        DO_BODY(p, catalog, catalog + nbody)
            BINARY_READ(p - catalog + 1, "mask", &Mask(p), sizeof(short), 1);
    }
    if (cballs_input_validate_bodies(cmd, filename, catalog, nbody) == FAILURE)
        goto fail;
    if (fclose(instr) != 0) { instr = NULL; goto close_fail; }
    bodytable[ifile] = catalog;
    gd->bodytable_allocated = TRUE;
    gd->nbodyTable[ifile] = cmd->nbody = nbody;
    gd->bytes_tot += (INTEGER)((size_t)nbody * sizeof(body));
    DO_COORD(k) gd->Box[k] = box[k];
    return SUCCESS;
close_fail:
    snprintf(cmd->error_message, _ERRORMSGSIZE_,
             "inputdata_binary: %s: close error", filename);
fail:
    if (instr != NULL) fclose(instr);
    free(catalog);
    return FAILURE;
#undef BINARY_READ
}

local int inputdata_bin(struct cmdline_data *cmd, struct global_data *gd,
                         string filename, int ifile)
{
    return inputdata_binary_columns(cmd, gd, filename, ifile, FALSE);
}

local int inputdata_bin_all(struct cmdline_data *cmd, struct global_data *gd,
                             string filename, int ifile)
{
    return inputdata_binary_columns(cmd, gd, filename, ifile, TRUE);
}

//B BEGIN:: Reading Takahasi simulations
//From Takahashi web page. Adapted to our needs

#include<math.h>
#include<stdio.h>
#include<stdlib.h>


void pix2ang(long pix, int nside, double *theta, double *phi);

local int Takahasi_region_selection(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, long npix,
            float *conv, float *shear1, float *shear2, float *rotat, int ifile);
local int Takahasi_region_selection_3d_all(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, long npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax,
            real *zmin, real *zmax, int ifile);
local int Takahasi_region_selection_3d(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, long npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax,
            real *zmin, real *zmax, int ifile);

#if NDIM == 2
local int Takahasi_region_selection_2d(struct cmdline_data* cmd, 
            struct  global_data* gd,
            int nside, long npix,
            float *conv, float *shear1, float *shear2, float *rotat,
            real dtheta_rot, real thetaL, real thetaR,
            real dphi_rot, real phiL, real phiR,
            real *xmin, real *xmax, real *ymin, real *ymax, int ifile);
#endif

/* Native Takahashi layout: fixed header followed by four float maps and
 * legacy inter-record markers. Check every transfer before selection/override. */
local int inputdata_takahashi(struct cmdline_data *cmd, struct global_data *gd,
                              string filename, int ifile)
{
    FILE *fp = NULL;
    int marker, nside, status = FAILURE;
    long npix, record;
    float *maps[4] = {NULL, NULL, NULL, NULL};
    const char *fields[] = {"convergence", "gamma1", "gamma2", "rotation"};
    const long boundaries[] = {536870908L, 1073741818L, 1610612728L,
                              2147483638L, 2684354547L, 3221225457L};
    gd->input_comment = "Takahashi input file";
    if (stropen_checked(filename, "rb", &fp, cmd->error_message,
                         _ERRORMSGSIZE_) == FAILURE) return FAILURE;
#define TAK_READ(row, field, ptr, size, count) \
    do { if (binary_input_read(cmd, fp, filename, row, field, ptr, size, count) \
                == FAILURE) goto cleanup; } while (0)
    TAK_READ(0, "record marker", &marker, sizeof(marker), 1);
    TAK_READ(0, "nside", &nside, sizeof(nside), 1);
    TAK_READ(0, "npix", &npix, sizeof(npix), 1);
    TAK_READ(0, "record marker", &record, sizeof(record), 1);
    if (nside <= 0 || npix <= 0
        || (uintmax_t)nside > (uintmax_t)LONG_MAX / 12 / (uintmax_t)nside
        || (uintmax_t)npix != 12 * (uintmax_t)nside * (uintmax_t)nside
        || (uintmax_t)npix > SIZE_MAX / sizeof(body)
        || (uintmax_t)npix > (uintmax_t)LONG_MAX / sizeof(body)) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "Takahashi input '%s': invalid nside/npix or catalog byte size", filename);
        goto cleanup;
    }
    if (scanopt(cmd->options, "header-info")) {
        verb_print(cmd->verbose, "Takahashi header: nside=%d npix=%ld\n", nside, npix);
        if (scanopt(cmd->options, "stop")) {
            gd->inputHeaderFlag = gd->stopflag = TRUE;
            status = SUCCESS;
            goto cleanup;
        }
    }
    for (int field = 0; field < 4; ++field) {
        if (field) TAK_READ(0, "record marker", &record, sizeof(record), 1);
        if (cballs_malloc_checked((void **)&maps[field], (size_t)npix, sizeof(float),
                                   "Takahashi map", cmd->error_message,
                                   _ERRORMSGSIZE_) == FAILURE) goto cleanup;
        /* Chunking preserves the historical record boundaries without one
         * fread call per pixel on survey-size maps. */
        long first = 0;
        for (int block = 0; block <= 6 && first < npix; ++block) {
            long end = block < 6 && boundaries[block] < npix
                ? boundaries[block] + 1 : npix;
            TAK_READ(first + 1, fields[field], maps[field] + first,
                     sizeof(float), (size_t)(end - first));
            if (block < 6 && end == boundaries[block] + 1)
                TAK_READ(end, "record marker", &record, sizeof(record), 1);
            first = end;
        }
        if (cballs_input_finite_array(cmd, filename, maps[field], (size_t)npix,
                                      sizeof(float), fields[field]) == FAILURE) goto cleanup;
    }
    status = Takahasi_region_selection(cmd, gd, nside, npix,
                                       maps[0], maps[1], maps[2], maps[3], ifile);
cleanup:
    if (fp && fclose(fp) != 0 && status == SUCCESS) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "Takahashi input '%s': close error", filename);
        status = FAILURE;
    }
    for (int i = 0; i < 4; ++i) free(maps[i]);
#undef TAK_READ
    return status;
}

local int Takahasi_region_selection(struct cmdline_data* cmd, 
                                    struct  global_data* gd,
                                    int nside, long npix,
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
        if (Takahasi_region_selection_3d(cmd, gd,
                                     nside, npix, conv, shear1, shear2, rotat,
                dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, ifile) == FAILURE)
            return FAILURE;
    } else {
        if (Takahasi_region_selection_3d_all(cmd, gd,
                                         nside, npix, conv, shear1, shear2, rotat,
            dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
            &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, ifile) == FAILURE)
            return FAILURE;
    }
#else   // ! TREEDIM
    real xmin, ymin;
    real xmax, ymax;

    if (Takahasi_region_selection_2d(cmd, gd,
                nside, npix, conv, shear1, shear2, rotat,
                dtheta_rot, thetaL, thetaR, dphi_rot, phiL, phiR,
                &xmin, &xmax, &ymin, &ymax, ifile) == FAILURE)
            return FAILURE;
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
                                           int nside, long npix,
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
    if (cmd->nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "Takahashi input: selection contains no bodies");
        return FAILURE;
    }
    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)cmd->nbody,
                               sizeof(body), "Takahashi catalog", cmd->error_message,
                               _ERRORMSGSIZE_) == FAILURE) return FAILURE;
    gd->bodytable_allocated = TRUE;
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

#if defined(NMultipoles) && defined(NONORMHIST)
    if (scanopt(cmd->options, "patch-with-all")) {
        verb_print(cmd->verbose,
            "\nsearchcalc_tc_kkk_omp: total number of pixels to be pivots: %ld\n",
            gd->pivotCount);
    }
#endif
    
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
                                       int nside, long npix,
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
    if (cballs_calloc_checked((void **)&bodytabtmp, (size_t)cmd->nbody,
                               sizeof(body), "Takahashi selection", cmd->error_message,
                               _ERRORMSGSIZE_) == FAILURE) return FAILURE;
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
    if (cmd->nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "Takahashi input: selection contains no bodies");
        free(bodytabtmp); return FAILURE;
    }
    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)cmd->nbody,
                               sizeof(body), "Takahashi catalog", cmd->error_message,
                               _ERRORMSGSIZE_) == FAILURE) { free(bodytabtmp); return FAILURE; }
    gd->bodytable_allocated = TRUE;

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
            int nside, long npix,
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
    if (cballs_calloc_checked((void **)&bodytabtmp, (size_t)cmd->nbody,
                               sizeof(body), "Takahashi selection", cmd->error_message,
                               _ERRORMSGSIZE_) == FAILURE) return FAILURE;
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
    if (cmd->nbody < 1) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "Takahashi input: selection contains no bodies");
        free(bodytabtmp); return FAILURE;
    }
    if (cballs_calloc_checked((void **)&bodytable[ifile], (size_t)cmd->nbody,
                               sizeof(body), "Takahashi catalog", cmd->error_message,
                               _ERRORMSGSIZE_) == FAILURE) { free(bodytabtmp); return FAILURE; }
    gd->bodytable_allocated = TRUE;
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
    int output_status = SUCCESS;
    int output_owner = TRUE;
#ifdef CBALLS_MPI_ENABLED
    output_owner = cballs_mpi_output_enabled(cmd);
#endif
    double cpustart = CPUTIME;
    if (output_owner && !strnull(cmd->outfile))
        output_status = outputdata(cmd, gd, btable[ifile], nbody[ifile]);
    if (output_owner && output_status == SUCCESS)
        gd->cputotalinout += CPUTIME - cpustart;

#ifdef CBALLS_MPI_ENABLED
    output_status = cballs_mpi_consensus(
        cmd, output_status, "MPI catalog output");
#endif
    return output_status;
}

local int outputdata(struct cmdline_data* cmd, struct  global_data* gd,
                     bodyptr btable, INTEGER nbody)
{
    switch(outfilefmt_int) {
        case OUTCOLUMNS:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            class_call_cballs(outputdata_ascii(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;
        case OUTCOLUMNSALL:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            class_call_cballs(outputdata_ascii_all(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;
        case OUTCOLUMNSBIN:
            verb_print(cmd->verbose, "\n\tbinary format output\n");
            class_call_cballs(outputdata_bin(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;
        case OUTCOLUMNSBINALL:
            verb_print(cmd->verbose, "\n\tbinary-all format output\n");
            class_call_cballs(outputdata_bin_all(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;
        case OUTNULL:
            verb_print(cmd->verbose, "\n\tcolumns-ascii format output\n");
            class_call_cballs(outputdata_ascii(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;

//B socket:
#ifdef ADDONS
#include "cballsio_include_03.h"
#endif
//E

        default:
            verb_print(cmd->verbose, 
                    "\n\toutput: Unknown output format...\n\tprinting in default format (columns-ascii)...\n");
            class_call_cballs(outputdata_ascii(cmd, gd, btable, nbody),
                                  errmsg, errmsg);
            break;
    }

    return SUCCESS;
}

local int outputdata_ascii(struct cmdline_data* cmd, struct  global_data* gd,
                             bodyptr bodytab, INTEGER nbody)
{
    string routineName = "outputdata_ascii";
    char namebuf[256];
    stream outstr;
    bodyptr p;

    if (format_checked(namebuf, sizeof(namebuf),
                       "output filename", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;

    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

#if NDIM == 3
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "# nbody NDIM Lx Ly Lz\n# %" INTEGER_FMT " %d ",
                         nbody, NDIM);
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "%lf %lf %lf\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "# nbody NDIM Lx Ly\n# %" INTEGER_FMT " %d ",
                         nbody, NDIM);
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "%lf %lf\n",gd->Box[0],gd->Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_vector_mar_checked(outstr, Pos(p), routineName, namebuf,
                                   cmd->error_message,
                                   _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
        if (out_real_mar_checked(outstr, Kappa(p), routineName, namebuf,
                                 cmd->error_message,
                                 _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
#ifdef DEBUG
        if (out_bool_mar_checked(outstr, HIT(p), routineName, namebuf,
                                 cmd->error_message,
                                 _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
#endif
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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

    if (format_checked(namebuf, sizeof(namebuf),
                       "output filename", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;

    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

#if NDIM == 3
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "# nbody NDIM Lx Ly Lz\n# %" INTEGER_FMT " %d ",
                         nbody, NDIM);
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "%lf %lf %lf\n",gd->Box[0],gd->Box[1],gd->Box[2]);
#else
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "# nbody NDIM Lx Ly\n# %" INTEGER_FMT " %d ",
                         nbody, NDIM);
    WRITE_OUTPUT_OR_FAIL(outstr, namebuf,
                         "%lf %lf\n",gd->Box[0],gd->Box[1]);
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_vector_mar_checked(outstr, Pos(p), routineName, namebuf,
                                 cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }

        if (out_real_mar_checked(outstr, Kappa(p), routineName, namebuf,
                                 cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }

            if (out_real_mar_checked(outstr, Weight(p), routineName, namebuf,
                                     cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
                if (outstr != NULL) fclose(outstr);
                return FAILURE;
            }

        if (out_short_mar_checked(outstr, Mask(p), routineName, namebuf,
                                 cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }

#ifdef DEBUG
        if (out_bool_mar_checked(outstr, HIT(p), routineName, namebuf,
                                 cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
#endif
        WRITE_OUTPUT_OR_FAIL(outstr, namebuf, "\n");
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    verb_print(cmd->verbose, "\t%s: data output to file %s\n",
               routineName, namebuf);

    return SUCCESS;
}

local int outputdata_bin(struct cmdline_data* cmd, struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody)
{
    string routineName = "outputdata_bin";
    char namebuf[256];
    stream outstr;
    bodyptr p;

    //B
    if (format_checked(namebuf, sizeof(namebuf),
                       "output filename", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;
    //E

    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

    if (out_int_bin_long_checked(outstr, nbody, routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_int_bin_checked(outstr, NDIM, routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_real_bin_checked(outstr, gd->Box[0], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_real_bin_checked(outstr, gd->Box[1], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
#if NDIM == 3
    if (out_real_bin_checked(outstr, gd->Box[2], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_vector_bin_checked(outstr, Pos(p), routineName, namebuf,
                                   cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }

    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_real_bin_checked(outstr, Kappa(p), routineName, namebuf,
                                   cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
    verb_print(cmd->verbose, "\tdata output to file %s\n", namebuf);

    return SUCCESS;
}

local int outputdata_bin_all(struct cmdline_data* cmd, struct  global_data* gd,
                         bodyptr bodytab, INTEGER nbody)
{
    string routineName = "outputdata_bin_all";
    char namebuf[256];
    stream outstr;
    bodyptr p;

    if (format_checked(namebuf, sizeof(namebuf),
                       "output filename", "%s", gd->fpfnameOutputFileName) != 0)
        return FAILURE;

    OPEN_OUTPUT_OR_FAIL(outstr, namebuf, "w!");

    if (out_int_bin_long_checked(outstr, nbody, routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_int_bin_checked(outstr, NDIM, routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_real_bin_checked(outstr, gd->Box[0], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
    if (out_real_bin_checked(outstr, gd->Box[1], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
#if NDIM == 3
    if (out_real_bin_checked(outstr, gd->Box[2], routineName, namebuf,
                               cmd->error_message, _ERRORMSGSIZE_) == FAILURE) {
        if (outstr != NULL) fclose(outstr);
        return FAILURE;
    }
#endif
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_vector_bin_checked(outstr, Pos(p), routineName, namebuf,
                                   cmd->error_message,
                                   _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_real_bin_checked(outstr, Kappa(p), routineName, namebuf,
                                   cmd->error_message,
                                 _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_real_bin_checked(outstr, Weight(p), routineName, namebuf,
                                   cmd->error_message,
                                 _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_short_bin_checked(outstr, Mask(p), routineName, namebuf,
                                   cmd->error_message,
                                  _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
#ifdef DEBUG
    DO_BODY(p, bodytab, bodytab+nbody) {
        if (out_bool_bin_checked(outstr, HIT(p), routineName, namebuf,
                                   cmd->error_message,
                                 _ERRORMSGSIZE_) == FAILURE) {
            if (outstr != NULL) fclose(outstr);
            return FAILURE;
        }
    }
#endif
    CLOSE_OUTPUT_OR_FAIL(outstr, namebuf);
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

//B I/O directories:
global int setFilesDirs_log(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "setFilesDirs_log";

    if (cmd->verbose_log>0) {           // gd->logfilePath is defined
        if (format_checked(gd->tmpDir, sizeof(gd->tmpDir),
                           "tmpDir", "%s/tmp", cmd->rootDir) != 0)
            return FAILURE;

        double cpustart = CPUTIME;

        if (mkdir_p(gd->tmpDir, 0777) != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: cannot create directory '%s': %s",
                     routineName, gd->tmpDir, strerror(errno));
            return FAILURE;
        }

        gd->cputotalinout += CPUTIME - cpustart;

        if (format_checked(gd->logfilePath, sizeof(gd->logfilePath),
                           "gd->logfilePath", "%s/cballs%s.log",
                           gd->tmpDir,cmd->suffixOutFiles) != 0)
            return FAILURE;
    }
    
    return SUCCESS;
}

global int setFilesDirs(struct cmdline_data* cmd, struct  global_data* gd)
{
    string routineName = "setFilesDirs";
    double cpustart = CPUTIME;
    int rc = FAILURE;

    if (gd->rootDirFlag == TRUE) {
        if (mkdir_p(cmd->rootDir, 0777) != 0) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "%s: cannot create directory '%s': %s",
                     routineName, cmd->rootDir, strerror(errno));
            goto fail;
        }
        gd->cputotalinout += CPUTIME - cpustart;
        
        if (format_checked(gd->fpfnameOutputFileName, sizeof(gd->fpfnameOutputFileName),
            "fpfnameOutputFileName", "%s/%s%s%s",
            cmd->rootDir,cmd->outfile,cmd->suffixOutFiles,EXTFILES) != 0)
            goto fail;
    
        if (format_checked(gd->fpfnamehistNNFileName, sizeof(gd->fpfnamehistNNFileName),
            "fpfnamehistNNFileName", "%s/%s%s%s",
            cmd->rootDir,cmd->histNNFileName,cmd->suffixOutFiles,EXTFILES) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistCFFileName, sizeof(gd->fpfnamehistCFFileName),
            "fpfnamehistCFFileName", "%s/%s%s%s",
            cmd->rootDir,"histCF",cmd->suffixOutFiles,EXTFILES) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistrBinsFileName, sizeof(gd->fpfnamehistrBinsFileName),
            "fpfnamehistrBinsFileName", "%s/%s%s%s",
            cmd->rootDir,"rbins",cmd->suffixOutFiles,EXTFILES) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistXi2pcfFileName, sizeof(gd->fpfnamehistXi2pcfFileName),
            "fpfnamehistXi2pcfFileName", "%s/%s",
            cmd->rootDir,cmd->histXi2pcfFileName) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistZetaGFileName, sizeof(gd->fpfnamehistZetaGFileName),
            "fpfnamehistZetaGFileName", "%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"G",cmd->suffixOutFiles) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistZetaGmFileName, sizeof(gd->fpfnamehistZetaGmFileName),
            "fpfnamehistZetaGmFileName", "%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"G",cmd->suffixOutFiles) != 0)
            goto fail;

        if (format_checked(gd->fpfnamehistZetaMFileName, sizeof(gd->fpfnamehistZetaMFileName),
            "fpfnamehistZetaMFileName", "%s/%s%s%s",
            cmd->rootDir,cmd->histZetaFileName,"M",cmd->suffixOutFiles) != 0)
            goto fail;

        if (format_checked(gd->fpfnamemhistZetaMFileName, sizeof(gd->fpfnamemhistZetaMFileName),
            "fpfnamemhistZetaMFileName", "%s/%s%s%s%s",
            cmd->rootDir,"m",cmd->histZetaFileName,"M",cmd->suffixOutFiles) != 0)
            goto fail;

        if (format_checked(gd->fpfnameCPUFileName, sizeof(gd->fpfnameCPUFileName),
            "fpfnameCPUFileName", "%s/cputime%s%s",
            cmd->rootDir,cmd->suffixOutFiles,EXTFILES) != 0)
            goto fail;

//B socket:
#ifdef ADDONS
#include "cballsio_include_09b.h"
#endif
//E
    } // ! rootDirFlag
    
    rc = SUCCESS;

fail:
    return rc;


}
//E

//B
local int EndRun_CloseLog(struct global_data *gd)
{
    if (gd->outlogFlagFree == TRUE && gd->outlog != NULL) {
        fclose(gd->outlog);
        gd->outlog = NULL;
        gd->outlogFlagFree = FALSE;
    }

    return SUCCESS;
}
//E

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
    string routineName = "EndRun";
    stream outstr;

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
            OPEN_OUTPUT_OR_FAIL(outstr, gd->fpfnameCPUFileName, "a");
            WRITE_OUTPUT_OR_FAIL(outstr, gd->fpfnameCPUFileName,
                                "%g %g %g\n",
                                (double) cmd->nbody, cpuTotal, gd->cpusearch);
            CLOSE_OUTPUT_OR_FAIL(outstr, gd->fpfnameCPUFileName);
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
    string routineName = "EndRun_FreeMemory";

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

    EndRun_CloseLog(gd);

    return SUCCESS;
}

global int EndRun_FreeMemory_tree(struct cmdline_data* cmd,
                                  struct global_data* gd)
{
    int ifile;

    for (ifile = 0; ifile < MAXITEMS; ifile++) {
        if (nodetablescanlev[ifile] != NULL) {
            free(nodetablescanlev[ifile]);
            nodetablescanlev[ifile] = NULL;
        }
        gd->nnodescanlevTable[ifile] = 0;

        if (nodetablescanlev_root[ifile] != NULL) {
            free(nodetablescanlev_root[ifile]);
            nodetablescanlev_root[ifile] = NULL;
        }
        gd->nnodescanlev_rootTable[ifile] = 0;

    #ifdef CBALLS_NEEDS_BALLS4_SCAN
        if (nodetablescanlevB4[ifile] != NULL) {
            free(nodetablescanlevB4[ifile]);
            nodetablescanlevB4[ifile] = NULL;
        }
        gd->nnodescanlevTableB4[ifile] = 0;
    #endif
    }

    if (!scanopt(cmd->searchMethod, "kdtree-omp")
        && !scanopt(cmd->searchMethod, "kdtree-box-omp")
        && !scanopt(cmd->searchMethod, "balltree-omp")
        && !scanopt(cmd->searchMethod, "balltree-mpi")
        && !scanopt(cmd->searchMethod, "balltree-2balls-omp")) {
        freeTree(cmd, gd);
    }

    gd->tree_allocated = FALSE;
    return SUCCESS;
}


global int EndRun_FreeMemory_bodytable(struct cmdline_data* cmd,
                                       struct global_data* gd)
{
    int ifile;

    for (ifile = 0; ifile < MAXITEMS; ifile++) {
        if (bodytable[ifile] != NULL) {
            free(bodytable[ifile]);
            bodytable[ifile] = NULL;
        }
        gd->nbodyTable[ifile] = 0;
    }

#if defined(DEBUG) && defined(BODYTABBF_ON)
    if (bodytabbf != NULL) {
        free(bodytabbf);
        bodytabbf = NULL;
    }
#endif
    
    gd->bodytable_allocated = FALSE;

    return SUCCESS;
}


global int EndRun_FreeMemory_histograms(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    cballs_scalar_window_free(gd);
    //B added by cBalls
#define FREE_DVECTOR_NULL(p,nl,nh) \
    do { if ((p) != NULL) { free_dvector((p),(nl),(nh)); (p) = NULL; } } while (0)

#define FREE_DMATRIX_NULL(p,nrl,nrh,ncl,nch) \
    do { if ((p) != NULL) { free_dmatrix((p),(nrl),(nrh),(ncl),(nch)); (p) = NULL; } } while (0)

#define FREE_DMATRIX3D_NULL(p,nrl,nrh,ncl,nch,ndl,ndh) \
    do { if ((p) != NULL) { free_dmatrix3D((p),(nrl),(nrh),(ncl),(nch),(ndl),(ndh)); (p) = NULL; } } while (0)
    //E

    FREE_DVECTOR_NULL(gd->histN2pcf, 1, cmd->sizeHistN);
    // 2pcf
#ifdef SMOOTHPIVOT
    FREE_DVECTOR_NULL(gd->histNNSubN2pcftotal,1,cmd->sizeHistN);
#endif
    FREE_DVECTOR_NULL(gd->histNNSubN2pcf, 1, cmd->sizeHistN);
    //E


//B socket:
#ifdef ADDONS
#include "cballsio_include_10.h"                    // this is empty and can
                                                    //be remove these 3 lines
#endif
//E

    
#ifdef TPCF
    FREE_DMATRIX3D_NULL(gd->histZetaGmIm,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaGmRe,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

        // Transpose of Zm(ti) X Ym(tj) = Zm(tj) X Ym(ti)
    FREE_DMATRIX3D_NULL(gd->histZetaMcossin,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaMsincos,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaMsin,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaMcos,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

        // (EE) edge_effects
    FREE_DMATRIX3D_NULL(gd->histZetaM_EE_Im,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaM_EE,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);

    FREE_DMATRIX3D_NULL(gd->histZetaM,
                        1, cmd->mChebyshev+1,
                        1, cmd->sizeHistN,
                        1, cmd->sizeHistN);
    
    FREE_DMATRIX_NULL(gd->histXisin, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN);
    FREE_DMATRIX_NULL(gd->histXicos, 1, cmd->mChebyshev+1, 1, cmd->sizeHistN);

#endif


    //B cross
    FREE_DVECTOR_NULL(gd->histXi2pcf13, 1, cmd->sizeHistN);
    FREE_DVECTOR_NULL(gd->histXi2pcf12, 1, cmd->sizeHistN);
    //E
    FREE_DVECTOR_NULL(gd->histXi2pcf, 1, cmd->sizeHistN);

    FREE_DVECTOR_NULL(gd->histNNN, 1, cmd->sizeHistN);
    // 2pcf
#ifdef SMOOTHPIVOT
    FREE_DVECTOR_NULL(gd->histNNSubXi2pcftotal, 1, cmd->sizeHistN);
#endif
    FREE_DVECTOR_NULL(gd->histNNSubXi2pcf, 1, cmd->sizeHistN);
    //
    FREE_DVECTOR_NULL(gd->histNNSub, 1, cmd->sizeHistN);
    FREE_DVECTOR_NULL(gd->histCF, 1, cmd->sizeHistN);
    FREE_DVECTOR_NULL(gd->histNN, 1, cmd->sizeHistN);

    //B Histogram arrays PXD versions
#ifdef PXD
    FREE_DVECTOR_NULL(gd->histZetaMFlatten, 1, 0); /* unused compatibility slot */
    FREE_DVECTOR_NULL(gd->rBins, 1, cmd->sizeHistN);
    //B offset at 0 in order to work with Cython
    FREE_DMATRIX_NULL(gd->matPXD, 0, cmd->sizeHistN-1, 0, cmd->sizeHistN-1);
    //E
    FREE_DVECTOR_NULL(gd->vecPXD, 1, cmd->sizeHistN);
#endif
    //E Histogram arrays PXD versions

    gd->histograms_allocated = FALSE;
    gd->histogram_results_ready = FALSE;
    gd->histogram_products = 0;

#undef FREE_DVECTOR_NULL
#undef FREE_DMATRIX_NULL
#undef FREE_DMATRIX3D_NULL

    return SUCCESS;
}

global int EndRun_FreeMemory_gd(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "EndRun_FreeMemory_gd";
    int ifile;

    //B Set gsl uniform random :: If not needed globally
    //      this line have to go to testdata
    #ifdef USEGSL
        if (gd->random_allocated == TRUE && gd->r != NULL) {
            gsl_rng_free(gd->r);        // allocated by random_init
            gd->r = NULL;
        }
        r_gsl = NULL;
    #endif
    //E

    gd->random_allocated = FALSE;

    gd->gd_allocated = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_gd_2(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    if (gd->deltaRV != NULL) {
        free_dvector(gd->deltaRV, 1, cmd->sizeHistN);
        gd->deltaRV = NULL;
    }
    if (gd->ddeltaRV != NULL) {
        free_dvector(gd->ddeltaRV, 1, cmd->sizeHistN - 1);
        gd->ddeltaRV = NULL;
    }

    gd->gd_allocated_2 = FALSE;

    return SUCCESS;
}

global int EndRun_FreeMemory_cmd(struct cmdline_data* cmd,
                             struct  global_data* gd)
{
    string routineName = "EndRun_FreeMemory_cmd";

#ifdef SAVERESTORE
if (gd->restorefileFlag == TRUE && cmd->restorefile != NULL) {
    free(cmd->restorefile);
    cmd->restorefile = NULL;
    gd->restorefileFlag = FALSE;
}

if (gd->statefileFlag == TRUE && cmd->statefile != NULL) {
    free(cmd->statefile);
    cmd->statefile = NULL;
    gd->statefileFlag = FALSE;
}
#endif

#ifdef IOLIB
    if (gd->columnsFlag == TRUE && cmd->columns != NULL) {
        free(cmd->columns);
        cmd->columns = NULL;
        gd->columnsFlag = FALSE;
    }
#endif
    
    //B
    if (gd->outfileFlag == TRUE && cmd->outfile != NULL) {
        free(cmd->outfile);
        cmd->outfile = NULL;
        gd->outfileFlag = FALSE;
    }

    if (gd->outfilefmtFlag == TRUE && cmd->outfilefmt != NULL) {
        free(cmd->outfilefmt);
        cmd->outfilefmt = NULL;
        gd->outfilefmtFlag = FALSE;
    }

    if (gd->histNNFileNameFlag == TRUE && cmd->histNNFileName != NULL) {
        free(cmd->histNNFileName);
        cmd->histNNFileName = NULL;
        gd->histNNFileNameFlag = FALSE;
    }

    if (gd->histXi2pcfFileNameFlag == TRUE && cmd->histXi2pcfFileName != NULL) {
        free(cmd->histXi2pcfFileName);
        cmd->histXi2pcfFileName = NULL;
        gd->histXi2pcfFileNameFlag = FALSE;
    }

    if (gd->histZetaFileNameFlag == TRUE && cmd->histZetaFileName != NULL) {
        free(cmd->histZetaFileName);
        cmd->histZetaFileName = NULL;
        gd->histZetaFileNameFlag = FALSE;
    }

    if (gd->suffixOutFilesFlag == TRUE && cmd->suffixOutFiles != NULL) {
        free(cmd->suffixOutFiles);
        cmd->suffixOutFiles = NULL;
        gd->suffixOutFilesFlag = FALSE;
    }

    if (gd->testmodelFlag == TRUE && cmd->testmodel != NULL) {
        free(cmd->testmodel);
        cmd->testmodel = NULL;
        gd->testmodelFlag = FALSE;
    }

    if (gd->preScriptFlag == TRUE && cmd->preScript != NULL) {
        free(cmd->preScript);
        cmd->preScript = NULL;
        gd->preScriptFlag = FALSE;
    }
    
    if (gd->posScriptFlag == TRUE && cmd->posScript != NULL) {
        free(cmd->posScript);
        cmd->posScript = NULL;
        gd->posScriptFlag = FALSE;
    }
    //E
    
    if (gd->optionsFlag == TRUE && cmd->options != NULL) {
        free(cmd->options);
        cmd->options = NULL;
        gd->optionsFlag = FALSE;
    }
    
    if (gd->rootDirFlagFree == TRUE && cmd->rootDir != NULL) {
        free(cmd->rootDir);
        cmd->rootDir = NULL;
        gd->rootDirFlagFree = FALSE;
    }
    
    if (gd->iCatalogsFlag == TRUE && cmd->iCatalogs != NULL) {
        free(cmd->iCatalogs);
        cmd->iCatalogs = NULL;
        gd->iCatalogsFlag = FALSE;
    }
    
    if (gd->infilefmtFlag == TRUE && cmd->infilefmt != NULL) {
        free(cmd->infilefmt);
        cmd->infilefmt = NULL;
        gd->infilefmtFlag = FALSE;
    }
    
    if (gd->infileFlag == TRUE && cmd->infile != NULL) {
        free(cmd->infile);
        cmd->infile = NULL;
        gd->infileFlag = FALSE;
    }
    
    if (gd->rsmoothFlagFree == TRUE && cmd->rsmooth != NULL) {
        free(cmd->rsmooth);
        cmd->rsmooth = NULL;
        gd->rsmoothFlagFree = FALSE;
    }
    
    if (gd->searchMethodFlag == TRUE && cmd->searchMethod != NULL) {
        free(cmd->searchMethod);
        cmd->searchMethod = NULL;
        gd->searchMethodFlag = FALSE;
    }
    
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
