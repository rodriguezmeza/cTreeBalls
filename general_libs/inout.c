/*==============================================================================
	MODULE: inout.c					[general_libs]
	Written by: M.A. Rodriguez-Meza
	Starting date:	January, 2005
	Purpose: Utilities for input and output data
	Language: C
	Use:
	Major revisions:
==============================================================================*/
//        1          2          3          4        ^ 5          6          7

// STATIC problem: gcc version 11
#include "globaldefs.h"

#include "stdinc.h"
#include "vectmath.h"
#include "inout.h"
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include "input_contracts.h"


// ------------[	inout normal definitions	 	]------------

void in_int(stream str, int *iptr)
{
    if (fscanf(str, "%d", iptr) != 1)
        error("in_int: input conversion error\n");
}

void in_int_long(stream str, INTEGER *iptr)
{
    if (fscanf(str, "%ld", iptr) != 1)
        error("in_int: input conversion error\n");
}

void in_short(stream str, short *iptr)
{
    int tmp;

    if (fscanf(str, "%d", &tmp) != 1) {
        error("in_short: input conversion error\n");
    }
    *iptr = tmp;
}

//B not finished yet...
void in_bool(stream str, bool *iptr)
{
    char tmp[10];

    if (fscanf(str, "%9s", tmp) != 1) {
        error("in_bool: input conversion error\n");
    }
//    *iptr = tmp;
    if (strchr("tTyY1", *tmp) != NULL) {
//        *iptr=TRUE;
        *iptr=1;
//        printf("\tin_bool: found: %d", *iptr);
    } else
        if (strchr("fFnN0", *tmp) != NULL)  {
//            *iptr=FALSE;
            *iptr=0;
//            printf("\tin_bool: found: %d", *iptr);
        } else {
            error("in_bool: %s not bool\n",*tmp);
        }

}
//E

void in_real(stream str, real *rptr)
{
    double tmp;

    if (fscanf(str, "%lf", &tmp) != 1)
        error("in_real: input conversion error\n");
    *rptr = tmp;
}

void in_real_double(stream str, double *rptr)
{
    in_real(str, rptr);
}

#if defined(THREEDIM)
void in_vector(stream str, vector vec)
{
    double tmpx, tmpy, tmpz;

    if (fscanf(str, "%lf%lf%lf", &tmpx, &tmpy, &tmpz) != 3)
        error("in_vector: input conversion error\n");
    vec[0] = tmpx;
    vec[1] = tmpy;
    vec[2] = tmpz;
}
#endif

#if defined(TWODIM)
void in_vector(stream str, vector vec)
{
    double tmpx, tmpy, tmpz;

    if (fscanf(str, "%lf%lf", &tmpx, &tmpy) != 2)
        error("in_vector: input conversion error\n");
    vec[0] = tmpx;
    vec[1] = tmpy;
}
#endif

#if defined(ONEDIM)
void in_vector(stream str, vector vec)
{
    double tmpx, tmpy, tmpz;

    if (fscanf(str, "%lf", &tmpx) != 1)
        error("in_vector: input conversion error\n");
    vec[0] = tmpx;
}
#endif

#define IFMT  " %d"
#define IFMTL  " %ld"
#define RFMT  " %18.10E"
#define RFMTLG  " %lg"

void out_int(stream str, int ival)
{
    if (fprintf(str, IFMT "\n", ival) < 0)
        error("out_int: fprintf failed\n");
}

void out_int_long(stream str, INTEGER ival)
{
    if (fprintf(str, IFMTL "\n", ival) < 0)
        error("out_int: fprintf failed\n");
}

void out_short(stream str, short ival)
{
    if (fprintf(str, IFMT "\n", ival) < 0)
        error("out_short: fprintf failed\n");
}

void out_real(stream str, real rval)
{
    if (fprintf(str, RFMT "\n", rval) < 0)
        error("out_real: fprintf failed\n");
}

void out_vector(stream str, vector vec)
{
#if defined(THREEDIM)
    if (fprintf(str, RFMT RFMT RFMT "\n", vec[0], vec[1], vec[2]) < 0)
        error("out_vector: fprintf failed\n");
#endif
#if defined(TWODIM)
    if (fprintf(str, RFMT RFMT "\n", vec[0], vec[1]) < 0)
        error("out_vector: fprintf failed\n");
#endif
#if defined(ONEDIM)
    if (fprintf(str, RFMT "\n", vec[0]) < 0)
        error("out_vector: fprintf failed\n");
#endif
}


// ------------[	inout mar definitions	 		]------------

void out_bool_mar(stream str, bool bval)
{
    if (fprintf(str," %s", bval ? "T" : "F" ) < 0)
        error("out_bool_mar: fprintf failed\n");
}

void out_int_mar(stream str, int ival)
{
    if (fprintf(str, IFMT, ival) < 0)
        error("out_int_mar: fprintf failed\n");
}

void out_short_mar(stream str, short ival)
{
    if (fprintf(str, " %1d", ival) < 0)
        error("out_short_mar: fprintf failed\n");
}

void out_real_mar(stream str, real rval)
{
    if (fprintf(str, RFMT , rval) < 0)
        error("out_real_mar: fprintf failed\n");
}

void out_vector_mar(stream str, vector vec)
{
#if defined(THREEDIM)
    if (fprintf(str, RFMT RFMT RFMT , vec[0], vec[1], vec[2]) < 0)
        error("out_vector_mar: fprintf failed\n");
#endif
#if defined(TWODIM)
    if (fprintf(str, RFMT RFMT , vec[0], vec[1]) < 0)
        error("out_vector_mar: fprintf failed\n");
#endif
#if defined(ONEDIM)
    if (fprintf(str, RFMT , vec[0]) < 0)
        error("out_vector_mar: fprintf failed\n");
#endif
}


// ------------[	inout binary definitions	 	]------------

void in_int_bin(stream str, int *iptr)
{
    if (fread((void *) iptr, sizeof(int), 1, str) != 1)
        error("in_int_bin: fread failed\n");
}

void in_int_bin_long(stream str, INTEGER *iptr)
{
    if (fread((void *) iptr, sizeof(INTEGER), 1, str) != 1)
        error("in_int_bin_long: fread failed\n");
}

void in_short_bin(stream str, short *iptr)
{
    if (fread((void *) iptr, sizeof(short), 1, str) != 1)
        error("in_short_bin: fread failed\n");
}

void in_real_bin(stream str, real *rptr)
{
    if (fread((void *) rptr, sizeof(real), 1, str) != 1)
        error("in_real_bin: fread failed\n");
}

void in_real_bin_double(stream str, double *rptr)
{
    in_real_bin(str, rptr);
}

void in_vector_bin(stream str, vector vec)
{
    real disk_vector[NDIM];
    int k;

    if (fread((void *) disk_vector, sizeof(real), NDIM, str) != NDIM)
        error("in_vector_bin: fread failed\n");
    DO_COORD(k)
        vec[k] = (cballs_storage_real)disk_vector[k];
}

void out_int_bin(stream str, int ival)
{
    if (fwrite((void *) &ival, sizeof(int), 1, str) != 1)
        error("out_int_bin: fwrite failed\n");
}

void out_int_bin_long(stream str, INTEGER ival)
{
    if (fwrite((void *) &ival, sizeof(INTEGER), 1, str) != 1)
        error("out_int_bin_long: fwrite failed\n");
}

void out_short_bin(stream str, short ival)
{
    if (fwrite((void *) &ival, sizeof(short), 1, str) != 1)
        error("out_short_bin: fwrite failed\n");
}

void out_real_bin(stream str, real rval)
{
    if (fwrite((void *) &rval, sizeof(real), 1, str) != 1)
        error("out_real_bin: fwrite failed\n");
}

void out_vector_bin(stream str, vector vec)
{
    real disk_vector[NDIM];
    int k;

    DO_COORD(k)
        disk_vector[k] = (real)vec[k];
    if (fwrite((void *) disk_vector, sizeof(real), NDIM, str) != NDIM)
        error("out_vector_bin: fwrite failed\n");
}

void out_bool_mar_bin(stream str, bool bval)
{
    if (fwrite((void *) &bval, sizeof(bool), 1, str) != 1)
        error("out_bool_mar_bin: fwrite failed\n");
}

void out_short_mar_bin(stream str, short bval)
{
    if (fwrite((void *) &bval, sizeof(short), 1, str) != 1)
        error("out_short_mar_bin: fwrite failed\n");
}


// ------------[	inout other definitions	 		]------------

void in_vector_ruben(stream str, vector vec)
{
    double tmpx, tmpy, tmpz;

    fscanf(str, "%lf%lf%lf%*s", &tmpx, &tmpy, &tmpz);
    vec[0] = tmpx;
    vec[1] = tmpy;
    vec[2] = tmpz;
}


void out_vector_ndim(stream str, double *vec, int ndim)
{
    int i;
    
    for (i=0; i<ndim-1; i++)
        if (fprintf(str, RFMT, vec[i]) < 0)
            error("out_vector_ndim: fprintf failed\n");
    if (fprintf(str, RFMT "\n", vec[ndim-1]) < 0)
            error("out_vector_ndim: fprintf failed\n");
}

void in_vector_ndim(stream str, double *vec, int ndim)
{
    int i, c;
	real tmp;

    for (i=0; i<ndim; i++) {
        if (fscanf(str, "%lg", &tmp) != 1)
            error("in_vector_ndim: fscanf failed %d\n",i);
		vec[i] = tmp;
	}
	while ((c = getc(str)) != EOF)		// Reading rest of the line ...
		if (c=='\n') break;
}

size_t gdgt_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;

  if((nread=fread(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("my_fread: I/O error (fread) has occured.\n");
      fflush(stdout);
      endrun(778);
    }
  return nread;
}

size_t gdgt_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("my_fwrite: I/O error (fwrite) on has occured.\n");
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


//B added by cBalls
int out_vector_checked(stream str, vector vec, string routineName,
                       string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

#if defined(THREEDIM)
    rc = fprintf(str, RFMT RFMT RFMT "\n", vec[0], vec[1], vec[2]);
#elif defined(TWODIM)
    rc = fprintf(str, RFMT RFMT "\n", vec[0], vec[1]);
#elif defined(ONEDIM)
    rc = fprintf(str, RFMT "\n", vec[0]);
#else
    rc = -1;
    errno = EINVAL;
#endif

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing vector to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}

int out_real_mar_checked(stream str, real rval, string routineName,
                  string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

    rc = fprintf(str, RFMT , rval);

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing real to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}

int out_real_checked(stream str, real rval, string routineName,
                     string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

    rc = fprintf(str, RFMT "\n", rval);

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing real to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}

int out_short_mar_checked(stream str, short ival, string routineName,
                  string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

    rc = fprintf(str, " %1d", ival);

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing short to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}

int out_bool_mar_checked(stream str, bool bval, string routineName,
                         string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

    rc = fprintf(str," %s", bval ? "T" : "F" );

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing bool to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;

}

int out_int_bin_checked(stream str, int ival, string routineName,
                         string filename, char *errmsg, size_t errmsg_size)
{
    size_t nwritten = fwrite((void *) &ival, sizeof(int), 1, str);

    if (nwritten != 1) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing int to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, 1,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}

int out_short_bin_checked(stream str, short ival, string routineName,
                        string filename, char *errmsg, size_t errmsg_size)
{
    size_t nwritten = fwrite((void *) &ival, sizeof(short), 1, str);

    if (nwritten != 1) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing short to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, 1,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}

int out_bool_bin_checked(stream str, bool ival, string routineName,
                        string filename, char *errmsg, size_t errmsg_size)
{
    size_t nwritten = fwrite((void *) &ival, sizeof(bool), 1, str);

    if (nwritten != 1) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing bool to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, 1,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}


int out_int_bin_long_checked(stream str, INTEGER ival, string routineName,
                      string filename, char *errmsg, size_t errmsg_size)
{
    size_t nwritten = fwrite((void *) &ival, sizeof(INTEGER), 1, str);

    if (nwritten != 1) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing INTEGER to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, 1,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}

int out_real_bin_checked(stream str, real rval, string routineName,
                          string filename, char *errmsg, size_t errmsg_size)
{
    size_t nwritten = fwrite((void *) &rval, sizeof(real), 1, str);

    if (nwritten != 1) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing real to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, 1,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}

int out_vector_mar_checked(stream str, vector vec, string routineName,
                           string filename, char *errmsg, size_t errmsg_size)
{
    int rc;

#if defined(THREEDIM)
    rc = fprintf(str, RFMT RFMT RFMT, vec[0], vec[1], vec[2]);
#elif defined(TWODIM)
    rc = fprintf(str, RFMT RFMT, vec[0], vec[1]);
#elif defined(ONEDIM)
    rc = fprintf(str, RFMT, vec[0]);
#else
    rc = -1;
    errno = EINVAL;
#endif

    if (rc < 0) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing vector to '%s': %s",
                 routineName, filename,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }

    return SUCCESS;
}

int out_vector_bin_checked(stream str, vector vec, string routineName,
                            string filename, char *errmsg, size_t errmsg_size)
{
    real disk_vector[NDIM];
    int k;
    size_t nwritten;

    DO_COORD(k)
        disk_vector[k] = (real)vec[k];
    nwritten = fwrite((void *) disk_vector, sizeof(real), NDIM, str);

    if (nwritten != NDIM) {
        int saved_errno = errno;
        snprintf(errmsg, errmsg_size,
                 "%s: error writing vector to '%s': wrote %zu of %d values: %s",
                 routineName, filename, nwritten, NDIM,
                 strerror(saved_errno ? saved_errno : EIO));
        return FAILURE;
    }
    
    return SUCCESS;
}
//E added by cBalls


// IN/OUT ROUTINES TO HANDLE STRINGS ...

void ReadInString(stream instr, char *path)
{
	int i, c;
	char word[200];

	i=0;
	while ((c = getc(instr)) != EOF) {
		if (c==' ' || c=='\n' || c=='\t') break;
		word[i++]=c;
	}
    word[i] = '\0';
	strcpy(path, word);
}

void ReadInLineString(stream instr, char *text)
{
	int i, c;
	char word[200];

	i=0;
	while ((c = getc(instr)) != EOF) {
		if (c=='\n') break;
		word[i++]=c;
	}
    word[i] = '\0';
	strcpy(text, word);
}


// InputData's (como el que está en nplot2d)
//B PARA IMPLEMENTAR LECTURA GENERAL DE ARCHIVOS DE DATOS
//  CON FORMATO DE COLUMNAS

#define IN 1
#define OUT 0
#define SI 1
#define NO 0

void inout_InputData(string filename, int col1, int col2, int *npts)
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    
    instr = stropen(filename, "r");
    
    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    
    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    
    rewind(instr);
    
    npoint=nrow;
    row = (realptr) allocate(ncol*sizeof(real));
    
    *npts = npoint;
    inout_xval = (real *) allocate(npoint * sizeof(real));
    inout_yval = (real *) allocate(npoint * sizeof(real));
    
    ip = 0;
    for (i=0; i<nl; i++) {
        if (lineQ[i]) {
            in_vector_ndim(instr, row, ncol);
            inout_xval[ip] = row[col1-1];
            inout_yval[ip] = row[col2-1];
            ++ip;
        } else {
            while ((c = getc(instr)) != EOF)		// Reading dummy line ...
                if (c=='\n') break;
        }
    }
    
    fclose(instr);
}

void InputData_check_file(string filename)
{
    stream instr;
    int ncol, nrow;
    int c, nl, nw, nc, state, salto, nwxc, i;
    short int *lineQ;
    
    instr = stropen(filename, "r");

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    printf("\n\tFile info: number of lines, words, and characters : %d %d %d\n",
           nl, nw, nc);
    
    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    printf("\tFile info: nrow, ncol, nvalues : %d %d %d\n", nrow, ncol, nw);
    
    fclose(instr);
}

/* Read selected numeric columns without legacy fatal conversion helpers.
 * getline preserves final lines without a newline and avoids fixed row limits.
 * Publish only fully read arrays; a failed read leaves runtime buffers alone. */
static int input_columns_checked(struct cmdline_data *cmd, string filename,
                                  int nselected, const int *columns, int *npts)
{
    FILE *instr = NULL;
    char *line = NULL;
    size_t line_capacity = 0, rows = 0, capacity = 0, physical_row = 0;
    real *values[5] = {NULL}, *grown = NULL;
    int expected_columns = -1, rc = FAILURE;
    ssize_t length;
    *npts = 0;
    for (int j = 0; j < nselected; ++j) {
        if (columns[j] < 1) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "column input '%s': column numbers must be positive", filename);
            goto done;
        }
    }
    if (stropen_checked(filename, "r", &instr, cmd->error_message,
                         _ERRORMSGSIZE_) == FAILURE)
        goto done;
    while ((length = getline(&line, &line_capacity, instr)) >= 0) {
        char *cursor = line, *end;
        real selected[5] = {0};
        int ncolumns = 0;
        ++physical_row;
        if (memchr(line, '\0', (size_t)length) != NULL) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "column input '%s': row %zu contains an embedded NUL", filename, physical_row);
            goto done;
        }
        while (1) {
            double parsed;
            while (isspace((unsigned char)*cursor)) ++cursor;
            if (!*cursor || *cursor == '#' || *cursor == '%') break;
            errno = 0;
            parsed = strtod(cursor, &end);
            if (end == cursor || (*end && !isspace((unsigned char)*end)
                                  && *end != '#' && *end != '%')) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "column input '%s': row %zu, column %d: invalid real number",
                         filename, physical_row, ncolumns + 1);
                goto done;
            }
            ++ncolumns;
            for (int j = 0; j < nselected; ++j) {
                if (columns[j] == ncolumns) {
                    char field[48];
                    snprintf(field, sizeof(field), "column %d", ncolumns);
                    if (cballs_input_finite(cmd, filename, physical_row, field, parsed) == FAILURE)
                        goto done;
                    if (errno == ERANGE) {
                        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                                 "column input '%s': row %zu, %s: out-of-range value",
                                 filename, physical_row, field);
                        goto done;
                    }
                    selected[j] = (real)parsed;
                    if (cballs_input_finite(cmd, filename, physical_row, field, selected[j]) == FAILURE)
                        goto done;
                }
            }
            cursor = end;
        }
        if (ncolumns == 0) continue;
        if (expected_columns < 0) expected_columns = ncolumns;
        if (ncolumns != expected_columns) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "column input '%s': row %zu: expected %d columns, found %d",
                     filename, physical_row, expected_columns, ncolumns);
            goto done;
        }
        for (int j = 0; j < nselected; ++j)
            if (columns[j] > ncolumns) {
                snprintf(cmd->error_message, _ERRORMSGSIZE_,
                         "column input '%s': row %zu: requested column %d exceeds %d columns",
                         filename, physical_row, columns[j], ncolumns);
                goto done;
            }
        if (rows == INT_MAX) {
            snprintf(cmd->error_message, _ERRORMSGSIZE_,
                     "column input '%s': row count exceeds INT_MAX", filename);
            goto done;
        }
        if (rows == capacity) {
            size_t next = capacity ? (capacity > INT_MAX / 2 ? INT_MAX : 2 * capacity) : 256;
            for (int j = 0; j < nselected; ++j) {
                if (cballs_malloc_checked((void **)&grown, next, sizeof(real),
                                           "ASCII column buffer", cmd->error_message,
                                           _ERRORMSGSIZE_) == FAILURE)
                    goto done;
                if (rows) memcpy(grown, values[j], rows * sizeof(real));
                free(values[j]); values[j] = grown; grown = NULL;
            }
            capacity = next;
        }
        for (int j = 0; j < nselected; ++j) values[j][rows] = selected[j];
        ++rows;
    }
    if (!feof(instr) || ferror(instr) || rows == 0) {
        snprintf(cmd->error_message, _ERRORMSGSIZE_,
                 "column input '%s': %s", filename,
                 rows == 0 ? "no data rows" : "read error");
        goto done;
    }
    if (fclose(instr) != 0) {
        instr = NULL;
        snprintf(cmd->error_message, _ERRORMSGSIZE_, "column input '%s': close error", filename);
        goto done;
    }
    instr = NULL;
    for (int j = 0; j < nselected; ++j) {
        cballs_runtime_current()->io_columns[j] = values[j];
        values[j] = NULL;
    }
    *npts = (int)rows;
    rc = SUCCESS;
done:
    if (instr) fclose(instr);
    free(line); free(grown);
    for (int j = 0; j < nselected; ++j) free(values[j]);
    return rc;
}

int InputData_2c(struct cmdline_data *cmd, string filename, int c1, int c2, int *npts)
{
    const int columns[] = {c1, c2};
    return input_columns_checked(cmd, filename, 2, columns, npts);
}
int InputData_3c(struct cmdline_data *cmd, string filename, int c1, int c2, int c3, int *npts)
{
    const int columns[] = {c1, c2, c3};
    return input_columns_checked(cmd, filename, 3, columns, npts);
}
int InputData_4c(struct cmdline_data *cmd, string filename, int c1, int c2, int c3, int c4, int *npts)
{
    const int columns[] = {c1, c2, c3, c4};
    return input_columns_checked(cmd, filename, 4, columns, npts);
}
int InputData_5c(struct cmdline_data *cmd, string filename, int c1, int c2, int c3, int c4, int c5, int *npts)
{
    const int columns[] = {c1, c2, c3, c4, c5};
    return input_columns_checked(cmd, filename, 5, columns, npts);
}

int inout_InputDataVector(
                          string filename, real *vec, int *npts,
                          short verbose, short verbose_log, FILE *outlog
                          )
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    int col1=1;
    
    instr = stropen(filename, "r");

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,
               "\nReading matrix elements from file %s... ",filename);

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,
                   "\n\tGeneral statistics : ");
        verb_log_print(verbose_log, outlog,
                   "number of lines, words, and characters : %d %d %d", nl, nw, nc);
    }

    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,"\n\tValid numbers statistics : ");
        verb_log_print(verbose_log, outlog,
                       "nrow, ncol, nvalues : %d %d %d", nrow, ncol, nw);
    }

    rewind(instr);
    
    npoint=nrow;
    row = (realptr) allocate(ncol*sizeof(real));

    *npts = npoint;
    
    ip = 1;
    for (i=0; i<nl; i++) {
        if (lineQ[i]) {
            in_vector_ndim(instr, row, ncol);
            vec[ip] = row[col1-1];
            ++ip;
        } else {
            while ((c = getc(instr)) != EOF)        // Reading dummy line ...
                if (c=='\n') break;
        }
    }
    fclose(instr);

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,"\n\t... done.\n");

    return SUCCESS;
}


int inout_InputDataMatrix_info(
                               string filename, int *nrow_out, int *ncol_out,
                               short verbose, short verbose_log, FILE *outlog
                               )
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    int col1=1;
    int col2=2;
    
    instr = stropen(filename, "r");

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,
               "\nReading matrix elements from file %s... ",filename);

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,
                   "\n\tGeneral statistics : ");
        verb_log_print(verbose_log, outlog,
                "number of lines, words, and characters : %d %d %d",
                nl, nw, nc);
    }

    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,"\n\tValid numbers statistics : ");
        verb_log_print(verbose_log, outlog,
                       "nrow, ncol, nvalues : %d %d %d", nrow, ncol, nw);
    }

    rewind(instr);

    if (nrow != ncol) {
        verb_print(verbose_log,
                   "\ninout_InputDataMatrix: warning!!! in input file: %s %d %d",
                   "nrow and ncol are different:", nrow, ncol);
    }

    *nrow_out = nrow;
    *ncol_out = ncol;

    fclose(instr);

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,"\n\t... done.\n");

    return SUCCESS;
}

int inout_InputDataMatrix(
                          string filename, real **mat, int *npts,
                          short verbose, short verbose_log, FILE *outlog
                          )
{
    stream instr;
    int ncol, nrow;
    real *row;
    int c, nl, nw, nc, state, salto, nwxc, i, npoint, ip;
    short int *lineQ;
    int col1=1;
    int col2=2;
    
    instr = stropen(filename, "r");

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,
               "\nReading matrix elements from file %s... ",filename);

    state = OUT;
    nl = nw = nc = 0;
    while ((c = getc(instr)) != EOF) {
        ++nc;
        if (c=='\n')
            ++nl;
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else if (state == OUT) {
            state = IN;
            ++nw;
        }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,
                   "\n\tGeneral statistics : ");
        verb_log_print(verbose_log, outlog,
                   "number of lines, words, and characters : %d %d %d", nl, nw, nc);
    }

    rewind(instr);
    
    lineQ = (short int *) allocate(nl * sizeof(short int));
    for (i=0; i<nl; i++) lineQ[i]=FALSE;
    
    nw = nrow = ncol = nwxc = 0;
    state = OUT;
    salto = NO;
    
    i=0;
    
    while ((c = getc(instr)) != EOF) {
        
        if(c=='%' || c=='#') {
            while ((c = getc(instr)) != EOF)
                if (c=='\n') break;
            ++i;
            continue;
        }
        
        if (c=='\n' && nw > 0) {
            if (salto==NO) {
                ++nrow;
                salto=SI;
                if (ncol != nwxc && nrow>1) {
                    printf("\nvalores diferentes : ");
                    error("(nrow, ncol before, ncol after) : %d %d %d\n\n",
                          nrow, ncol, nwxc);
                }
                ncol = nwxc;
                lineQ[i]=TRUE;
                ++i;
                nwxc=0;
            } else {
                ++i;
            }
        }
        
        if (c==' ' || c=='\n' || c=='\t')
            state = OUT;
        else
            if (state == OUT) {
                state = IN;
                ++nw; ++nwxc;
                salto=NO;
            }
    }
    if (verbose_log>=3) {
        verb_log_print(verbose_log, outlog,"\n\tValid numbers statistics : ");
        verb_log_print(verbose_log, outlog,
                       "nrow, ncol, nvalues : %d %d %d", nrow, ncol, nw);
    }

    rewind(instr);

    if (nrow != ncol) {
        error("\ninout_InputDataMatrix: in input file: %s %d %d",
              "nrow and ncol are different:", nrow, ncol);
    }

    npoint=nrow;
    row = (realptr) allocate(ncol*sizeof(real));

    *npts = npoint;
    
    ip = 1;
    for (i=0; i<nl; i++) {
        if (lineQ[i]) {
            in_vector_ndim(instr, row, ncol);
            for (int j=0; j<ncol; j++) {
                mat[ip][j+1] = row[j];
            }
            ++ip;
        } else {
            while ((c = getc(instr)) != EOF)        // Reading dummy line ...
                if (c=='\n') break;
        }
    }
    fclose(instr);

    if (verbose_log>=3)
        verb_log_print(verbose_log, outlog,"\n\t... done.\n");

    return SUCCESS;
}

//E additions


#undef IN
#undef OUT
#undef SI
#undef NO

//E PARA IMPLEMENTAR LECTURA GENERAL DE ARCHIVOS DE DATOS CON FORMATO DE COLUMNAS

// Extract input rootDirPath and preFileName
int extractInputRootDir(const char *infilenames,
                        char *rootDirPath, size_t rootDirPathSize,
                        char *preFileName, size_t preFileNameSize, int ifile,
                        short verbose, short verbose_log, FILE *outlog
                        )
{
    const char *slash;
    const char *prefix;
    size_t root_length;
    size_t prefix_length;

    if (infilenames == NULL || rootDirPath == NULL || preFileName == NULL
        || rootDirPathSize == 0 || preFileNameSize == 0)
        return FAILURE;

    slash = strrchr(infilenames, '/');
    if (slash == NULL) {
        prefix = infilenames;
        root_length = 1;
        if (rootDirPathSize < 2)
            return FAILURE;
        rootDirPath[0] = '.';
        rootDirPath[1] = '\0';
    } else {
        prefix = slash + 1;
        root_length = (size_t)(slash - infilenames);
        if (root_length + 1 > rootDirPathSize)
            return FAILURE;
        memcpy(rootDirPath, infilenames, root_length);
        rootDirPath[root_length] = '\0';
    }

    prefix_length = strlen(prefix);
    if (prefix_length + 1 > preFileNameSize)
        return FAILURE;
    memcpy(preFileName, prefix, prefix_length + 1);

    verb_print_q(2, verbose,
                 "extractInputRootDir: file %d root '%s', prefix '%s'\n",
                 ifile, rootDirPath, preFileName);
    return SUCCESS;
}


#ifdef ADDONS
//#include "inout_02.h"
#endif

void error_open_file_kd(char *fname)
{
// Open error handler
  fprintf(stderr,"error_open_file: Couldn't open file %s \n",fname);
  exit(1);
}

