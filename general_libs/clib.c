//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

#include <math.h>
#include <time.h>
#include <sys/time.h>                               // To get time of the day

#include "stdinc.h"
#ifdef GETPARAM
#include "getparam.h"
#endif

                                                    // To remove:
#include <unistd.h>                                 // warning: implicit
                                                    //  declaration of function
                                                    //  ‘dup’
                                                    //  [-Wimplicit-function
                                                    //  -declaration]
                                                    //  fds = dup(fileno(inflag ?
                                                    //  stdin : stdout));
#include <sys/timeb.h>
#include <string.h>	

#include <sys/stat.h>                               // it has also mkdir()
#include <stdarg.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <setjmp.h>
#include <stdatomic.h>


//B Memory section

typedef struct cballs_allocation_context {
    jmp_buf recovery;
    char *errmsg;
    size_t errmsg_size;
    struct cballs_allocation_context *previous;
} cballs_allocation_context;

static _Thread_local cballs_allocation_context *active_allocation_context = NULL;
static _Atomic long allocation_failure_after = -1;

static int cballs_should_fail_allocation(void)
{
    long remaining = atomic_load(&allocation_failure_after);

    while (remaining >= 0) {
        long replacement = remaining == 0 ? -1 : remaining - 1;
        if (atomic_compare_exchange_weak(&allocation_failure_after,
                                         &remaining, replacement))
            return remaining == 0;
    }

    return FALSE;
}

static void *cballs_raw_malloc(size_t size)
{
    if (cballs_should_fail_allocation())
        return NULL;
    return malloc(size == 0 ? 1 : size);
}

static void *cballs_raw_calloc(size_t count, size_t item_size)
{
    if (cballs_should_fail_allocation())
        return NULL;
    if (count == 0 || item_size == 0)
        return calloc(1, 1);
    return calloc(count, item_size);
}

void cballs_allocation_failure(size_t bytes, const char *label)
{
    cballs_allocation_context *context = active_allocation_context;

    if (context != NULL) {
        if (context->errmsg != NULL && context->errmsg_size > 0) {
            snprintf(context->errmsg, context->errmsg_size,
                     "memory allocation failed for %s (%zu bytes)",
                     label != NULL ? label : "allocation", bytes);
        }
        longjmp(context->recovery, 1);
    }

#ifdef GETPARAM
    fprintf(stderr, "allocate in %s: not enough memory (%zu bytes)\n",
            getargv0(), bytes);
#else
    fprintf(stderr, "allocate: not enough memory (%zu bytes)\n", bytes);
#endif
    fflush(stderr);
    exit(1);
}

static int cballs_allocation_bytes(size_t count, size_t item_size,
                                   size_t *bytes,
                                   const char *label,
                                   char *errmsg, size_t errmsg_size)
{
    if (bytes == NULL || (item_size != 0 && count > SIZE_MAX / item_size)) {
        if (errmsg != NULL && errmsg_size > 0) {
            snprintf(errmsg, errmsg_size,
                     "memory allocation size overflow for %s (%zu x %zu)",
                     label != NULL ? label : "allocation", count, item_size);
        }
        return FAILURE;
    }

    *bytes = count * item_size;
    return SUCCESS;
}

int cballs_allocation_guard(cballs_allocation_callback callback,
                            void *argument,
                            char *errmsg, size_t errmsg_size)
{
    cballs_allocation_context context;
    int status;

    if (callback == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "invalid null callback for guarded allocation call");
        return FAILURE;
    }

    context.errmsg = errmsg;
    context.errmsg_size = errmsg_size;
    context.previous = active_allocation_context;
    active_allocation_context = &context;

    if (setjmp(context.recovery) != 0) {
        active_allocation_context = context.previous;
        return FAILURE;
    }

    status = callback(argument);
    active_allocation_context = context.previous;
    return status;
}

int cballs_malloc_checked(void **pointer, size_t count, size_t item_size,
                          const char *label,
                          char *errmsg, size_t errmsg_size)
{
    size_t bytes;

    if (pointer == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "missing destination for allocation of %s",
                     label != NULL ? label : "allocation");
        return FAILURE;
    }
    *pointer = NULL;

    if (cballs_allocation_bytes(count, item_size, &bytes, label,
                                errmsg, errmsg_size) == FAILURE)
        return FAILURE;

    *pointer = cballs_raw_malloc(bytes);
    if (*pointer == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "memory allocation failed for %s (%zu bytes)",
                     label != NULL ? label : "allocation", bytes);
        return FAILURE;
    }

    return SUCCESS;
}

int cballs_calloc_checked(void **pointer, size_t count, size_t item_size,
                          const char *label,
                          char *errmsg, size_t errmsg_size)
{
    size_t bytes;

    if (pointer == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "missing destination for allocation of %s",
                     label != NULL ? label : "allocation");
        return FAILURE;
    }
    *pointer = NULL;

    if (cballs_allocation_bytes(count, item_size, &bytes, label,
                                errmsg, errmsg_size) == FAILURE)
        return FAILURE;

    *pointer = cballs_raw_calloc(count, item_size);
    if (*pointer == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "memory allocation failed for %s (%zu bytes)",
                     label != NULL ? label : "allocation", bytes);
        return FAILURE;
    }

    return SUCCESS;
}

void cballs_test_fail_allocation_after(long successful_allocations)
{
    atomic_store(&allocation_failure_after,
                 successful_allocations < 0 ? -1 : successful_allocations);
}

void cballs_test_reset_allocation_failure(void)
{
    atomic_store(&allocation_failure_after, -1);
}

local int numeric_parse_error(char *errmsg, size_t errmsg_size,
                              const char *label, const char *text,
                              const char *reason)
{
    if (errmsg != NULL && errmsg_size > 0) {
        snprintf(errmsg, errmsg_size, "invalid numeric value for %s: '%s' (%s)",
                 label != NULL ? label : "value",
                 text != NULL ? text : "(null)", reason);
    }
    return FAILURE;
}

int parse_long_checked(const char *text, long *value,
                       char *errmsg, size_t errmsg_size,
                       const char *label)
{
    char *end;
    long parsed;

    if (text == NULL || value == NULL)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "missing input or destination");

    errno = 0;
    parsed = strtol(text, &end, 10);
    if (end == text)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "expected an integer");
    while (isspace((unsigned char)*end))
        end++;
    if (*end != '\0')
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "trailing characters");
    if (errno == ERANGE)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "out of range");

    *value = parsed;
    return SUCCESS;
}

int parse_int_checked(const char *text, int *value,
                      char *errmsg, size_t errmsg_size,
                      const char *label)
{
    long parsed;

    if (value == NULL)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "missing destination");
    if (parse_long_checked(text, &parsed, errmsg, errmsg_size, label) == FAILURE)
        return FAILURE;
    if (parsed < INT_MIN || parsed > INT_MAX)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "out of int range");

    *value = (int)parsed;
    return SUCCESS;
}

int parse_double_checked(const char *text, double *value,
                         char *errmsg, size_t errmsg_size,
                         const char *label)
{
    char *end;
    double parsed;

    if (text == NULL || value == NULL)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "missing input or destination");

    errno = 0;
    parsed = strtod(text, &end);
    if (end == text)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "expected a finite number");
    while (isspace((unsigned char)*end))
        end++;
    if (*end != '\0')
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "trailing characters");
    if (errno == ERANGE || !isfinite(parsed))
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "not finite or out of range");

    *value = parsed;
    return SUCCESS;
}

local int ascii_string_equal_ignore_case(const char *left, const char *right)
{
    unsigned char left_char;
    unsigned char right_char;

    if (left == NULL || right == NULL)
        return FALSE;

    while (*left != '\0' && *right != '\0') {
        left_char = (unsigned char)*left++;
        right_char = (unsigned char)*right++;
        if (left_char >= 'A' && left_char <= 'Z')
            left_char = (unsigned char)(left_char - 'A' + 'a');
        if (right_char >= 'A' && right_char <= 'Z')
            right_char = (unsigned char)(right_char - 'A' + 'a');
        if (left_char != right_char)
            return FALSE;
    }

    return *left == '\0' && *right == '\0';
}

int parse_bool_checked(const char *text, bool *value,
                       char *errmsg, size_t errmsg_size,
                       const char *label)
{
    if (text == NULL || value == NULL)
        return numeric_parse_error(errmsg, errmsg_size, label, text,
                                   "missing input or destination");

    if (strcmp(text, "1") == 0
        || ascii_string_equal_ignore_case(text, "true")
        || ascii_string_equal_ignore_case(text, "t")
        || ascii_string_equal_ignore_case(text, "yes")
        || ascii_string_equal_ignore_case(text, "y")) {
        *value = TRUE;
        return SUCCESS;
    }

    if (strcmp(text, "0") == 0
        || ascii_string_equal_ignore_case(text, "false")
        || ascii_string_equal_ignore_case(text, "f")
        || ascii_string_equal_ignore_case(text, "no")
        || ascii_string_equal_ignore_case(text, "n")) {
        *value = FALSE;
        return SUCCESS;
    }

    if (errmsg != NULL && errmsg_size > 0) {
        snprintf(errmsg, errmsg_size,
                 "invalid boolean value for %s: '%s' "
                 "(expected true/false, yes/no, or 1/0)",
                 label != NULL ? label : "value", text);
    }
    return FAILURE;
}

int read_parameter_line_checked(FILE *input, char *line, size_t line_size,
                                int *has_line, unsigned long *line_number,
                                const char *filename,
                                char *errmsg, size_t errmsg_size)
{
    int next;
    size_t length = 0;

    if (input == NULL || line == NULL || line_size < 2
        || has_line == NULL || line_number == NULL) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "invalid destination passed to parameter-file reader");
        return FAILURE;
    }

    *has_line = FALSE;
    next = fgetc(input);
    if (next == EOF) {
        if (ferror(input)) {
            if (errmsg != NULL && errmsg_size > 0)
                snprintf(errmsg, errmsg_size,
                         "could not read parameter file '%s'",
                         filename != NULL ? filename : "(unknown)");
            return FAILURE;
        }
        return SUCCESS;
    }

    if (*line_number == ULONG_MAX) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "parameter file '%s' has too many physical lines",
                     filename != NULL ? filename : "(unknown)");
        return FAILURE;
    }
    (*line_number)++;
    while (next != '\n' && next != EOF) {
        if (next == '\0') {
            do {
                next = fgetc(input);
            } while (next != '\n' && next != EOF);
            if (errmsg != NULL && errmsg_size > 0) {
                snprintf(errmsg, errmsg_size,
                         "%s:%lu: parameter line contains an embedded NUL",
                         filename != NULL ? filename : "parameter input",
                         *line_number);
            }
            return FAILURE;
        }
        if (length + 1 >= line_size) {
            do {
                next = fgetc(input);
            } while (next != '\n' && next != EOF);
            if (errmsg != NULL && errmsg_size > 0) {
                snprintf(errmsg, errmsg_size,
                         "%s:%lu: parameter line exceeds %zu characters",
                         filename != NULL ? filename : "parameter input",
                         *line_number, line_size - 1);
            }
            return FAILURE;
        }
        line[length++] = (char)next;
        next = fgetc(input);
    }

    if (next == EOF && ferror(input)) {
        if (errmsg != NULL && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "could not finish reading parameter file '%s'",
                     filename != NULL ? filename : "(unknown)");
        return FAILURE;
    }

    if (next == '\n' && length + 1 < line_size)
        line[length++] = '\n';
    line[length] = '\0';
    *has_line = TRUE;

    return SUCCESS;
}

local void parameter_syntax_error(char *errmsg, size_t errmsg_size,
                                  const char *filename,
                                  unsigned long line_number,
                                  const char *reason)
{
    if (errmsg == NULL || errmsg_size == 0)
        return;

    if (line_number > 0) {
        snprintf(errmsg, errmsg_size, "%s:%lu: %s",
                 filename != NULL ? filename : "parameter input",
                 line_number, reason);
    } else {
        snprintf(errmsg, errmsg_size, "%s: %s",
                 filename != NULL ? filename : "parameter input", reason);
    }
}

int parse_parameter_line_checked(const char *line,
                                 char *name, size_t name_size,
                                 char *value, size_t value_size,
                                 int *is_data,
                                 const char *filename,
                                 unsigned long line_number,
                                 char *errmsg, size_t errmsg_size)
{
    const char *cursor;
    const char *equal = NULL;
    const char *name_begin;
    const char *name_end;
    const char *value_begin;
    const char *value_end;
    const char *scan;
    char quote = '\0';
    int escaped = FALSE;
    size_t length;

    if (line == NULL || name == NULL || name_size == 0
        || value == NULL || value_size == 0 || is_data == NULL) {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "invalid parameter-line destination");
        return FAILURE;
    }

    *is_data = FALSE;
    name[0] = '\0';
    value[0] = '\0';

    cursor = line;
    while (isspace((unsigned char)*cursor))
        cursor++;
    if (*cursor == '\0' || *cursor == '#' || *cursor == '%')
        return SUCCESS;

    for (scan = cursor; *scan != '\0'; scan++) {
        if (quote != '\0') {
            if (escaped) {
                escaped = FALSE;
            } else if (*scan == '\\') {
                escaped = TRUE;
            } else if (*scan == quote) {
                quote = '\0';
            }
        } else if (*scan == '\'' || *scan == '"') {
            quote = *scan;
        } else if (*scan == '#' || *scan == '%') {
            return SUCCESS;
        } else if (*scan == '=') {
            equal = scan;
            break;
        }
    }

    if (equal == NULL) {
        if (quote != '\0') {
            parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                                   "unterminated quote before '='");
            return FAILURE;
        }
        return SUCCESS;
    }

    name_begin = cursor;
    name_end = equal;
    while (name_end > name_begin
           && isspace((unsigned char)name_end[-1]))
        name_end--;

    if (name_end > name_begin
        && (*name_begin == '\'' || *name_begin == '"')) {
        if (name_end - name_begin < 2 || name_end[-1] != *name_begin) {
            parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                                   "parameter name has unmatched quotes");
            return FAILURE;
        }
        name_begin++;
        name_end--;
    } else if (name_end > name_begin
               && (name_end[-1] == '\'' || name_end[-1] == '"')) {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "parameter name has unmatched quotes");
        return FAILURE;
    }

    if (name_end == name_begin) {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "no parameter name before '='");
        return FAILURE;
    }

    length = (size_t)(name_end - name_begin);
    if (length >= name_size) {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "parameter name is too long");
        return FAILURE;
    }
    memcpy(name, name_begin, length);
    name[length] = '\0';

    value_begin = equal + 1;
    while (isspace((unsigned char)*value_begin))
        value_begin++;

    value_end = line + strlen(line);
    quote = '\0';
    escaped = FALSE;
    for (scan = value_begin; scan < value_end; scan++) {
        if (quote != '\0') {
            if (escaped) {
                escaped = FALSE;
            } else if (*scan == '\\') {
                escaped = TRUE;
            } else if (*scan == quote) {
                quote = '\0';
            }
        } else if (*scan == '\'' || *scan == '"') {
            quote = *scan;
        } else if (*scan == '#' || *scan == '%') {
            value_end = scan;
            break;
        }
    }

    if (quote != '\0') {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "parameter value has an unterminated quote");
        return FAILURE;
    }

    while (value_end > value_begin
           && isspace((unsigned char)value_end[-1]))
        value_end--;
    if (value_end == value_begin)
        return SUCCESS;

    length = (size_t)(value_end - value_begin);
    if (length >= value_size) {
        parameter_syntax_error(errmsg, errmsg_size, filename, line_number,
                               "parameter value is too long");
        return FAILURE;
    }
    memcpy(value, value_begin, length);
    value[length] = '\0';
    *is_data = TRUE;

    return SUCCESS;
}

void *allocate_array(int nb)
{
    void *mem;

    if (nb < 0)
        cballs_allocation_failure(0, "allocate_array");
    mem = cballs_raw_calloc((size_t)nb, 1);
    if (mem == NULL)
        cballs_allocation_failure((size_t)nb, "allocate_array");
    return (mem);
}

void *allocate(long int nb)
{
    void *mem;

    if (nb < 0)
        cballs_allocation_failure(0, "allocate");
    mem = cballs_raw_calloc((size_t)nb, 1);
    if (mem == NULL)
        cballs_allocation_failure((size_t)nb, "allocate");
    return (mem);
}


void FreeVecINormal(int *v)
{
	free(v);
}

//B Following NR
#define NOFFSET_END 1
#define FREE_ARG char*

int *ivector(long nl, long nh)
{
    int *v;

    size_t bytes = (size_t)(nh-nl+1+NOFFSET_END)*sizeof(int);

    v=(int *)cballs_raw_malloc(bytes);
    if (!v) cballs_allocation_failure(bytes, "ivector");
    return v-nl+NOFFSET_END;
}

void free_ivector(int *v, long nl, long nh)
{
    if (v == NULL)
        return;
    free((FREE_ARG) (v+nl-NOFFSET_END));
}

double *dvector(long nl, long nh)
{
    double *v;

    size_t bytes = (size_t)(nh-nl+1+NOFFSET_END)*sizeof(double);

    v=(double *)cballs_raw_malloc(bytes);
    if (!v) cballs_allocation_failure(bytes, "dvector");

    return v-nl+NOFFSET_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,
            ncol=nch-ncl+1;
    double **m, **row_base;
    size_t pointer_bytes = (size_t)(nrow+NOFFSET_END)*sizeof(double*);
    size_t data_bytes = (size_t)(nrow*ncol+NOFFSET_END)*sizeof(double);

    m=(double **) cballs_raw_malloc(pointer_bytes);
    if (!m) cballs_allocation_failure(pointer_bytes, "dmatrix row pointers");
    row_base = m;
    m += NOFFSET_END;
    m -= nrl;

    m[nrl]=(double *) cballs_raw_malloc(data_bytes);
    if (!m[nrl]) {
        free(row_base);
        cballs_allocation_failure(data_bytes, "dmatrix data");
    }
    m[nrl] += NOFFSET_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}


double ***dmatrix3D(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    double ***t, ***plane_base;
    double **row_base;
    size_t plane_bytes = (size_t)(nrow+NOFFSET_END)*sizeof(double**);
    size_t row_bytes = (size_t)(nrow*ncol+NOFFSET_END)*sizeof(double*);
    size_t data_bytes =
        (size_t)(nrow*ncol*ndep+NOFFSET_END)*sizeof(double);

    /* allocate pointers to pointers to rows */
    t=(double ***) cballs_raw_malloc(plane_bytes);
    if (!t) cballs_allocation_failure(plane_bytes, "dmatrix3D plane pointers");
    plane_base = t;
    t += NOFFSET_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(double **) cballs_raw_malloc(row_bytes);
    if (!t[nrl]) {
        free(plane_base);
        cballs_allocation_failure(row_bytes, "dmatrix3D row pointers");
    }
    row_base = t[nrl];
    t[nrl] += NOFFSET_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(double *) cballs_raw_malloc(data_bytes);
    if (!t[nrl][ncl]) {
        free(row_base);
        free(plane_base);
        cballs_allocation_failure(data_bytes, "dmatrix3D data");
    }
    t[nrl][ncl] += NOFFSET_END;
    t[nrl][ncl] -= ndl;

    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

void free_dmatrix3D(double ***t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh)
/* free double f3tensor allocated by f3tensor() */
{
    if (t == NULL)
        return;
    free((FREE_ARG) (t[nrl][ncl]+ndl-NOFFSET_END));
    free((FREE_ARG) (t[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (t+nrl-NOFFSET_END));
}


void free_dvector(double *v, long nl, long nh)
{
    if (v == NULL)
        return;
    if (nl>0)
        free((FREE_ARG) (v+nl-NOFFSET_END));
    else
        free((FREE_ARG) (v));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    if (m == NULL)
        return;
    free((FREE_ARG) (m[nrl]+ncl-NOFFSET_END));
    free((FREE_ARG) (m+nrl-NOFFSET_END));
}


#undef NOFFSET_END
#undef FREE_ARG
//E

void error_mem_out_kd(void)
{
// Memory shortage handler
  fprintf(stderr,"error_mem_out_kd: Out of memory!!\n");
  exit(1);
}


//E Memory section


//B Time section

double cputime(void)
{
    time_t ltime;
    time(&ltime);
    return(((double)ltime)/((double) 60.0));
}

double cpu_time(void)
{
    double value;
      value = (double) clock () / (double) CLOCKS_PER_SEC;
      return value;
}

// Gives time of the day
#define REALTIMECLIB   (gettimeofday(&current_time, NULL))
#define GETREALTIME(time) \
    {REALTIMECLIB; \
    (time) = current_time.tv_sec;}

long rcpu_time(void)
{
    struct timeval current_time;
    long cpurealfinal;
    GETREALTIME(cpurealfinal);
    return cpurealfinal;
}

#undef REALTIMECLIB
#undef GETREALTIME

double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

}

double timediff(double t0, double t1)
{
  double dt;
  
  dt=t1-t0;

    if(dt<0) {                                      // overflow has occured (for
                                                    //  systems with 32bit tick
                                                    //  counter)
#ifdef WALLCLOCK
        dt = 0;
#else
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}

//  Timing variables
#ifdef OPENMPCODE
#include <omp.h>
static double relbeg_kd,relend_kd,absbeg_kd,absend_kd;
#else //OPENMPCODE
#include <time.h>
static time_t relbeg_kd,relend_kd,absbeg_kd,absend_kd;
#endif //OPENMPCODE

void timer_kd(int i)
{
// Timing routine
// timer(0) -> initialize relative clock
// timer(1) -> read relative clock
// timer(2) -> read relative clock and initialize it afterwards
// timer(4) -> initialize absolute clock
// timer(5) -> read absolute clock
#ifdef OPENMPCODE
  if(i==0)
    relbeg_kd=omp_get_wtime();
  else if(i==1) {
    relend_kd=omp_get_wtime();
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend_kd-relbeg_kd));
  }
  else if(i==2) {
    relend_kd=omp_get_wtime();
    printf("    Relative time ellapsed %.1lf ms\n",1000*(relend_kd-relbeg_kd));
    relbeg_kd=omp_get_wtime();
  }
  else if(i==4)
    absbeg_kd=omp_get_wtime();
  else if(i==5) {
    absend_kd=omp_get_wtime();
    printf("    Total time ellapsed %.1lf ms\n",1000*(absend_kd-absbeg_kd));
  }
#else //OPENMPCODE
  int diff;
  
  if(i==0)
    relbeg_kd=time(NULL);
  else if(i==1) {
    relend_kd=time(NULL);
    diff=(int)(difftime(relend_kd,relbeg_kd));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
  }
  else if(i==2) {
    relend_kd=time(NULL);
    diff=(int)(difftime(relend_kd,relbeg_kd));
    printf("    Relative time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
    relbeg_kd=time(NULL);
  }
  else if(i==4)
    absbeg_kd=time(NULL);
  else if(i==5) {
    absend_kd=time(NULL);
    diff=(int)(difftime(absend_kd,absbeg_kd));
    printf("    Total time ellapsed %02d:%02d:%02d \n",
       diff/3600,(diff/60)%60,diff%60);
  }
#endif //OPENMPCODE
}

//E Time section


void error(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);
    exit(1);
}

void verb_print_q(int iq, int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose == iq) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

// DEBUG WARNING!!
//  check lines: 91--106 routine InputData in cballsio.c
void verb_print(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose > 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_print_zero(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose >= 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_print_warning(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose >= 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_print_debug(int verbose, string fmt, ...)
{
    va_list ap;

    if (verbose > 0) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
}

void verb_log_print(int verbose, stream sout,  string fmt, ...)
{
    va_list ap;

    if (verbose > 0 && sout != NULL) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }
}

void verb_print_no_info(int verbose, int verbose_log, stream sout,
                  string fmt, ...)
{
    va_list ap;

    if (verbose >= VERBOSENOINFO) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }

    if (verbose_log >= VERBOSEMININFO && sout != NULL) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }

}

void verb_print_min_info(int verbose, int verbose_log, stream sout,
                  string fmt, ...)
{
    va_list ap;

    if (verbose >= VERBOSEMININFO) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
    if (verbose_log >= VERBOSEMININFO && sout != NULL) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }

}

void verb_print_normal_info(int verbose, int verbose_log, stream sout,
                            string fmt, ...)
{
    va_list ap;

    if (verbose >= VERBOSENORMALINFO) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
    if (verbose_log >= VERBOSENORMALINFO && sout != NULL) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }

}

void verb_print_debug_info(int verbose, int verbose_log, stream sout,
                  string fmt, ...)
{
    va_list ap;

    if (verbose >= VERBOSEDEBUGINFO) {
        va_start(ap, fmt);
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
        va_end(ap);
    }
    if (verbose_log >= VERBOSEDEBUGINFO && sout != NULL) {
        va_start(ap, fmt);
        vfprintf(sout, fmt, ap);
        fflush(sout);
        va_end(ap);
    }

}

void endrun(int ierr)
{
  if(ierr)
    {
      fprintf(stdout,
		  "endrun called with an error level of %d\n\n\n", ierr);
      exit(1);
    }
  exit(0);
}


void eprintf(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);
}

bool scanopt(string opt, string key)
{
    char *op, *kp;
    
    op = (char *) opt;
    while (*op != '\0') {
        kp = key;
        while ((*op != ',' ? *op : '\0') == *kp) {
            if (*kp++ == '\0')
                return (TRUE);
            op++;
        }
        while (*op != '\0' && *op++ != ',')
            
            continue;
    }
    return (FALSE);
}


stream stropen(string name, string mode)
{
    bool inflag;
    int fds;
    stream res;
    struct stat buf;

    inflag = streq(mode, "r");
    if (name[0] == '-') {
        if (streq(name, "-")) {
            fds = dup(fileno(inflag ? stdin : stdout));
            if (fds == -1)
                error("stropen: cannot dup %s\n", inflag ? "stdin" : "stdout");
        } else {
            char errmsg[128];
            if (parse_int_checked(&name[1], &fds,
                                  errmsg, sizeof(errmsg),
                                  "file descriptor") == FAILURE)
                error("stropen: %s\n", errmsg);
        }
        res = fdopen(fds, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen: cannot open f.d. %d for %s\n", fds, inflag ? "input" : "output");
    } else {
        if (streq(mode, "w") && stat(name, &buf) == 0)
            error("stropen: file \"%s\" already exists\n", name);
        res = fopen(name, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen: cannot open file \"%s\" for %s\n", name, inflag ? "input" : "output");
    }
    return (res);
}

//B from cBalls
//stream stropen(string name, string mode)
int stropen_checked(string name, string mode, stream *out,
                    char *errmsg, size_t errmsg_size)
{
    string routineName = "stropen_checked";
    bool inflag;
    int fds;
    stream res;
    struct stat buf;

    //B
    if (out == NULL)
        return FAILURE;
    *out = NULL;

    if (name == NULL || mode == NULL) {
        if (errmsg && errmsg_size > 0)
            snprintf(errmsg, errmsg_size,
                     "%s: NULL name or mode\n", routineName);
        return FAILURE;
    }
    //E
    
    inflag = streq(mode, "r");
    if (name[0] == '-') {
        if (streq(name, "-")) {
            fds = dup(fileno(inflag ? stdin : stdout));
            if (fds == -1) {
                snprintf(errmsg, errmsg_size,
                         "%s: cannot dup %s\n",
                         routineName, inflag ? "stdin" : "stdout");
                return FAILURE;
            }
        } else if (parse_int_checked(&name[1], &fds,
                                     errmsg, errmsg_size,
                                     "file descriptor") == FAILURE) {
            return FAILURE;
        }
        res = fdopen(fds, streq(mode, "w!") ? "w" : mode);
        if (res == NULL) {
            snprintf(errmsg, errmsg_size,
                     "%s: cannot open f.d. %d for %s\n",
                     routineName, fds, inflag ? "input" : "output");
            return FAILURE;
        }
    } else {
        if (streq(mode, "w") && stat(name, &buf) == 0) {
            snprintf(errmsg, errmsg_size,
                     "%s: file \"%s\" already exists\n", routineName, name);
            return FAILURE;
        }
        res = fopen(name, streq(mode, "w!") ? "w" : mode);
        if (res == NULL) {
            snprintf(errmsg, errmsg_size,
                     "%s: cannot open file \"%s\" for %s\n", routineName, name, inflag ? "input" : "output");
            return FAILURE;
        }
    }

    *out = res;

    return SUCCESS;
}
//E


//B From G

double dmax(double x,double y)
{
  if(x>y)
    return x;
  else
    return y;
}

double dmin(double x,double y)
{
  if(x<y)
    return x;
  else
    return y;
}

int imax(int x,int y)
{
  if(x>y)
    return x;
  else
    return y;
}

int imin(int x,int y)
{
  if(x<y)
    return x;
  else
    return y;
}

//E


//B common section
#ifndef CLASSLIB

void class_protect_sprintf(char *dest, size_t dest_size, const char *tpl, ...) {
  va_list args;

  if (dest == NULL || dest_size == 0)
    return;

  va_start(args, tpl);
  vsnprintf(dest, dest_size, tpl, args);
  va_end(args);

  dest[dest_size - 1] = '\0';
}

void class_protect_fprintf(FILE *stream, const char *tpl, ...) {
  va_list args;
  char dest[6000];

  va_start(args, tpl);
  vsnprintf(dest, sizeof(dest), tpl, args);
  va_end(args);

  fprintf(stream, "%s", dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}

int get_number_of_titles(char * titlestring){
  int i;
  int number_of_titles=0;

  for (i=0; i<strlen(titlestring); i++){
    if (titlestring[i] == '\t')
      number_of_titles++;
  }
  return number_of_titles;
}

/*
 * Finds wether or not a file exists.
 *
 * param fname  Input: File name
 * return boolean
 */
int file_exists(const char *fname){
  FILE *file = fopen(fname, "r");
  if (file != NULL){
    fclose(file);
    return TRUE;
  }

  return FALSE;

}

/*
 * Finds whether two doubles are equal or which one is bigger
 *
 * param a Input: first number
 * param b Input: second number
 * return -1, 1 or 0
 */
int compare_doubles(const void *a,
                    const void *b){
  double *x = (double *) a;
  double *y = (double *) b;
  if (*x < *y)
    return -1;
  else if
    (*x > *y) return 1;
  return 0;
}

/*
 * This function detects if a string begins with a character,
 * ignoring whitespaces during its search
 *
 * returns the result, NOT the SUCCESS or _FAILURE_ codes.
 * (This is done such that it can be used inside of an if statement)
 *
 * param thestring  Input: string to test
 * param beginchar  Input: the character by which the string begins or not
 * return boolean
 */

int string_begins_with(char* thestring, char beginchar){

  // Define temporary variables
  int int_temp=0;
  int strlength = strlen((thestring));
  int result = FALSE;

  // Check through the beginning of the string to see if the beginchar is met
  for(int_temp=0;int_temp<strlength;++int_temp){
    // Skip over whitespaces (very important)
    if(thestring[int_temp]==' ' || thestring[int_temp]=='\t'){continue;}
    // If the beginchar is met, everything is good
    else if(thestring[int_temp]==beginchar){result=TRUE;}
    // If something else is met, cancel
    else{break;}
  }

  return result;
}

#endif
//E common section

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

// to stop using this system()...
// sprintf(outputDir,cmd->rootDir);
// sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi", outputDir, outputDir);
// system(buf);
int mkdir_p(const char *path, mode_t mode) {
    char tmp[MAXLENGTHOFFILES];
    char *p = NULL;
    size_t len;

    int nwritten;

    nwritten = snprintf(tmp, sizeof(tmp), "%s", path);
    if (nwritten < 0 || (size_t)nwritten >= sizeof(tmp)) {
        errno = ENAMETOOLONG;
        return -1;
    }
    
    len = strlen(tmp);
    if (len == 0) {
        return 0;
    }

    // Remove a trailing slash to avoid double-processing the final directory
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = '\0';
    }

    // Traverse the path looking for directory delimiters
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0'; // Temporarily truncate the string
            
            // Attempt to create the parent directory segment
            if (mkdir(tmp, mode) != 0) {
                if (errno != EEXIST) {
                    return -1; // Fail if error is something other than "already exists"
                }
            }
            
            *p = '/'; // Restore the delimiter
        }
    }

    // Create the final target directory
    if (mkdir(tmp, mode) != 0 && errno != EEXIST) {
        return -1;
    }

    return 0;
}

//B to test
int main_mkdir(void) {
    // 0755: Read, write, execute for owner; read and execute for group/others
    if (mkdir_p("nested/path/to/my/folder", 0755) == 0) {
        printf("Directory tree created successfully.\n");
    } else {
        perror("Error creating directories");
    }
    return 0;
}
//E


//B to solve unsafe strcpy / sprintf usage
#include <errno.h>
#include <stdio.h>
#include <string.h>

int copy_checked(char *dst, size_t dst_size, const char *src,
                        const char *label)
{
    int n;

    if (dst == NULL || src == NULL || dst_size == 0)
        return -1;

    n = snprintf(dst, dst_size, "%s", src);

    if (n < 0 || (size_t)n >= dst_size) {
        fprintf(stderr, "%s too long: max %zu bytes\n", label, dst_size - 1);
        errno = ENAMETOOLONG;
        return -1;
    }

    return 0;
}
//E

//B added by cBalls
int format_checked(char *dst, size_t dst_size,
                   const char *label, const char *fmt, ...)
{
    va_list args;
    int n;

    if (dst == NULL || fmt == NULL || dst_size == 0)
        return -1;

    va_start(args, fmt);
    n = vsnprintf(dst, dst_size, fmt, args);
    va_end(args);

    if (n < 0 || (size_t)n >= dst_size) {
        dst[0] = '\0';
        fprintf(stderr, "%s too long: max %zu bytes\n", label, dst_size - 1);
        errno = ENAMETOOLONG;
        return -1;
    }

    return 0;
}
//E
