//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"
#include <stdarg.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef __COMMON__
#define __COMMON__

//#define _VERSION_ "v3.2.1"
#define _VERSION_ "v1.0.0"         // It must match version in cmdline_defs.h

/* @cond INCLUDE_WITH_DOXYGEN */

#if !defined(TRUE)
#define TRUE 1 /**< integer associated to true statement */
#endif
#if !defined(FALSE)
#define FALSE 0 /**< integer associated to false statement */
#endif

#define SUCCESS 0 /**< integer returned after successful call of a function */
#define FAILURE 1 /**< integer returned after failure in a function */

#define _ERRORMSGSIZE_ 2048 /**< generic error messages are cut beyond this number of characters */
typedef char ErrorMsg[_ERRORMSGSIZE_]; /**< Generic error messages (there is such a field in each structure) */

#define _FILENAMESIZE_ 256 /**< size of the string read in each line of the file (extra characters not taken into account) */
typedef char FileName[_FILENAMESIZE_];

#define _SUFFIXNAMESIZE_ 4 /**< maximum size of the short string appended to file names to account for initial conditions, etc. */

#define _PI_ 3.1415926535897932384626433832795e0 /**< The number pi */

#define _PIHALF_ 1.57079632679489661923132169164e0 /**< pi divided by 2 */

#define _TWOPI_ 6.283185307179586476925286766559e0 /**< 2 times pi */

#define _SQRT2_ 1.41421356237309504880168872421e0 /** < square root of 2. */

#define _SQRT6_ 2.4494897427831780981972840747059e0 /**< square root of 6. */

#define _SQRT_PI_ 1.77245385090551602729816748334e0 /**< square root of pi. */

#define _E_ 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427427466391932003059921817413596629043572900334295260595630738132328627943490763233829880753195251019011573834187930702154089149934884167509244761460668082264800168477411853742345442437107539077744992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977209310141692836819025515108657463772111252389784425056953696 /**< exponential of one */

#define _MAX_IT_ 10000/**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ 250 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _QUADRATURE_MAX_BG_ 800 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _TOLVAR_ 100. /**< The minimum allowed variation is the machine precision times this number */

#define _HUGE_ 1.e99

#define _EPSILON_ 1.e-10

#define _OUTPUTPRECISION_ 12 /**< Number of significant digits in some output files */

#define _COLUMNWIDTH_ 24 /**< Must be at least _OUTPUTPRECISION_+8 for guaranteed fixed width columns */

#define _MAXTITLESTRINGLENGTH_ 8000 /**< Maximum number of characters in title strings */

#define _DELIMITER_ "\t" /**< character used for delimiting titles in the title strings */

#ifndef __CLASSDIR__
#define __CLASSDIR__ "." /**< The directory of CLASS. This is set to the absolute path to the CLASS directory so this is just a failsafe. */
#endif

#define MIN(a,b) (((a)<(b)) ? (a) : (b) ) /**< the usual "min" function */
#define MAX(a,b) (((a)<(b)) ? (b) : (a) ) /**< the usual "max" function */
#define SIGN(a) (((a)>0) ? 1. : -1. )
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* @endcond */

/* needed because of weird openmp bug on macosx lion... */

void class_protect_sprintf(char* dest, char* tpl,...);
void class_protect_fprintf(FILE* dest, char* tpl,...);
void* class_protect_memcpy(void* dest, void* from, size_t sz);

/* some general functions */

int get_number_of_titles(char * titlestring);
int file_exists(const char *fname);
int compare_doubles(const void * a,
                    const void * b);
int string_begins_with(char* thestring, char beginchar);

/* general CLASS macros */

#define class_build_error_string(dest,tmpl,...) {                               \
  ErrorMsg FMsg;                                                                \
  class_protect_sprintf(FMsg,tmpl,__VA_ARGS__);                                 \
  class_protect_sprintf(dest,"%s(L:%d) :%s",__func__,__LINE__,FMsg);            \
}

// Error reporting macros

// Call
#define class_call_message(err_out,extra,err_mess)                              \
  class_build_error_string(err_out,"error in %s;\n=>%s",extra,err_mess);

/* macro for calling function and returning error if it failed */
#define class_call_except(function, error_message_from_function, error_message_output,list_of_commands) { \
  if (function == FAILURE) {                                                        \
    class_call_message(error_message_output,#function,error_message_from_function); \
    list_of_commands;                                                               \
    return FAILURE;                                                                 \
  }                                                                                 \
}

/* macro for trying to call function */
#define class_call_try(function, error_message_from_function, error_message_output,list_of_commands) { \
  if (function == FAILURE) {                                                         \
    class_call_message(error_message_output,#function,error_message_from_function);  \
    list_of_commands;                                                                \
  }                                                                                  \
}

/* macro for calling function and returning error if it failed */
#define class_call(function, error_message_from_function, error_message_output)      \
  class_call_except(function, error_message_from_function,error_message_output,)

/* same in parallel region */
#define class_call_parallel(function, error_message_from_function, error_message_output) { \
  if (abort == FALSE) {                                                                    \
    if (function == FAILURE) {                                                             \
      class_call_message(error_message_output,#function,error_message_from_function);      \
      abort=TRUE;                                                                          \
    }                                                                                      \
  }                                                                                        \
}




// Alloc
#define class_alloc_message(err_out,extra,sz)                                       \
  class_build_error_string(err_out,"could not allocate %s with size %d",extra,sz);

/* macro for allocating memory and returning error if it failed */
#define class_alloc(pointer, size, error_message_output)  {                         \
  pointer=malloc(size);                                                             \
  if (pointer == NULL) {                                                            \
    int size_int;                                                                   \
    size_int = size;                                                                \
    class_alloc_message(error_message_output,#pointer, size_int);                   \
    return FAILURE;                                                                 \
  }                                                                                 \
}

/* same inside parallel structure */
#define class_alloc_parallel(pointer, size, error_message_output)  {                \
  pointer=NULL;                                                                     \
  if (abort == FALSE) {                                                             \
    pointer=malloc(size);                                                           \
    if (pointer == NULL) {                                                          \
      int size_int;                                                                 \
      size_int = size;                                                              \
      class_alloc_message(error_message_output,#pointer, size_int);                 \
      abort=TRUE;                                                                   \
    }                                                                               \
  }                                                                                 \
}

/* macro for allocating memory, initializing it with zeros/ and returning error if it failed */
#define class_calloc(pointer, init,size, error_message_output)  {                   \
  pointer=calloc(init,size);                                                        \
  if (pointer == NULL) {                                                            \
    int size_int;                                                                   \
    size_int = size;                                                                \
    class_alloc_message(error_message_output,#pointer, size_int);                   \
    return FAILURE;                                                                 \
  }                                                                                 \
}

/* macro for re-allocating memory, returning error if it failed */
#define class_realloc(pointer, newname, size, error_message_output)  {              \
    pointer=realloc(newname,size);                                                  \
  if (pointer == NULL) {                                                            \
    int size_int;                                                                   \
    size_int = size;                                                                \
    class_alloc_message(error_message_output,#pointer, size_int);                   \
    return FAILURE;                                                                 \
  }                                                                                 \
}

// Testing

#define class_test_message(err_out,extra,args...) {                                 \
  ErrorMsg Optional_arguments;                                                      \
  class_protect_sprintf(Optional_arguments,args);                                   \
  class_build_error_string(err_out,"condition (%s) is true; %s",extra,Optional_arguments); \
}

/* macro for testing condition and returning error if condition is true;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_test_except(condition, error_message_output,list_of_commands, args...) { \
  if (condition) {                                                                     \
    class_test_message(error_message_output,#condition, args);                         \
    list_of_commands;                                                                  \
    return FAILURE;                                                                    \
  }                                                                                    \
}

#define class_test(condition, error_message_output, args...) {                         \
  if (condition) {                                                                     \
    class_test_message(error_message_output,#condition, args);                         \
    return FAILURE;                                                                    \
  }                                                                                    \
}

#define class_test_parallel(condition, error_message_output, args...) {                \
  if (abort == FALSE) {                                                                \
    if (condition) {                                                                   \
      class_test_message(error_message_output,#condition, args);                       \
      abort=TRUE;                                                                      \
    }                                                                                  \
  }                                                                     \
}

/* macro for returning error message;
   args is a variable list of optional arguments, e.g.: args="x=%d",x
   args cannot be empty, if there is nothing to pass use args="" */
#define class_stop(error_message_output,args...) {                                      \
  ErrorMsg Optional_arguments;                                                          \
  class_protect_sprintf(Optional_arguments,args);                                       \
  class_build_error_string(error_message_output,"error; %s",Optional_arguments);        \
  return FAILURE;                                                                       \
}

// IO
/* macro for opening file and returning error if it failed */
#define class_open(pointer, filename,	mode, error_output) {                           \
  pointer=fopen(filename,mode);                                                         \
  if (pointer == NULL) {                                                                \
    class_build_error_string(error_output,"could not open %s with name %s and mode %s",#pointer,filename,#mode); \
    return FAILURE;                                                                     \
  }                                                                                     \
}

/* macro for defining indices (usually one, sometimes a block) */
#define class_define_index(index,                                       \
                           condition,                                   \
                           running_index,                               \
                           number_of_indices) {                         \
    if (condition) {                                                    \
      index = running_index;                                            \
      running_index += number_of_indices;                               \
    }                                                                   \
  }

/* macros for writing formatted output */
#define class_fprintf_double(file,                                      \
                             output,                                    \
                             condition){                                \
    if (condition == TRUE)                                              \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,output);    \
  }

#define class_fprintf_double_or_default(file,                           \
                                        output,                         \
                                        condition,                      \
                                        defaultvalue){                  \
    if (condition == TRUE)                                              \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,output);    \
    else                                                                \
      fprintf(file,"%*.*e ",_COLUMNWIDTH_,_OUTPUTPRECISION_,defaultvalue); \
}

#define class_fprintf_int(file,                                         \
                          output,                                       \
                          condition){                                   \
    if (condition == TRUE)                                              \
      fprintf(file,"%*d%*s ",                                           \
              MAX(0,_COLUMNWIDTH_-_OUTPUTPRECISION_-5),                 \
              output, _OUTPUTPRECISION_+5," ");                         \
  }

#define class_fprintf_columntitle(file,                                 \
                                  title,                                \
                                  condition,                            \
                                  colnum){                              \
    if (condition == TRUE)                                              \
      fprintf(file,"%*s%2d:%-*s ",                                      \
              MAX(0,MIN(_COLUMNWIDTH_-_OUTPUTPRECISION_-6-3,_COLUMNWIDTH_-((int) strlen(title))-3)), \
              "",colnum++,_OUTPUTPRECISION_+6,title);                   \
  }

#define class_store_columntitle(titlestring,                            \
				title,					\
				condition){				\
    if (condition == TRUE){                                           \
      strcat(titlestring,title);                                        \
      strcat(titlestring,_DELIMITER_);                                  \
    }                                                                   \
  }
//,_MAXTITLESTRINGLENGTH_-strlen(titlestring)-1);

#define class_store_double(storage,					\
			   value,					\
			   condition,                                   \
                           dataindex){                                  \
    if (condition == TRUE)                                            \
      storage[dataindex++] = value;                                     \
  }

#define class_store_double_or_default(storage,                          \
                                      value,                            \
                                      condition,                        \
                                      dataindex,                        \
                                      defaultvalue){                    \
    if (condition == TRUE)                                            \
      storage[dataindex++] = value;                                     \
    else                                                                \
      storage[dataindex++] = defaultvalue;                              \
}


#endif
