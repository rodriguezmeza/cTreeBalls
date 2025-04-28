
#ifndef __INPUT__
#define __INPUT__

#include "globaldefs.h"
#include "common.h"
#include "parser.h"

#define _N_FILEROOT_ 100

#define class_read_double(name,destination)                                     \
  do {                                                                          \
    double param_temp; int flag_temp;                                           \
    class_call(parser_read_double(pfc,name,&param_temp,&flag_temp,errmsg),      \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      destination = param_temp;                                                 \
    }                                                                           \
  } while(0);


#define class_read_int(name,destination)                                        \
  do {                                                                          \
    int int_temp,flag_temp;                                                     \
    class_call(parser_read_int(pfc,name,&int_temp,&flag_temp,errmsg),           \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      destination = int_temp;                                                   \
    }                                                                           \
  } while(0);

#define class_read_string(name,destination)                                     \
  do {                                                                          \
    char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;                \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      strcpy(destination,string_temp);                                          \
    }                                                                           \
  } while(0);

#define class_read_flag(name,destination)                                       \
  do {                                                                          \
char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;            \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      if( string_begins_with(string_temp,'y')                                   \
         || string_begins_with(string_temp,'Y') ){                              \
        destination = TRUE;                                                   \
      }                                                                         \
      else if( string_begins_with(string_temp,'n')                              \
         || string_begins_with(string_temp,'N') ){                              \
        destination = FALSE;                                                  \
      }                                                                         \
      else {                                                                    \
        class_stop(errmsg,"incomprehensible input '%s' for the field '%s'.",    \
                   string_temp, name);                                          \
      }                                                                         \
    }                                                                           \
  } while(0);

#define class_read_flag_or_deprecated(name,oldname,destination)                 \
  do {                                                                          \
    char string_temp[_ARGUMENT_LENGTH_MAX_]; int flag_temp;                     \
    class_call(parser_read_string(pfc,name,&string_temp,&flag_temp,errmsg),     \
               errmsg,                                                          \
               errmsg);                                                         \
    /* Compatibility code BEGIN */                                              \
    if(flag_temp == FALSE){                                                   \
      class_call(parser_read_string(pfc,oldname,&string_temp,&flag_temp,errmsg),\
                 errmsg,                                                        \
                 errmsg);                                                       \
    }                                                                           \
    /* Compatibility code END */                                                \
    if (flag_temp == TRUE){                                                   \
      if( string_begins_with(string_temp,'y')                                   \
         || string_begins_with(string_temp,'Y') ){                              \
        destination = TRUE;                                                   \
      }                                                                         \
      else if( string_begins_with(string_temp,'n')                              \
         || string_begins_with(string_temp,'N') ){                              \
        destination = FALSE;                                                  \
      }                                                                         \
      else {                                                                    \
        class_stop(errmsg,"incomprehensible input '%s' for the field '%s'.",    \
                   string_temp, name);                                          \
      }                                                                         \
    }                                                                           \
  } while(0);

#define class_read_double_one_of_two(name1,name2,destination)                   \
  do {                                                                          \
    int flag_temp1,flag_temp2;                                                  \
    double param_temp1,param_temp2;                                             \
    class_call(parser_read_double(pfc,name1,&param_temp1,&flag_temp1,errmsg),   \
               errmsg,                                                          \
               errmsg);                                                         \
    class_call(parser_read_double(pfc,name2,&param_temp2,&flag_temp2,errmsg),   \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test((flag_temp1 == TRUE) && (flag_temp2 == TRUE),                \
               errmsg,                                                          \
               "You can only enter one of '%s' or '%s'.",                       \
               name1,name2);                                                    \
    if (flag_temp1 == TRUE){                                                  \
      destination = param_temp1;                                                \
    }                                                                           \
    if (flag_temp2 == TRUE){                                                  \
      destination = param_temp2;                                                \
    }                                                                           \
  } while(0);

#define class_at_least_two_of_three(a,b,c)                                      \
  (((a) == TRUE) && ((b) == TRUE)) ||                                       \
  (((a) == TRUE) && ((c) == TRUE)) ||                                       \
  (((b) == TRUE) && ((c) == TRUE))

#define class_none_of_three(a,b,c)                                              \
  ((a) == FALSE) && ((b) == FALSE) && ((c) == FALSE)

#define class_read_list_of_doubles_or_default(name,destination,val_default,siz) \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_doubles(pfc,name,                            \
                                           &entries_read_temp,                  \
                                           &(destination),                      \
                                           &flag_temp,                          \
                                           errmsg),                             \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      class_test(entries_read_temp != siz, errmsg,                              \
                 "Number of entries of '%s' (%d) does not match expected number (%d).", \
                 name, entries_read_temp, siz);                                 \
    }else{                                                                      \
      class_alloc(destination,siz*sizeof(double),errmsg);                       \
      for(n=0; n<siz; n++){destination[n] = val_default;}                       \
    }                                                                           \
  } while(0);

#define class_read_list_of_integers_or_default(name,destination,val_default,siz)\
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_integers(pfc,name,                           \
                                            &entries_read_temp,                 \
                                            &(destination),                     \
                                            &flag_temp,                         \
                                            errmsg),                            \
               errmsg,                                                          \
               errmsg);                                                         \
    if (flag_temp == TRUE){                                                   \
      class_test(entries_read_temp != siz, errmsg,                              \
                 "Number of entries of '%s' (%d) does not match expected number (%d).", \
                 name, entries_read_temp, siz);                                 \
    }else{                                                                      \
      class_alloc(destination,siz*sizeof(int),errmsg);                          \
      for(n=0; n<siz; n++){destination[n] = val_default;}                       \
    }                                                                           \
  } while(0);

#define class_read_list_of_doubles(name,destination,siz)                        \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_doubles(pfc,name,                            \
                                           &entries_read_temp,                  \
                                           &(destination),                      \
                                           &flag_temp,                          \
                                           errmsg),                             \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test(flag_temp == FALSE, errmsg,                                    \
               "Entry '%s' is required but not found!", name)                   \
    class_test(entries_read_temp != siz, errmsg,                                \
               "Number of entries of '%s' (%d) does not match expected number (%d).", \
               name,entries_read_temp, siz);                                    \
  } while(0);

#define class_read_list_of_integers(name,destination,siz)                       \
  do {                                                                          \
    int flag_temp,entries_read_temp;                                            \
    class_call(parser_read_list_of_integers(pfc,name,                           \
                                            &entries_read_temp,                 \
                                            &(destination),                     \
                                            &flag_temp,                         \
                                            errmsg),                            \
               errmsg,                                                          \
               errmsg);                                                         \
    class_test(flag_temp == FALSE, errmsg,                                    \
               "Entry '%s' is required but not found!", name)                   \
    class_test(entries_read_temp != siz, errmsg,                                \
               "Number of entries of '%s' (%d) does not match expected number (%d).", \
               name,entries_read_temp, siz);                                    \
  } while(0);


/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int input_find_file(struct  cmdline_data* cmd, char *fname,
                    struct file_content * fc,
                    ErrorMsg errmsg);

  int input_set_root(char* input_file,
                     struct file_content** ppfc_input,
                     struct file_content* pfc_setroot,
                     ErrorMsg errmsg);

  int input_read_from_file(struct cmdline_data *cmd, struct file_content * pfc,
                           ErrorMsg errmsg);


  int input_read_parameters(struct cmdline_data *cmd, struct file_content * pfc,
                            ErrorMsg errmsg);

    int input_read_parameters_general_free(struct file_content * pfc,
                                           ErrorMsg errmsg);

  int input_read_parameters_general(struct cmdline_data *cmd, struct file_content * pfc,
                                    ErrorMsg errmsg);

  int input_default_params(struct cmdline_data *);

  int tpcf_version( char * version);

#ifdef __cplusplus
}
#endif

#endif // ! __INPUT__
