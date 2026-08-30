// added by cBalls
//=============================================================================*/
//        1          2          3          4        ^ 5          6          7

#ifndef _cballs_class_call_h
#define _cballs_class_call_h

#include <stdio.h>

#define class_call_cballs(function, error_message_from_function, error_message_output)  \
    do {                                                                                \
        int _cballs_status = (function);                                                \
        if (_cballs_status == FAILURE) {                                                \
            ErrorMsg _cballs_local_error;                                               \
            snprintf(_cballs_local_error, sizeof(_cballs_local_error), "%s",            \
                     (error_message_from_function));                                    \
            snprintf((error_message_output), _ERRORMSGSIZE_,                            \
                     "%s(L:%d): error in %.512s;\n=>%.1024s",                          \
                     __func__, __LINE__, #function, _cballs_local_error);               \
            return FAILURE;                                                             \
        }                                                                               \
    } while (0)

#endif // ! _cballs_class_call_h
