
cfitsiolib folder

cfitsio version may cause small differences in the results when using fits input format.

Chose cfitsio version by using a linux link:

$ rm cfitsio
$ ln -s cfitsio_ver_4.4.1 cfitsio

Almost only .c and .h files kept from de original cfitsio_4.4.1

Almost only .c and .h files kept from de original cfitsio_4.6.3
(
In: cfitsiolib/cfitsio/group.c ::
//B to correct: error: implicit declaration of function 'getcwd'; did you mean 'getw'?
#include <unistd.h> // Required for getcwd()
//E
)

