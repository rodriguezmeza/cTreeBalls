
#ifndef _cballsio_gadget_00_h
#define _cballsio_gadget_00_h

// This tag number must be diferent than others use in I/O
#define INGADGET                5


global int inputdata_gadget(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           string filename, int ifile);

#endif	// ! _cballsio_gadget_00_h
