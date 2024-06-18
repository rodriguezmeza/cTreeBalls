
#ifndef _tpc_io_gadget_00_h
#define _tpc_io_gadget_00_h

#define INGADGET                5

//B Some gadget definitions. Needed for reading gadget format
//global particle_data_ptr P;
global particle_data_ptr PP;
//global int NumPart;
global INTEGER NumPart;
global int GridScale;
global int NROW,NGRID;
global float BoxSize;
global float Scale;
//E

//B CTELIB
global float l_box_tpcf;
global float l_box_half_tpcf;
global int input_format_tpcf;
global lint n_objects_tpcf;
//E

local int inputdata_gadget(struct cmdline_data* cmd, struct  global_data* gd);

#endif	// ! _tpc_io_gadget_00_h
