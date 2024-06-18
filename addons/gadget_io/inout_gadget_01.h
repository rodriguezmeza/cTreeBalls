// Use:
//ADDONS:
//#include "inout_gadget_01.h"

#ifndef _inout_gadget_01_h
#define _inout_gadget_01_h

//B Some gadget definitions
typedef struct {
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8]; // fills to 256 Bytes
} io_header_1;

typedef struct {
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  int    Type;

  float  Rho, U, Temp, Ne;
} particle_data, *particle_data_ptr;

int allocateMemory(int**IdPointer, particle_data **Ppointer);
int load_snapshot(char *shortfile,int empty, particle_data **Ppointer);
void PtoMesh(float *densityz, particle_data *Pz);// Does CiC assignment.
//E

//B cute_box
#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
//typedef int lint;
typedef INTEGER lint;
#endif //_LONGIDS

typedef struct {
  lint np;          //#objects in the catalog
  double *pos;
} Catalog_tpcf;         //Catalog (double precision)

typedef struct {
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8];
  // fills to 256 Bytes
} gad_header_tpcf;

typedef struct {
  char label[4];
  int size;
} gad_title_tpcf;

//static float l_box_tpcf;
//static float l_box_half_tpcf;
//static int input_format_tpcf;
//static lint n_objects_tpcf;

lint linecount_tpcf(FILE *f);
void error_mem_out_tpcf(void);
void error_read_line_tpcf(char *fname,lint nlin);
void error_open_file_tpcf(char *fname);

//Catalog read_gadget(char *prefix,lint *np,int input)
int read_catalog_tpcf(char *fname, lint *np, particle_data **Ppointer);

//E

#endif	// ! _inout_gadget_01_h
