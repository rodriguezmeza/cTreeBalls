
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"

//B Some gadget definitions
typedef struct {
    double  Pos[3];
    double  Vel[3];
    double  Mass;
  int    Type;

    double  Rho, U, Temp, Ne;
} particle_data, *particle_data_ptr;
//E

void PtoMesh(struct cmdline_data* cmd,
             struct  global_data* gd,
             float *densityz, particle_data *Pz);   // Does CiC assignment.

//B cute_box
#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef INTEGER lint;
#endif //_LONGIDS

typedef struct {
  lint np;                                          //#objects in the catalog
  double *pos;
} Catalog_tpcf;                                     //Catalog (double precision)

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


lint linecount_tpcf(FILE *f);
void error_mem_out_tpcf(void);
void error_read_line_tpcf(char *fname,lint nlin);
void error_open_file_tpcf(char *fname);

static Catalog_tpcf read_gadget(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                char *prefix,lint *np,int input);
int read_catalog_tpcf(struct cmdline_data* cmd,
                      struct  global_data* gd, int ifile,
                      char *fname, lint *np, particle_data **Ppointer);
//E cute_box


//B Some gadget definitions. Needed for reading gadget format
local INTEGER NumPart;
local int GridScale;
local int NROW,NGRID;
local double BoxSize;
local double Scale;

local double l_box_tpcf;
local double l_box_half_tpcf;
local int input_format_tpcf;
local lint n_objects_tpcf;
//E

global int inputdata_gadget(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           string filename, int ifile)
{
    float * density;
    particle_data_ptr Ppointer[1];
    bodyptr p;
    INTEGER ip;
    int k;
    real mass=1.0;
    real weight=1.0;

    gd->input_comment = "Gadget input file";

    read_catalog_tpcf(cmd, gd, ifile, cmd->infile, &NumPart, Ppointer);
    particle_data_ptr PP;
    PP = (particle_data_ptr) allocate(cmd->nbody * sizeof(particle_data));

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        ip = p-bodytable[ifile];
        DO_COORD(k)
            PP[ip].Pos[k] = Pos(p)[k];
    }

    GridScale=3;                                    // How many times finer
                                                    //  (in 1-dim) is the PM grid
                                                    //  compared to the
                                                    //  particle grid.
    NROW=(int) (0.5+powf(((float )NumPart),1./3.));
    NGRID=GridScale*NROW;
    Scale=2.*M_PI/BoxSize;
    verb_print_q(2, cmd->verbose,
                 "input_gadget: NROW,NGRID,Scale: %i %i %g\n",
                 NROW,NGRID,Scale);

    density=malloc(NGRID*NGRID*NGRID*sizeof(float)); // PM density grid.
    PtoMesh(cmd, gd, density, PP);

    DO_COORD(k)
    gd->Box[k] = BoxSize;

    cmd->nbody = NumPart;
    verb_print_debug(1, "\ninputdata_gadget:: cmd->nbody = %ld\n",cmd->nbody);
    fprintf(stdout,"inputdata_gadget:: Pos0 %g\n",PP[1].Pos[0]);
    fflush(stdout);

    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        ip = p-bodytable[ifile];
        Kappa(p) = 1.0 + rcos(2.0*PI*Pos(p)[0]/(gd->Box[0]/20.0))
                        * rsin(2.0*PI*Pos(p)[1]/(gd->Box[1]/20.0));
        if (scanopt(cmd->options, "kappa-constant"))
            Kappa(p) = 2.0;
        Type(p) = BODY;
        Mass(p) = mass;
        Weight(p) = weight;
        Id(p) = p-bodytable[ifile]+1;
    }

    real kavg=0.0;
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        kavg += Kappa(p);
    }
    verb_print(cmd->verbose,
               "inputdata_gadget: average of kappa (%ld particles) = %le\n",
               cmd->nbody, kavg/((real)cmd->nbody) );

    return SUCCESS;
}


void WRt(float *d,int i,int j,int k,float f)
     {
     d[k + NGRID * (j + NGRID * i)]=f;
     }

float REd(float *d,int i,int j,int k)
     {
     return d[k + NGRID * (j + NGRID * i)];
     }

void PtoMesh(struct cmdline_data* cmd,
             struct  global_data* gd,
             float *densityz, particle_data *Pz)    // Does CiC assignment.
{
    int i;

     printf("Calculating PtoMesh\n");

     float scaleBox=((float) NGRID)/((float) BoxSize);
     float WPAR=pow(((float) NGRID)/((float)NROW),3);

     for(i=0; i<NGRID*NGRID*NGRID; i++)
            densityz[i] = -1.0;
            
     float X,Y,Z;
     float T1,T2,T3;
     float D1,D2,D3;
     float D2W,T2W;
     int iI,J,K,I1,J1,K1;

    for(i=0; i<NumPart; i++)
     {
            X=Pz[i].Pos[0]*scaleBox;
            Y=Pz[i].Pos[1]*scaleBox;
            Z=Pz[i].Pos[2]*scaleBox;
            
            iI=(int) X;
            J=(int) Y;
            K=(int) Z;
            D1=X-((float) iI);
            D2=Y-((float) J);
            D3=Z-((float) K);
            T1=1.-D1;
            T2=1.-D2;
            T3=1.-D3;


            T2W =T2*WPAR;
            D2W =D2*WPAR;

            if(iI >= NGRID)iI=0;
            if(J >= NGRID)J=0;
            if(K >= NGRID)K=0;
            
            I1=iI+1;
               if(I1 >= NGRID)I1=0;
            J1=J+1;
               if(J1 >= NGRID)J1=0;
            K1=K+1;
               if(K1 >= NGRID)K1=0;

              WRt(densityz,iI ,J ,K ,REd(densityz,iI ,J ,K ) +T3*T1*T2W);
              WRt(densityz,I1,J ,K , REd(densityz,I1,J ,K ) +T3*D1*T2W);
              WRt(densityz,iI ,J1,K ,REd(densityz,iI ,J1,K ) +T3*T1*D2W);
              WRt(densityz,I1,J1,K , REd(densityz,I1,J1,K ) +T3*D1*D2W);
    
              WRt(densityz,iI ,J ,K1 ,REd(densityz,iI ,J ,K1 ) +D3*T1*T2W);
              WRt(densityz,I1,J ,K1 , REd(densityz,I1,J ,K1 ) +D3*D1*T2W);
              WRt(densityz,iI ,J1,K1 ,REd(densityz,iI ,J1,K1 ) +D3*T1*D2W);
              WRt(densityz,I1,J1,K1 , REd(densityz,I1,J1,K1 ) +D3*D1*D2W);
              
    }

    verb_print_q(2, cmd->verbose, "CIC density assignment finished\n");

      return;
}

//
//B Cute_box :: I/O from cute_box's David Alonso.
//  Adapted to cBalls
//

void error_open_file_tpcf(char *fname)
{
  // Open error handler
  fprintf(stderr,"CUTE: Couldn't open file %s \n",fname);
  exit(1);
}

static size_t my_fread(void *p,size_t size,
               size_t nmemb,FILE *stream)
{
  // Self-checked binary reading routine
  size_t nread;

  if((nread=fread(p,size,nmemb,stream))!=nmemb) {
    fprintf(stderr,"CUTE: error reading binary file \n");
    exit(1);
  }
  return nread;
}

void error_mem_out_tpcf(void)
{
  // Memory shortage handler
  fprintf(stderr,"CUTE: Out of memory!!\n");
  exit(1);
}

static void gad_check_block(int b1,int b2)
{
  // Checks that a block from a snapshot is
  // consistent from its begin/end values
  if(b1!=b2) {
    fprintf(stderr,"CUTE: Corrupted block!\n");
    exit(1);
  }
}


static int gad_seek_block(FILE *snap,char name[])
{
  // Seeks block from title
  gad_title_tpcf tit;
  int block1,block2;

  rewind(snap);

  while(1>0) {
    if(!(fread(&block1,sizeof(int),1,snap))||
       feof(snap)||ferror(snap)) {
      fprintf(stderr,"CUTE: Block %s not found!!\n",name);
      exit(1);
    }
    my_fread(&tit,sizeof(gad_title_tpcf),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);
    if(strncmp(tit.label,name,3)!=0)
      fseek(snap,tit.size,SEEK_CUR);
    else
      break;
  }

  return 0;
}

static int check_num_files(char *prefix)
{
  FILE *fil;

  fil=fopen(prefix,"r");
  if(fil!=NULL) {
    fclose(fil);
    return 1;
  }
  else {
    int nfils=0;
    while(nfils>=0) {
      char fname[256];
      sprintf(fname,"%s.%d",prefix,nfils);
      fil=fopen(fname,"r");
      if(fil!=NULL) {
    fclose(fil);
    nfils++;
      }
      else {
    if(nfils==0) {
      fprintf(stderr,"check_num_files: can't find file %s or %s.x\n",prefix,prefix);
      return -1;
    }
    else if(nfils==1) {
      fprintf(stderr,"check_num_files: only file %s found. Weird.\n",fname);
      return -1;
    }
    else {
      return nfils;
    }
      }
    }
  }
  
  fprintf(stderr,"check_num_files: this shouldn't have happened \n");
  return -1;
}

static Catalog_tpcf read_snapshot_single(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                         char *fname,lint *np,int input)
{
    // Creates catalog from a single snapshot file
    Catalog_tpcf cat;
    lint ii;
    gad_header_tpcf head;
    int block1,block2;
    
    FILE *snap=fopen(fname,"r");
    if(snap==NULL) error_open_file_tpcf(fname);
    
    //Read header
    if(input==2)
        gad_seek_block(snap,"HEAD");
    my_fread(&block1,sizeof(int),1,snap);
    my_fread(&head,sizeof(gad_header_tpcf),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);
    
    if(head.num_files!=1) {
        fprintf(stderr,"CUTE: Multi-file input not expected \n");
        exit(1);
    }

    if (cmd->verbose>=2) {
        printf("  The cosmological model is:\n");
    printf("   - Omega_M = %.3lf\n",head.Omega0);
    printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
    printf("   - h = %.3lf\n",head.HubbleParam);
    printf("  This file contains: \n");
    for(ii=0;ii<6;ii++) {
        printf("   - %d particles of type %d with mass",
               head.npart[ii],(int)ii);
        printf(" %.3lE (%d in total)\n",
               head.mass[ii],head.npartTotal[ii]);
    }
    printf("  The box size is %.3lf\n",head.BoxSize);
    printf("  Redshift z = %.3lf \n",head.redshift);
    } // ! verbose

  l_box_tpcf=head.BoxSize;
  l_box_half_tpcf=l_box_tpcf*0.5;
  cat.np=0;
  for(ii=0;ii<6;ii++) {
    cat.np+=head.npart[ii];
    if(head.npart[ii]!=head.npartTotal[ii]) {
      fprintf(stderr,"CUTE: error reading snapshot \n");
      exit(1);
    }
  }
  *np=cat.np;

  cat.pos=(double *)malloc(3*cat.np*sizeof(double));
  if(cat.pos==NULL)
    error_mem_out_tpcf();

  if(input==2)
    gad_seek_block(snap,"POS");
  my_fread(&block1,sizeof(int),1,snap);
  for(ii=0;ii<cat.np;ii++) {
    float pos[3];
    my_fread(pos,sizeof(float),3,snap);
    cat.pos[3*ii]=(double)(pos[0]);
    cat.pos[3*ii+1]=(double)(pos[1]);
    cat.pos[3*ii+2]=(double)(pos[2]);
  }
  my_fread(&block2,sizeof(int),1,snap);
  gad_check_block(block1,block2);
  fclose(snap);

  return cat;
}

static Catalog_tpcf read_gadget(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                char *prefix,lint *np,int input)
{
  Catalog_tpcf cat;
  int nfils=check_num_files(prefix);
  if(nfils<=0) exit(1);

    verb_print_q(2, cmd->verbose, "  Reading from GADGET snapshot format \n");
  
  if(nfils==1) {
    printf("  Reading single snapshot file\n");
    cat=read_snapshot_single(cmd, gd, prefix,np,input);
    return cat;
  }
  else {
    lint ii;
    char fname[256];
    gad_header_tpcf head;
    int block1,block2;
    FILE *snap;

      verb_print_q(2, cmd->verbose, "  Reading %d snapshot files \n",nfils);

    sprintf(fname,"%s.0",prefix);
    snap=fopen(fname,"r");
    if(snap==NULL) error_open_file_tpcf(fname);
  
    //Read header
    if(input==2)
      gad_seek_block(snap,"HEAD");
    my_fread(&block1,sizeof(int),1,snap);
    my_fread(&head,sizeof(gad_header_tpcf),1,snap);
    my_fread(&block2,sizeof(int),1,snap);
    gad_check_block(block1,block2);

    if(head.num_files!=nfils) {
      fprintf(stderr,
          "CUTE: Header and existing files do not match %d != %d.\n",
          nfils,head.num_files);
      fprintf(stderr,"      There may be some files missing\n");
      exit(1);
    }

      if (cmd->verbose>=2) {
          printf("  The cosmological model is:\n");
          printf("   - Omega_M = %.3lf\n",head.Omega0);
          printf("   - Omega_L = %.3lf\n",head.OmegaLambda);
          printf("   - h = %.3lf\n",head.HubbleParam);
          printf("  This file contains: \n");
          for(ii=0;ii<6;ii++) {
              printf("   - %d particles of type %d with mass %.3lE\n",
                     head.npartTotal[ii],(int)ii,head.mass[ii]);
          }
          printf("  The box size is %.3lf\n",head.BoxSize);
          printf("  Redshift z = %.3lf \n",head.redshift);
      }

      BoxSize=head.BoxSize;

    l_box_tpcf=head.BoxSize;
    l_box_half_tpcf=l_box_tpcf*0.5;
    cat.np=0;
    for(ii=0;ii<6;ii++)
      cat.np+=head.npartTotal[ii];
    *np=cat.np;

    cat.pos=(double *)malloc(3*cat.np*sizeof(double));
    if(cat.pos==NULL)
      error_mem_out_tpcf();
    fclose(snap);

    lint np_read=0;
      for(ii=0;ii<nfils;ii++) {
          lint np_new;
          lint jj;

          sprintf(fname,"%s.%d",prefix,(int)ii);
      snap=fopen(fname,"r");
        if(snap==NULL) error_open_file_tpcf(fname);
        verb_print_q(2, cmd->verbose,"  Reading file  %s \n",fname);

      //Read header
      if(input==2)
    gad_seek_block(snap,"HEAD");
      my_fread(&block1,sizeof(int),1,snap);
      my_fread(&head,sizeof(gad_header_tpcf),1,snap);
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);

      np_new=0;
      for(jj=0;jj<6;jj++)
    np_new+=head.npart[jj];
      printf("  %ld parts in file %ld \n",(long)np_new,(long)ii);

      if(np_read+np_new>cat.np) {
    fprintf(stderr,
        "CUTE: files seem to contain too many particles\n");
    fprintf(stderr,"      file %s, %ld > %ld \n",
        fname,(long)(np_read+np_new),(long)(cat.np));
    exit(1);
      }

      if(input==2)
    gad_seek_block(snap,"POS");
      my_fread(&block1,sizeof(int),1,snap);
      for(jj=np_read;jj<np_read+np_new;jj++) {
    float pos[3];
    my_fread(pos,sizeof(float),3,snap);
    cat.pos[3*jj]=(double)(pos[0]);
    cat.pos[3*jj+1]=(double)(pos[1]);
    cat.pos[3*jj+2]=(double)(pos[2]);
      }
      my_fread(&block2,sizeof(int),1,snap);
      gad_check_block(block1,block2);
      fclose(snap);

      np_read+=np_new;
    }
      
    if(np_read!=cat.np) {
      fprintf(stderr,
          "CUTE: #particles read disagrees with header: %ld != %ld\n",
          (long)np_read,(long)(cat.np));
      exit(1);
    }

    return cat;
  }
}


static double wrap_double(double x)
{
  //////
  // Returns x mod(l_box)
  if(x<0)
    return wrap_double(x+l_box_tpcf);
  else if(x>=l_box_tpcf)
    return wrap_double(x-l_box_tpcf);
  else
    return x;
}

lint linecount_tpcf(FILE *f)
{
  //////
  // Returns #lines in f
  lint i0=0;
  char ch[1024];
  while((fgets(ch,sizeof(ch),f))!=NULL) {
    i0++;
  }
  return i0;
}

void error_read_line_tpcf(char *fname,lint nlin)
{
  // Reading error handler
  fprintf(stderr,"CUTE: Couldn't read file %s, line %d \n",
      fname,(int)nlin);
  exit(1);
}

static Catalog_tpcf read_ascii(struct cmdline_data* cmd,
                               struct  global_data* gd,
                               char *fname,lint *np)
{
  // Reads catalog from ascii file with
  // default format.
  Catalog_tpcf cat;
  lint ii,n_lin;
  FILE *fd;

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file_tpcf(fname);
  if(n_objects_tpcf==-1)
    n_lin=linecount_tpcf(fd);
  else
    n_lin=n_objects_tpcf;
  rewind(fd);
      
    verb_print_q(2, cmd->verbose,"  %ld objects will be read \n",(long)n_lin);
      
  //Allocate catalog memory
  *np=n_lin;
  cat.np=n_lin;
  cat.pos=(double *)malloc(3*cat.np*sizeof(double));
  if(cat.pos==NULL)
    error_mem_out_tpcf();
  
  rewind(fd);
  //Read galaxies in mask
  for(ii=0;ii<n_lin;ii++) {
    double xx,yy,zz;
    char s0[1024];
    int sr;
    if(fgets(s0,sizeof(s0),fd)==NULL)
      error_read_line_tpcf(fname,ii+1);
    sr=sscanf(s0,"%lf %lf %lf",&xx,&yy,&zz);
    if(sr!=3)
      error_read_line_tpcf(fname,ii+1);
    cat.pos[3*ii]=xx;
    cat.pos[3*ii+1]=yy;
    cat.pos[3*ii+2]=zz;
  }
  fclose(fd);

  return cat;
}

int read_catalog_tpcf(struct cmdline_data* cmd,
                      struct  global_data* gd, int ifile,
                      char *fname,lint *np, particle_data **Ppointer)
{
  // Creates catalog from file fname
  lint ii;
  double x_mean=0,y_mean=0,z_mean=0;
  Catalog_tpcf cat;

    verb_print_q(2, cmd->verbose,"Reading catalog from file %s\n",fname);

    input_format_tpcf=1;

    if(input_format_tpcf) {
        cat=read_gadget(cmd, gd, fname,np,input_format_tpcf);
    } else {
        cat=read_ascii(cmd, gd, fname,np);
    }

    gd->nbodyTable[ifile] = cmd->nbody = *np;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    bodyptr p;

    //Correct particles out of bounds and calculate CoM
  for(ii=0;ii<cat.np;ii++) {
    double xx,yy,zz;
      p= bodytable[ifile] + ii;

    xx=cat.pos[3*ii];
    yy=cat.pos[3*ii+1];
    zz=cat.pos[3*ii+2];
    if((xx<0)||(xx>=l_box_tpcf)) xx=wrap_double(xx);
    if((yy<0)||(yy>=l_box_tpcf)) yy=wrap_double(yy);
    if((zz<0)||(zz>=l_box_tpcf)) zz=wrap_double(yy);
        Pos(p)[0] = xx;
        Pos(p)[1] = yy;
        Pos(p)[2] = zz;
    cat.pos[3*ii]=xx;
    cat.pos[3*ii+1]=yy;
    cat.pos[3*ii+2]=zz;
    x_mean+=xx/cat.np;
    y_mean+=yy/cat.np;
    z_mean+=zz/cat.np;
  }

    verb_print_q(2, cmd->verbose,
                 "  The center of mass is (%.3lf,%.3lf,%.3lf) \n\n",
                 x_mean,y_mean,z_mean);
  
    return 0;
}

//
//E !Cute_box
//
