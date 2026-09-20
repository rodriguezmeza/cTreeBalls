
//=============================================================================
//        1          2          3          4        ^ 5          6          7

#include "globaldefs.h"
#include <stdint.h>

#ifdef CLASSLIB
#define GADGET_SET_ERROR(cmd, ...) \
    snprintf((cmd)->error_message, _ERRORMSGSIZE_, __VA_ARGS__)
#else
#define GADGET_SET_ERROR(cmd, ...) error(__VA_ARGS__)
#endif

//B cute_box
#ifdef _LONGIDS
typedef long lint;
#else //_LONGIDS
typedef INTEGER lint;
#endif //_LONGIDS

typedef struct {
  lint np;                                          //#objects in the catalog
  double *pos;
  double box_size;
} Catalog_tpcf;                                     //Catalog (double precision)

static Catalog_tpcf catalog_failure_tpcf(void)
{
    Catalog_tpcf cat;
    cat.np = -1;
    cat.pos = NULL;
    cat.box_size = 0.0;
    return cat;
}



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



static Catalog_tpcf read_gadget(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                char *prefix,lint *np,int input);
int read_catalog_tpcf(struct cmdline_data* cmd,
                      struct  global_data* gd, int ifile,
                      char *fname, lint *np);
//E cute_box


global int inputdata_gadget(struct cmdline_data* cmd,
                           struct  global_data* gd,
                           string filename, int ifile)
{
    bodyptr p;
    lint numpart;
    const int synthetic = scanopt(cmd->options, "gadget-kappa-synthetic");
    const int constant = scanopt(cmd->options, "kappa-constant");
    const int constant_one = scanopt(cmd->options, "kappa-constant-one");

    gd->input_comment = "Gadget positions; explicitly assigned scalar field";
    if (synthetic && (constant || constant_one)) {
        GADGET_SET_ERROR(cmd,
            "inputdata_gadget: gadget-kappa-synthetic conflicts with kappa-constant / kappa-constant-one");
        return FAILURE;
    }
    if (read_catalog_tpcf(cmd, gd, ifile, filename, &numpart) == FAILURE)
        return FAILURE;

    /* Gadget input reads positions only: no density or convergence field.
       The historical sinusoid is available only as an explicit test option. */
    verb_print(cmd->verbose, "inputdata_gadget: scalar field = %s\n",
        synthetic ? "synthetic sinusoid (gadget-kappa-synthetic)" :
        (constant && !constant_one ? "constant 2 (kappa-constant)" : "constant 1"));
    DO_BODY(p, bodytable[ifile], bodytable[ifile]+cmd->nbody) {
        Kappa(p) = constant && !constant_one ? 2.0 : 1.0;
        if (synthetic)
            Kappa(p) = 1.0 + rcos(40.0*PI*(Pos(p)[0]/gd->Box[0]))
                            * rsin(40.0*PI*(Pos(p)[1]/gd->Box[1]));
        Type(p) = BODY;
        Mass(p) = 1.0;
        Weight(p) = 1.0;
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


static int my_fread(struct cmdline_data *cmd, void *p, size_t size,
                    size_t nmemb, FILE *stream)
{
    string routineName = "my_fread";
    if (fread(p, size, nmemb, stream) != nmemb) {
        GADGET_SET_ERROR(cmd, "%s: error reading binary file\n", routineName);
        return FAILURE;
    }
    return SUCCESS;
}

static int gad_check_block(struct cmdline_data *cmd, int b1, int b2)
{
    string routineName = "gad_check_block";
    if (b1 != b2) {
        GADGET_SET_ERROR(cmd, "%s: Corrupted block!\n", routineName);
        return FAILURE;
    }
    return SUCCESS;
}


//static int gad_seek_block(FILE *snap,char name[])
static int gad_seek_block(struct cmdline_data *cmd, FILE *snap, char name[])
{
  // Seeks block from title
    string routineName = "gad_seek_block";
  gad_title_tpcf tit;
  int block1,block2;

  rewind(snap);

  while(1>0) {
    if(!(fread(&block1,sizeof(int),1,snap))||
       feof(snap)||ferror(snap)) {
        GADGET_SET_ERROR(cmd, "%s: Block %s not found!!\n", routineName, name);
        goto fail;
        
    }
      if (my_fread(cmd, &tit,sizeof(gad_title_tpcf),1,snap) == FAILURE)
          goto fail;
    
      if (my_fread(cmd, &block2,sizeof(int),1,snap) == FAILURE)
          goto fail;

      if (gad_check_block(cmd, block1, block2) == FAILURE)
          goto fail;
      
    if(strncmp(tit.label,name,3)!=0)
      fseek(snap,tit.size,SEEK_CUR);
    else
      break;
  }
    
    return SUCCESS;

fail:
    return FAILURE;
}
    
//}

static int check_num_files(char *prefix)
{
    string routineName = "check_num_files";
  FILE *fil;

  fil=fopen(prefix,"rb");
  if(fil!=NULL) {
    fclose(fil);
    return 1;
  }
  else {
    int nfils=0;
    while(nfils>=0) {
      char fname[256];
//      sprintf(fname,"%s.%d",prefix,nfils);
//        if (format_checked(fname, sizeof(fname),
//            "fname", "!%s/%s", "%s.%d",prefix,nfils) != 0)
//            return FAILURE;
        
        if (format_checked(fname, sizeof(fname), "fname", "%s.%d", prefix, nfils) != 0)
            return -1;
        
      fil=fopen(fname,"rb");
      if(fil!=NULL) {
    fclose(fil);
    nfils++;
      }
      else {
    if(nfils==0) {
      fprintf(stderr,"%s: can't find file %s or %s.x\n",
              routineName, prefix,prefix);
      return -1;
    }
    else if(nfils==1) {
      fprintf(stderr,"%s: only file %s found. Weird.\n",
              routineName, fname);
      return -1;
    }
    else {
      return nfils;
    }
      }
    }
  }
  
  fprintf(stderr,"%s: this shouldn't have happened \n", routineName);
  return -1;
}

static Catalog_tpcf read_snapshot_single(struct cmdline_data* cmd,
                                         struct  global_data* gd,
                                         char *fname,lint *np,int input)
{
    // Creates catalog from a single snapshot file
    string routineName = "read_snapshot_single";
    lint ii;
    gad_header_tpcf head;
    int block1,block2;
    
    Catalog_tpcf cat = catalog_failure_tpcf();
    FILE *snap = NULL;
    
     snap=fopen(fname,"rb");
    
    if (snap == NULL) {
        GADGET_SET_ERROR(cmd, "%s: Couldn't open file %s\n",
                         routineName, fname);
        goto fail;
    }
    
    //Read header
    if (input == 2) {
        if (gad_seek_block(cmd, snap, "HEAD") == FAILURE)
            goto fail;
    }
    if (my_fread(cmd, &block1,sizeof(int),1,snap) == FAILURE)
        goto fail;

    if (my_fread(cmd, &head,sizeof(gad_header_tpcf),1,snap) == FAILURE)
        goto fail;

    if (my_fread(cmd, &block2, sizeof(int), 1, snap) == FAILURE)
        goto fail;

    if (gad_check_block(cmd, block1, block2) == FAILURE)
        goto fail;
    

    if(head.num_files!=1) {
        GADGET_SET_ERROR(cmd,
            "%s: Multi-file input not expected \n", routineName);
        goto fail;
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

  if (!isfinite(head.BoxSize) || head.BoxSize <= 0.0) {
      GADGET_SET_ERROR(cmd, "%s: BoxSize must be finite and positive", routineName);
      goto fail;
  }
  cat.box_size = head.BoxSize;
  cat.np=0;
  for(ii=0;ii<6;ii++) {
    if (head.npart[ii] < 0) {
        GADGET_SET_ERROR(cmd, "%s: negative particle count", routineName);
        goto fail;
    }
    cat.np+=head.npart[ii];
    if(head.npart[ii]!=head.npartTotal[ii]) {
        GADGET_SET_ERROR(cmd,
                         "%s: error reading snapshot \n", routineName);
        goto fail;
    }
  }
  if (cat.np <= 0 || (size_t)cat.np > SIZE_MAX / (3*sizeof(double))) {
      GADGET_SET_ERROR(cmd, "%s: invalid catalog size", routineName);
      goto fail;
  }
  *np=cat.np;

  cat.pos=(double *)malloc(3*cat.np*sizeof(double));

    if (cat.pos==NULL) {
        GADGET_SET_ERROR(cmd, "%s: Out of memory!!\n", routineName);
        goto fail;
    }
    


  if(input==2)
      if (gad_seek_block(cmd, snap, "POS") == FAILURE)
          goto fail;
    if (my_fread(cmd, &block1, sizeof(int), 1, snap) == FAILURE)
        goto fail;

  for(ii=0;ii<cat.np;ii++) {
    float pos[3];
      if (my_fread(cmd, pos, sizeof(float), 3, snap) == FAILURE)
          goto fail;

    cat.pos[3*ii]=(double)(pos[0]);
    cat.pos[3*ii+1]=(double)(pos[1]);
    cat.pos[3*ii+2]=(double)(pos[2]);
  }
    if (my_fread(cmd, &block2,sizeof(int),1,snap) == FAILURE)
        goto fail;

    if (gad_check_block(cmd, block1, block2) == FAILURE)
        goto fail;

    fclose(snap);
    snap = NULL;
    return cat;

fail:
    if (snap != NULL) {
        fclose(snap);
        snap = NULL;
    }
    free(cat.pos);
    *np = 0;
    return catalog_failure_tpcf();
    
}

static Catalog_tpcf read_gadget(struct cmdline_data* cmd,
                                struct  global_data* gd,
                                char *prefix,lint *np,int input)
{
    string routineName = "read_gadget";
    Catalog_tpcf cat = catalog_failure_tpcf();
    FILE *snap = NULL;
    
  int nfils=check_num_files(prefix);
    if(nfils<=0) {
        GADGET_SET_ERROR(cmd,
            "%s: nfils <= 0\n", routineName);
        goto fail;
    }

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

      verb_print_q(2, cmd->verbose, "  Reading %d snapshot files \n",nfils);

//    sprintf(fname,"%s.0",prefix);
//      if (format_checked(fname, sizeof(fname),
//          "fname", "%s.0",prefix) != 0)
//          return FAILURE;
      if (format_checked(fname, sizeof(fname), "fname", "%s.0", prefix) != 0) {
          GADGET_SET_ERROR(cmd, "%s: filename too long for %s.0\n", routineName, prefix);
          goto fail;
      }
      
    snap=fopen(fname,"rb");
      if (snap == NULL) {
          GADGET_SET_ERROR(cmd, "%s: Couldn't open file %s\n",
                           routineName, fname);
          goto fail;
      }
      

    //Read header
    if(input==2)
        if (gad_seek_block(cmd, snap, "HEAD") == FAILURE)
            goto fail;
      if (my_fread(cmd, &block1,sizeof(int),1,snap) == FAILURE)
          goto fail;
      
      if (my_fread(cmd, &head,sizeof(gad_header_tpcf),1,snap) == FAILURE)
          goto fail;
      
      if (my_fread(cmd, &block2,sizeof(int),1,snap) == FAILURE)
          goto fail;

      if (gad_check_block(cmd, block1, block2) == FAILURE)
          goto fail;

    if(head.num_files!=nfils) {
        GADGET_SET_ERROR(cmd,
            "%s: Header and existing files do not match %d != %d.\n %s",
            routineName, nfils,head.num_files,
            "      There may be some files missing\n");
        goto fail;
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

    if (!isfinite(head.BoxSize) || head.BoxSize <= 0.0) {
        GADGET_SET_ERROR(cmd, "%s: BoxSize must be finite and positive", routineName);
        goto fail;
    }
    cat.box_size = head.BoxSize;
    cat.np=0;
    for(ii=0;ii<6;ii++) {
        if (head.npartTotal[ii] < 0) {
            GADGET_SET_ERROR(cmd, "%s: negative total particle count", routineName);
            goto fail;
        }
        cat.np+=head.npartTotal[ii];
    }
    if (cat.np <= 0 || (size_t)cat.np > SIZE_MAX / (3*sizeof(double))) {
        GADGET_SET_ERROR(cmd, "%s: invalid catalog size", routineName);
        goto fail;
    }
    *np=cat.np;

    cat.pos=(double *)malloc(3*cat.np*sizeof(double));
      if (cat.pos==NULL) {
          GADGET_SET_ERROR(cmd, "%s: Out of memory!!\n", routineName);
          goto fail;
      }
      
    fclose(snap);
      snap = NULL;

    lint np_read=0;
      for(ii=0;ii<nfils;ii++) {
          lint np_new;
          lint jj;

//          sprintf(fname,"%s.%d",prefix,(int)ii);
//          if (format_checked(fname, sizeof(fname),
//              "fname", "%s.%d",prefix,(int)ii) != 0)
//              return FAILURE;
          if (format_checked(fname, sizeof(fname), "fname", "%s.%d", prefix, (int)ii) != 0) {
              GADGET_SET_ERROR(cmd, "%s: filename too long for %s.%d\n", routineName, prefix, (int)ii);
              goto fail;
          }
          
      snap=fopen(fname,"rb");
          if (snap == NULL) {
              GADGET_SET_ERROR(cmd, "%s: Couldn't open file %s\n",
                               routineName, fname);
              goto fail;
          }
          
        verb_print_q(2, cmd->verbose,"  Reading file  %s \n",fname);

      //Read header
      if (input == 2) {
          if (gad_seek_block(cmd, snap, "HEAD") == FAILURE)
              goto fail;
      }
      if (my_fread(cmd, &block1,sizeof(int),1,snap) == FAILURE)
              goto fail;

          if (my_fread(cmd, &head,sizeof(gad_header_tpcf),1,snap) == FAILURE)
              goto fail;

            if (my_fread(cmd, &block2, sizeof(int), 1, snap) == FAILURE)
              goto fail;

          if (gad_check_block(cmd, block1, block2) == FAILURE)
              goto fail;
          

      if (!isfinite(head.BoxSize) || head.BoxSize != cat.box_size) {
          GADGET_SET_ERROR(cmd, "%s: inconsistent BoxSize in %s", routineName, fname);
          goto fail;
      }
      if (head.num_files != nfils) {
          GADGET_SET_ERROR(cmd, "%s: inconsistent num_files in %s", routineName, fname);
          goto fail;
      }
      np_new=0;
      for(jj=0;jj<6;jj++) {
          if (head.npart[jj] < 0) {
              GADGET_SET_ERROR(cmd, "%s: negative particle count in %s", routineName, fname);
              goto fail;
          }
          np_new+=head.npart[jj];
      }
      printf("  %ld parts in file %ld \n",(long)np_new,(long)ii);

      if(np_read+np_new>cat.np) {
          GADGET_SET_ERROR(cmd,
    "%s: files seem to contain too many particles\n     file %s, %ld > %ld \n",
            routineName, fname,(long)(np_read+np_new),(long)(cat.np));
          goto fail;
      }

      if (input == 2) {
          if (gad_seek_block(cmd, snap, "POS") == FAILURE)
              goto fail;
      }
      if (my_fread(cmd, &block1,sizeof(int),1,snap) == FAILURE)
              goto fail;

      for(jj=np_read;jj<np_read+np_new;jj++) {
    float pos[3];
          if (my_fread(cmd, pos,sizeof(float),3,snap) == FAILURE)
              goto fail;

    cat.pos[3*jj]=(double)(pos[0]);
    cat.pos[3*jj+1]=(double)(pos[1]);
    cat.pos[3*jj+2]=(double)(pos[2]);
      }
          if (my_fread(cmd, &block2,sizeof(int),1,snap) == FAILURE)
              goto fail;

          if (gad_check_block(cmd, block1, block2) == FAILURE)
              goto fail;
          
      fclose(snap);
          snap = NULL;

      np_read+=np_new;
    }
      
    if(np_read!=cat.np) {
        GADGET_SET_ERROR(cmd,
            "%s: #particles read disagrees with header: %ld != %ld\n",
            routineName, (long)np_read,(long)(cat.np));
        goto fail;
    }

    return cat;
      
  }

fail:
    if (snap != NULL) {
        fclose(snap);
        snap = NULL;
    }
    free(cat.pos);
    *np = 0;
    return catalog_failure_tpcf();

}


static double wrap_double(double x, double box_size)
{
    /* Constant-time wrapping also handles positions many boxes away.
       Rounding a tiny negative remainder can produce box_size exactly. */
    double wrapped = fmod(x, box_size);
    if (wrapped < 0.0) wrapped += box_size;
    return wrapped < box_size ? wrapped : 0.0;
}

int read_catalog_tpcf(struct cmdline_data* cmd,
                      struct  global_data* gd, int ifile,
                      char *fname,lint *np)
{
  // Creates catalog from file fname
  lint ii;
  double x_mean=0,y_mean=0,z_mean=0;
  Catalog_tpcf cat;

    verb_print_q(2, cmd->verbose,"Reading catalog from file %s\n",fname);

    cat = read_gadget(cmd, gd, fname, np, 1);
    if (cat.np < 0 || cat.pos == NULL)
        return FAILURE;
    for (ii = 0; ii < 3*cat.np; ii++) {
        if (!isfinite(cat.pos[ii])) {
            free(cat.pos);
            GADGET_SET_ERROR(cmd,
                "read_catalog_tpcf: non-finite coordinate at particle %ld", (long)(ii/3+1));
            return FAILURE;
        }
    }
    for (int k = 0; k < NDIM; k++) gd->Box[k] = cat.box_size;
    /* Startup later copies lengthBox back into gd->Box. Keep the snapshot
       header authoritative for both initialization and periodic searches. */
    cmd->lengthBox = cat.box_size;

    gd->nbodyTable[ifile] = cmd->nbody = *np;
    bodytable[ifile] = (bodyptr) allocate(cmd->nbody * sizeof(body));
    gd->bodytable_allocated = TRUE;
    bodyptr p;

    //Correct particles out of bounds and calculate CoM
  for(ii=0;ii<cat.np;ii++) {
    double xx,yy,zz;
      p= bodytable[ifile] + ii;

    xx=cat.pos[3*ii];
    yy=cat.pos[3*ii+1];
    zz=cat.pos[3*ii+2];
    xx=wrap_double(xx, cat.box_size);
    yy=wrap_double(yy, cat.box_size);
    zz=wrap_double(zz, cat.box_size);
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
  
    free(cat.pos);
    return SUCCESS;
}
