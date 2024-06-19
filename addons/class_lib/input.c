/*==============================================================================
 * input.c module.
 *
 * Julien Lesgourgues, 27.08.2010
 *
 * Adapted to be used in cTreeBalls by Mario A. Rodriguez-Meza
==============================================================================*/
//        1          2          3          4          5          6          7


#include "globaldefs.h"
#include "input.h"

int input_find_file(char *fname,
                    struct file_content * fc,
                    ErrorMsg errmsg){

  struct file_content fc_input;
  struct file_content fc_precision;
  struct file_content * pfc_input;
  struct file_content fc_setroot;

  int i;
  char extension[5];
  char input_file[_ARGUMENT_LENGTH_MAX_];
  char precision_file[_ARGUMENT_LENGTH_MAX_];

  pfc_input = &fc_input;

  fc->size = 0;
  fc_input.size = 0;
  fc_precision.size = 0;
  input_file[0]='\0';
  precision_file[0]='\0';

    strcpy(input_file,fname);

    if (!strnull(fname)) {
        class_call(parser_read_file(fname,&fc_input,errmsg), errmsg, errmsg);
        class_call(input_set_root(input_file,&pfc_input,&fc_setroot,errmsg),
               errmsg, errmsg);
  }

  if ((input_file[0]!='\0') || (precision_file[0]!='\0')){
      class_call(parser_cat(pfc_input, 
                             &fc_precision,
                            fc, errmsg), errmsg, errmsg);
  }

  class_call(parser_free(pfc_input), errmsg, errmsg);
  class_call(parser_free(&fc_precision), errmsg, errmsg);

  return _SUCCESS_;
}

int input_set_root(char* input_file,
                   struct file_content** ppfc_input,
                   struct file_content * pfc_setroot,
                   ErrorMsg errmsg) {

  int flag1, filenum, iextens;
  int index_root_in_fc_input = -1;
  int overwrite_root;
  int found_filenum;

  FileArg outfname;

  char tmp_file[_ARGUMENT_LENGTH_MAX_+26];          // 26 is enough to extend
                                                    //  the file name [...] with
                                                    //  the characters
                                                // "output/[...]%02d_parameters.ini"
                                                    //  (as done below)
  struct file_content fc_root;                      // Temporary structure with
                                                    //  only the root name

  FileArg string1;                                  //Is ignored

//B Dummy so far...
  int n_extensions = 7;                             // Keep this as the length
                                                    //  of the below list
  char* output_extensions[7] = {"cl.dat","pk.dat","tk.dat","parameters.ini","background.dat","thermodynamics.dat","perturbations_k0.dat"};
//E

  struct file_content * pfc = *ppfc_input;

  class_call(parser_read_string(pfc,"rootDir",&string1,&flag1,errmsg),
             errmsg, errmsg);

    overwrite_root = _TRUE_;
    class_read_flag("overwrite_root",overwrite_root);

    if (flag1 == _FALSE_){
        memcpy(outfname, "Output/", 7);
        memcpy(outfname+7, input_file, strlen(input_file)-4);
        outfname[7+strlen(input_file)-4] = '\0';
    } else {
        for (index_root_in_fc_input=0;index_root_in_fc_input<pfc->size;
             ++index_root_in_fc_input) {
            if(strcmp(pfc->name[index_root_in_fc_input],"rootDir") == 0){
                strcpy(outfname,pfc->value[index_root_in_fc_input]);
                break;
            }
        }
    }

    if(overwrite_root == _FALSE_) {                     // Segment no use so far...
        found_filenum = _TRUE_;
        for (filenum = 0; filenum < _N_FILEROOT_ && found_filenum; filenum++){
            found_filenum = _FALSE_;
            for(iextens = 0; iextens < n_extensions; ++iextens){
                sprintf(tmp_file,"%s%02d_%s", outfname, filenum, 
                        output_extensions[iextens]);
                if (file_exists(tmp_file) == _TRUE_){
                    found_filenum = _TRUE_;
                }
            }
            if(found_filenum == _FALSE_){
                break;
            }
        }
        if(flag1 == _FALSE_){
            class_call(parser_init(&fc_root, 1, pfc->filename, errmsg),
                       errmsg,errmsg);
            sprintf(fc_root.name[0],"rootDir");
            sprintf(fc_root.value[0],"%s%02d_",outfname,filenum);
            fc_root.read[0] = _FALSE_;
            class_call(parser_cat(pfc, &fc_root, pfc_setroot, errmsg),
                       errmsg, errmsg);
            class_call(parser_free(pfc), errmsg, errmsg);
            class_call(parser_free(&fc_root), errmsg, errmsg);
            (*ppfc_input) = pfc_setroot;
        } else {
            sprintf(pfc->value[index_root_in_fc_input],"%s%02d_",outfname,filenum);
            (*ppfc_input) = pfc;
        }
    } else { // ! overwrite_root
        if(flag1 == _FALSE_){
            class_call(parser_init(&fc_root, 1, pfc->filename, errmsg),
                       errmsg,errmsg);
            sprintf(fc_root.name[0],"rootDir");
            sprintf(fc_root.value[0],"%s_",outfname);
            fc_root.read[0] = _FALSE_;
            class_call(parser_cat(pfc, &fc_root, pfc_setroot, errmsg),
                       errmsg, errmsg);
            class_call(parser_free(pfc), errmsg, errmsg);
            class_call(parser_free(&fc_root), errmsg, errmsg);
            (*ppfc_input) = pfc_setroot;
        } else {
            sprintf(pfc->value[index_root_in_fc_input],"%s",outfname);
            (*ppfc_input) = pfc;
        }
    }  // ! overwrite_root

  return _SUCCESS_;
}


int input_read_from_file(struct cmdline_data *cmd, struct file_content * pfc,
                         ErrorMsg errmsg) {

    int input_verbose = 0;

    class_read_int("verbose",input_verbose);
    verb_print(input_verbose, "Reading input parameters...\n");

    class_call(input_read_parameters(cmd, pfc, errmsg),errmsg,errmsg);

  return _SUCCESS_;
}


int input_read_parameters(struct cmdline_data *cmd, struct file_content * pfc,
                          ErrorMsg errmsg)
{
  int input_verbose=0;

  class_call(input_default_params(cmd),errmsg,errmsg);

  class_read_int("input_verbose",input_verbose);
  class_call(input_read_parameters_general(cmd, pfc,errmsg),errmsg,errmsg);

  return _SUCCESS_;
}


int input_read_parameters_general_free(struct file_content * pfc,
                                       ErrorMsg errmsg) {
    
}

int input_read_parameters_general(struct cmdline_data *cmd,
                                  struct file_content * pfc, ErrorMsg errmsg){

    int flag;
    int flag1,flag2;
    int param;
    int index;
    size_t slen;
  double param1,param2;
  char string1[_ARGUMENT_LENGTH_MAX_];

// All malloc have to be freed at the end of the run (EndRun)

    class_call(parser_read_string(pfc,"searchMethod",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"searchMethod") == 0){
              slen = strlen(pfc->value[index]);
              cmd->searchMethod = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->searchMethod,pfc->value[index],slen+1);
            break;
          }
        }
    }

  class_call(parser_read_double(pfc,"theta",&param1,&flag1,errmsg),errmsg,errmsg);
  if (flag1 == _TRUE_){
    cmd->theta = param1;
  }

    class_call(parser_read_string(pfc,"nsmooth",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"nsmooth") == 0){
              slen = strlen(pfc->value[index]);
              cmd->nsmooth = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->nsmooth,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"rsmooth",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"rsmooth") == 0){
              slen = strlen(pfc->value[index]);
              cmd->rsmooth = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->rsmooth,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"options",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"options") == 0){
              slen = strlen(pfc->value[index]);
              cmd->options = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->options,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"infile",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"infile") == 0){
              cmd->infile = strdup(pfc->value[index]);  // Only check that
                                                        //  is part of the current
                                                        //  ISO C or it is part
                                                        //  of POSIX...
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"infileformat",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"infileformat") == 0){
              slen = strlen(pfc->value[index]);
              cmd->infilefmt = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->infilefmt,pfc->value[index],slen+1);
            break;
          }
        }
    }

    class_call(parser_read_string(pfc,"iCatalogs",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"iCatalogs") == 0){
              slen = strlen(pfc->value[index]);
              cmd->iCatalogs = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->iCatalogs,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"rootDir",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"rootDir") == 0){
              slen = strlen(pfc->value[index]);
              cmd->rootDir = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->rootDir,pfc->value[index],slen+1);
            break;
          }
        }
    }

    class_call(parser_read_double(pfc,"rangeN",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->rangeN = param1;
    }
    class_call(parser_read_double(pfc,"rminHist",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->rminHist = param1;
    }
    class_call(parser_read_int(pfc,"sizeHistN",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->sizeHistN = param;

    class_call(parser_read_int(pfc,"mChebyshev",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->mchebyshev = param;

    class_call(parser_read_int(pfc,"stepState",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->stepState = param;

    class_call(parser_read_int(pfc,"verbose",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->verbose = param;
    class_call(parser_read_int(pfc,"verbose_log",&param,&flag,errmsg),
               errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->verbose_log = param;

    class_call(parser_read_int(pfc,"numberThreads",&param,&flag,errmsg),
               errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->numthreads = param;

    class_call(parser_read_string(pfc,"outfile",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"outfile") == 0){
              slen = strlen(pfc->value[index]);
              cmd->outfile = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->outfile,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"outfileformat",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"outfileformat") == 0){
              slen = strlen(pfc->value[index]);
              cmd->outfilefmt = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->outfilefmt,pfc->value[index],slen+1);
            break;
          }
        }
    }

    class_call(parser_read_double(pfc,"thetaL",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->thetaL = param1;
    }
    class_call(parser_read_double(pfc,"thetaR",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->thetaR = param1;
    }
    class_call(parser_read_double(pfc,"phiL",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->phiL = param1;
    }
    class_call(parser_read_double(pfc,"phiR",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->phiR = param1;
    }

    class_call(parser_read_string(pfc,"histNFileName",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histNFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histNFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histNFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"histXi2pcfFileName",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histXi2pcfFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histXi2pcfFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histXi2pcfFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"histZetaMFileName",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histZetaMFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histZetaMFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histZetaMFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"mhistZetaFileName",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"mhistZetaFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->mhistZetaFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->mhistZetaFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"suffixOutFiles",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"suffixOutFiles") == 0){
              slen = strlen(pfc->value[index]);
              cmd->suffixOutFiles = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->suffixOutFiles,pfc->value[index],slen+1);
            break;
          }
        }
    }

    class_call(parser_read_string(pfc,"script",&string1,&flag1,errmsg),
               errmsg,errmsg);
    char *script1;
    char *script2;
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"script") == 0){
              slen = strlen(pfc->value[index]);
              cmd->script = (char*) malloc((slen-2)*sizeof(char));
              script1 = (char*) malloc(slen*sizeof(char));
              memcpy(script1,pfc->value[index],slen);
              script2 = strchr(script1, '"');
              memcpy(cmd->script,script2+1,slen-2);
            break;
          }
        }
    }

    class_call(parser_read_int(pfc,"seed",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->seed = param;
    class_call(parser_read_string(pfc,"testmodel",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"testmodel") == 0){
              slen = strlen(pfc->value[index]);
              cmd->testmodel = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->testmodel,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_int(pfc,"nbody",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->nbody = param;
    class_call(parser_read_double(pfc,"lengthBox",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_){
      cmd->lengthBox = param1;
    }

    class_call(parser_read_string(pfc,"computeTPCF",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == _TRUE_) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->computeTPCF=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->computeTPCF=0;
    }

    class_call(parser_read_string(pfc,"useLogHist",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == _TRUE_) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->useLogHist=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->useLogHist=0;
    }
    class_call(parser_read_int(pfc,"logHistBinsPD",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->logHistBinsPD = param;

    class_call(parser_read_string(pfc,"usePeriodic",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == _TRUE_) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->usePeriodic=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->usePeriodic=0;
    }

#ifdef ADDONS
#ifdef BALLS
    class_call(parser_read_int(pfc,"scanLevel",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->scanLevel = param;
    class_call(parser_read_int(pfc,"scanLevelRoot",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->scanLevelRoot = param;

    class_call(parser_read_string(pfc,"scanLevelMin",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == _TRUE_) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"scanLevelMin") == 0){
              cmd->scanLevelMin = strdup(pfc->value[index]);
            break;
          }
        }
    }

    class_call(parser_read_int(pfc,"ntosave",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->ntosave = param;
#endif
#ifdef TREE3PCFDIRECTOMP
    class_call(parser_read_int(pfc,"sizeHistTheta",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == _TRUE_)
      cmd->sizeHistTheta = param;
#endif
#endif


  return _SUCCESS_;
}


//B cTreeBalls default values
int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

#ifdef MPICODE
    if (ThisTask==0) {                                  // Input all parameters on
                                                        //  proccess 0
#endif
    cmd->searchMethod = "tree-omp-sincos";
    cmd->theta = 1.0;
    cmd->mchebyshev = 8;
//B Parameters to set a region in the sky, for example for Takahasi data set.
    cmd->thetaL = 1.279928;
    cmd->thetaR = 1.861664;
    cmd->phiL = 1.280107;
    cmd->phiR = 1.861869;
//E
    cmd->sizeHistN = 40;
    cmd->rangeN = 100.0;
    cmd->rminHist = 1.0e-3;
    cmd->infile = "";
    cmd->infilefmt = "columns-ascii";
    cmd->iCatalogs = "1";
    cmd->rootDir = "output";
    cmd->outfile = "";
    cmd->outfilefmt = "columns-ascii";

    cmd->histNFileName = "histN";
    cmd->histXi2pcfFileName = "histXi2pcf";
    cmd->histZetaMFileName = "histZetaM";
    cmd->mhistZetaFileName = "mhistZeta";
    cmd->suffixOutFiles = "";

    cmd->stepState = 10000;

    cmd->verbose = 1;
    cmd->verbose_log = 1;

#ifdef OPENMPCODE
    cmd->numthreads = 4;
#endif

    cmd->script = "";
    cmd->options = "";

    cmd->seed=123;                                          // to always have
                                                            //  defaults Check in gsl
    cmd->nsmooth = "1";
    cmd->testmodel = "simple-cubic-random";
    cmd->nbody = 16384;
    cmd->lengthBox = 10000.0;

    cmd->rsmooth = "";

#ifdef ADDONS
#ifdef BALLS
//#include "startrun_include_02.h"
        cmd->scanLevel = 6;
        cmd->scanLevelRoot = 3;
        cmd->scanLevelMin = "-0";
        cmd->ntosave = 1000;
#endif
#ifdef TREE3PCFDIRECTOMP
        cmd->sizeHistTheta = 40;
#endif
#endif

    cmd->computeTPCF = 1;
    cmd->useLogHist = 1;
    cmd->logHistBinsPD = 5;
    cmd->usePeriodic = 0;

        
#ifdef MPICODE
    }
    MPI_Bcast(&cmd, sizeof(cmdline_data), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  return _SUCCESS_;
}
//E
