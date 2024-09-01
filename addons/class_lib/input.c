/*==============================================================================
 * input.c module.
 *
 * Julien Lesgourgues, 27.08.2010
 *
 * Adapted to be used in cTreeBalls by Mario A. Rodriguez-Meza
==============================================================================*/
//        1          2          3          4        ^ 5          6          7


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

    stream outstr;
    if (strnull(fname)) {
        verb_print(1, "If you intend to use a parameter file call with <ParameterFileName>\n");
        return FAILURE;
    } else {
        outstr = stropen(fname, "r");
        fclose(outstr);
    }

    strcpy(input_file,fname);

    if (!strnull(input_file)) {
        class_call(parser_read_file(input_file,&fc_input,errmsg), errmsg, errmsg);
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

  return SUCCESS;
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

//B To behave as not class_lib parameter file
    overwrite_root = TRUE;
    class_read_flag("overwrite_root",overwrite_root);
    overwrite_root = TRUE;
//E

    if (flag1 == FALSE){
        memcpy(outfname, "Output", 7);
//        memcpy(outfname, "Output/", 7);
//B To behave as not class_lib parameter file
//        memcpy(outfname+7, input_file, strlen(input_file)-4);
//        outfname[7+strlen(input_file)-4] = '\0';
//E
    } else {
        for (index_root_in_fc_input=0;index_root_in_fc_input<pfc->size;
             ++index_root_in_fc_input) {
            if(strcmp(pfc->name[index_root_in_fc_input],"rootDir") == 0){
                strcpy(outfname,pfc->value[index_root_in_fc_input]);
                break;
            }
        }
    }

    if(overwrite_root == FALSE) {                   // Segment no use so far...
        found_filenum = TRUE;
        for (filenum = 0; filenum < _N_FILEROOT_ && found_filenum; filenum++){
            found_filenum = FALSE;
            for(iextens = 0; iextens < n_extensions; ++iextens){
                sprintf(tmp_file,"%s%02d_%s", outfname, filenum, 
                        output_extensions[iextens]);
                if (file_exists(tmp_file) == TRUE){
                    found_filenum = TRUE;
                }
            }
            if(found_filenum == FALSE){
                break;
            }
        }
        if(flag1 == FALSE){
            class_call(parser_init(&fc_root, 1, pfc->filename, errmsg),
                       errmsg,errmsg);
            sprintf(fc_root.name[0],"rootDir");
            sprintf(fc_root.value[0],"%s%02d_",outfname,filenum);
            fc_root.read[0] = FALSE;
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
        if(flag1 == FALSE){
            class_call(parser_init(&fc_root, 1, pfc->filename, errmsg),
                       errmsg,errmsg);
            sprintf(fc_root.name[0],"rootDir");
//B To behave as not class_lib parameter file
//            sprintf(fc_root.value[0],"%s_",outfname);
            sprintf(fc_root.value[0],"%s",outfname);
//E
            fc_root.read[0] = FALSE;
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

  return SUCCESS;
}


int input_read_from_file(struct cmdline_data *cmd, struct file_content * pfc,
                         ErrorMsg errmsg) {

    int input_verbose = 0;

    class_read_int("verbose",input_verbose);
    verb_print(input_verbose, "Reading input parameters...\n");

    class_call(input_read_parameters(cmd, pfc, errmsg),errmsg,errmsg);

    return SUCCESS;
}


int input_read_parameters(struct cmdline_data *cmd, struct file_content * pfc,
                          ErrorMsg errmsg)
{
  int input_verbose=0;

  class_call(input_default_params(cmd),errmsg,errmsg);

  class_read_int("input_verbose",input_verbose);
  class_call(input_read_parameters_general(cmd, pfc,errmsg),errmsg,errmsg);

  return SUCCESS;
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

    if (flag1 == TRUE) {
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
  if (flag1 == TRUE){
    cmd->theta = param1;
  }

    class_call(parser_read_string(pfc,"nsmooth",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE){
      cmd->rangeN = param1;
    }
    class_call(parser_read_double(pfc,"rminHist",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->rminHist = param1;
    }
    class_call(parser_read_int(pfc,"sizeHistN",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistN = param;

    class_call(parser_read_int(pfc,"mChebyshev",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->mChebyshev = param;

    class_call(parser_read_int(pfc,"sizeHistTheta",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistTheta = param;

    class_call(parser_read_int(pfc,"stepState",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->stepState = param;

    class_call(parser_read_int(pfc,"verbose",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->verbose = param;
    class_call(parser_read_int(pfc,"verbose_log",&param,&flag,errmsg),
               errmsg,errmsg);
    if (flag == TRUE)
      cmd->verbose_log = param;

#ifdef OPENMPCODE
    class_call(parser_read_int(pfc,"numberThreads",&param,&flag,errmsg),
               errmsg,errmsg);
    if (flag == TRUE)
      cmd->numthreads = param;
#endif

    class_call(parser_read_string(pfc,"outfile",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"outfile") == 0){
              slen = strlen(pfc->value[index]);
              cmd->outfile = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->outfile,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"outfileformat",&string1,&flag1,
                                  errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE){
      cmd->thetaL = param1;
    }
    class_call(parser_read_double(pfc,"thetaR",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->thetaR = param1;
    }
    class_call(parser_read_double(pfc,"phiL",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->phiL = param1;
    }
    class_call(parser_read_double(pfc,"phiR",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->phiR = param1;
    }

    class_call(parser_read_string(pfc,"histNNFileName",&string1,&flag1,
                                  errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histNNFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histNNFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histNNFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"histXi2pcfFileName",&string1,&flag1,
                                  errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histXi2pcfFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histXi2pcfFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histXi2pcfFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"histZetaFileName",&string1,&flag1,
                                  errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"histZetaFileName") == 0){
              slen = strlen(pfc->value[index]);
              cmd->histZetaFileName = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->histZetaFileName,pfc->value[index],slen+1);
            break;
          }
        }
    }

    class_call(parser_read_string(pfc,"suffixOutFiles",&string1,&flag1,
                                  errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
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
    if (flag1 == TRUE) {
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
    if (flag == TRUE)
      cmd->seed = param;
    class_call(parser_read_string(pfc,"testmodel",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
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
    if (flag == TRUE)
      cmd->nbody = param;
    class_call(parser_read_double(pfc,"lengthBox",&param1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->lengthBox = param1;
    }

    class_call(parser_read_string(pfc,"computeTPCF",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->computeTPCF=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->computeTPCF=0;
    }

    class_call(parser_read_string(pfc,"useLogHist",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->useLogHist=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->useLogHist=0;
    }
    class_call(parser_read_int(pfc,"logHistBinsPD",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->logHistBinsPD = param;

    class_call(parser_read_string(pfc,"usePeriodic",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->usePeriodic=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->usePeriodic=0;
    }

/*
#ifdef ADDONS
#ifdef BALLS
    class_call(parser_read_int(pfc,"scanLevel",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->scanLevel = param;
    class_call(parser_read_int(pfc,"scanLevelRoot",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->scanLevelRoot = param;

    class_call(parser_read_string(pfc,"scanLevelMin",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"scanLevelMin") == 0){
              cmd->scanLevelMin = strdup(pfc->value[index]);
            break;
          }
        }
    }

    class_call(parser_read_int(pfc,"ntosave",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->ntosave = param;
#endif

#ifdef TREE3PCFDIRECTOMP
#endif

#ifdef IOLIB
    class_call(parser_read_string(pfc,"columns",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"columns") == 0){
              slen = strlen(pfc->value[index]);
              cmd->columns = (char*) malloc(slen*sizeof(char));
            memcpy(cmd->columns,pfc->value[index],slen+1);
            break;
          }
        }
    }
#endif
#endif
*/
#ifdef ADDONS
#include "class_lib_include_01.h"
#endif

    class_call(parser_read_int(pfc,"sizeHistTheta",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistTheta = param;


  return SUCCESS;
}


//B cTreeBalls default values
int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

    cmd->searchMethod = "tree-omp-sincos";
    cmd->theta = 1.0;
    cmd->mChebyshev = 7;
    cmd->sizeHistTheta = 32;
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

    cmd->histNNFileName = "histNN";
    cmd->histXi2pcfFileName = "histXi2pcf";
    cmd->histZetaFileName = "histZeta";
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
/*
#ifdef ADDONS
#ifdef BALLS
        cmd->scanLevel = 6;
        cmd->scanLevelRoot = 3;
        cmd->scanLevelMin = "-0";
        cmd->ntosave = 1000;
#endif
#ifdef IOLIB
        // pos, kappa, gamma1, gamma2, weight, seven places maximum
        // use in multi-columns-ascii and cfitsio
        // default is 3D: three pos and one convergence (kappa)
        cmd->columns = "1,2,3,4";
#endif
#ifdef TREE3PCFDIRECTOMP
#endif
#endif
*/
#ifdef ADDONS
#include "class_lib_include_02.h"
#endif


    cmd->sizeHistTheta = 32;

    cmd->computeTPCF = 1;
    cmd->useLogHist = 1;
    cmd->logHistBinsPD = 5;
    cmd->usePeriodic = 0;

  return SUCCESS;
}
//E
