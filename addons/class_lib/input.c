/*==============================================================================
 * input.c module.
 *
 * Julien Lesgourgues, 27.08.2010
 *
 * Adapted to be used in cTreeBalls by Mario A. Rodriguez-Meza
==============================================================================*/
//        1          2          3          4        ^ 5          6          7

//
// lines where there is a "//B socket:" string are places to include module files
//  that can be found in addons/addons_include folder
//

#include "globaldefs.h"
#include "input.h"

local void testParameterFile(struct  cmdline_data*,
                             struct  global_data*,
                             char *);


int input_find_file(struct  cmdline_data* cmd, struct  global_data* gd,
                    char *fname,
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

    //B test if cmd->paramfile exist...
    if (!strnull(fname)) {
        testParameterFile(cmd, gd, cmd->paramfile);
    }
    //E

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

//  char tmp_file[_ARGUMENT_LENGTH_MAX_+26];          // 26 is enough to extend
                                                    //  the file name [...] with
                                                    //  the characters
                                                // "output/[...]%02d_parameters.ini"
                                                    //  (as done below)
  struct file_content fc_root;                      // Temporary structure with
                                                    //  only the root name

  FileArg string1;                                  //Is ignored

//B Dummy so far...
//  int n_extensions = 7;                             // Keep this as the length
                                                    //  of the below list
//  char* output_extensions[7] = {"cl.dat","pk.dat","tk.dat","parameters.ini","background.dat","thermodynamics.dat","perturbations_k0.dat"};
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
/*
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
    } else { // ! overwrite_root */
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
//    }  // ! overwrite_root

  return SUCCESS;
}


int input_read_from_file(struct cmdline_data *cmd,
                         struct file_content * pfc,
                         ErrorMsg errmsg)
{

    int input_verbose = 0;

    class_read_int("verbose",input_verbose);
//    if (!strnull(cmd->paramfile)) {
//        testParameterFile(cmd, cmd->paramfile);
//    }
    verb_print(input_verbose, "\nReading input parameters...\n");

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
                                  struct file_content * pfc, ErrorMsg errmsg)
{

    int flag;
    int flag1,flag2;
    int param;
    int index;
    size_t slen;
    double param1,param2;
    char string1[_ARGUMENT_LENGTH_MAX_];

// All malloc have to be freed at the end of the run (EndRun)

    //B Parameters related to the searching method
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

    class_call(parser_read_int(pfc,"mChebyshev",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->mChebyshev = param;

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

    class_call(parser_read_double(pfc,"theta",&param1,&flag1,errmsg),errmsg,errmsg);
    if (flag1 == TRUE){
      cmd->theta = param1;
    }

    class_call(parser_read_string(pfc,"computeTPCF",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->computeTPCF=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->computeTPCF=0;
    }

    //B correction 2025-04-06
    // Move this to addon that compute shear correlations
//    class_call(parser_read_string(pfc,"computeShearCF",&string1,&flag1,errmsg),
//             errmsg,errmsg);
//    if (flag1 == TRUE) {
//        if (strchr("tTyY1", *string1) != NULL)
//            cmd->computeShearCF=1;
//        if (strchr("fFnN0", *string1) != NULL)
//            cmd->computeShearCF=0;
//    }
    //E

    class_call(parser_read_string(pfc,"usePeriodic",&string1,&flag1,errmsg),
             errmsg,errmsg);
    if (flag1 == TRUE) {
        if (strchr("tTyY1", *string1) != NULL)
            cmd->usePeriodic=1;
        if (strchr("fFnN0", *string1) != NULL)
            cmd->usePeriodic=0;
    }
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
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
    // Output parameters
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
    // Parameters to set a region in the sky, for example for Takahasi data set
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
    //E

    //B Parameters to control histograms and their output files
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
    //
    class_call(parser_read_int(pfc,"sizeHistN",&param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistN = param;
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

    class_call(parser_read_int(pfc,"sizeHistPhi",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistPhi = param;
    //
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
    //E

    //B Set of parameters needed to construct a test model
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
    //E

    //B Miscellaneous parameters
    class_call(parser_read_string(pfc,"preScript",&string1,&flag1,errmsg),
               errmsg,errmsg);
    char *script1;
    char *script2;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"preScript") == 0){
              slen = strlen(pfc->value[index]);
              cmd->preScript = (char*) malloc((slen-2)*sizeof(char));
              script1 = (char*) malloc(slen*sizeof(char));
              memcpy(script1,pfc->value[index],slen);
              script2 = strchr(script1, '"');
              memcpy(cmd->preScript,script2+1,slen-2);
              free(script1);
            break;
          }
        }
    }
    class_call(parser_read_string(pfc,"posScript",&string1,&flag1,errmsg),
               errmsg,errmsg);
//    char *script1;
//    char *script2;
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"posScript") == 0){
              slen = strlen(pfc->value[index]);
              cmd->posScript = (char*) malloc((slen-2)*sizeof(char));
              script1 = (char*) malloc(slen*sizeof(char));
              memcpy(script1,pfc->value[index],slen);
              script2 = strchr(script1, '"');
              memcpy(cmd->posScript,script2+1,slen-2);
              free(script1);
            break;
          }
        }
    }

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

    class_call(parser_read_string(pfc,"options",&string1,&flag1,errmsg),
               errmsg,errmsg);
    if (flag1 == TRUE) {
        for (index=0;index<pfc->size;++index){
          if (strcmp(pfc->name[index],"options") == 0){
//B The use of memcpy to copy gives bad string on string options
//      when we use parser_free(&fc) in StartRun
//      memcpy is faster than strcpy...
//              slen = strlen(pfc->value[index]);
//              cmd->options = (char*) malloc(slen*sizeof(char));
//              memcpy(cmd->options,pfc->value[index],slen+1);
              cmd->options = (char*) malloc(MAXLENGTHOFSTRSCMD);
              strcpy(cmd->options,pfc->value[index]);
//E
            break;
          }
        }
    }
    //E

//B socket:
#ifdef ADDONS
#include "class_lib_include_01.h"
#endif
//

    class_call(parser_read_int(pfc,"sizeHistPhi",
                               &param,&flag,errmsg),errmsg,errmsg);
    if (flag == TRUE)
      cmd->sizeHistPhi = param;


  return SUCCESS;
}


//B cTreeBalls default values
#ifndef CMDLINE_DEFS_UNITSPHERE

int input_default_params(struct cmdline_data *cmd)
{
// Every item in cmdline_defs.h must have an item here::

    //B Parameters related to the searching method
    cmd->searchMethod = "tree-omp-sincos";
    cmd->mChebyshev = 7;
    cmd->nsmooth = "1";
    cmd->rsmooth = "";
    cmd->theta = 1.0;
    cmd->computeTPCF = 1;
    //B correction 2025-04-06
    // Move this to addon that computes shear correlations
//    cmd->computeShearCF = 0;
    //E
    cmd->usePeriodic = 0;
    //E

    //B Parameters about the I/O file(s)
    // Input catalog parameters
    cmd->infile = "";
    cmd->infilefmt = "columns-ascii";
    cmd->iCatalogs = "1";
    // Output parameters
    cmd->rootDir = "Output";
    cmd->outfile = "";
    cmd->outfilefmt = "columns-ascii";
    // Parameters to set a region in the sky, for example for Takahasi data set
    cmd->thetaL = 1.279928;
    cmd->thetaR = 1.861664;
    cmd->phiL = 1.280107;
    cmd->phiR = 1.861869;
    //E

    //B Parameters to control histograms and their output files
    cmd->useLogHist = 1;
    cmd->logHistBinsPD = 5;
    //
    cmd->sizeHistN = 40;
    cmd->rangeN = 100.0;
    cmd->rminHist = 1.0e-3;
    cmd->sizeHistPhi = 32;
    //
    cmd->histNNFileName = "histNN";
    cmd->histXi2pcfFileName = "histXi2pcf";
    cmd->histZetaFileName = "histZeta";
    cmd->suffixOutFiles = "";
    //E

    //B Set of parameters needed to construct a test model
    cmd->seed=123;                                          // to always have
                                                            //  defaults Check in gsl
    cmd->testmodel = "simple-cubic-random";
    cmd->nbody = 16384;
    cmd->lengthBox = 10000.0;
    //E

    //B Miscellaneous parameters
//    cmd->script = "";
    cmd->preScript = "";
    cmd->posScript = "";
    cmd->stepState = 10000;
    cmd->verbose = 0;
    cmd->verbose_log = 0;
#ifdef OPENMPCODE
    cmd->numthreads = 4;
#endif
    cmd->options = "";
    //E

//B socket:
#ifdef ADDONS
#include "class_lib_include_02.h"
#endif
//

  return SUCCESS;
}

#else
#include "input_default_params.h"
#endif

//E


//B parameter reading/testing from a file
local void testParameterFile(struct  cmdline_data* cmd,
                             struct  global_data* gd,
                             char *fname)
{
// Every item in cmdline_defs.h must have an item here::
#define DOUBLE 1
#define STRING 2
#define INT 3
#define LONG 6
#define BOOLEAN 4
#define MAXTAGS 300
#define MAXCHARBUF 1024

    FILE *fd;

  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

//    printf("\nTesting input parameters file: %s...\n", fname);
//    if (strnull(fname)) return;
    int input_verbose = 2;
    verb_print(input_verbose, "\nparsing input parameters...\n");
//    verb_print_normal_info(input_verbose, 0, gd->outlog,
//                           "\nparsing input parameters...\n");

    nt=0;

    //B Parameters related to the searching method
    SPName(cmd->searchMethod,"searchMethod",MAXLENGTHOFSTRSCMD);
    IPName(cmd->mChebyshev,"mChebyshev");
    SPName(cmd->nsmooth,"nsmooth",MAXLENGTHOFSTRSCMD);
    SPName(cmd->rsmooth,"rsmooth",MAXLENGTHOFSTRSCMD);
    RPName(cmd->theta,"theta");
    BPName(cmd->computeTPCF,"computeTPCF");
    //B correction 2025-04-06
    // Move this to addon that compute shear correlations
//    BPName(cmd->computeShearCF,"computeShearCF");
    //E
    BPName(cmd->usePeriodic,"usePeriodic");
    //E

    //B Parameters to control the I/O file(s)
    // Input catalog parameters
    SPName(cmd->infile,"infile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->infilefmt,"infileformat",MAXLENGTHOFSTRSCMD);
    SPName(cmd->iCatalogs,"iCatalogs",MAXLENGTHOFSTRSCMD);
    // Output parameters
    SPName(cmd->rootDir,"rootDir",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfile,"outfile",MAXLENGTHOFSTRSCMD);
    SPName(cmd->outfilefmt,"outfileformat",MAXLENGTHOFSTRSCMD);
    //B Parameters to set a region in the sky, for example for Takahasi data set.
    RPName(cmd->thetaL,"thetaL");
    RPName(cmd->thetaR,"thetaR");
    RPName(cmd->phiL,"phiL");
    RPName(cmd->phiR,"phiR");
    //E

    //B Parameters to control histograms and their output files
    BPName(cmd->useLogHist,"useLogHist");
    IPName(cmd->logHistBinsPD,"logHistBinsPD");
    //
    IPName(cmd->sizeHistN,"sizeHistN");
    RPName(cmd->rangeN,"rangeN");
    RPName(cmd->rminHist,"rminHist");
    IPName(cmd->sizeHistPhi,"sizeHistPhi");
    //
    SPName(cmd->histNNFileName,"histNNFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histXi2pcfFileName,"histXi2pcfFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->histZetaFileName,"histZetaFileName",MAXLENGTHOFSTRSCMD);
    SPName(cmd->suffixOutFiles,"suffixOutFiles",MAXLENGTHOFSTRSCMD);
    //E

    //B Set of parameters needed to construct a test model
    IPName(cmd->seed,"seed");                       // to always have
                                                    //  defaults
    SPName(cmd->testmodel,"testmodel",MAXLENGTHOFSTRSCMD);
#ifdef LONGINT
    LPName(cmd->nbody,"nbody");
#else
    IPName(cmd->nbody,"nbody");
#endif
    RPName(cmd->lengthBox,"lengthBox");
    //E

    //B Miscellaneous parameters
    SPName(cmd->preScript,"preScript",MAXLENGTHOFSTRSCMD);
    SPName(cmd->posScript,"posScript",MAXLENGTHOFSTRSCMD);
#ifdef LONGINT
    LPName(cmd->stepState,"stepState");
#else
    IPName(cmd->stepState,"stepState");
#endif
    IPName(cmd->verbose,"verbose");
    IPName(cmd->verbose_log,"verbose_log");
#ifdef OPENMPCODE
    IPName(cmd->numthreads,"numberThreads");
#endif
    SPName(cmd->options,"options",MAXLENGTHOFSTRSCMD);
    //E
    
//B socket:
#ifdef ADDONS
#include "startrun_include_03.h"
#endif
//E

    size_t slen;
    char *script1;
    char *script2;
    char *script3;
    char *script4;
    char buf4[MAXCHARBUF];
    char buf5[MAXCHARBUF];

//B
#ifndef _LINE_LENGTH_MAX_
#define _LINE_LENGTH_MAX_ 1024
#endif
#define _ARGUMENT_LENGTH_MAX_ 1024
        char line[_LINE_LENGTH_MAX_];
        char name[_ARGUMENT_LENGTH_MAX_];
        char value[_ARGUMENT_LENGTH_MAX_];
        char * phash;
        char * pequal;
        char * left;
        char * right;
//E

    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
//B
            fgets(line,MAXCHARBUF,fd);

            pequal=strchr(line,'=');
            if (pequal == NULL)
                continue;
            phash=strchr(line,'#');
            if ((phash != NULL) && (phash-pequal<2))
                continue;
            phash=strchr(line,'%');
            if ((phash != NULL) && (phash-pequal<2))
                continue;

            left=line;
            while (left[0]==' ') {
              left++;
            }
            if(left[0]=='\'' || left[0]=='\"'){
              left++;
            }
            right=pequal-1;
            while (right[0]==' ') {
              right--;
            }
            if(right[0]=='\'' || right[0]=='\"'){
              right--;
            }

            if (right-left < 0) {
                fprintf(stdout,
        "Error in file %s: there is no variable name before '=' in line: '%s'\n",
                    fname, line);
                errorFlag=1;
                continue;
            }

            strncpy(name,left,right-left+1);
            name[right-left+1]='\0';

            left = pequal+1;
            while (left[0]==' ') {
              left++;
            }

            if (phash == NULL)
              right = line+strlen(line)-1;
            else
              right = phash-1;

            while (right[0]<=' ') {
              right--;
            }

            if (right-left < 0)
                continue;

            strncpy(value,left,right-left+1);
            value[right-left+1]='\0';
//E

            for(i=0,j=-1;i<nt;i++)
                if(strcmp(name,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(value);
                        break;
                    case STRING:
                        if (strcmp(name,"preScript") == 0){ // To remove both '"'
                            int index;
                            size_t slen;
                            slen = strlen(value);
                            //B
                            script1 = (char*) malloc(1*sizeof(char));
                            memcpy(script1,value,1);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                      value);
                            free(script1);
                            //E
                            //B
                            script1 = (char*) malloc((slen-1)*sizeof(char));
                            memcpy(script1,value+1,slen-1);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                      value);
                            free(script1);
                            //E
                            cmd->preScript = (char*) malloc((slen-2)*sizeof(char));
                            script1 = (char*) malloc(slen*sizeof(char));
                            memcpy(script1,value,slen);
                            script2 = strchr(script1, '"');
                            if (script2 == NULL)
                                error("preScript parameter needs enclosing script with \"\"!!\n\n");
                            memcpy(cmd->preScript,script2+1,slen-2);
                            free(script1);
                        } else {
                            if (strcmp(name,"posScript") == 0){ // To remove both '"'
                                int index;
                                size_t slen;
                                slen = strlen(value);
                                //B
                                script1 = (char*) malloc(1*sizeof(char));
                                memcpy(script1,value,1);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          value);
                                free(script1);
                                //E
                                //B
                                script1 = (char*) malloc((slen-1)*sizeof(char));
                                memcpy(script1,value+1,slen-1);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          value);
                                free(script1);
                                //E
                                cmd->posScript = (char*) malloc((slen-2)*sizeof(char));
                                script1 = (char*) malloc(slen*sizeof(char));
                                memcpy(script1,value,slen);
                                script2 = strchr(script1, '"');
                                if (script2 == NULL)
                                    error("posScript parameter needs enclosing script with \"\"!! (%s)\n\n",
                                          script1);
                                memcpy(cmd->posScript,script2+1,slen-2);
                                free(script1);
                            } else
                                strcpy(addr[j],value);
                        }
                        break;
                    case INT:
                        *((int*)addr[j])=atoi(value);
                        break;
                    case LONG:
                        *((long*)addr[j])=atol(value);
                        break;
                    case BOOLEAN:
                        if (strchr("tTyY1", *value) != NULL) {
                            *((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *value) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",name,value);
                            }
                        break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s\n",
                        fname, name,
                        "not allowed or multiple defined...\n");
                errorFlag=1;
            }
        } // ! while loop
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=2;
        exit(0);
    }

    if (errorFlag==1)
        error("\ntestParameterFile: going out\n");

    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            if (cmd->verbose>2)
                fprintf(stdout,
                "Warning! I miss a value for tag '%s' in parameter file '%s'.\n",
                    tag[i],fname);
            switch(id[i]) {
                case DOUBLE:
                    *((double*)addr[i])=GetdParam(tag[i]);
                    break;
                case STRING:
                    strcpy(addr[i],GetParam(tag[i]));
                    break;
                case INT:
                    *((int*)addr[i])=GetiParam(tag[i]);
                    break;
                case LONG:
                    *((long*)addr[i])=GetlParam(tag[i]);
                    break;
                case BOOLEAN:
                    *((bool*)addr[i])=GetbParam(tag[i]);
                    break;
            }
            errorFlag=3;
        }
    }

#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
#undef MAXCHARBUF
}
//E
