#ifndef CBALLS_CLI_CONTRACTS_H
#define CBALLS_CLI_CONTRACTS_H
#include "input.h"
#include "parser.h"

/* The active CLASSLIB CLI uses the same checked parser as Python. In particular,
 * malformed rank-local parameters must not enter legacy Get*Param exit paths. */
static int cballs_cli_assign(struct file_content *fc, const char *text, char *error)
{
    const char *equal = strchr(text, '=');
    if (!equal || equal == text || (size_t)(equal-text) >= sizeof(FileArg)
        || strlen(equal+1) >= sizeof(FileArg)) {
        snprintf(error, _ERRORMSGSIZE_, "invalid or overlong CLI parameter: %.200s", text);
        return FAILURE;
    }
    char name[sizeof(FileArg)];
    memcpy(name, text, equal-text); name[equal-text] = '\0';
    const char *canonical = name, *previous = NULL;
    for (int i=0; defv[i]; i++) {
        if (defv[i][0] == ':' && !strcmp(defv[i]+1,name) && previous) {
            size_t length = strcspn(previous,"=");
            memcpy(name,previous,length); name[length]='\0'; break;
        }
        if (defv[i][0] != ';' && defv[i][0] != ':' && strchr(defv[i],'=')) previous=defv[i];
    }
    for (int i=0; i<fc->size; i++) if (!strcmp(fc->name[i], canonical)) {
        strcpy(fc->value[i],equal+1); fc->read[i]=FALSE; return SUCCESS;
    }
    snprintf(error, _ERRORMSGSIZE_, "unknown CLI parameter: %.200s", name);
    return FAILURE;
}
static int cballs_cli_input(struct cmdline_data *cmd, struct global_data *gd,
                            int argc, char **argv)
{
    struct file_content fc={0}, file={0};
    int count=0, status=FAILURE, index=0;
    const char *parameter_file="";
    for (int i=0; defv[i]; i++)
        if (defv[i][0] != ';' && defv[i][0] != ':' && strchr(defv[i],'=')) count++;
    if (parser_init(&fc,count,"command line",cmd->error_message)==FAILURE) goto done;
    for (int i=0; defv[i]; i++) {
        const char *equal=strchr(defv[i],'=');
        if (defv[i][0]==';' || defv[i][0]==':' || !equal) continue;
        size_t length=equal-defv[i];
        if(length>=sizeof(FileArg) || strlen(equal+1)>=sizeof(FileArg)) goto done;
        memcpy(fc.name[index],defv[i],length); fc.name[index][length]='\0';
        strcpy(fc.value[index],equal+1); fc.read[index++]=FALSE;
    }
    for (int i=1; i<argc; i++) {
        if (!strncmp(argv[i],"paramfile=",10)) parameter_file=argv[i]+10;
        else if (!strchr(argv[i],'=') && argv[i][0]!='-' && i==1) parameter_file=argv[i];
    }
    if (*parameter_file) {
        if (parser_read_file((char *)parameter_file,&file,cmd->error_message)==FAILURE) goto done;
        for (int i=0; i<file.size; i++) {
            char assignment[2*sizeof(FileArg)+2];
            snprintf(assignment,sizeof(assignment),"%s=%s",file.name[i],file.value[i]);
            if (cballs_cli_assign(&fc,assignment,cmd->error_message)==FAILURE) goto done;
        }
    }
    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i],parameter_file) || !strncmp(argv[i],"paramfile=",10)) continue;
        if (cballs_cli_assign(&fc,argv[i],cmd->error_message)==FAILURE) goto done;
    }
    for (int j=0; j<fc.size; j++)
        if ((!strcmp(fc.name[j],"preScript") || !strcmp(fc.name[j],"posScript"))
            && fc.value[j][0] != '"') {
            char script[sizeof(FileArg)];
            if (strlen(fc.value[j])+3 > sizeof(script)) goto done;
            snprintf(script,sizeof(script),"\"%s\"",fc.value[j]);
            strcpy(fc.value[j],script);
        }
    char detail[_ERRORMSGSIZE_];
    status=input_read_from_file_guarded(cmd,gd,&fc,detail);
    if (status == FAILURE) snprintf(cmd->error_message,_ERRORMSGSIZE_,"%s",detail);
 done:
    parser_free(&file); parser_free(&fc);
    return status;
}
static int cballs_cli_start(struct cmdline_data *cmd, struct global_data *gd,
                            int argc, char **argv)
{
    const char *method="MPI startup";
    for (int i=1; i<argc; i++) if (!strncmp(argv[i],"searchMethod=",13)) method=argv[i]+13;
    if (cballs_mpi_bootstrap(cmd,gd,method)==FAILURE) return FAILURE;
    int status=cballs_cli_input(cmd,gd,argc,argv);
    status=cballs_mpi_context_consensus(cmd,status,"MPI CLI input preflight");
    if(status==FAILURE) return FAILURE;
    gd->headline0=argv[0]; gd->headline1=HEAD1; gd->headline2=HEAD2; gd->headline3=HEAD3;
    if(cballs_start_run_common_guarded(cmd,gd)==FAILURE) return FAILURE;
    if(cballs_print_parameter_file_guarded(cmd,gd,"parameters_null-cballs")==FAILURE) return FAILURE;
    return cballs_set_number_threads_guarded(cmd);
}
#endif
