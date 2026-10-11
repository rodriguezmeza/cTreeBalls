#include "globaldefs.h"
#include "parser.h"
#include "getparam.h"
#include <sys/wait.h>
#include <unistd.h>
#define CHECK(c) do { if (!(c)) {fprintf(stderr,"configuration check failed line %d\n",__LINE__);return 1;} } while(0)
int main(void)
{
    struct file_content fc={0}; ErrorMsg message={0}; int found=9, iv=77, pos=88;
    double dv=3.25; bool bv=FALSE; FileArg text="unchanged";
    CHECK(parser_init(&fc,2,"reader-test",message)==SUCCESS);
    strcpy(fc.name[0],"x"); strcpy(fc.name[1],"x");
    strcpy(fc.value[0],"1"); strcpy(fc.value[1],"2");
    CHECK(parser_read_int(&fc,"x",&iv,&found,message)==FAILURE && iv==77 && !found && !fc.read[0]);
    CHECK(parser_read_double_and_position(&fc,"x",&dv,&pos,&found,message)==FAILURE && dv==3.25 && pos==88 && !found);
    CHECK(parser_read_string(&fc,"x",&text,&found,message)==FAILURE && !strcmp(text,"unchanged"));
    CHECK(parser_read_bool(&fc,"x",&bv,&found,message)==FAILURE && !bv && !found);
    strcpy(fc.name[1],"other");
    CHECK(parser_read_int(&fc,"absent",&iv,&found,message)==SUCCESS && !found && iv==77);
    strcpy(fc.value[0],"YES");
    CHECK(parser_read_bool(&fc,"x",&bv,&found,message)==SUCCESS && found && bv);
    fc.read[0]=FALSE; strcpy(fc.value[0],"truejunk"); bv=FALSE;
    CHECK(parser_read_bool(&fc,"x",&bv,&found,message)==FAILURE && !found && !bv && !fc.read[0]);
    strcpy(fc.value[0],"17");
    CHECK(parser_read_int(&fc,"x",&iv,&found,message)==SUCCESS && found && iv==17);
    strcpy(fc.value[0],"3.5");
    CHECK(parser_read_double_and_position(&fc,"x",&dv,&pos,&found,message)==SUCCESS && dv==3.5 && pos==0);
    CHECK(parser_free(&fc)==SUCCESS);
    /* Exercise the standalone GetParam adapter, including its process error path. */
    char *defaults[]={"flag=false","Version=test",NULL};
    char *valid[]={"reader-test","flag=YeS",NULL};
    InitParam(valid,defaults); CHECK(GetbParam("flag")==TRUE);
    pid_t child=fork(); CHECK(child>=0);
    if(child==0){char *bad[]={"reader-test","flag=truejunk",NULL};InitParam(bad,defaults);GetbParam("flag");_exit(0);}
    int status; CHECK(waitpid(child,&status,0)==child && WIFEXITED(status) && WEXITSTATUS(status)!=0);
    puts("shared configuration readers: PASS"); return 0;
}
