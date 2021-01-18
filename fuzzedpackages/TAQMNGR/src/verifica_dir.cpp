#include <taq.h>


int VerificaDir(const char *ndir)
{
	struct stat myStat;
	if(stat(ndir,&myStat) == 0) return 0;
	else return 1;
}
