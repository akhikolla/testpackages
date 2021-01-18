#include <taq.h>

int create_directory(const char *name)
{
	#ifdef _WIN32
		if(CreateDirectory(name, NULL)) return 0;
		else return 1;
	#else
		return mkdir(name, S_IRWXU | S_IRWXG | S_IRWXO);
	#endif
}
	
