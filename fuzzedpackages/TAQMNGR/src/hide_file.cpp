#include <taq.h>

void hide_file(const char *name)
{
	#ifdef _WIN32
		SetFileAttributes(name, FILE_ATTRIBUTE_HIDDEN);
		return;
	#else
		return;
	#endif
}
