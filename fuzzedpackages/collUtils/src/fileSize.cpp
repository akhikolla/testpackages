#include "fileSize.h"

long fileSize(std::string fn)
{
    struct stat fsStatBuf;
    stat(fn.c_str(), &fsStatBuf);
    long fileLen = fsStatBuf.st_size;
    return fileLen;
}
