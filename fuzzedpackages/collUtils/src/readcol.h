#ifndef READCOL_H
#define READCOL_H

#include <vector>
#include <string>

std::vector< std::string > readcol(std::string fileName, long colNum, long nSkip, long maxRowNum);

#endif // READCOL_H
