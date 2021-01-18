// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTError_H
#define DTError_H

/*
 A central point for error messages.

 Set a breakpoint inside the source code to catch errors, such as out of bounds access to arrays.

*/

#include <string>
#include <vector>
#include <unistd.h>

extern void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t m);
extern void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t j,ssize_t m,ssize_t n);
extern void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t j,ssize_t k,ssize_t m,ssize_t n,ssize_t o);

extern void DTWarningMessage(std::string fcn,std::string msg);
extern void DTErrorMessage(std::string fcn,std::string msg);
extern void DTErrorMessage(std::string msg);

extern ssize_t DTHowManyErrors();

extern std::vector<std::string> &DTErrorList(void);

#endif
