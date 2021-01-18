
// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTError.h"
#include "DTUtilities.h"

#include <iostream>

static vector<std::string> errorList;

std::vector<std::string> &DTErrorList(void)
{
    return errorList;
}

ssize_t DTHowManyErrors()
{
    return (ssize_t)errorList.size();
}

void DTErrorMessage(std::string fcn,std::string msg)
{
    std::string theErr = fcn + ": " + msg;
    DTErrorMessage(theErr);
}

void DTErrorMessage(std::string msg)
{
    // Set a breakpoint here, and then trace it back.
    if (msg.length()==0) return;
    errorList.push_back(msg);
#ifndef DG_NOSTDErrOut
    std::cerr << msg << std::endl;
    std::cerr.flush();
#endif
}

void DTWarningMessage(std::string fcn,std::string msg)
{
    // DTSource does not use this, so you don't need to set a breakpoint unless you use it.
    std::string theErr = fcn + ": " + msg;
    errorList.push_back(theErr);
#ifndef DG_NOSTDErrOut
    std::cerr << theErr << std::endl;
    std::cerr.flush();
#endif
}

void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t m)
{
    std::string toReturn = type + "(" + DTSize2String(i) + ") is not valid, needs to be lie in [0," + DTSize2String(m-1) + "].";
    DTErrorMessage(toReturn);
}

void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t j,ssize_t m,ssize_t n)
{
    std::string toReturn = type + "(" + DTSize2String(i) + "," + DTSize2String(j) + ") is not valid, needs to lie in [0,"  + DTSize2String(m-1) + "]x[0," + DTSize2String(n-1) + "].";
    DTErrorMessage(toReturn);
}

void DTErrorOutOfRange(std::string type,ssize_t i,ssize_t j,ssize_t k,ssize_t m,ssize_t n,ssize_t o)
{
    std::string toReturn = type + "(" + DTSize2String(i) + "," + DTSize2String(j) + "," + DTSize2String(k) + ") is not valid, needs to lie in [0," + DTSize2String(m-1) + "]x[0," + DTSize2String(n-1) + "]x[0," + DTSize2String(o-1) + "].";
    DTErrorMessage(toReturn);
}

