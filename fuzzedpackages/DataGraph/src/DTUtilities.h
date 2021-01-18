// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTUtilities_H
#define DTUtilities_H

#include <string>
#include <unistd.h>

using namespace std;

extern std::string DTInt2String(int n);
extern std::string DTSize2String(ssize_t n);
extern std::string DTFloat2StringShort(double); // The %lg action

struct DTRange {
    DTRange() : start(0), length(0) {}
    DTRange(ssize_t st,ssize_t ln) : start(st), length(ln) {}
    ssize_t start;
    ssize_t length;

    ssize_t end(void) const {return start+length;}
};

extern DTRange Intersection(const DTRange &,const DTRange &);

#endif
