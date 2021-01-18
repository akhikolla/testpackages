// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTUtilities.h"

#include <cstdio>

string DTInt2String(int n)
{
    char temp[20];
    sprintf(temp,"%d",n);
    return string(temp);
}

string DTSize2String(ssize_t n)
{
    char temp[30];
    sprintf(temp,"%ld",(long int)n);
    return string(temp);
}

string DTFloat2StringShort(double v)
{
    char temp[20];
    sprintf(temp,"%g",v);
    return string(temp);
}

DTRange Intersection(const DTRange &a,const DTRange &b)
{
    if (a.start<b.start) {
        if (a.end()<b.start) return DTRange();
        if (a.end()<b.end())
            return DTRange(b.start,a.end()-b.start);
        else
            return b;
    }
    else {
        if (b.end()<a.start) return DTRange();
        if (a.end()<b.end())
            return a;
        else
            return DTRange(a.start,b.end()-a.start);
    }
}
