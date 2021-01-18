// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTDoubleArrayOperators_Header
#define DTDoubleArrayOperators_Header

#include "DTDoubleArray.h"

// Array operator Array
extern DTMutableDoubleArray operator+(const DTDoubleArray &A,const DTDoubleArray &B);
extern DTMutableDoubleArray operator-(const DTDoubleArray &A,const DTDoubleArray &B);
extern DTMutableDoubleArray operator*(const DTDoubleArray &A,const DTDoubleArray &B);
extern DTMutableDoubleArray operator/(const DTDoubleArray &A,const DTDoubleArray &B);

// Array operator Number
extern DTMutableDoubleArray operator+(const DTDoubleArray &A,double b);
extern DTMutableDoubleArray operator-(const DTDoubleArray &A,double b);
extern DTMutableDoubleArray operator*(const DTDoubleArray &A,double b);
extern DTMutableDoubleArray operator/(const DTDoubleArray &A,double b);

// Number operator Array
extern DTMutableDoubleArray operator+(double a,const DTDoubleArray &B);
extern DTMutableDoubleArray operator-(double a,const DTDoubleArray &B);
extern DTMutableDoubleArray operator*(double a,const DTDoubleArray &B);
extern DTMutableDoubleArray operator/(double a,const DTDoubleArray &B);

// Negating
extern DTMutableDoubleArray operator-(const DTDoubleArray &A);

#endif
