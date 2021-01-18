// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTFloatArrayOperators_Header
#define DTFloatArrayOperators_Header

#include "DTFloatArray.h"

// Array operator Array
extern DTMutableFloatArray operator+(const DTFloatArray &A,const DTFloatArray &B);
extern DTMutableFloatArray operator-(const DTFloatArray &A,const DTFloatArray &B);
extern DTMutableFloatArray operator*(const DTFloatArray &A,const DTFloatArray &B);
extern DTMutableFloatArray operator/(const DTFloatArray &A,const DTFloatArray &B);

// Array operator Number
extern DTMutableFloatArray operator+(const DTFloatArray &A,float b);
extern DTMutableFloatArray operator-(const DTFloatArray &A,float b);
extern DTMutableFloatArray operator*(const DTFloatArray &A,float b);
extern DTMutableFloatArray operator/(const DTFloatArray &A,float b);

// Number operator Array
extern DTMutableFloatArray operator+(float a,const DTFloatArray &B);
extern DTMutableFloatArray operator-(float a,const DTFloatArray &B);
extern DTMutableFloatArray operator*(float a,const DTFloatArray &B);
extern DTMutableFloatArray operator/(float a,const DTFloatArray &B);

// Negating
extern DTMutableFloatArray operator-(const DTFloatArray &A);

#endif
