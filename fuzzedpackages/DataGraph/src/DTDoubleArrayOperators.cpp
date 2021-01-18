// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTDoubleArrayOperators.h"

#include "DTArrayTemplates.h"

DTMutableDoubleArray operator+(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return DTAddArrays<DTDoubleArray,DTMutableDoubleArray,double>("DoubleArray+DoubleArray",A,B);
}

DTMutableDoubleArray operator-(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return DTSubtractArrays<DTDoubleArray,DTMutableDoubleArray,double>("DoubleArray-DoubleArray",A,B);
}

DTMutableDoubleArray operator*(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return DTMultiplyArrays<DTDoubleArray,DTMutableDoubleArray,double>("DoubleArray*DoubleArray",A,B);
}

DTMutableDoubleArray operator/(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return DTDivideArrays<DTDoubleArray,DTMutableDoubleArray,double>("DoubleArray/DoubleArray",A,B);
}

DTMutableDoubleArray operator+(const DTDoubleArray &A,double b)
{
    return DTArrayPlusNumber<DTDoubleArray,DTMutableDoubleArray,double>(A,b);
}

DTMutableDoubleArray operator-(const DTDoubleArray &A,double b)
{
    return DTArrayPlusNumber<DTDoubleArray,DTMutableDoubleArray,double>(A,-b);
}

DTMutableDoubleArray operator*(const DTDoubleArray &A,double b)
{
    return DTArrayTimesNumber<DTDoubleArray,DTMutableDoubleArray,double>(A,b);
}

DTMutableDoubleArray operator/(const DTDoubleArray &A,double b)
{
    return DTArrayTimesNumber<DTDoubleArray,DTMutableDoubleArray,double>(A,1.0/b);
}

DTMutableDoubleArray operator+(double a,const DTDoubleArray &B)
{
    return DTArrayPlusNumber<DTDoubleArray,DTMutableDoubleArray,double>(B,a);
}

DTMutableDoubleArray operator-(double a,const DTDoubleArray &B)
{
    return DTNumberMinusArray<DTDoubleArray,DTMutableDoubleArray,double>(a,B);
}

DTMutableDoubleArray operator*(double a,const DTDoubleArray &B)
{
    return DTArrayTimesNumber<DTDoubleArray,DTMutableDoubleArray,double>(B,a);
}

DTMutableDoubleArray operator/(double a,const DTDoubleArray &B)
{
    return DTNumberDividedByArray<DTDoubleArray,DTMutableDoubleArray,double>(a,B);
}

DTMutableDoubleArray operator-(const DTDoubleArray &A)
{
    return DTNegateArray<DTDoubleArray,DTMutableDoubleArray,double>(A);
}


