// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTFloatArrayOperators.h"

#include "DTArrayTemplates.h"

DTMutableFloatArray operator+(const DTFloatArray &A,const DTFloatArray &B)
{
    return DTAddArrays<DTFloatArray,DTMutableFloatArray,float>("FloatArray+FloatArray",A,B);
}

DTMutableFloatArray operator-(const DTFloatArray &A,const DTFloatArray &B)
{
    return DTSubtractArrays<DTFloatArray,DTMutableFloatArray,float>("FloatArray-FloatArray",A,B);
}

DTMutableFloatArray operator*(const DTFloatArray &A,const DTFloatArray &B)
{
    return DTMultiplyArrays<DTFloatArray,DTMutableFloatArray,float>("FloatArray*FloatArray",A,B);
}

DTMutableFloatArray operator/(const DTFloatArray &A,const DTFloatArray &B)
{
    return DTDivideArrays<DTFloatArray,DTMutableFloatArray,float>("FloatArray/FloatArray",A,B);
}

DTMutableFloatArray operator+(const DTFloatArray &A,float b)
{
    return DTArrayPlusNumber<DTFloatArray,DTMutableFloatArray,float>(A,b);
}

DTMutableFloatArray operator-(const DTFloatArray &A,float b)
{
    return DTArrayPlusNumber<DTFloatArray,DTMutableFloatArray,float>(A,-b);
}

DTMutableFloatArray operator*(const DTFloatArray &A,float b)
{
    return DTArrayTimesNumber<DTFloatArray,DTMutableFloatArray,float>(A,b);
}

DTMutableFloatArray operator/(const DTFloatArray &A,float b)
{
    return DTArrayTimesNumber<DTFloatArray,DTMutableFloatArray,float>(A,1.0f/b);
}

DTMutableFloatArray operator+(float a,const DTFloatArray &B)
{
    return DTArrayPlusNumber<DTFloatArray,DTMutableFloatArray,float>(B,a);
}

DTMutableFloatArray operator-(float a,const DTFloatArray &B)
{
    return DTNumberMinusArray<DTFloatArray,DTMutableFloatArray,float>(a,B);
}

DTMutableFloatArray operator*(float a,const DTFloatArray &B)
{
    return DTArrayTimesNumber<DTFloatArray,DTMutableFloatArray,float>(B,a);
}

DTMutableFloatArray operator/(float a,const DTFloatArray &B)
{
    return DTNumberDividedByArray<DTFloatArray,DTMutableFloatArray,float>(a,B);
}

DTMutableFloatArray operator-(const DTFloatArray &A)
{
    return DTNegateArray<DTFloatArray,DTMutableFloatArray,float>(A);
}

