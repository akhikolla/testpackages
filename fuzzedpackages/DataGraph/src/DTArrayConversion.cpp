// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTArrayConversion.h"

#include "DTIntArray.h"
#include "DTFloatArray.h"
#include "DTDoubleArray.h"
#include "DTShortIntArray.h"
#include "DTUShortIntArray.h"
#include "DTUCharArray.h"
#include "DTCharArray.h"

#include "DTError.h"

template <class TA,class SA,class T,class S> void DTConversionBetweenPointers(const TA &A,const SA &B,
                                                                              const T *from,S *to,ssize_t length)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("ConvertArray(A,B)","Incompatible array sizes.");
        return;
    }
    for (ssize_t i=0;i<length;i++) to[i] = (S)from[i];
}

DTMutableDoubleArray ConvertToDouble(const DTFloatArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableDoubleArray ConvertToDouble(const DTIntArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableDoubleArray ConvertToDouble(const DTShortIntArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableDoubleArray ConvertToDouble(const DTUShortIntArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableDoubleArray ConvertToDouble(const DTCharArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableDoubleArray ConvertToDouble(const DTUCharArray &A)
{
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableDoubleArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTDoubleArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTIntArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTShortIntArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTUShortIntArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTCharArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableFloatArray ConvertToFloat(const DTUCharArray &A)
{
    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableFloatArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTDoubleArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTFloatArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTShortIntArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTUShortIntArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTCharArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableIntArray ConvertToInt(const DTUCharArray &A)
{
    DTMutableIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTDoubleArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTFloatArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTIntArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTShortIntArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTUShortIntArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUCharArray ConvertToUChar(const DTCharArray &A)
{
    DTMutableUCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableUCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTDoubleArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTFloatArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTIntArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTShortIntArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTUShortIntArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableCharArray ConvertToChar(const DTUCharArray &A)
{
    DTMutableCharArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableCharArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTDoubleArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTFloatArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTIntArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTUShortIntArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUShortIntArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTCharArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableShortIntArray ConvertToShortInt(const DTUCharArray &A)
{
    DTMutableShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTDoubleArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTDoubleArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTFloatArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTFloatArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTIntArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTIntArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTShortIntArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTShortIntArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTCharArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTCharArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

DTMutableUShortIntArray ConvertToUShortInt(const DTUCharArray &A)
{
    DTMutableUShortIntArray toReturn(A.m(),A.n(),A.o());
    ConvertArray(A,toReturn);
    return toReturn;
}

void ConvertArray(const DTUCharArray &A,DTMutableUShortIntArray &B)
{
    DTConversionBetweenPointers(A,B,A.Pointer(),B.Pointer(),B.Length());
}

