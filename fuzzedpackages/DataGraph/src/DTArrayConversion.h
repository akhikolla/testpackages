// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTArrayConversion_Header
#define DTArrayConversion_Header

// Convert between array types.  Two forms
//
// newDoubleArray = ConvertToDouble(fromarray);
// newFloatArray = ConvertToFloat(fromarray);
// ConvertToInt(...)
// ConvertToChar(...)
// ConvertToUChar(...)
// ConvertToShortInt()
// ConvertToUShortInt()
// 
// and
//
// ConvertArray(fromArray,toArray);
//
//    In this case, both arrays have to have the same size,
//    and the toArray has to be mutable.

// Functions in DTSource, with the exception of the IO routines, will not automatically
// convert an array from one type to another.  This is done to avoid the common "problem"
// where an automatic conversion will not be evident from looking the source code
// but can cost a significant amount of execution time.
// This way the compiler will complain when you pass a float array to a function
// that expects a double, and even though this conversion is loss-less it is very likely
// that the conversion is an unwanted side effect.  If you really want to do that, you can
// always use 
//              foo(ConvertToDouble(theArray))

#include <stdio.h>
#include <unistd.h>

class DTDoubleArray;
class DTMutableDoubleArray;
class DTFloatArray;
class DTMutableFloatArray;
class DTIntArray;
class DTMutableIntArray;
class DTCharArray;
class DTMutableCharArray;
class DTUCharArray;
class DTMutableUCharArray;
class DTShortIntArray;
class DTMutableShortIntArray;
class DTUShortIntArray;
class DTMutableUShortIntArray;

// Conversion to DTMutableDoubleArray
extern DTMutableDoubleArray ConvertToDouble(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableDoubleArray &);
extern DTMutableDoubleArray ConvertToDouble(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableDoubleArray &);
extern DTMutableDoubleArray ConvertToDouble(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableDoubleArray &);
extern DTMutableDoubleArray ConvertToDouble(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableDoubleArray &);
extern DTMutableDoubleArray ConvertToDouble(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableDoubleArray &);
extern DTMutableDoubleArray ConvertToDouble(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableDoubleArray &);

// Conversion to DTMutableFloatArray
extern DTMutableFloatArray ConvertToFloat(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableFloatArray &);
extern DTMutableFloatArray ConvertToFloat(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableFloatArray &);
extern DTMutableFloatArray ConvertToFloat(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableFloatArray &);
extern DTMutableFloatArray ConvertToFloat(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableFloatArray &);
extern DTMutableFloatArray ConvertToFloat(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableFloatArray &);
extern DTMutableFloatArray ConvertToFloat(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableFloatArray &);

// Conversion to DTMutableIntArray
extern DTMutableIntArray ConvertToInt(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableIntArray &);
extern DTMutableIntArray ConvertToInt(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableIntArray &);
extern DTMutableIntArray ConvertToInt(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableIntArray &);
extern DTMutableIntArray ConvertToInt(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableIntArray &);
extern DTMutableIntArray ConvertToInt(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableIntArray &);
extern DTMutableIntArray ConvertToInt(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableIntArray &);

// Conversion to DTMutableUCharArray - unsigned char
extern DTMutableUCharArray ConvertToUChar(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableUCharArray &);
extern DTMutableUCharArray ConvertToUChar(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableUCharArray &);
extern DTMutableUCharArray ConvertToUChar(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableUCharArray &);
extern DTMutableUCharArray ConvertToUChar(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableUCharArray &);
extern DTMutableUCharArray ConvertToUChar(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableUCharArray &);
extern DTMutableUCharArray ConvertToUChar(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableUCharArray &);

// Conversion to DTMutableCharArray - char
extern DTMutableCharArray ConvertToChar(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableCharArray &);
extern DTMutableCharArray ConvertToChar(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableCharArray &);
extern DTMutableCharArray ConvertToChar(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableCharArray &);
extern DTMutableCharArray ConvertToChar(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableCharArray &);
extern DTMutableCharArray ConvertToChar(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableCharArray &);
extern DTMutableCharArray ConvertToChar(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableCharArray &);

// Conversion to DTMutableShortIntArray - short int
extern DTMutableShortIntArray ConvertToShortInt(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableShortIntArray &);
extern DTMutableShortIntArray ConvertToShortInt(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableShortIntArray &);

// Conversion to DTMutableUShortIntArray - unsigned short int
extern DTMutableUShortIntArray ConvertToUShortInt(const DTDoubleArray &A);
extern void ConvertArray(const DTDoubleArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTFloatArray &A);
extern void ConvertArray(const DTFloatArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTIntArray &A);
extern void ConvertArray(const DTIntArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTShortIntArray &A);
extern void ConvertArray(const DTShortIntArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTUShortIntArray &A);
extern void ConvertArray(const DTUShortIntArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTCharArray &A);
extern void ConvertArray(const DTCharArray &A,DTMutableUShortIntArray &);
extern DTMutableUShortIntArray ConvertToUShortInt(const DTUCharArray &A);
extern void ConvertArray(const DTUCharArray &A,DTMutableUShortIntArray &);

template <class T,class S> void DTConversionBetweenRowPointers(const T *from,S *to,size_t length)
{
    // It is the callers responsibility to make sure that the size is valid.
    for (size_t i=0;i<length;i++) to[i] = (S)from[i];
}

#endif
