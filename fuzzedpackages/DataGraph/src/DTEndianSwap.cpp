// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTEndianSwap.h"

#include "DTUCharArray.h"
#include "DTShortIntArray.h"
#include "DTUShortIntArray.h"
#include "DTIntArray.h"
#include "DTDoubleArray.h"
#include "DTFloatArray.h"

#include "DTError.h"

void DTSwap2Bytes(unsigned char *data,size_t length)
{
    if (length==0 || length%2!=0)
        return;
    
    size_t i;
    unsigned char temp;
    for (i=0;i<length;i+=2) {
        temp = data[i];
        data[i] = data[i+1];
        data[i+1] = temp;
    }
}

void DTSwap4Bytes(unsigned char *data,size_t length)
{
    if (length==0 || length%4!=0)
        return;

    size_t i;
    unsigned char temp1,temp2;
    for (i=0;i<length;i+=4) {
        temp1 = data[i];
        temp2 = data[i+1];
        data[i] = data[i+3];
        data[i+1] = data[i+2];
        data[i+3] = temp1;
        data[i+2] = temp2;
    }
}

void DTSwap8Bytes(unsigned char *data,size_t length)
{
    if (length==0 || length%8!=0)
        return;

    size_t i;
    unsigned char temp1,temp2,temp3,temp4;
    for (i=0;i<length;i+=8) {
        temp1 = data[i];
        temp2 = data[i+1];
        temp3 = data[i+2];
        temp4 = data[i+3];
        data[i] = data[i+7];
        data[i+1] = data[i+6];
        data[i+2] = data[i+5];
        data[i+3] = data[i+4];
        data[i+7] = temp1;
        data[i+6] = temp2;
        data[i+5] = temp3;
        data[i+4] = temp4;
    }
}

void DTSwap10Bytes(unsigned char *data,size_t length)
{
    if (length==0 || length%10!=0)
        return;
    
    size_t i;
    unsigned char temp1,temp2,temp3,temp4,temp5;
    for (i=0;i<length;i+=10) {
        temp1 = data[i];
        temp2 = data[i+1];
        temp3 = data[i+2];
        temp4 = data[i+3];
        temp5 = data[i+4];
        data[i] = data[i+9];
        data[i+1] = data[i+8];
        data[i+2] = data[i+7];
        data[i+3] = data[i+6];
        data[i+4] = data[i+5];
        data[i+9] = temp1;
        data[i+8] = temp2;
        data[i+7] = temp3;
        data[i+6] = temp4;
        data[i+5] = temp5;
    }
}

void Swap2Bytes(DTMutableUCharArray &arr)
{
    if (arr.m()%2!=0) {
        DTErrorMessage("Swap2Bytes(UCharArray)","First array dimension needs to be even.");
        return;
    }

    DTSwap2Bytes((unsigned char *)arr.Pointer(),arr.Length());
}

void Swap4Bytes(DTMutableUCharArray &arr)
{
    if (arr.m()%4!=0) {
        DTErrorMessage("Swap4Bytes(UCharArray)","First array dimension needs to be divisible by 4.");
        return;
    }

    DTSwap4Bytes((unsigned char *)arr.Pointer(),arr.Length());
}

void Swap8Bytes(DTMutableUCharArray &arr)
{
    if (arr.m()%8!=0) {
        DTErrorMessage("Swap8Bytes(UCharArray)","First array dimension needs to be divisible by 8.");
        return;
    }

    DTSwap8Bytes((unsigned char *)arr.Pointer(),arr.Length());
}

void SwapEndian(DTMutableShortIntArray &arr)
{
    DTSwap2Bytes((unsigned char *)arr.Pointer(),arr.Length()*2);
}

void SwapEndian(DTMutableUShortIntArray &arr)
{
    DTSwap2Bytes((unsigned char *)arr.Pointer(),arr.Length()*2);
}

void SwapEndian(DTMutableIntArray &arr)
{
    DTSwap4Bytes((unsigned char *)arr.Pointer(),arr.Length()*4);
}

void SwapEndian(DTMutableDoubleArray &arr)
{
    DTSwap8Bytes((unsigned char *)arr.Pointer(),arr.Length()*8);
}

void SwapEndian(DTMutableFloatArray &arr)
{
    DTSwap4Bytes((unsigned char *)arr.Pointer(),arr.Length()*4);
}

