// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTEndianSwap_Header
#define DTEndianSwap_Header

// Very basic little-big endian conversion routines. Used by the IO routines.
#include <cstring>

class DTMutableUCharArray;
class DTMutableUShortIntArray;
class DTMutableShortIntArray;
class DTMutableIntArray;
class DTMutableFloatArray;
class DTMutableDoubleArray;

extern void Swap2Bytes(DTMutableUCharArray &arr);
extern void Swap4Bytes(DTMutableUCharArray &arr);
extern void Swap8Bytes(DTMutableUCharArray &arr);

extern void SwapEndian(DTMutableShortIntArray &arr);
extern void SwapEndian(DTMutableUShortIntArray &arr);
extern void SwapEndian(DTMutableIntArray &arr);
extern void SwapEndian(DTMutableFloatArray &arr);
extern void SwapEndian(DTMutableDoubleArray &arr);

// Raw calls
extern void DTSwap2Bytes(unsigned char *data,size_t numberOfBytes);
extern void DTSwap4Bytes(unsigned char *data,size_t numberOfBytes);
extern void DTSwap8Bytes(unsigned char *data,size_t numberOfBytes);
extern void DTSwap10Bytes(unsigned char *data,size_t numberOfBytes);

#endif
