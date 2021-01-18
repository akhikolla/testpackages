// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTRegion1D.h"

#include "DTDataStorage.h"

#include "DTIntArray.h"
#include "DTFloatArray.h"
#include "DTDoubleArray.h"
#include "DTError.h"

#include <math.h>
#include <limits>

void DTRegion1D::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    std::cerr << "[" << minV << ", " << maxV << "]\n" << flush;
#endif
}

bool operator==(const DTRegion1D &A,const DTRegion1D &B)
{
    return (A.isSet==B.isSet && A.minV==B.minV && A.maxV==B.maxV);
}

bool operator!=(const DTRegion1D &A,const DTRegion1D &B)
{
    return (A.isSet!=B.isSet || A.minV!=B.minV || A.maxV!=B.maxV);
}

DTRegion1D ValueRange(const DTIntArray &vals)
{
    if (vals.IsEmpty())
        return DTRegion1D();

    size_t len = vals.Length();
    int minV = 2147483647;
    int maxV = -2147483647;

    int v;
    size_t i;

    const int *valD = vals.Pointer();

    for (i=0;i<len;i++) {
        v = valD[i];
        minV = (v < minV ? v : minV);
        maxV = (maxV < v ? v : maxV);
    }

    if (minV>maxV)
        return DTRegion1D();
    else
        return DTRegion1D(minV,maxV);
}

DTRegion1D ValueRange(const DTFloatArray &vals)
{
    if (vals.IsEmpty())
        return DTRegion1D();

#if defined(WIN32) && !defined(INFINITY)
#define INFINITY std::numeric_limits<float>::infinity();
#endif
    size_t len = vals.Length();
    float minV = INFINITY;
    float maxV = -INFINITY;

    float v;
    size_t i;

    const float *valD = vals.Pointer();

    for (i=0;i<len;i++) {
        v = valD[i];
        minV = (v < minV ? v : minV);
        maxV = (maxV < v ? v : maxV);
    }

    if (minV>maxV)
        return DTRegion1D();
    else
        return DTRegion1D(minV,maxV);
}

DTRegion1D ValueRange(const DTDoubleArray &vals)
{
    if (vals.IsEmpty())
        return DTRegion1D();

#if defined(WIN32) && !defined(INFINITY)
#define INFINITY std::numeric_limits<float>::infinity();
#endif
    size_t len = vals.Length();
    double minV = INFINITY;
    double maxV = -INFINITY;

    double v;
    size_t i;

    const double *valD = vals.Pointer();

    for (i=0;i<len;i++) {
        v = valD[i];
        minV = (v < minV ? v : minV);
        maxV = (maxV < v ? v : maxV);
    }

    if (minV>maxV)
        return DTRegion1D();
    else
        return DTRegion1D(minV,maxV);
}

DTRegion1D Union(const DTRegion1D &A,const DTRegion1D &B)
{
    if (A.isSet==false) return B;
    if (B.isSet==false) return A;

    return DTRegion1D(A.minV<B.minV ? A.minV : B.minV,A.maxV>B.maxV ? A.maxV : B.maxV); // min(,),max(,) - but portable.
}

void Read(const DTDataStorage &input,const std::string &name,DTRegion1D &toReturn)
{
    DTDoubleArray theArr = input.ReadDoubleArray(name);
    if (theArr.Length()==0)
        toReturn = DTRegion1D();
    else if (theArr.Length()==2)
        toReturn = DTRegion1D(theArr(0),theArr(1));
    else {
        DTErrorMessage("ReadFromArray(DTRegion1D)","Invalid length of array.");
        toReturn = DTRegion1D();
    }
}

void Write(DTDataStorage &output,const std::string &name,const DTRegion1D &theVar)
{
    if (theVar.isSet==false)
        output.Save(DTDoubleArray(),name); // Empty array means not set/invalid.
    else {
        DTMutableDoubleArray theArr(2);
        theArr(0) = theVar.minV;
        theArr(1) = theVar.maxV;
        output.Save(theArr,name);
    }
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTRegion1D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"1D Region");
    output.Flush();
}

