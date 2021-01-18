// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTPoint3D.h"

#include "DTDoubleArray.h"
#include "DTError.h"

void DTPoint3D::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    std::cerr << "(" << x << ", " << y << ", " << z << ")\n" << std::flush;
#endif
}

void Read(const DTDataStorage &input,const std::string &name,DTPoint3D &toReturn)
{
    DTDoubleArray theArr = input.ReadDoubleArray(name);
    if (theArr.Length()==0) {
        // There isn't a way to indicate that this point is empty.
        toReturn = DTPoint3D(NAN,NAN,NAN);
    }
    else if (theArr.Length()!=3) {
        DTErrorMessage("ReadFromArray(DTPoint3D)","Invalid length of array.");
        toReturn = DTPoint3D();
    }
    else {
        toReturn = DTPoint3D(theArr(0),theArr(1),theArr(2));
    }
}

void Write(DTDataStorage &output,const std::string &name,const DTPoint3D &theVar)
{
    DTMutableDoubleArray theArr(3);
    theArr(0) = theVar.x;
    theArr(1) = theVar.y;
    theArr(2) = theVar.z;
    output.Save(theArr,name);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTPoint3D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"3D Point");
    output.Flush();
}

