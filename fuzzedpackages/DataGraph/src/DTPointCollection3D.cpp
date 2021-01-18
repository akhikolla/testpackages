// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTPointCollection3D.h"

#include "DTError.h"
#include "DTDataStorage.h"
#include "DTRegion3D.h"
#include "DTArrayConversion.h"
#include "DTDoubleArrayOperators.h"
#include "DTFloatArrayOperators.h"
#include "DTUtilities.h"
#include "DTTransform3D.h"

#include <cstring>

DTPointCollection3D::DTPointCollection3D(const DTFloatArray &input)
{
    if (input.IsEmpty()) return;

    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray)","Invalid size of array.");
        return;
    }

    floatData = input;
}

DTPointCollection3D::DTPointCollection3D(const DTFloatArray &input,const DTIntArray &ptN)
{
    if (input.IsEmpty()) return;

    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray,DTIntArray)","Invalid size of array.");
        return;
    }

    floatData = input;

    if (ptN.Length()!=input.n() || ptN.m()!=ptN.Length()) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray,DTIntArray)","Point number array has the wrong size.");
    }
    else {
        pointNumbers = ptN;
    }
}

DTPointCollection3D::DTPointCollection3D(const DTDoubleArray &input)
{
    if (input.IsEmpty()) return;
    
    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray)","Invalid size of array.");
        return;
    }
    
    doubleData = input;
}

DTPointCollection3D::DTPointCollection3D(const DTDoubleArray &input,const DTIntArray &ptN)
{
    if (input.IsEmpty()) return;
    
    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray,DTIntArray)","Invalid size of array.");
        return;
    }
    
    doubleData = input;
    
    if (ptN.Length()!=input.n() || ptN.m()!=ptN.Length()) {
        DTErrorMessage("DTPointCollection3D(DTFloatArray,DTIntArray)","Point number array has the wrong size.");
    }
    else {
        pointNumbers = ptN;
    }
}

ssize_t DTPointCollection3D::NumberOfPoints(void) const
{
    if (doubleData.NotEmpty()) {
        return doubleData.n();
    }
    else {
        return floatData.n();
    }
}

DTDoubleArray DTPointCollection3D::DoubleData(void) const
{
    if (floatData.NotEmpty()) {
        DTErrorMessage("DTPointCollection3D::DoubleData","Array saved as float."); // Call ConvertToDouble(points) to convert it
        return DTDoubleArray();
    }
    else {
        return doubleData;
    }
}

DTFloatArray DTPointCollection3D::FloatData(void) const
{
    if (doubleData.NotEmpty()) {
        DTErrorMessage("DTPointCollection3D::FloatData","Array saved as double."); // Call ConvertToFloat(points) to convert it
        return DTFloatArray();
    }
    else {
        return floatData;
    }
}

DTMutablePointCollection3D DTPointCollection3D::Copy(void) const
{
    if (doubleData.NotEmpty()) {
        if (pointNumbers.NotEmpty())
            return DTMutablePointCollection3D(DoubleData().Copy(),pointNumbers);
        else
            return DTMutablePointCollection3D(DoubleData().Copy());
    }
    else {
        if (pointNumbers.NotEmpty())
            return DTMutablePointCollection3D(FloatData().Copy(),pointNumbers);
        else
            return DTMutablePointCollection3D(FloatData().Copy());
    }
}

DTFloatArray DTPointCollection3D::Data(void) const
{
    DTErrorMessage("PointCollection::Data","Obsolete.  Use FloatData or DoubleData.");
    return FloatData();
}

DTPoint3D DTPointCollection3D::operator()(ssize_t i) const
{
    if (doubleData.NotEmpty()) {
        return DTPoint3D(doubleData(0,i),doubleData(1,i),doubleData(2,i));
    }
    else {
        return DTPoint3D(floatData(0,i),floatData(1,i),floatData(2,i));
    }
}

void DTPointCollection3D::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    if (IsEmpty())
        std::cerr << "No points";
    else if (NumberOfPoints()==1)
        std::cerr << "1 point";
    else
        std::cerr << NumberOfPoints() << " points";
    if (pointNumbers.Length()!=0) {
        std::cerr << ", point numbers specified";
    }
    if (DoublePrecision()) {
        std::cerr << " (double)";
    }
    else if (FloatPrecision()) {
        std::cerr << " (float)";
    }
    
    std::cerr << std::endl << std::flush;
#endif
}

void DTPointCollection3D::pall(void) const
{
#ifndef DG_NOSTDErrOut
	ssize_t i,howMany = NumberOfPoints();
    DTPoint3D P;
    if (pointNumbers.Length()==0) {
        for (i=0;i<howMany;i++) {
            P = operator()(i);
            std::cerr << i << " : " << P.x << ", " << P.y << ", " << P.z << std::endl;
        }
    }
    else {
        for (i=0;i<howMany;i++) {
            P = operator()(i);
            std::cerr << i << " : " << P.x << ", " << P.y << ", " << P.z  << " # = " << pointNumbers(i) << std::endl;
        }
    }
#endif
}

#pragma mark Mutable version

DTMutablePointCollection3D::DTMutablePointCollection3D(const DTMutableDoubleArray &input)
{
    if (input.IsEmpty()) return;

    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(Array)","Invalid size of array.");
        return;
    }

    doubleData = input;
    mutableDoubleData = input;
}

DTMutablePointCollection3D::DTMutablePointCollection3D(const DTMutableDoubleArray &input,const DTIntArray &ptN)
{
    if (input.IsEmpty()) return;

    if (input.m()!=2 || input.o()!=1 || ptN.m()!=ptN.Length() || ptN.Length()!=input.n()) {
        DTErrorMessage("DTPointCollection3D(Array,Point Numbers)","Invalid size of array.");
        return;
    }

    doubleData = input;
    mutableDoubleData = input;
    pointNumbers = ptN;
}

DTMutablePointCollection3D::DTMutablePointCollection3D(const DTMutableFloatArray &input)
{
    if (input.IsEmpty()) return;

    if (input.m()!=3 || input.o()!=1) {
        DTErrorMessage("DTPointCollection3D(Array)","Invalid size of array.");
        return;
    }

    floatData = input;
    mutableFloatData = input;
}

DTMutablePointCollection3D::DTMutablePointCollection3D(const DTMutableFloatArray &input,const DTIntArray &ptN)
{
    if (input.IsEmpty()) return;

    if (input.m()!=3 || input.o()!=1 || ptN.m()!=ptN.Length() || ptN.Length()!=input.n()) {
        DTErrorMessage("DTPointCollection3D(Array,Point Numbers)","Invalid size of array.");
        return;
    }

    floatData = input;
    mutableFloatData = input;
    pointNumbers = ptN;
}

void DTMutablePointCollection3D::operator-=(DTPoint3D P)
{
    if (DoublePrecision()) {
        ssize_t i,howMany = mutableDoubleData.n();
        for (i=0;i<howMany;i++) {
            mutableDoubleData(0,i) -= P.x;
            mutableDoubleData(1,i) -= P.y;
            mutableDoubleData(2,i) -= P.z;
        }
    }
    else {
        ssize_t i,howMany = mutableFloatData.n();
        float x = float(P.x);
        float y = float(P.y);
        float z = float(P.z);
        for (i=0;i<howMany;i++) {
            mutableFloatData(0,i) -= x;
            mutableFloatData(1,i) -= y;
            mutableFloatData(2,i) -= z;
        }
    }
}

void DTMutablePointCollection3D::operator+=(const DTPointCollection3D &P)
{
    if (NumberOfPoints()!=P.NumberOfPoints()) {
        DTErrorMessage("MutablePointCollection3D+=PointCollection3D","Sizes don't match");
    }
    else if (P.DoublePrecision()!=DoublePrecision()) {
        DTErrorMessage("MutablePointCollection3D+=PointCollection3D","Both have to be double or float");
    }
    else {
        if (P.DoublePrecision()) {
            mutableDoubleData += P.DoubleData();
        }
        else {
            mutableFloatData += P.FloatData();
        }
    }
}

void DTMutablePointCollection3D::operator-=(const DTPointCollection3D &P)
{
    if (NumberOfPoints()!=P.NumberOfPoints()) {
        DTErrorMessage("MutablePointCollection3D-=PointCollection3D","Sizes don't match");
    }
    else if (P.DoublePrecision()!=DoublePrecision()) {
        DTErrorMessage("MutablePointCollection3D-=PointCollection3D","Both have to be double or float");
    }
    else {
        if (P.DoublePrecision()) {
            mutableDoubleData -= P.DoubleData();
        }
        else {
            mutableFloatData -= P.FloatData();
        }
    }
}

void DTMutablePointCollection3D::Overwrite(const DTPointCollection3D &P,const DTRange &r)
{
    if (r.length!=P.NumberOfPoints()) {
        DTErrorMessage("PointCollection::Overwite(PointCollection3D,Range)", "Range has a different length than the point collection");
        return;
    }
    if (r.start+r.length>=NumberOfPoints()) {
        DTErrorMessage("PointCollection::Overwite(PointCollection3D,Range)", "Range is out of bounds");
        return;
    }
    else if (P.DoublePrecision()!=DoublePrecision()) {
        DTErrorMessage("PointCollection::Overwite(PointCollection3D,Range)","Both have to be double or float");
    }
    if (DoublePrecision()) {
        std::memcpy(mutableDoubleData.Pointer()+2*r.start, P.DoubleData().Pointer(), sizeof(double)*3*r.length);
    }
    else {
        std::memcpy(mutableFloatData.Pointer()+2*r.start, P.FloatData().Pointer(), sizeof(float)*3*r.length);
    }
}

bool operator==(const DTPointCollection3D &A,const DTPointCollection3D &B)
{
    if (A.DoublePrecision()!=B.DoublePrecision()) return false;
    if (A.DoublePrecision()) {
        if (A.DoubleData()!=B.DoubleData()) return false;
    }
    else {
        if (A.FloatData()!=B.FloatData()) return false;
    }
    return (A.PointNumbers()==B.PointNumbers());
}

bool operator!=(const DTPointCollection3D &A,const DTPointCollection3D &B)
{
    if (A.DoublePrecision()!=B.DoublePrecision()) return true;
    if (A.DoublePrecision()) {
        if (A.DoubleData()!=B.DoubleData()) return true;
    }
    else {
        if (A.FloatData()!=B.FloatData()) return true;
    }
    return (A.PointNumbers()!=B.PointNumbers());
}

#pragma mark Functions

void Read(const DTDataStorage &input,const std::string &name,DTPointCollection3D &toReturn)
{
    std::string theName = input.ResolveName(name);
    
    // The main array stores the loop structure
    if (input.SavedAsDouble(theName)) {
        DTDoubleArray points = input.ReadDoubleArray(theName);
        if (input.Contains(name+"_ptN")) {
            DTIntArray pointNumbers = input.ReadIntArray(theName+"_ptN");
            toReturn = DTPointCollection3D(points,pointNumbers);
        }
        else {
            toReturn = DTPointCollection3D(points);
        }
    }
    else {
        DTFloatArray points = input.ReadFloatArray(theName);
        if (input.Contains(name+"_ptN")) {
            DTIntArray pointNumbers = input.ReadIntArray(theName+"_ptN");
            toReturn = DTPointCollection3D(points,pointNumbers);
        }
        else {
            toReturn = DTPointCollection3D(points);
        }
    }
}

DTPointCollection3D ConvertToFloat(const DTPointCollection3D &p)
{
    if (p.FloatPrecision()) {
        return p;
    }
    else {
        if (p.PointNumbers().NotEmpty()) {
            return DTPointCollection3D(ConvertToFloat(p.DoubleData()),p.PointNumbers());
        }
        else {
            return DTPointCollection3D(ConvertToFloat(p.DoubleData()));
        }
    }
}

DTPointCollection3D ConvertToDouble(const DTPointCollection3D &p)
{
    if (p.DoublePrecision()) {
        return p;
    }
    else {
        if (p.PointNumbers().NotEmpty()) {
            return DTPointCollection3D(ConvertToDouble(p.FloatData()),p.PointNumbers());
        }
        else {
            return DTPointCollection3D(ConvertToDouble(p.FloatData()));
        }
    }
}

DTMutablePointCollection3D operator*(const DTTransform3D &T,const DTPointCollection3D &P)
{
    if (P.IsEmpty()) return DTMutablePointCollection3D();
    if (P.DoublePrecision()) {
        if (P.PointNumbersSpecified()) {
            return DTMutablePointCollection3D(TransformPoints(T,P.DoubleData()),P.PointNumbers());
        }
        else {
            return DTMutablePointCollection3D(TransformPoints(T,P.DoubleData()));
        }
    }
    else {
        if (P.PointNumbersSpecified()) {
            return DTMutablePointCollection3D(TransformPoints(T,P.FloatData()),P.PointNumbers());
        }
        else {
            return DTMutablePointCollection3D(TransformPoints(T,P.FloatData()));
        }
    }
}

DTMutablePointCollection3D operator-(const DTPointCollection3D &A,const DTPointCollection3D &B)
{
    if (A.NumberOfPoints()!=B.NumberOfPoints()) {
        DTErrorMessage("PointCollection3D-PointCollection3D","Incompatible lengths");
        return DTMutablePointCollection3D();
    }
    else if (A.DoublePrecision()!=B.DoublePrecision()) {
        DTErrorMessage("PointCollection::Overwite(PointCollection3D,Range)","Both have to be double or float");
        return DTMutablePointCollection3D();
    }

    if (A.DoublePrecision()) {
        return DTMutablePointCollection3D(A.DoubleData()-B.DoubleData());
    }
    else {
        return DTMutablePointCollection3D(A.FloatData()-B.FloatData());
    }
}

DTMutablePointCollection3D operator-(const DTPointCollection3D &A,const DTPoint3D &b)
{
    ssize_t i,howMany = A.NumberOfPoints();
    if (A.DoublePrecision()) {
        DTMutableDoubleArray toReturn(3,howMany);
        DTDoubleArray Aa = A.DoubleData();
        for (i=0;i<howMany;i++) {
            toReturn(0,i) = Aa(0,i) - b.x;
            toReturn(1,i) = Aa(1,i) - b.y;
            toReturn(2,i) = Aa(2,i) - b.z;
        }
        return DTMutablePointCollection3D(toReturn);
    }
    else {
        DTMutableFloatArray toReturn(3,howMany);
        DTFloatArray Aa = A.FloatData();
        float bx = float(b.x);
        float by = float(b.y);
        float bz = float(b.z);
        for (i=0;i<howMany;i++) {
            toReturn(0,i) = Aa(0,i) - bx;
            toReturn(1,i) = Aa(1,i) - by;
            toReturn(2,i) = Aa(2,i) - bz;
        }
        return DTMutablePointCollection3D(toReturn);
    }
}

DTRegion3D BoundingBox(const DTPointCollection3D &coll)
{
    if (coll.FloatPrecision()) {
        return BoundingBox3D(coll.FloatData());
    }
    else {
        return BoundingBox3D(coll.DoubleData());
    }
}

DTMutablePointCollection3D ExtractIndices(const DTPointCollection3D &A,const DTRange &r)
{
    if (r.start+r.length>A.NumberOfPoints()) {
        DTErrorMessage("ExtractIndices(PointCollection,Range)","Range is out of bounds");
        return DTMutablePointCollection3D();
    }
    if (A.PointNumbersSpecified())
        return DTMutablePointCollection3D(ExtractColumns(A.DoubleData(),r),ExtractIndices(A.PointNumbers(),r));
    else
        return DTMutablePointCollection3D(ExtractColumns(A.DoubleData(),r));
}

#pragma mark I/O Routines

void Write(DTDataStorage &output,const std::string &name,const DTPointCollection3D &theVar)
{
    Write(output,name+"_bbox3D",BoundingBox(theVar));
    if (theVar.PointNumbers().NotEmpty()) Write(output,name+"_ptN",theVar.PointNumbers());
    if (theVar.DoublePrecision())
        Write(output,name,theVar.DoubleData());
    else
        Write(output,name,theVar.FloatData());
}

void WriteFloat(DTDataStorage &output,const std::string &name,const DTPointCollection3D &theVar)
{
    Write(output,name+"_bbox3D",BoundingBox(theVar));
    if (theVar.PointNumbers().NotEmpty()) Write(output,name+"_ptN",theVar.PointNumbers());
    if (theVar.DoublePrecision())
        Write(output,name,ConvertToFloat(theVar.DoubleData()));
    else
        Write(output,name,theVar.FloatData());
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"3D Point Collection");
    output.Flush();
}

void Read(const DTDataStorage &input,const std::string &name,DTPointCollection3D &toReturn,DTPointCollection3D_SaveInfo &SaveInfo)
{
    std::string theName = input.ResolveName(name);
    if (SaveInfo.name==theName) {
        toReturn = SaveInfo.points; // Already read this grid in.
        return;
    }
    
    Read(input,theName,toReturn);

    SaveInfo.points = toReturn;
    SaveInfo.name = theName;
}

void Write(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite,DTPointCollection3D_SaveInfo &SaveInfo)
{
    if (SaveInfo.name.length() && SaveInfo.points==toWrite) {
        // Just save the reference.
        Write(output,name,SaveInfo.name);
    }
    else {
        Write(output,name,toWrite);
        SaveInfo.points = toWrite;
        SaveInfo.name = name;
    }
}

void WriteFloat(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite,DTPointCollection3D_SaveInfo &SaveInfo)
{
    if (SaveInfo.name.length() && SaveInfo.points==toWrite) {
        // Just save the reference.
        Write(output,name,SaveInfo.name);
    }
    else {
        WriteFloat(output,name,toWrite);
        SaveInfo.points = toWrite;
        SaveInfo.name = name;
    }
}

