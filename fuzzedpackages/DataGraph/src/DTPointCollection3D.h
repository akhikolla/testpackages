// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTPointCollection3D_H
#define DTPointCollection3D_H

#include "DTFloatArray.h"
#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTPoint3D.h"

#include <string>
using namespace std;

class DTDataStorage;
struct DTRegion3D;
struct DTRange;
class DTTransform3D;
class DTMutablePointCollection3D;

class DTPointCollection3D {
public:
    DTPointCollection3D() : floatData(), doubleData(), pointNumbers() {}
    explicit DTPointCollection3D(const DTFloatArray &input);
    explicit DTPointCollection3D(const DTDoubleArray &input);
    DTPointCollection3D(const DTFloatArray &input,const DTIntArray &ptNumbers);
    DTPointCollection3D(const DTDoubleArray &input,const DTIntArray &ptNumbers);
    
    bool IsEmpty(void) const {return (floatData.IsEmpty() && doubleData.IsEmpty());}
    ssize_t NumberOfPoints(void) const;
    ssize_t MemoryUsed(void) const {return floatData.MemoryUsed()+doubleData.MemoryUsed()+pointNumbers.MemoryUsed();}

    DTMutablePointCollection3D Copy(void) const;
    
    bool FloatPrecision(void) const {return (floatData.NotEmpty());}
    bool DoublePrecision(void) const {return (doubleData.NotEmpty());}
    
    DTDoubleArray DoubleData(void) const;
    DTFloatArray FloatData(void) const;
    
    DTFloatArray Data(void) const; // Obsolete (9/20/14)
    DTIntArray PointNumbers(void) const {return pointNumbers;}
    bool PointNumbersSpecified(void) const {return pointNumbers.NotEmpty();}

    DTPoint3D operator()(ssize_t) const;

    void pinfo(void) const;
	void pall(void) const;

protected:
    DTFloatArray floatData;
    DTDoubleArray doubleData;
    DTIntArray pointNumbers; // Empty or N entries
};

class DTMutablePointCollection3D : public DTPointCollection3D {
public:
    DTMutablePointCollection3D() : DTPointCollection3D() {}
    explicit DTMutablePointCollection3D(const DTMutableDoubleArray &input);
    explicit DTMutablePointCollection3D(const DTMutableFloatArray &input);
    DTMutablePointCollection3D(const DTMutableDoubleArray &input,const DTIntArray &);
    DTMutablePointCollection3D(const DTMutableFloatArray &input,const DTIntArray &);

    DTMutableDoubleArray DoubleData(void) {return mutableDoubleData;}
    DTDoubleArray DoubleData(void) const {return doubleData;}

    void Overwrite(const DTPointCollection3D &,const DTRange &);

    void operator-=(DTPoint3D); // Subtract a point from every location
    void operator+=(const DTPointCollection3D &);
    void operator-=(const DTPointCollection3D &);

private:
    DTMutableFloatArray mutableFloatData;
    DTMutableDoubleArray mutableDoubleData;
};


extern bool operator==(const DTPointCollection3D &,const DTPointCollection3D &);
extern bool operator!=(const DTPointCollection3D &,const DTPointCollection3D &);

extern DTMutablePointCollection3D operator*(const DTTransform3D &,const DTPointCollection3D &);

extern DTMutablePointCollection3D operator-(const DTPointCollection3D &,const DTPointCollection3D &);
extern DTMutablePointCollection3D operator-(const DTPointCollection3D &,const DTPoint3D &);
extern DTMutablePointCollection3D ExtractIndices(const DTPointCollection3D &,const DTRange &);

extern DTRegion3D BoundingBox(const DTPointCollection3D &);

extern DTPointCollection3D ConvertToFloat(const DTPointCollection3D &);
extern DTPointCollection3D ConvertToDouble(const DTPointCollection3D &);

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTPointCollection3D &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTPointCollection3D &theVar);
extern void WriteFloat(DTDataStorage &output,const std::string &name,const DTPointCollection3D &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite);

// A utility class to enable grid sharing.
struct DTPointCollection3D_SaveInfo {
    DTPointCollection3D points;
    std::string name;
};

// Pass in as the third argument to Read and Write to enable the the read/write routines
// to avoid having to read or write identical grids.  This is used when reading DTTrianglarMesh2D and DTTriangularVectorField2D
extern void Read(const DTDataStorage &input,const std::string &name,DTPointCollection3D &toReturn,DTPointCollection3D_SaveInfo &);
extern void Write(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite,DTPointCollection3D_SaveInfo &);
extern void WriteFloat(DTDataStorage &output,const std::string &name,const DTPointCollection3D &toWrite,DTPointCollection3D_SaveInfo &);

#endif
