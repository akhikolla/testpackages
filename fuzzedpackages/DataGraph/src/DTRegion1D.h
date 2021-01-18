// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTRegion1D_H
#define DTRegion1D_H

class DTDataStorage;
class DTIntArray;
class DTFloatArray;
class DTDoubleArray;

#include <string>
using namespace std;

struct DTRegion1D {
    DTRegion1D() : isSet(false), minV(0.0), maxV(0.0) {}
    DTRegion1D(double mn,double mx) : isSet(true), minV(mn<mx ? mn : mx), maxV(mn<mx ? mx : mn) {}

    bool isSet;
    double minV,maxV;

    void pinfo(void) const;
};

bool operator==(const DTRegion1D &,const DTRegion1D &);
bool operator!=(const DTRegion1D &,const DTRegion1D &);

extern DTRegion1D ValueRange(const DTIntArray &values);
extern DTRegion1D ValueRange(const DTFloatArray &values);
extern DTRegion1D ValueRange(const DTDoubleArray &values);

extern DTRegion1D Union(const DTRegion1D &,const DTRegion1D &);

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTRegion1D &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTRegion1D &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTRegion1D &toWrite); // One time value, self documenting.

#endif
