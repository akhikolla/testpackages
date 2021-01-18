// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTPoint3D_H
#define DTPoint3D_H

#include "DTDataStorage.h"
#include <math.h>

struct DTPoint3D {
    DTPoint3D() : x(0.0), y(0.0), z(0.0) {}
    DTPoint3D(double xv,double yv,double zv) : x(xv), y(yv), z(zv) {}

    void pinfo(void) const;
    
    double x,y,z;
};

inline DTPoint3D operator+(const DTPoint3D &A,const DTPoint3D &B) {return DTPoint3D(A.x+B.x,A.y+B.y,A.z+B.z);}
inline DTPoint3D operator-(const DTPoint3D &A,const DTPoint3D &B) {return DTPoint3D(A.x-B.x,A.y-B.y,A.z-B.z);}
inline DTPoint3D operator*(const DTPoint3D &A,double v) {return DTPoint3D(A.x*v,A.y*v,A.z*v);}
inline DTPoint3D operator*(double v,const DTPoint3D &A) {return DTPoint3D(A.x*v,A.y*v,A.z*v);}
inline DTPoint3D operator/(const DTPoint3D &A,double v) {return DTPoint3D(A.x/v,A.y/v,A.z/v);}
inline DTPoint3D operator-(const DTPoint3D &A) {return DTPoint3D(-A.x,-A.y,-A.z);}

inline bool operator==(const DTPoint3D &A,const DTPoint3D &B) {return (A.x==B.x && A.y==B.y && A.z==B.z);}
inline bool operator!=(const DTPoint3D &A,const DTPoint3D &B) {return (A.x!=B.x || A.y!=B.y || A.z!=B.z);}

inline double Norm(const DTPoint3D &A) {return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);}
inline DTPoint3D Cross(const DTPoint3D &A,const DTPoint3D &B) {return DTPoint3D(A.y*B.z-A.z*B.y,B.x*A.z-A.x*B.z,A.x*B.y-A.y*B.x);}
inline double Dot(const DTPoint3D &A,const DTPoint3D &B) {return (A.x*B.x+A.y*B.y+A.z*B.z);}
inline double Distance(const DTPoint3D &A,const DTPoint3D &B) {return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));}

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTPoint3D &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTPoint3D &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTPoint3D &toWrite); // One time value, self documenting.

#endif
