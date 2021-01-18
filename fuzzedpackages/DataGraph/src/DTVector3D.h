// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTVector3D_H
#define DTVector3D_H

#include <math.h>
#include "DTPoint3D.h"
#include <iostream>

struct DTVector3D {
    DTVector3D() : x(0.0), y(0.0), z(0.0) {}
    DTVector3D(double xv,double yv,double zv) : x(xv), y(yv), z(zv) {}
    explicit DTVector3D(const DTPoint3D &p) : x(p.x), y(p.y), z(p.z) {}
    
    void operator*=(double s) {x*=s; y*=s; z*=s;}
    
    void pinfo(void) const {std::cerr << "(" << x << ", " << y << ", " << z << ")\n" << std::flush;}

    double x,y,z;
};

inline bool operator==(const DTVector3D &A,const DTVector3D &B) {return (A.x==B.x && A.y==B.y && A.z==B.z);}
inline bool operator!=(const DTVector3D &A,const DTVector3D &B) {return (A.x!=B.x || A.y!=B.y || A.z!=B.z);}

inline DTVector3D operator+(const DTVector3D &A,const DTVector3D &B) {return DTVector3D(A.x+B.x,A.y+B.y,A.z+B.z);}
inline DTVector3D operator-(const DTVector3D &A,const DTVector3D &B) {return DTVector3D(A.x-B.x,A.y-B.y,A.z-B.z);}
inline DTVector3D operator*(const DTVector3D &A,double v) {return DTVector3D(A.x*v,A.y*v,A.z*v);}
inline DTVector3D operator*(double v,const DTVector3D &A) {return DTVector3D(A.x*v,A.y*v,A.z*v);}
inline DTVector3D operator/(const DTVector3D &A,double v) {return DTVector3D(A.x/v,A.y/v,A.z/v);}
inline DTVector3D operator-(const DTVector3D &A) {return DTVector3D(-A.x,-A.y,-A.z);}

inline DTPoint3D operator-(const DTPoint3D &P,const DTVector3D &v) {return DTPoint3D(P.x-v.x,P.y-v.y,P.z-v.z);}
inline DTPoint3D operator+(const DTPoint3D &P,const DTVector3D &v) {return DTPoint3D(P.x+v.x,P.y+v.y,P.z+v.z);}

inline double Norm(const DTVector3D &A) {return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);}
inline DTVector3D Cross(const DTVector3D &A,const DTVector3D &B) {return DTVector3D(A.y*B.z-A.z*B.y,B.x*A.z-A.x*B.z,A.x*B.y-A.y*B.x);}
inline double Dot(const DTVector3D &A,const DTVector3D &B) {return (A.x*B.x+A.y*B.y+A.z*B.z);}

#endif
