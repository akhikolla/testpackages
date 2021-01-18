// Part of DTSource. Copyright 2006. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTTransform3D_H
#define DTTransform3D_H

// [x]   [a b c sx][x]
// [y] = [d e f sy][y]
// [z]   [g h i sz][z]
// [1]   [0 0 0  1][1]

#include "DTDataStorage.h"
#include "DTPoint3D.h"
#include "DTVector3D.h"

class DTDoubleArray;
class DTFloatArray;
class DTMutableDoubleArray;
class DTMutableFloatArray;

// Used to map x -> A*x + b where A is a linear transformation and b is a shift vector.  Apply it to a point or vector
// as if it was a matrix transformation, i.e. if T is a DTTransform3D and x is a DTPoint3D or DTVector3D you can use the
// syntax T*x to apply the transformation.
// If you have two transformations S followed by T, you use T*S*x, which is equal to (T*S)*.
// You can take the inverse of a transformation by either using the function convention Inverse(T) or the method call
// T.Inverse()
// If you want to apply the transformation to a large number of points use
// the method TransformPoints(T,array).  This works for both DTDoubleArray and DTFloatArray, but requires the size to be a 3xN
// that is array.m()==3 and array.o()==1.

class DTTransform3D {
public:
    DTTransform3D();
    DTTransform3D(double t11,double t12,double t13,double s1,
                  double t21,double t22,double t23,double s2,
                  double t31,double t32,double t33,double s3);
    
    void pinfo(void) const;
    
    DTPoint3D operator*(const DTPoint3D &P) const {double scale = 1.0/(1.0+P.x*P1+P.y*P2+P.z*P3); return DTPoint3D((T11*P.x+T12*P.y+T13*P.z+S1)*scale,
                                                                                                                   (T21*P.x+T22*P.y+T23*P.z+S2)*scale,
                                                                                                                   (T31*P.x+T32*P.y+T33*P.z+S3)*scale);}
    DTVector3D operator*(const DTVector3D &P) const {double scale = 1.0/(1.0+P.x*P1+P.y*P2+P.z*P3); return DTVector3D((T11*P.x+T12*P.y+T13*P.z)*scale,
                                                                                                                      (T21*P.x+T22*P.y+T23*P.z)*scale,
                                                                                                                      (T31*P.x+T32*P.y+T33*P.z)*scale);}
    DTTransform3D operator*(const DTTransform3D &) const;
    
    bool IsAShift(void) const {return (T11==1.0 && T12==0.0 && T13==0.0 &&
                                       T21==0.0 && T22==1.0 && T23==0.0 &&
                                       T31==0.0 && T32==0.0 && T33==1.0 &&
                                       P1==0.0 && P2==0.0 && P3==0.0);}
        
    bool IsOrthogonal(void) const;
    DTTransform3D Inverse(void) const;
    DTPoint3D ShiftVector(void) const {return DTPoint3D(S1,S2,S3);}
    DTTransform3D RotationOnly(void) const {return DTTransform3D(T11,T12,T13,0.0,T21,T22,T23,0.0,T31,T32,T33,0.0);}
    
    // Access the underlying 4x4 matrix
    // [A s]
    // [0 1]
    double operator()(int i,int j) const;
    
    static DTTransform3D RotateAroundX(double);
    static DTTransform3D RotateAroundY(double);
    static DTTransform3D RotateAroundZ(double);
    static DTTransform3D RotateAroundVector(const DTVector3D &,double);

    static DTTransform3D Shift(const DTPoint3D &);
    static DTTransform3D Scale(double,double,double);
    static DTTransform3D Scale(double);

private:
    double T11,T12,T13,S1;
    double T21,T22,T23,S2;
    double T31,T32,T33,S3;
    double P1,P2,P3;
};

extern DTTransform3D TransformFromCoordinateToGrid(const DTPoint3D &origin,double dx,double dy,double dz);
extern DTTransform3D TransformFromGridToCoordinate(const DTPoint3D &origin,double dx,double dy,double dz);
extern DTTransform3D Inverse(const DTTransform3D &);

extern DTMutableDoubleArray TransformPoints(const DTTransform3D &,const DTDoubleArray &); // 3xN array or empty.
extern DTMutableFloatArray TransformPoints(const DTTransform3D &,const DTFloatArray &); // 3xN array or empty.

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTTransform3D &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTTransform3D &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTTransform3D &toWrite); // One time value, self documenting.

#endif
