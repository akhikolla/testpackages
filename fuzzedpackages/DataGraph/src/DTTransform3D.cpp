// Part of DTSource. Copyright 2006. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTTransform3D.h"
#include "DTVector3D.h"
#include "DTError.h"
#include "DTDoubleArray.h"
#include "DTFloatArray.h"

DTTransform3D::DTTransform3D()
{
    T11 = 1.0; T12 = 0.0; T13 = 0.0; S1 = 0.0;
    T21 = 0.0; T22 = 1.0; T23 = 0.0; S2 = 0.0;
    T31 = 0.0; T32 = 0.0; T33 = 1.0; S3 = 0.0;
    P1 = 0.0; P2 = 0.0; P3 = 0.0;
}

DTTransform3D::DTTransform3D(double t11,double t12,double t13,double s1,
                             double t21,double t22,double t23,double s2,
                             double t31,double t32,double t33,double s3)
{
    T11 = t11; T12 = t12; T13 = t13; S1 = s1;
    T21 = t21; T22 = t22; T23 = t23; S2 = s2;
    T31 = t31; T32 = t32; T33 = t33; S3 = s3;
    P1 = 0.0; P2 = 0.0; P3 = 0.0;
}

void DTTransform3D::pinfo() const
{
#ifndef DG_NOSTDErrOut
    std::cerr << T11 << ", " << T12 << ", " << T13 << ", " << S1 << "\n";
    std::cerr << T21 << ", " << T22 << ", " << T23 << ", " << S2 << "\n";
    std::cerr << T31 << ", " << T32 << ", " << T33 << ", " << S3 << "\n";
    if (P1 || P2 || P3) {
        std::cerr << P1 << ", " << P2 << ", " << P3 << ", 1\n";
    }
#endif
}

DTTransform3D DTTransform3D::operator*(const DTTransform3D &T) const
{
    return DTTransform3D(T11*T.T11+T12*T.T21+T13*T.T31,
                         T11*T.T12+T12*T.T22+T13*T.T32,
                         T11*T.T13+T12*T.T23+T13*T.T33,
                         T11*T.S1 +T12*T.S2+T13*T.S3+S1,
                         T21*T.T11+T22*T.T21+T23*T.T31,
                         T21*T.T12+T22*T.T22+T23*T.T32,
                         T21*T.T13+T22*T.T23+T23*T.T33,
                         T21*T.S1 +T22*T.S2+T23*T.S3+S2,
                         T31*T.T11+T32*T.T21+T33*T.T31,
                         T31*T.T12+T32*T.T22+T33*T.T32,
                         T31*T.T13+T32*T.T23+T33*T.T33,
                         T31*T.S1 +T32*T.S2+T33*T.S3+S3);
}

bool DTTransform3D::IsOrthogonal() const
{
    double ip12 = T11*T12 + T21*T22 + T31*T32;
    double ip13 = T11*T13 + T21*T23 + T31*T33;
    double ip23 = T12*T13 + T22*T23 + T32*T33;
    
    if (fabs(ip12)<1e-15 && fabs(ip13)<1e-15 && fabs(ip23)<1e-15) {
        double n1 = T11*T11 + T21*T21 + T31*T31;
        double n2 = T12*T12 + T22*T22 + T32*T32;
        double n3 = T13*T13 + T23*T23 + T33*T33;
        return (fabs(n1-1.0)<1e-14 && fabs(n2-1.0)<1e-14 && fabs(n3-1.0)<1e-14);
    }
    else {
        return false;
    }
}

DTTransform3D DTTransform3D::Inverse() const
{
    // Use Cramers rule to find the inverse.
    double overallDet = T11*(T22*T33 - T23*T32) - T12*(T21*T33 - T31*T23) + T13*(T21*T32 - T31*T22);
	
	double t11 =  (T22*T33 - T32*T23)/overallDet;
	double t21 = -(T21*T33 - T31*T23)/overallDet;
	double t31 =  (T21*T32 - T31*T22)/overallDet;
    
	double t12 = -(T12*T33 - T32*T13)/overallDet;
	double t22 =  (T11*T33 - T31*T13)/overallDet;
	double t32 = -(T11*T32 - T31*T12)/overallDet;
	
	double t13 =  (T12*T23 - T22*T13)/overallDet;
	double t23 = -(T11*T23 - T21*T13)/overallDet;
	double t33 =  (T11*T22 - T21*T12)/overallDet;
	
	// The mapping is really 
	
	// x -> A*x + b = y
	// x = A^-1 * (y - b) = A^-1 * y - A^-1*b
	// so the new column is - A^-1 * old column
	
	double s1 = - (t11*S1 + t12*S2 + t13*S3);
	double s2 = - (t21*S1 + t22*S2 + t23*S3);
	double s3 = - (t31*S1 + t32*S2 + t33*S3);
	
	return DTTransform3D(t11,t12,t13,s1,
                         t21,t22,t23,s2,
                         t31,t32,t33,s3);
}

double DTTransform3D::operator()(int i,int j) const
{
    if (i<0 || j<0 || i>3 || j>3) {
        DTErrorMessage("DTTransform::operator(i,j)","Index out of bounds");
        return 0.0;
    }
    
    double toReturn = 0.0;
    switch (i+j*4) {
        case 0:  toReturn = T11; break;
        case 1:  toReturn = T21; break;
        case 2:  toReturn = T31; break;
        case 3:  toReturn = 0.0; break;
        case 4:  toReturn = T12; break;
        case 5:  toReturn = T22; break;
        case 6:  toReturn = T32; break;
        case 7:  toReturn = 0.0; break;
        case 8:  toReturn = T13; break;
        case 9:  toReturn = T23; break;
        case 10: toReturn = T33; break;
        case 11: toReturn = 0.0; break;
        case 12: toReturn =  S1; break;
        case 13: toReturn =  S2; break;
        case 14: toReturn =  S3; break;
        case 15: toReturn = 1.0; break;
        default: toReturn = 0.0; break;
    }
    
    return toReturn;
}

DTTransform3D DTTransform3D::RotateAroundX(double rad)
{
    return DTTransform3D(1.0,     0.0,      0.0,0.0,
                         0.0,cos(rad),-sin(rad),0.0,
                         0.0,sin(rad), cos(rad),0.0);
}

DTTransform3D DTTransform3D::RotateAroundY(double rad)
{
    return DTTransform3D(cos(rad) ,0.0, sin(rad),0.0,
                         0.0      ,1.0,      0.0,0.0,
                         -sin(rad),0.0, cos(rad),0.0);
    
}

DTTransform3D DTTransform3D::RotateAroundZ(double rad)
{
    return DTTransform3D(cos(rad),-sin(rad),0.0,0.0,
                         sin(rad), cos(rad),0.0,0.0,
                         0.0,      0.0     ,1.0,0.0);
}

DTTransform3D DTTransform3D::RotateAroundVector(const DTVector3D &vector,double angle)
{
    double length = Norm(vector);
    if (length==0 || isfinite(length)==0) {
        DTErrorMessage("DTTransform3D::RotateAroundVector","Invalid vector");
        return DTTransform3D();
    }
    if (angle==0) return DTTransform3D();
    double nx = vector.x/length;
    double ny = vector.y/length;
    double nz = vector.z/length;
    
    double c = cos(angle);
    double s = sin(angle);
    double oneMc = 1-c;
    
    return DTTransform3D(nx*nx*oneMc +    c,nx*ny*oneMc - nz*s,nx*nz*oneMc + ny*s,0.0,
                         nx*ny*oneMc + nz*s,ny*ny*oneMc +    c,ny*nz*oneMc - nx*s,0.0,
                         nx*nz*oneMc - ny*s,ny*nz*oneMc + nx*s,nz*nz*oneMc +    c,0.0);

}

DTTransform3D DTTransform3D::Shift(const DTPoint3D &p)
{
    return DTTransform3D(1.0,0.0,0.0,p.x,
                         0.0,1.0,0.0,p.y,
                         0.0,0.0,1.0,p.z);
}

DTTransform3D DTTransform3D::Scale(double xs,double ys,double zs)
{
    return DTTransform3D(xs ,0.0,0.0,0.0,
                         0.0, ys,0.0,0.0,
                         0.0,0.0, zs,0.0);
}

DTTransform3D DTTransform3D::Scale(double s)
{
    return DTTransform3D(s ,0.0,0.0,0.0,
                         0.0, s,0.0,0.0,
                         0.0,0.0, s,0.0);
}

DTTransform3D TransformFromCoordinateToGrid(const DTPoint3D &origin,double dx,double dy,double dz)
{
    // i -> i*dx + xzero
    return DTTransform3D(dx,  0,  0, origin.x,
                         0,  dy,  0, origin.y,
                         0,   0, dz, origin.z);
}

DTTransform3D TransformFromGridToCoordinate(const DTPoint3D &origin,double dx,double dy,double dz)
{
    // x -> (x-xzero)/dx
    return DTTransform3D(1.0/dx,      0,      0, -origin.x/dx,
                         0,      1.0/dy,      0, -origin.y/dx,
                         0,           0, 1.0/dz,-origin.z/dy);
}

DTTransform3D Inverse(const DTTransform3D &T)
{
    return T.Inverse();
}

void Read(const DTDataStorage &input,const std::string &name,DTTransform3D &toReturn)
{
    DTDoubleArray theArray = input.ReadDoubleArray(name);
    if (theArray.Length()==0) {
        toReturn = DTTransform3D();
    }
    else if (theArray.m()!=3 || theArray.n()!=4 || theArray.Length()!=12) {
        DTErrorMessage("ReadFromArray(Transform3D)","Invalid size for array.");
        toReturn = DTTransform3D();
    }
    else {
        toReturn = DTTransform3D(theArray(0,0),theArray(0,1),theArray(0,2),theArray(0,3),
                                 theArray(1,0),theArray(1,1),theArray(1,2),theArray(1,3),
                                 theArray(2,0),theArray(2,1),theArray(2,2),theArray(2,3));
    }
}

DTMutableDoubleArray TransformPoints(const DTTransform3D &T,const DTDoubleArray &A)
{
    if (A.IsEmpty()) return DTMutableDoubleArray();
    if (A.m()!=3 || A.o()!=1) {
        DTErrorMessage("TransformPoints(Transform3D,DoubleArray","Invalid array size");
        return DTMutableDoubleArray();
    }
    
    size_t i,howMany = A.n();
    DTMutableDoubleArray toReturn(3,howMany);
    double x,y,z;
    
    double T11 = T(0,0);
    double T21 = T(1,0);
    double T31 = T(2,0);
    double T12 = T(0,1);
    double T22 = T(1,1);
    double T32 = T(2,1);
    double T13 = T(0,2);
    double T23 = T(1,2);
    double T33 = T(2,2);
    double S1 = T(0,3);
    double S2 = T(1,3);
    double S3 = T(2,3);
    
    for (i=0;i<howMany;i++) {
        x = A(0,i);
        y = A(1,i);
        z = A(2,i);
        toReturn(0,i) = T11*x+T12*y+T13*z + S1;
        toReturn(1,i) = T21*x+T22*y+T23*z + S2;
        toReturn(2,i) = T31*x+T32*y+T33*z + S3;
    }
    
    return toReturn;
}

DTMutableFloatArray TransformPoints(const DTTransform3D &T,const DTFloatArray &A)
{
    if (A.IsEmpty()) return DTMutableFloatArray();
    if (A.m()!=3 || A.o()!=1) {
        DTErrorMessage("TransformPoints(Transform3D,FloatArray","Invalid array size");
        return DTMutableFloatArray();
    }
    
    size_t i,howMany = A.n();
    DTMutableFloatArray toReturn(3,howMany);
    double x,y,z;
    
    double T11 = T(0,0);
    double T21 = T(1,0);
    double T31 = T(2,0);
    double T12 = T(0,1);
    double T22 = T(1,1);
    double T32 = T(2,1);
    double T13 = T(0,2);
    double T23 = T(1,2);
    double T33 = T(2,2);
    double S1 = T(0,3);
    double S2 = T(1,3);
    double S3 = T(2,3);
    
    for (i=0;i<howMany;i++) {
        x = A(0,i);
        y = A(1,i);
        z = A(2,i);
        toReturn(0,i) = float(T11*x+T12*y+T13*z + S1);
        toReturn(1,i) = float(T21*x+T22*y+T23*z + S2);
        toReturn(2,i) = float(T31*x+T32*y+T33*z + S3);
    }
    
    return toReturn;
}

void Write(DTDataStorage &output,const std::string &name,const DTTransform3D &theVar)
{
    DTMutableDoubleArray theArr(3,4);
    for (int j=0;j<4;j++) {
        for (int i=0;i<3;i++) {
            theArr(i,j) = theVar(i,j);
        }
    }
    output.Save(theArr,name);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTTransform3D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"3D Transform");
    output.Flush();
}

