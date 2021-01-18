// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTRegion3D.h"

#include "DTDoubleArray.h"
#include "DTDataStorage.h"
#include "DTFloatArray.h"
#include "DTError.h"

#include <math.h>
#include <limits>

void DTRegion3D::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    std::cerr << "[" << xmin << ", " << xmax << "] x [" << ymin << ", " << ymax << "] x [" << zmin << ", " << zmax << "]\n" << std::flush;
#endif
}

bool operator==(const DTRegion3D &A,const DTRegion3D &B)
{
    return (A.isSet==B.isSet && A.xmin==B.xmin && A.xmax==B.xmax && A.ymin==B.ymin && A.ymax==B.ymax && A.zmin==B.zmin && A.zmax==B.zmax);
}

bool operator!=(const DTRegion3D &A,const DTRegion3D &B)
{
    return (A.isSet!=B.isSet || A.xmin!=B.xmin || A.xmax!=B.xmax || A.ymin!=B.ymin || A.ymax!=B.ymax || A.zmin!=B.zmin || A.zmax!=B.zmax);
}

DTRegion3D BoundingBox3D(const DTFloatArray &points)
{
    if (points.IsEmpty())
        return DTRegion3D();
    if (points.m()%3!=0) {
        DTErrorMessage("BoundingBox(Array)","The array had an incorrect size.");
        return DTRegion3D();
    }

    ssize_t len = points.Length();

#if defined(WIN32) && !defined(INFINITY)
#define INFINITY std::numeric_limits<float>::infinity();
#endif
    float xmin,ymin,zmin,xmax,ymax,zmax;
    xmin = ymin = zmin = INFINITY;
    xmax = ymax = zmax = -INFINITY;

    float x,y,z;
    ssize_t i;

    const float *pointD = points.Pointer();

    for (i=0;i<len;i+=3) {
        x = pointD[i];
        y = pointD[i+1];
        z = pointD[i+2];
        if (isfinite(x) & isfinite(y) & isfinite(z)) {
            // Non-finite points are not included in the bounding box.
            xmin = (x < xmin ? x : xmin);
            xmax = (xmax < x ? x : xmax);
            ymin = (y < ymin ? y : ymin);
            ymax = (ymax < y ? y : ymax);
            zmin = (z < zmin ? z : zmin);
            zmax = (zmax < z ? z : zmax);
        }
    }

    if (xmin>xmax || ymin>ymax || zmin>zmax)
        return DTRegion3D();
    else
        return DTRegion3D(DTPoint3D(xmin,ymin,zmin),DTPoint3D(xmax,ymax,zmax));
}

DTRegion3D BoundingBox3D(const DTDoubleArray &points)
{
    if (points.IsEmpty())
        return DTRegion3D();
    if (points.m()%3!=0) {
        DTErrorMessage("BoundingBox(Array)","The array had an incorrect size.");
        return DTRegion3D();
    }
    
    ssize_t len = points.Length();
    
#if defined(WIN32) && !defined(INFINITY)
#define INFINITY std::numeric_limits<float>::infinity();
#endif
    double xmin,ymin,zmin,xmax,ymax,zmax;
    xmin = ymin = zmin = INFINITY;
    xmax = ymax = zmax = -INFINITY;
    
    double x,y,z;
    size_t i;
    
    const double *pointD = points.Pointer();
    
    for (i=0;i<len;i+=3) {
        x = pointD[i];
        y = pointD[i+1];
        z = pointD[i+2];
        if (isfinite(x) & isfinite(y) & isfinite(z)) {
            // Non-finite points are not included in the bounding box.
            xmin = (x < xmin ? x : xmin);
            xmax = (xmax < x ? x : xmax);
            ymin = (y < ymin ? y : ymin);
            ymax = (ymax < y ? y : ymax);
            zmin = (z < zmin ? z : zmin);
            zmax = (zmax < z ? z : zmax);
        }
    }
    
    if (xmin>xmax || ymin>ymax || zmin>zmax)
        return DTRegion3D();
    else
        return DTRegion3D(DTPoint3D(xmin,ymin,zmin),DTPoint3D(xmax,ymax,zmax));
}

DTRegion3D Offset(const DTRegion3D &A,double v)
{
	if (A.isSet==false) return A;
	
	DTRegion3D toReturn = A;
	toReturn.xmin -= v;
	toReturn.ymin -= v;
	toReturn.zmin -= v;
	toReturn.xmax += v;
	toReturn.ymax += v;
	toReturn.zmax += v;
	
	if (toReturn.xmin>toReturn.xmax || toReturn.ymin>toReturn.ymax || toReturn.zmin>toReturn.zmax)
		return DTRegion3D();
	
	return toReturn;
}

DTRegion3D Union(const DTRegion3D &A,const DTRegion3D &B)
{
    if (A.isSet==false) return B;
    if (B.isSet==false) return A;

    DTPoint3D ll(std::min(A.xmin,B.xmin),std::min(A.ymin,B.ymin),std::min(A.zmin,B.zmin));
    DTPoint3D ur(std::max(A.xmax,B.xmax),std::max(A.ymax,B.ymax),std::max(A.zmax,B.zmax));
        
    return DTRegion3D(ll,ur);
}

DTRegion3D Intersection(const DTRegion3D &A,const DTRegion3D &B)
{
	double xmin = std::max(A.xmin,B.xmin);
	double xmax = std::min(A.xmax,B.xmax);
	double ymin = std::max(A.ymin,B.ymin);
	double ymax = std::min(A.ymax,B.ymax);
	double zmin = std::max(A.zmin,B.zmin);
	double zmax = std::min(A.zmax,B.zmax);
	
	if (xmax<xmin || ymax<ymin || zmax<zmin) return DTRegion3D(DTPoint3D(0,0,0),DTPoint3D(0,0,0));
	return DTRegion3D(DTPoint3D(xmin,ymin,zmin),DTPoint3D(xmax,ymax,zmax));
}

bool BoxesIntersect(const DTRegion3D &A,const DTRegion3D &B)
{
	if (A.isSet==false || B.isSet==false) return false;
	double xmin = std::max(A.xmin,B.xmin);
	double xmax = std::min(A.xmax,B.xmax);
	double ymin = std::max(A.ymin,B.ymin);
	double ymax = std::min(A.ymax,B.ymax);
	double zmin = std::max(A.zmin,B.zmin);
	double zmax = std::min(A.zmax,B.zmax);
	
	if (xmax<=xmin || ymax<=ymin || zmax<=zmin) return false;
	return true;
}

void Read(const DTDataStorage &input,const std::string &name,DTRegion3D &region)
{
    DTDoubleArray onDisk = input.ReadDoubleArray(name);
    region = DTRegion3D();
    if (onDisk.Length()==6) {
        region.isSet = true;
        region.xmin = onDisk(0);
        region.xmax = onDisk(1);
        region.ymin = onDisk(2);
        region.ymax = onDisk(3);
        region.zmin = onDisk(4);
        region.zmax = onDisk(5);
    }
}

void Write(DTDataStorage &output,const std::string &name,const DTRegion3D &region)
{
    DTMutableDoubleArray toSave;
    
    if (region.isSet) {
        toSave = DTMutableDoubleArray(6);
        toSave(0) = region.xmin;
        toSave(1) = region.xmax;
        toSave(2) = region.ymin;
        toSave(3) = region.ymax;
        toSave(4) = region.zmin;
        toSave(5) = region.zmax;
    }

    output.Save(toSave,name);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTRegion3D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"3D Region");
    output.Flush();
}

