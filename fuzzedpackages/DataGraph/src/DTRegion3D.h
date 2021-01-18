// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTRegion3D_H
#define DTRegion3D_H

#include "DTDataStorage.h"
#include "DTPoint3D.h"

class DTFloatArray;
class DTDoubleArray;

struct DTRegion3D {
    DTRegion3D() : isSet(false), xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0) {}
    DTRegion3D(DTPoint3D ll,DTPoint3D ur) : isSet(true), xmin(ll.x), xmax(ur.x), ymin(ll.y), ymax(ur.y), zmin(ll.z), zmax(ur.z) {}

    DTPoint3D Minimum(void) const {return DTPoint3D(xmin,ymin,zmin);}
    DTPoint3D Maximum(void) const {return DTPoint3D(xmax,ymax,zmax);}
    
    bool IsEmpty(void) const {return (xmin==xmax && ymin==ymax && zmin==zmax);}
	
    void pinfo(void) const;
	
	bool Contains(const DTRegion3D &a) const {return (xmin<=a.xmin && xmax>=a.xmax && ymin<=a.ymin && ymax>=a.ymax && zmin<=a.zmin && zmax>=a.zmax);}
    
    bool isSet;
    double xmin,xmax,ymin,ymax,zmin,zmax;
};

extern bool operator==(const DTRegion3D &,const DTRegion3D &);
extern bool operator!=(const DTRegion3D &,const DTRegion3D &);

extern DTRegion3D BoundingBox3D(const DTFloatArray &points);
extern DTRegion3D BoundingBox3D(const DTDoubleArray &points);

extern DTRegion3D Offset(const DTRegion3D &,double);
extern DTRegion3D Union(const DTRegion3D &,const DTRegion3D &);
extern DTRegion3D Intersection(const DTRegion3D &,const DTRegion3D &);

extern bool BoxesIntersect(const DTRegion3D &,const DTRegion3D &);

extern void Read(const DTDataStorage &input,const std::string &name,DTRegion3D &region);
extern void Write(DTDataStorage &output,const std::string &name,const DTRegion3D &region);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTRegion3D &toWrite); // One time value, self documenting.

#endif
