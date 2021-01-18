// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTSurface3D_H
#define DTSurface3D_H

#include "DTFloatArray.h"
#include "DTIntArray.h"
#include "DTPointCollection3D.h"

#include "DTDataStorage.h"
struct DTRegion3D;
struct DTPoint3D;
class DTTransform3D;

// This is a triangulated surface.
// The vertices are stored as a 3xN array of xyz values.
// The triangles are stored as a 3xM integer array of offsets.
// The normals are stored as a 3xN array of (nx,ny,nz) values, but are optional.

// Surfaces coming from DataTank are saved in this format, but if a surface is intended
// to only be read by DataTank, you can create the DTSurface3D object from a 9xM array of floating points
// where the points have to be repeated.

class DTSurface3D {
public:
    DTSurface3D() : points(), connections(), normals(), triangles() {}
    DTSurface3D(const DTIntArray &conn,const DTFloatArray &pts);
    DTSurface3D(const DTIntArray &conn,const DTPointCollection3D &pts);
    DTSurface3D(const DTIntArray &conn,const DTPointCollection3D &pts,const DTFloatArray &nrm);
    explicit DTSurface3D(const DTFloatArray &); // 9xN array
    // Next constructor is only used by Add... functions.  Do not use, since it does not check input.
    DTSurface3D(const DTIntArray &,const DTPointCollection3D &,const DTFloatArray &,const DTIntArray &);

    bool IsEmpty(void) const {return (connections.IsEmpty() && triangles.IsEmpty());}
    ssize_t MemoryUsed(void) const {return points.MemoryUsed()+connections.MemoryUsed()+normals.MemoryUsed()+nextTriangles.MemoryUsed()+triangles.MemoryUsed();}

    bool SavedAsOffsets(void) const {return triangles.IsEmpty();}
    bool SavedAsTriangles(void) const {return triangles.NotEmpty();}
    
    // Surface saved as points, connections, normals.  This is how they are saved from DataTank
    DTPointCollection3D Points(void) const;
    DTIntArray Connections(void) const;
    
    ssize_t NumberOfPoints(void) const {return points.NumberOfPoints();}
    ssize_t NumberOfTriangles(void) const {return connections.n();}

    bool NormalsDefined(void) const {return normals.NotEmpty();}
    DTFloatArray Normals(void) const;
    
    // Gives the index of the neigboring triangle.  If the points in the triangle are
    // A = offset(0,i), B = offset(1,i), C = offset(2,i) then
    // the next triangle offsets refer to the edges AB, BC and CA in that order.
    // nextTriangle(0,A) = 3*(triangle that AB connects to) + (which edge in that triangle)
    // So to get the next triangle number you use nextTriangle(0,A)/3
    bool NextTrianglesDefined(void) const {return (!nextTriangles.IsEmpty());}
    DTIntArray NextTriangles() const;
    
    // Raw 9xN array.
    DTFloatArray Triangles(void) const {return triangles;}

    void pinfo(void) const;
    void pt(int i) const; // Print a single triangle.

private:
    DTPointCollection3D points;
    DTIntArray connections;
    
    // Optional data.
    DTFloatArray normals;
    DTIntArray nextTriangles;

    DTFloatArray triangles;
};

extern bool operator==(const DTSurface3D &,const DTSurface3D &);
extern bool operator!=(const DTSurface3D &,const DTSurface3D &);

extern DTSurface3D ConvertToFloat(const DTSurface3D &);
extern DTSurface3D ConvertToDouble(const DTSurface3D &);

extern DTSurface3D OffsetForm(const DTSurface3D &);
extern DTSurface3D FlipOrientation(const DTSurface3D &);
extern DTSurface3D operator*(const DTTransform3D &,const DTSurface3D &);
extern DTSurface3D operator+(const DTSurface3D &,const DTSurface3D &);
extern DTSurface3D AddNextTriangleInformation(const DTSurface3D &);
extern DTSurface3D AddNextTriangleInformation(const DTSurface3D &tri,const DTIntArray &nextTri);
extern DTSurface3D Combine(const DTList<DTSurface3D> &);

extern DTRegion3D BoundingBox(const DTSurface3D &);

extern DTSurface3D ExtractTriangles(const DTSurface3D &,const DTIntArray &); // Same points, sub set of the triangles.  Meant for debugging.

// Create standard surfaces
extern DTSurface3D Sphere(const DTPoint3D &,float r,int divisions);

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTSurface3D &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite);
extern void WriteFloat(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite); // One time value, self documenting.

// A utility class to enable grid sharing.
struct DTSurface3D_SaveInfo {
    DTSurface3D surface;
    std::string name;
};

// Pass in as the third argument to Read and Write to enable the the read/write routines
// to avoid having to read or write identical grids.  This is used when reading DTTrianglarMesh2D and DTTriangularVectorField2D
extern void Read(const DTDataStorage &input,const std::string &name,DTSurface3D &toReturn,DTSurface3D_SaveInfo &);
extern void Write(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite,DTSurface3D_SaveInfo &);
extern void WriteFloat(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite,DTSurface3D_SaveInfo &);

#endif
