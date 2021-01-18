// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTSurface3D.h"

#include "DTError.h"
#include "DTRegion1D.h"
#include "DTRegion3D.h"
#include "DTList.h"
#include "DTTransform3D.h"

#include <algorithm>
using std::sort;

DTSurface3D::DTSurface3D(const DTIntArray &conn,const DTFloatArray &pts)
{
    // Check just the dimensions.
    // If one is empty, everything has to be empty.
    if (pts.IsEmpty() || conn.IsEmpty()) {
        if (!pts.IsEmpty() || !conn.IsEmpty()) {
            DTErrorMessage("DTSurface3D(conn,pts)","Invalid array sizes (one is empty).");
        }
        return;
    }

    if (pts.m()!=3 || conn.m()!=3 || conn.o()>1 || pts.o()>1) {
        DTErrorMessage("DTSurface3D(conn,pts)","Invalid array sizes.");
        return;
    }

    DTRegion1D offRange = ValueRange(conn);
    if (offRange.minV<0 || offRange.maxV>=pts.n()) {
        DTErrorMessage("DTSurface3D(conn,pts)",
                       "Offset array refers to points out of range.");
        return;
    }

    points = DTPointCollection3D(pts);
    connections = conn;
}

DTSurface3D::DTSurface3D(const DTIntArray &conn,const DTPointCollection3D &pts)
{
    // Check just the dimensions.
    // If one is empty, everything has to be empty.
    if (pts.IsEmpty() || conn.IsEmpty()) {
        if (!pts.IsEmpty() || !conn.IsEmpty()) {
            DTErrorMessage("DTSurface3D(conn,pts)","Invalid array sizes (one is empty).");
        }
        return;
    }

    if (conn.m()!=3 || conn.o()>1) {
        DTErrorMessage("DTSurface3D(conn,pts)","Invalid array sizes.");
        return;
    }

    DTRegion1D offRange = ValueRange(conn);
    if (offRange.minV<0 || offRange.maxV>=pts.NumberOfPoints()) {
        DTErrorMessage("DTSurface3D(conn,pts)",
                       "Offset array refers to points out of range.");
        return;
    }

    points = pts;
    connections = conn;
}

DTSurface3D::DTSurface3D(const DTIntArray &conn,const DTPointCollection3D &pts,const DTFloatArray &nrm)
{
    // Check just the dimensions.
    // If one is empty, everything has to be empty.
    if (pts.IsEmpty() || conn.IsEmpty() || nrm.IsEmpty()) {
        if (!pts.IsEmpty() || !conn.IsEmpty() || !nrm.IsEmpty()) {
            DTErrorMessage("DTSurface3D(conn,pts,nrm)","points, connections or normals are empty.");
        }
        return;
    }
    
    if (conn.m()!=3 || pts.NumberOfPoints()!=nrm.n() || nrm.m()*nrm.n()!=nrm.Length() || conn.o()>1) {
        DTErrorMessage("DTSurface3D(conn,pts,nrm)","Invalid array sizes.");
        return;
    }

    DTRegion1D offRange = ValueRange(conn);
    if (offRange.minV<0 || offRange.maxV>=pts.NumberOfPoints()) {
        DTErrorMessage("DTSurface3D(conn,pts,nrm)",
                       "Offset array refers to points out of range.");
        return;
    }

    points = DTPointCollection3D(pts);
    connections = conn;
    normals = nrm;
}

DTSurface3D::DTSurface3D(const DTFloatArray &tris)
{
    // Check if valid.
    if (tris.IsEmpty()) return;

    if (tris.m()!=9 || tris.o()!=1) {
        DTErrorMessage("DTSurface3D(triangles)","Invalid array size.  Needs to be a 9xN array.");
        return;
    }

    triangles = tris;
}

DTSurface3D::DTSurface3D(const DTIntArray &conn,const DTPointCollection3D &pts,const DTFloatArray &nrm,const DTIntArray &next)
{
    connections = conn;
    points = pts;
    normals = nrm;
    nextTriangles = next;
}

DTSurface3D ConvertToFloat(const DTSurface3D &s)
{
    if (s.Points().FloatPrecision())
        return s;
    else {
        if (s.NormalsDefined()) {
            return DTSurface3D(s.Connections(),ConvertToFloat(s.Points()),s.Normals());
        }
        else {
            return DTSurface3D(s.Connections(),ConvertToFloat(s.Points()));
        }
    }
}

DTSurface3D ConvertToDouble(const DTSurface3D &s)
{
    if (s.Points().DoublePrecision())
        return s;
    else {
        if (s.NormalsDefined()) {
            return DTSurface3D(s.Connections(),ConvertToDouble(s.Points()),s.Normals());
        }
        else {
            return DTSurface3D(s.Connections(),ConvertToDouble(s.Points()));
        }
    }
}

DTPointCollection3D DTSurface3D::Points(void) const
{
    return points;
}

DTIntArray DTSurface3D::Connections(void) const
{
    return connections;
}

DTFloatArray DTSurface3D::Normals(void) const
{
    if (normals.IsEmpty() && points.NumberOfPoints()!=0) {
        DTErrorMessage("DTSurface3D::Normals","No normals saved for surface.");
    }
    return normals;
}

DTIntArray DTSurface3D::NextTriangles() const
{
    if (nextTriangles.IsEmpty() && !connections.IsEmpty()) {
        DTErrorMessage("DTSurface3D::NextTriangles","Not specified");
    }
    return nextTriangles;
}

void DTSurface3D::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    if (points.IsEmpty()) {
        if (triangles.IsEmpty()) {
            std::cerr << "No triangles.";
        }
        else {
            std::cerr << "Only storing triangles - " << triangles.n() << " triangles.";
        }
    }
    else {
        std::cerr << connections.n() << " triangles, " << points.NumberOfPoints() << " points.";
        if (normals.NotEmpty())
            std::cerr << " Includes normals.";
    }
    
    if (points.DoublePrecision())
        std::cerr << " (double)";
    else if (points.FloatPrecision())
        std::cerr << " (float)";
    
    std::cerr << std::endl << std::flush;
#endif
}

void DTSurface3D::pt(int i) const
{
#ifndef DG_NOSTDErrOut
    if (i<0 || i>=connections.n())
        std::cerr << "Triangle offset out of bounds.\n";
    else {
        DTFloatArray pts = points.Data();
        std::cerr << "offsets = (" << connections(0,i) << ", " << connections(1,i) << ", " << connections(2,i) << ")\n";
        std::cerr << "    " << connections(0,i) << " = (" 
            << pts(0,connections(0,i)) << ", "
            << pts(1,connections(0,i)) << ", "
            << pts(2,connections(0,i)) << ")\n";
        std::cerr << "    " << connections(1,i) << " = (" 
            << pts(0,connections(1,i)) << ", "
            << pts(1,connections(1,i)) << ", "
            << pts(2,connections(1,i)) << ")\n";
        std::cerr << "    " << connections(2,i) << " = ("
            << pts(0,connections(2,i)) << ", "
            << pts(1,connections(2,i)) << ", "
            << pts(2,connections(2,i)) << ")\n";
    }
#endif
}

struct DTSurface3DPointAndIndex {
    // Dictionary sort.  That is first compare the x coordinate, and if they are the same, compare the y.
    int operator<(const DTSurface3DPointAndIndex &C) const {return (z<C.z || (z==C.z && y<C.y) || (z==C.z && y==C.y && x<C.x));}
    int operator==(const DTSurface3DPointAndIndex &C) const {return (x==C.x && y==C.y && z==C.z);}
    
    float x,y,z;
    int index;
};

DTSurface3D OffsetForm(const DTSurface3D &s)
{
    if (s.SavedAsOffsets())
        return s;
    
    DTFloatArray Tri = s.Triangles();
    
    // Remove all triangles that contain Nan or Inf.
    const float *TriD = Tri.Pointer();
    const ssize_t howManyEntries = Tri.Length();
    int i;
    for (i=0;i<howManyEntries;i++) {
        if (isfinite(TriD[i])==0) break;
    }
    
    if (i<howManyEntries) {
        // Remove the Nan triangles.  Begin by counting how many triangles
        // will be in the new triangle array.
        int howManyTri = 0;
        for (i=0;i<howManyEntries;i+=9) {
            if (isfinite(TriD[i  ]) && isfinite(TriD[i+1]) && isfinite(TriD[i+2]) &&
                isfinite(TriD[i+3]) && isfinite(TriD[i+4]) && isfinite(TriD[i+5]) &&
                isfinite(TriD[i+6]) && isfinite(TriD[i+7]) && isfinite(TriD[i+8]))
                howManyTri++;
        }
        DTMutableFloatArray NewTri(9,howManyTri);
        float *NewTriD = NewTri.Pointer();
        int pos = 0;
        for (i=0;i<howManyEntries;i+=9) {
            if (isfinite(TriD[i  ]) && isfinite(TriD[i+1]) && isfinite(TriD[i+2]) &&
                isfinite(TriD[i+3]) && isfinite(TriD[i+4]) && isfinite(TriD[i+5]) &&
                isfinite(TriD[i+6]) && isfinite(TriD[i+7]) && isfinite(TriD[i+8])) {
                NewTriD[pos  ] = TriD[i  ];
                NewTriD[pos+1] = TriD[i+1];
                NewTriD[pos+2] = TriD[i+2];
                NewTriD[pos+3] = TriD[i+3];
                NewTriD[pos+4] = TriD[i+4];
                NewTriD[pos+5] = TriD[i+5];
                NewTriD[pos+6] = TriD[i+6];
                NewTriD[pos+7] = TriD[i+7];
                NewTriD[pos+8] = TriD[i+8];
                pos+=9;
            }
        }
        Tri = NewTri;
    }
    
    if (Tri.IsEmpty())
        return DTSurface3D();
    
    ssize_t NumberOfTriangles = Tri.n();
    ssize_t NumberOfPoints = NumberOfTriangles*3;
    const float *TriData = Tri.Pointer();
    
    DTMutableList<DTSurface3DPointAndIndex> points(NumberOfPoints);
    DTSurface3DPointAndIndex *pointsD = points.Pointer();
    
    DTSurface3DPointAndIndex p;
    int i3 = 0;
    
    for (i=0;i<NumberOfPoints;i++) {
        p.x = TriData[i3  ];
        p.y = TriData[i3+1];
        p.z = TriData[i3+2];
        p.index = i;
        pointsD[i] = p;
        i3+=3;
    }
    
    // Now sort the points.
    sort(pointsD,pointsD+NumberOfPoints);
    
    // Now count the number of different points.
    int CompareTo,Compare;
    CompareTo = 0;
    Compare = 1;
    int NumberOfDifferentPoints = 1;
    while (Compare<NumberOfPoints) {
        while (Compare<NumberOfPoints && pointsD[CompareTo]==pointsD[Compare]) Compare++;
        if (Compare<NumberOfPoints) {
            // Found a new one.
            CompareTo = Compare;
            Compare++;
            NumberOfDifferentPoints++;
        }
    }
    
    // Now create the point array, and triangle array and fill in both.
    DTMutableFloatArray thePoints(3,NumberOfDifferentPoints);
    DTMutableIntArray triangles(3,NumberOfTriangles);
    triangles = -1;
    
    CompareTo = 0;
    Compare = 1;
    int posInPoints = 0;
    p = points(0);
    thePoints(0,posInPoints) = p.x;
    thePoints(1,posInPoints) = p.y;
    thePoints(2,posInPoints) = p.z;
    triangles(p.index) = posInPoints;
    while (Compare<NumberOfPoints) {
        while (Compare<NumberOfPoints && pointsD[CompareTo]==pointsD[Compare]) {
            p = pointsD[Compare];
            triangles(p.index) = posInPoints;
            Compare++;
        }
        if (Compare<NumberOfPoints) {
            // Found a new one.
            posInPoints++;
            CompareTo = Compare;
            p = pointsD[Compare];
            thePoints(0,posInPoints) = p.x;
            thePoints(1,posInPoints) = p.y;
            thePoints(2,posInPoints) = p.z;
            triangles(p.index) = posInPoints;
            Compare++;
        }
    }
    
    return DTSurface3D(triangles,thePoints);    
}

bool operator==(const DTSurface3D &A,const DTSurface3D &B)
{
    // Do not check optional data.
    if (A.Points()!=B.Points()) return false;
    if (A.Connections()!=B.Connections()) return false;
    if (A.Triangles()!=B.Triangles()) return false;
    if (A.NormalsDefined() && B.NormalsDefined()) {
        return A.Normals()==B.Normals();
    }
    if (A.NormalsDefined()==false && B.NormalsDefined()==false) {
        return true;
    }
    else {
        return false;
    }
}

bool operator!=(const DTSurface3D &A,const DTSurface3D &B)
{
    if (A.Points()!=B.Points()) return true;
    if (A.Connections()!=B.Connections()) return true;
    if (A.Triangles()!=B.Triangles()) return true;
    if (A.NormalsDefined() && B.NormalsDefined()) {
        return A.Normals()!=B.Normals();
    }
    if (A.NormalsDefined()==false && B.NormalsDefined()==false) {
        return false;
    }
    else {
        return true;
    }
}

DTSurface3D FlipOrientation(const DTSurface3D &S)
{
    DTIntArray connections = S.Connections();
    ssize_t numberOfTriangles = connections.n();
    DTMutableIntArray newConnections(3,numberOfTriangles);
    int triN;
    for (triN=0;triN<numberOfTriangles;triN++) {
        newConnections(0,triN) = connections(0,triN);
        newConnections(1,triN) = connections(2,triN);
        newConnections(2,triN) = connections(1,triN);
    }
    
    DTSurface3D toReturn;
    if (S.NormalsDefined()) {
        DTMutableFloatArray normals = S.Normals().Copy();
        ssize_t howMany = normals.Length();
        for (int i=0;i<howMany;i++) normals(i)*=-1.0;
        toReturn = DTSurface3D(newConnections,S.Points(),normals);
    }
    else {
        toReturn = DTSurface3D(newConnections,S.Points());
    }
    
    return toReturn;
}

DTSurface3D AddNextTriangleInformation(const DTSurface3D &tri,const DTIntArray &nextTri)
{
    if (nextTri.NotEmpty()) {
        if (nextTri.o()>1 || nextTri.m()!=3 || nextTri.n()!=tri.NumberOfTriangles()) {
            DTErrorMessage("AddNextTriangleInformation(DTSurface3D,DTIntArray)","Invalid array size.");
            return DTSurface3D();
        }
        DTRegion1D valRegion = ValueRange(nextTri);
        if (valRegion.minV<-1 || valRegion.maxV>=tri.NumberOfTriangles()*3) {
            DTErrorMessage("AddNextTriangleInformation(DTSurface3D,DTIntArray)","Invalid offset");
            return DTSurface3D();
        }
    }
    
    if (tri.NormalsDefined())  {
        return DTSurface3D(tri.Connections(),tri.Points(),tri.Normals(),nextTri);
    }
    else {
        return DTSurface3D(tri.Connections(),tri.Points(),DTFloatArray(),nextTri);
    }
}

struct DTTriangleSideLongSurface {
    
    bool operator<(const DTTriangleSideLongSurface &A) const {return (label<A.label);}
    
    size_t label;
    int positionInOffset;
};

DTSurface3D AddNextTriangleInformation(const DTSurface3D &tri)
{
    if (tri.NextTrianglesDefined()) return tri;
    
    DTIntArray toReturn;
    DTIntArray offsets = tri.Connections();
    ssize_t Len = offsets.n();
    int i;
    int Ai,Bi,Ci,ind1,ind2,temp;
    
    // Can save the offset as a single unsigned long long.
    DTMutableList<DTTriangleSideLongSurface> offsetInfo(3*Len);
    DTTriangleSideLongSurface side;
    
    for (i=0;i<Len;i++) {
        Ai = offsets(0,i);
        Bi = offsets(1,i);
        Ci = offsets(2,i);
        
        // Numbering is AB,BC,CA
        
        // The AB side.
        ind1 = Ai;
        ind2 = Bi;
        // Sort the indices with bubble sort.
        if (ind2>ind1) {temp = ind1; ind1 = ind2; ind2 = temp;}
        side.label = ind2;
        side.label *= Len;
        side.label += ind1;
        side.positionInOffset = 3*i;
        offsetInfo(3*i) = side;
        
        // BC
        ind1 = Bi;
        ind2 = Ci;
        if (ind2>ind1) {temp = ind1; ind1 = ind2; ind2 = temp;}
        side.label = ind2;
        side.label *= Len;
        side.label += ind1;
        side.positionInOffset = 3*i+1;
        offsetInfo(3*i+1) = side;
        
        // CA
        ind1 = Ci;
        ind2 = Ai;
        if (ind2>ind1) {temp = ind1; ind1 = ind2; ind2 = temp;}
        side.label = ind2;
        side.label *= Len;
        side.label += ind1;
        side.positionInOffset = 3*i+2;
        offsetInfo(3*i+2) = side;
    }
    
    // Sort the list
    sort(offsetInfo.Pointer(),offsetInfo.Pointer()+offsetInfo.Length());
    
    // Each side should be repeated.  Go through the sorted list and connect the pairs.
    
    DTMutableIntArray NextTri(3,Len);
    NextTri = -1;
    
    if (Len) {
        int pos = 0;
        int pos1,pos2;
        
        while (pos<3*Len-1) {
            if (offsetInfo(pos).label==offsetInfo(pos+1).label) {
                pos1 = offsetInfo(pos).positionInOffset;
                pos2 = offsetInfo(pos+1).positionInOffset;
                NextTri(pos1%3,pos1/3) = pos2;
                NextTri(pos2%3,pos2/3) = pos1;
                pos+=2;
            }
            else {
                // Not connected
                pos++;
            }
        }
    }
    
    return AddNextTriangleInformation(tri,NextTri);
}

DTSurface3D operator*(const DTTransform3D &T,const DTSurface3D &s)
{
    DTPointCollection3D rotatedPoints = T*s.Points();
    if (s.NormalsDefined()) {
        DTTransform3D mapNormals = T.RotationOnly().Inverse();
        DTMutableFloatArray rotatedNormals = TransformPoints(mapNormals,s.Normals());
        if (mapNormals.IsOrthogonal()==false) {
            // Need to normalize the vectors
            ssize_t i,howMany = rotatedNormals.n();
            float nx,ny,nz,length;
            for (i=0;i<howMany;i++) {
                nx = rotatedNormals(0,i);
                ny = rotatedNormals(1,i);
                nz = rotatedNormals(2,i);
                length = sqrtf(nx*nx+ny*ny+nz*nz);
                rotatedNormals(0,i) /= length;
                rotatedNormals(1,i) /= length;
                rotatedNormals(2,i) /= length;
            }
        }
        return DTSurface3D(s.Connections(),rotatedPoints,rotatedNormals);
    }
    else {
        return DTSurface3D(s.Connections(),rotatedPoints);
    }
}

DTSurface3D operator+(const DTSurface3D &A,const DTSurface3D &B)
{
	DTMutableList<DTSurface3D> L(2);
	L(0) = A;
	L(1) = B;
	return Combine(L);
}

DTSurface3D Combine(const DTList<DTSurface3D> &L)
{
	int surfN;
    ssize_t howManySurfaces = L.Length();
	int howManyPoints = 0;
	int howManyTriangles = 0;
	DTSurface3D S;
	bool hasNormals = true;
	for (surfN=0;surfN<howManySurfaces;surfN++) {
		S = L(surfN);
		howManyPoints += S.NumberOfPoints();
		howManyTriangles += S.NumberOfTriangles();
		if (S.NormalsDefined()==false) hasNormals = false;
	}
	
	DTMutableFloatArray newPoints(3,howManyPoints);
	DTMutableFloatArray newNormals(3,hasNormals ? howManyPoints : 0);
	DTMutableIntArray newTriangles(3,howManyTriangles);
	int posInPoints = 0;
	int posInTriangles = 0;
	for (surfN=0;surfN<howManySurfaces;surfN++) {
		S = L(surfN);
		MemoryCopy(newPoints,3*posInPoints,S.Points().FloatData());
		CopyValuesIntoAndAdd(newTriangles,3*posInTriangles,S.Connections(),posInPoints);
		if (hasNormals)
			MemoryCopy(newNormals,3*posInPoints,S.Normals());
		posInPoints += S.NumberOfPoints();
		posInTriangles += S.NumberOfTriangles();
	}
	
	if (hasNormals)
		return DTSurface3D(newTriangles,DTPointCollection3D(newPoints),newNormals);
	else
		return DTSurface3D(newTriangles,DTPointCollection3D(newPoints));
}

DTSurface3D ExtractTriangles(const DTSurface3D &surf,const DTIntArray &which)
{
	// Same points, sub set of the triangles.
	ssize_t i,howMany = which.Length();
    ssize_t index;
	DTIntArray connections = surf.Connections();
	ssize_t howManyConnections = connections.n();
	DTMutableIntArray newConnections(3,howMany);
	int pos = 0;
	for (i=0;i<howMany;i++) {
		index = which(i);
		if (index<0 || index>=howManyConnections) {
			DTErrorMessage("ExtractTriangles(DTSurface3D,DTIntArray","Index out of bounds");
		}
		else {
			newConnections(0,pos) = connections(0,index);
			newConnections(1,pos) = connections(1,index);
			newConnections(2,pos) = connections(2,index);
			pos++;
		}
	}
	if (pos!=howMany) {
		newConnections = TruncateSize(newConnections,3*pos);
	}
	if (surf.NormalsDefined()) {
		return DTSurface3D(newConnections,surf.Points(),surf.Normals());
	}
	else {
		return DTSurface3D(newConnections,surf.Points());
	}
}

DTRegion3D BoundingBox(const DTSurface3D &coll)
{
    if (coll.Triangles().IsEmpty()) {
        return BoundingBox(coll.Points());
    }
    else {
        return BoundingBox3D(coll.Triangles());
    }
}

void Read(const DTDataStorage &input,const std::string &name,DTSurface3D &toReturn)
{
    std::string theName = input.ResolveName(name);

    if (input.Contains(theName+"_P")) {
        DTIntArray conn = input.ReadIntArray(theName);
        DTPointCollection3D points;
        Read(input,name+"_P",points);

        if (input.Contains(theName+"_N")) {
            DTFloatArray nrms = input.ReadFloatArray(theName+"_N");
            toReturn = DTSurface3D(conn,points,nrms);
        }
        else {
            toReturn = DTSurface3D(conn,points);
        }
    }
    else {
        DTFloatArray arr = input.ReadFloatArray(theName);
        toReturn = DTSurface3D(arr);
    }
}

DTSurface3D Sphere(const DTPoint3D &P,float r,int divisions)
{
	int i,j;
	float x = float(P.x);
	float y = float(P.y);
	float z = float(P.z);
	
	int divisionInXYPlane = divisions;
	int divisionInZ = divisions;
	
	int posInPoints = 0;
	DTMutableFloatArray c(divisionInXYPlane);
	DTMutableFloatArray s(divisionInXYPlane);
	double dtheta = 2*M_PI/divisionInXYPlane;
	for (i=0;i<divisionInXYPlane;i++) {
		c(i) = cos(i*dtheta);
		s(i) = sin(i*dtheta);
	}

	// Compute points and normals
	DTMutableFloatArray points(3,2+(divisionInZ-2)*divisionInXYPlane);
	DTMutableFloatArray normals(3,2+(divisionInZ-2)*divisionInXYPlane);
	
	// Bottom point
	points(0,posInPoints) = x;
	points(1,posInPoints) = y;
	points(2,posInPoints) = z-r;
	normals(0,posInPoints) = 0;
	normals(1,posInPoints) = 0;
	normals(2,posInPoints) = -1;
	posInPoints++;
	// Center chunk
	float rc,rs,cs,sn;
	float dphi = M_PI/(divisionInZ-1.0);
	for (j=1;j<divisionInZ-1;j++) {
		cs = cos(j*dphi-M_PI_2);
		rc = r*cs;
		sn = sin(j*dphi-M_PI_2);
		rs = r*sn;
		for (i=0;i<divisionInXYPlane;i++) {
			points(0,posInPoints) = x+c(i)*rc;
			points(1,posInPoints) = y+s(i)*rc;
			points(2,posInPoints) = z+rs;
			normals(0,posInPoints) = c(i)*cs;
			normals(1,posInPoints) = s(i)*cs;
			normals(2,posInPoints) = sn;
			posInPoints++;
		}
	}
	// Top point
	points(0,posInPoints) = x;
	points(1,posInPoints) = y;
	points(2,posInPoints) = z+r;
	normals(0,posInPoints) = 0;
	normals(1,posInPoints) = 0;
	normals(2,posInPoints) = 1;
	posInPoints++;
	
	// Now triangle offsets
	DTMutableIntArray connections(3,2*divisionInXYPlane*(divisionInZ-2));

	// Triangles at bottom
	int posInTriangles = 0;
	for (i=0;i<divisionInXYPlane-1;i++) {
		connections(0,posInTriangles) = 0;
		connections(1,posInTriangles) = i+2;
		connections(2,posInTriangles) = i+1;
		posInTriangles++;
	}
	connections(0,posInTriangles) = 0;
	connections(1,posInTriangles) = 1;
	connections(2,posInTriangles) = divisionInXYPlane;
	posInTriangles++;
	
	// Bulk of the the triangles
	int offsetAtStartOfCircle = 1;
	for (j=1;j<divisionInZ-2;j++) {
		for (i=0;i<divisionInXYPlane-1;i++) {
			connections(0,posInTriangles) = i+offsetAtStartOfCircle;
			connections(1,posInTriangles) = i+1+offsetAtStartOfCircle;
			connections(2,posInTriangles) = i+1+offsetAtStartOfCircle+divisionInXYPlane;
			posInTriangles++;
			connections(0,posInTriangles) = i+1+offsetAtStartOfCircle+divisionInXYPlane;
			connections(1,posInTriangles) = i+offsetAtStartOfCircle+divisionInXYPlane;
			connections(2,posInTriangles) = i+offsetAtStartOfCircle;
			posInTriangles++;
		}
		connections(0,posInTriangles) = divisionInXYPlane-1+offsetAtStartOfCircle;
		connections(1,posInTriangles) = offsetAtStartOfCircle;
		connections(2,posInTriangles) = offsetAtStartOfCircle+divisionInXYPlane;
		posInTriangles++;
		connections(0,posInTriangles) = offsetAtStartOfCircle+divisionInXYPlane;
		connections(1,posInTriangles) = divisionInXYPlane-1+offsetAtStartOfCircle+divisionInXYPlane;
		connections(2,posInTriangles) = divisionInXYPlane-1+offsetAtStartOfCircle;
		posInTriangles++;
		
		offsetAtStartOfCircle+=divisionInXYPlane;
	}
	
	// Top
	int endOffset = offsetAtStartOfCircle+divisionInXYPlane;
	for (i=0;i<divisionInXYPlane-1;i++) {
		connections(0,posInTriangles) = i+offsetAtStartOfCircle;
		connections(1,posInTriangles) = i+1+offsetAtStartOfCircle;
		connections(2,posInTriangles) = endOffset;
		posInTriangles++;
	}
	connections(0,posInTriangles) = offsetAtStartOfCircle+(divisionInXYPlane-1);
	connections(1,posInTriangles) = offsetAtStartOfCircle;
	connections(2,posInTriangles) = endOffset;
	posInTriangles++;
	
	return DTSurface3D(connections,DTPointCollection3D(points),normals);
}

void Write(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite)
{
    Write(output,name+"_bbox3D",BoundingBox(toWrite));

    if (toWrite.Triangles().IsEmpty()) {
        Write(output,name+"_P",toWrite.Points());
        if (toWrite.NormalsDefined()) {
            Write(output,name+"_N",toWrite.Normals());
        }
        Write(output,name,toWrite.Connections());
    }
    else {
        Write(output,name,toWrite.Triangles());
    }
}

void WriteFloat(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite)
{
    Write(output,name+"_bbox3D",BoundingBox(toWrite));

    if (toWrite.Triangles().IsEmpty()) {
        WriteFloat(output,name+"_P",toWrite.Points());
        if (toWrite.NormalsDefined()) {
            Write(output,name+"_N",toWrite.Normals());
        }
        Write(output,name,toWrite.Connections());
    }
    else {
        Write(output,name,toWrite.Triangles());
    }
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"3D Surface");
    output.Flush();
}

void Read(const DTDataStorage &input,const std::string &name,DTSurface3D &toReturn,DTSurface3D_SaveInfo &SaveInfo)
{
    std::string theName = input.ResolveName(name);
    if (SaveInfo.name==theName) {
        toReturn = SaveInfo.surface; // Already read this grid in.
        return;
    }
    
    Read(input,theName,toReturn);
    
    SaveInfo.surface = toReturn;
    SaveInfo.name = theName;
}

void Write(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite,DTSurface3D_SaveInfo &SaveInfo)
{
    if (SaveInfo.name.length() && SaveInfo.surface==toWrite) {
        // Just save the reference.
        Write(output,name,SaveInfo.name);
    }
    else {
        Write(output,name,toWrite);
        SaveInfo.surface = toWrite;
        SaveInfo.name = name;
    }
}

void WriteFloat(DTDataStorage &output,const std::string &name,const DTSurface3D &toWrite,DTSurface3D_SaveInfo &SaveInfo)
{
    if (SaveInfo.name.length() && SaveInfo.surface==toWrite) {
        // Just save the reference.
        Write(output,name,SaveInfo.name);
    }
    else {
        WriteFloat(output,name,toWrite);
        SaveInfo.surface = toWrite;
        SaveInfo.name = name;
    }
}

