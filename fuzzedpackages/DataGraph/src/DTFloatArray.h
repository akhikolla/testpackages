// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTFloatArray_Header
#define DTFloatArray_Header

#include <iostream>
#include <unistd.h>
#include "DTLock.h"

// By default, range check is turned on.
#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

// An array of float numbers.  See comments inside DTDoubleArray for more information.
class DTFloatArrayStorage {
public:
    DTFloatArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov);
    ~DTFloatArrayStorage();

    DTLock accessLock;
    ssize_t m,n,o,mn,length;
    int referenceCount;
    float *Data;
    
private:
    DTFloatArrayStorage(const DTFloatArrayStorage &);
    DTFloatArrayStorage &operator=(const DTFloatArrayStorage &);
};

class DTMutableFloatArray;
class DTIndex;
class DTFloatArrayRegion;
class DTIntArray;
struct DTRange;

class DTFloatArray {

public:
    DTFloatArray() : Storage(new DTFloatArrayStorage(0,0,0)), invalidEntry(0.0) {}
    virtual ~DTFloatArray() {Storage->accessLock.Lock(); int refCnt = (--Storage->referenceCount); Storage->accessLock.Unlock(); if (refCnt==0) delete Storage;}
    DTFloatArray(const DTFloatArray &A) : Storage(A.Storage), invalidEntry(0.0) {Storage->accessLock.Lock(); Storage->referenceCount++;Storage->accessLock.Unlock(); }
    DTFloatArray &operator=(const DTFloatArray &A);

protected:
    // If you get a notice that this is protected, change DTFloatArray to DTMutableFloatArray
    explicit DTFloatArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : Storage(new DTFloatArrayStorage(mv,nv,ov)), invalidEntry(0.0) {}

public:
    DTMutableFloatArray Copy() const;

    // Size information.
    ssize_t m() const {return Storage->m;}
    ssize_t n() const {return Storage->n;}
    ssize_t o() const {return Storage->o;}
    ssize_t Length() const {return Storage->length;}
    bool IsEmpty() const {return (Storage->length==0);}
    bool NotEmpty() const {return (Storage->length!=0);}
    ssize_t MemoryUsed(void) const {return Length()*sizeof(float);}

    // Low level access
    int ReferenceCount() const {Storage->accessLock.Lock(); int refCnt = Storage->referenceCount; Storage->accessLock.Unlock(); return refCnt;}
    const float *Pointer() const {return Storage->Data;}

    // Allow A(i) and A(i,j), but check each access.
    float operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    float operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    float operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Debug functions, since gdb can't call the () operator.
    float e(int i) const;
    float e(int i,int j) const;
    float e(int i,int j,int k) const;
    void pinfo(void) const;
    void pi(int i) const; // (i,:)
    void pj(int j) const; // (:,j)
    void pall(void) const;  // Uses the same layout as DataTank in the variable monitor.
    void prange(int s,ssize_t e) const; // Offsets [s,e], both included.
    void pcrange(int s,int e) const; // Print (:,[s:e]), otherwise like pall().
 	void psigns(void) const; // -,0,+

    ssize_t Find(float v) const;
    
    // Support for subregions
    const DTFloatArrayRegion operator()(DTIndex) const;
    const DTFloatArrayRegion operator()(DTIndex,DTIndex) const;
    const DTFloatArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
protected:
    DTFloatArrayStorage *Storage;
    float invalidEntry;

    // Error messages for index access.
    void PrintErrorMessage(ssize_t i) const;
    void PrintErrorMessage(ssize_t i,ssize_t j) const;
    void PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const;
};

class DTMutableFloatArray : public DTFloatArray
{
public:
    DTMutableFloatArray() : DTFloatArray() {}
    explicit DTMutableFloatArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : DTFloatArray(mv,nv,ov) {}
    DTMutableFloatArray(const DTMutableFloatArray &A) : DTFloatArray(A) {}

    DTMutableFloatArray &operator=(const DTMutableFloatArray &A) {DTFloatArray::operator=(A); return *this;}

    // Assignment
    DTMutableFloatArray &operator=(float a);

    // Raw access
    float *Pointer() {return Storage->Data;}
    const float *Pointer() const {return Storage->Data;}

    // High level access
    float operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    float &operator()(ssize_t i)
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
        return Storage->Data[i];}
    float operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    float &operator()(ssize_t i,ssize_t j)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
        return Storage->Data[i+j*Storage->m];}
    float operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    float &operator()(ssize_t i,ssize_t j,ssize_t k)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
        return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Support for subregions
    const DTFloatArrayRegion operator()(DTIndex) const;
    const DTFloatArrayRegion operator()(DTIndex,DTIndex) const;
    const DTFloatArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
    DTFloatArrayRegion operator()(DTIndex);
    DTFloatArrayRegion operator()(DTIndex,DTIndex);
    DTFloatArrayRegion operator()(DTIndex,DTIndex,DTIndex);
	
	void operator*=(float v);
	void operator/=(float v);
	void operator+=(float v);
	void operator-=(float v);
	
	void operator+=(const DTFloatArray &);
	void operator-=(const DTFloatArray &);
	void operator*=(const DTFloatArray &);
	void operator/=(const DTFloatArray &);
};

bool operator==(const DTFloatArray &A,const DTFloatArray &B);
bool operator!=(const DTFloatArray &A,const DTFloatArray &B);

// Misc
extern DTMutableFloatArray Transpose(const DTFloatArray &A);
extern DTMutableFloatArray Reshape(const DTFloatArray &A,ssize_t m,ssize_t n=1,ssize_t o=1);
extern DTMutableFloatArray Sort(const DTFloatArray &A);
extern DTMutableFloatArray FlipJ(const DTFloatArray &A);

extern void Swap(DTMutableFloatArray &,DTMutableFloatArray &);
extern void Swap(DTFloatArray &,DTFloatArray &);
extern void CopyValues(DTMutableFloatArray &into,const DTFloatArray &from);
extern void CopyIntoColumns(DTMutableFloatArray &into,const DTRange &intoRange,const DTFloatArray &from,const DTRange &fromRange);  // Copies into into(:,j). Needs list.Length()==into.m().

// Was called CopyValuesInto, but moved to unify the function name.
extern void MemoryCopy(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from);
extern void MemoryCopy(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from,const DTRange &range);
extern void MemoryCopyColumns(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from,const DTRange &range);
extern void MemoryMove(DTMutableFloatArray &into,ssize_t intoLocation,const DTRange &range);
extern void MemoryMoveColumns(DTMutableFloatArray &into,ssize_t intoLocation,const DTRange &range);

extern DTMutableFloatArray Region(const DTFloatArray &,const DTRange &iRange,const DTRange &jRange,const DTRange &kRange);
extern DTMutableFloatArray Region(const DTFloatArray &,const DTRange &iRange,const DTRange &jRange);
extern DTMutableFloatArray ExtractIndices(const DTFloatArray &,const DTIntArray &indices);
extern DTMutableFloatArray ExtractIndices(const DTFloatArray &A,const DTRange &r);
extern DTMutableFloatArray ExtractColumns(const DTFloatArray &,const DTIntArray &indices);
extern DTMutableFloatArray ExtractColumns(const DTFloatArray &,const DTRange &);


extern float Minimum(const DTFloatArray &);
extern float Maximum(const DTFloatArray &);
extern float Mean(const DTFloatArray &);
extern DTMutableFloatArray Minimum(const DTFloatArray &,const DTFloatArray &);
extern DTMutableFloatArray Maximum(const DTFloatArray &,const DTFloatArray &);

// Changing the size of an array
extern DTMutableFloatArray TruncateSize(const DTFloatArray &A,ssize_t length);
extern DTMutableFloatArray IncreaseSize(const DTFloatArray &A,ssize_t addLength);

extern DTMutableFloatArray CombineColumns(const DTFloatArray &First,const DTFloatArray &Second);
extern DTMutableFloatArray CombineColumns(const DTFloatArray &First,const DTFloatArray &Second,ssize_t fromSecond);

#endif
