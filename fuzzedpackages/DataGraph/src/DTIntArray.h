// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTIntArray_Header
#define DTIntArray_Header

#include <iostream>
#include <unistd.h>
#include "DTLock.h"

// By default, range check is turned on.
#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

// An array of int numbers.  See comments inside DTDoubleArray for more information.
class DTIntArrayStorage {
public:
    DTIntArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov);
    ~DTIntArrayStorage();

    DTLock accessLock;
    ssize_t m,n,o,mn,length;
    int referenceCount;
    int *Data;
    
private:
    DTIntArrayStorage(const DTIntArrayStorage &);
    DTIntArrayStorage &operator=(const DTIntArrayStorage &);
};

class DTMutableIntArray;
class DTIndex;
class DTIntArrayRegion;
struct DTRange;

class DTIntArray {

public:
    DTIntArray() : Storage(new DTIntArrayStorage(0,0,0)), invalidEntry(0) {}
    virtual ~DTIntArray() {Storage->accessLock.Lock(); int refCnt = (--Storage->referenceCount); Storage->accessLock.Unlock(); if (refCnt==0) delete Storage;}
    DTIntArray(const DTIntArray &A) : Storage(A.Storage), invalidEntry(0) {Storage->accessLock.Lock(); Storage->referenceCount++;Storage->accessLock.Unlock(); }
    DTIntArray &operator=(const DTIntArray &A);

protected:
    // If you get a notice that this is protected, change DTIntArray to DTMutableIntArray
    explicit DTIntArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : Storage(new DTIntArrayStorage(mv,nv,ov)), invalidEntry(0) {}

public:
    DTMutableIntArray Copy() const;

    // Size information.
    ssize_t m() const {return Storage->m;}
    ssize_t n() const {return Storage->n;}
    ssize_t o() const {return Storage->o;}
    ssize_t Length() const {return Storage->length;}
    bool IsEmpty() const {return (Storage->length==0);}
    bool NotEmpty() const {return (Storage->length!=0);}
    ssize_t MemoryUsed(void) const {return Length()*sizeof(int);}

    // Low level access
    int ReferenceCount() const {Storage->accessLock.Lock(); int refCnt = Storage->referenceCount; Storage->accessLock.Unlock(); return refCnt;}
    const int *Pointer() const {return Storage->Data;}

    // Allow A(i) and A(i,j), but check each access.
    int operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    int operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    int operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Debug functions, since gdb can't call the () operator.
    int e(int i) const;
    int e(int i,int j) const;
    int e(int i,int j,int k) const;
    void pinfo(void) const;
    void pi(int i) const; // (i,:)
    void pj(int j) const; // (:,j)
    void pall(void) const;  // Uses the same layout as DataTank in the variable monitor.
    void prange(int s,int e) const; // Offsets [s,e], both included.
    void pcrange(int s,int e) const; // Print (:,[s:e]), otherwise like pall().
	void psigns(void) const; // -,0,+

    ssize_t Find(int v) const;
    
    // Support for subregions
    const DTIntArrayRegion operator()(DTIndex) const;
    const DTIntArrayRegion operator()(DTIndex,DTIndex) const;
    const DTIntArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
protected:
    DTIntArrayStorage *Storage;
    int invalidEntry;

    // Error messages for index access.
    void PrintErrorMessage(ssize_t i) const;
    void PrintErrorMessage(ssize_t i,ssize_t j) const;
    void PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const;
};

class DTMutableIntArray : public DTIntArray
{
public:
    DTMutableIntArray() : DTIntArray() {}
    explicit DTMutableIntArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : DTIntArray(mv,nv,ov) {}
    DTMutableIntArray(const DTMutableIntArray &A) : DTIntArray(A) {}

    DTMutableIntArray &operator=(const DTMutableIntArray &A) {DTIntArray::operator=(A); return *this;}

    // Assignment
    DTMutableIntArray &operator=(int a);

    // Raw access
    int *Pointer() {return Storage->Data;}
    const int *Pointer() const {return Storage->Data;}

    // High level access
    int operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    int &operator()(ssize_t i)
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
        return Storage->Data[i];}
    int operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    int &operator()(ssize_t i,ssize_t j)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
        return Storage->Data[i+j*Storage->m];}
    int operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    int &operator()(ssize_t i,ssize_t j,ssize_t k)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
        return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    
    // Support for subregions
    const DTIntArrayRegion operator()(DTIndex) const;
    const DTIntArrayRegion operator()(DTIndex,DTIndex) const;
    const DTIntArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
    DTIntArrayRegion operator()(DTIndex);
    DTIntArrayRegion operator()(DTIndex,DTIndex);
    DTIntArrayRegion operator()(DTIndex,DTIndex,DTIndex);
    
    void operator-=(int v);
    void operator+=(int v);
};

bool operator==(const DTIntArray &A,const DTIntArray &B);
bool operator!=(const DTIntArray &A,const DTIntArray &B);

// Misc
extern DTMutableIntArray Transpose(const DTIntArray &A);
extern DTMutableIntArray Reshape(const DTIntArray &A,ssize_t m,ssize_t n=1,ssize_t o=1);
extern DTMutableIntArray Sort(const DTIntArray &A);
extern DTMutableIntArray FlipJ(const DTIntArray &A);

extern void Swap(DTMutableIntArray &,DTMutableIntArray &);
extern void Swap(DTIntArray &,DTIntArray &);
extern void CopyValues(DTMutableIntArray &into,const DTIntArray &from);
extern void CopyValuesIntoAndAdd(DTMutableIntArray &into,ssize_t offset,const DTIntArray &from,ssize_t add);
extern void CopyIntoColumn(DTMutableIntArray &into,const DTIntArray &list,ssize_t j);  // Copies into into(:,j). Needs list.Length()==into.m().
extern void CopyIntoColumns(DTMutableIntArray &into,const DTRange &intoRange,const DTIntArray &from,const DTRange &fromRange);  // Copies into into(:,j). Needs list.Length()==into.m().

// Memory copies and moves that make sure you are not overwriting bounds
extern void MemoryCopy(DTMutableIntArray &into,ssize_t intoLocation,const DTIntArray &from,ssize_t fromLocation,ssize_t numberOfEntries);
extern void MemoryCopyColumns(DTMutableIntArray &into,ssize_t intoLocation,const DTIntArray &from,ssize_t fromLocation,ssize_t numberOfColumns);
extern void MemoryMove(DTMutableIntArray &into,ssize_t intoLocation,ssize_t fromLocation,ssize_t numberOfEntries);
extern void MemoryMoveColumns(DTMutableIntArray &into,ssize_t intoLocation,ssize_t fromLocation,ssize_t numberOfColumns);

extern DTMutableIntArray ExtractColumns(const DTIntArray &,const DTIntArray &indices);
extern DTMutableIntArray ExtractColumns(const DTIntArray &,const DTRange &);
extern DTMutableIntArray ExtractIndices(const DTIntArray &,const DTIntArray &indices);
extern DTMutableIntArray ExtractIndices(const DTIntArray &,const DTRange &); // A single list, range in offsets

extern int Minimum(const DTIntArray &);
extern int Maximum(const DTIntArray &);
extern void Range(const DTIntArray &,int &minV,int &maxV);

extern ssize_t FindEntry(const DTIntArray &,int); // Linear search for first entry (offset) -1 if not found
extern ssize_t FindEntryInSorted(const DTIntArray &,int); // Same as above, just expects list to be increasing.

extern DTMutableIntArray CombineRows(const DTIntArray &,const DTIntArray &);
extern DTMutableIntArray CombineColumns(const DTIntArray &,const DTIntArray &);

// Changing the size of an array
extern DTMutableIntArray TruncateSize(const DTIntArray &A,ssize_t length);
extern DTMutableIntArray IncreaseSize(const DTIntArray &A,ssize_t addLength);
extern DTMutableIntArray IncreaseSize(const DTIntArray &A); // Doubles the size.

#endif
