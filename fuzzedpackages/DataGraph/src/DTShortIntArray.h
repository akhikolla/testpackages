// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTShortIntArray_Header
#define DTShortIntArray_Header

#include <iostream>
#include <unistd.h>
#include "DTLock.h"

// By default, range check is turned on.
#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

// An array of short int numbers.  See comments inside DTDoubleArray for more information.
class DTShortIntArrayStorage {
public:
    DTShortIntArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov);
    ~DTShortIntArrayStorage();

    DTLock accessLock;
    ssize_t m,n,o,mn,length;
    int referenceCount;
    short int *Data;
    
private:
    DTShortIntArrayStorage(const DTShortIntArrayStorage &);
    DTShortIntArrayStorage &operator=(const DTShortIntArrayStorage &);
};

class DTMutableShortIntArray;
class DTIndex;
class DTShortIntArrayRegion;
struct DTRange;

class DTShortIntArray {

public:
    DTShortIntArray() : Storage(new DTShortIntArrayStorage(0,0,0)), invalidEntry(0) {}
    virtual ~DTShortIntArray() {Storage->accessLock.Lock(); int refCnt = (--Storage->referenceCount); Storage->accessLock.Unlock(); if (refCnt==0) delete Storage;}
    DTShortIntArray(const DTShortIntArray &A) : Storage(A.Storage), invalidEntry(0) {Storage->accessLock.Lock(); Storage->referenceCount++;Storage->accessLock.Unlock(); }
    DTShortIntArray &operator=(const DTShortIntArray &A);

protected:
    // If you get a notice that this is protected, change DTShortIntArray to DTMutableShortIntArray
    explicit DTShortIntArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : Storage(new DTShortIntArrayStorage(mv,nv,ov)), invalidEntry(0) {}

public:
    DTMutableShortIntArray Copy() const;

    // Size information.
    ssize_t m() const {return Storage->m;}
    ssize_t n() const {return Storage->n;}
    ssize_t o() const {return Storage->o;}
    ssize_t Length() const {return Storage->length;}
    bool IsEmpty() const {return (Storage->length==0);}
    bool NotEmpty() const {return (Storage->length!=0);}
    ssize_t MemoryUsed(void) const {return Length()*sizeof(short int);}

    // Low level access
    int ReferenceCount() const {Storage->accessLock.Lock(); int refCnt = Storage->referenceCount; Storage->accessLock.Unlock(); return refCnt;}
    const short int *Pointer() const {return Storage->Data;}

    // Allow A(i) and A(i,j), but check each access.
    short int operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    short int operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    short int operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Debug functions, since gdb can't call the () operator.
    short int e(int i) const;
    short int e(int i,int j) const;
    short int e(int i,int j,int k) const;
    void pinfo(void) const;
    void pi(int i) const; // (i,:)
    void pj(int j) const; // (:,j)
    void pall(void) const;  // Uses the same layout as DataTank in the variable monitor.
    
    // Support for subregions
    const DTShortIntArrayRegion operator()(DTIndex) const;
    const DTShortIntArrayRegion operator()(DTIndex,DTIndex) const;
    const DTShortIntArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;

protected:
    DTShortIntArrayStorage *Storage;
    short int invalidEntry;

    // Error messages for index access.
    void PrintErrorMessage(ssize_t i) const;
    void PrintErrorMessage(ssize_t i,ssize_t j) const;
    void PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const;
};

class DTMutableShortIntArray : public DTShortIntArray
{
public:
    DTMutableShortIntArray() : DTShortIntArray() {}
    explicit DTMutableShortIntArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : DTShortIntArray(mv,nv,ov) {}
    DTMutableShortIntArray(const DTMutableShortIntArray &A) : DTShortIntArray(A) {}

    DTMutableShortIntArray &operator=(const DTMutableShortIntArray &A) {DTShortIntArray::operator=(A); return *this;}

    // Assignment
    DTMutableShortIntArray &operator=(short int a);

    // Raw access
    short int *Pointer() {return Storage->Data;}
    const short int *Pointer() const {return Storage->Data;}

    // High level access
    short int operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    short int &operator()(ssize_t i)
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
        return Storage->Data[i];}
    short int operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    short int &operator()(ssize_t i,ssize_t j)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
        return Storage->Data[i+j*Storage->m];}
    short int operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    short int &operator()(ssize_t i,ssize_t j,ssize_t k)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
        return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    
    // Support for subregions
    const DTShortIntArrayRegion operator()(DTIndex) const;
    const DTShortIntArrayRegion operator()(DTIndex,DTIndex) const;
    const DTShortIntArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
    DTShortIntArrayRegion operator()(DTIndex);
    DTShortIntArrayRegion operator()(DTIndex,DTIndex);
    DTShortIntArrayRegion operator()(DTIndex,DTIndex,DTIndex);
};

// Misc
extern DTMutableShortIntArray Transpose(const DTShortIntArray &A);
extern bool operator==(const DTShortIntArray &A,const DTShortIntArray &B);
extern bool operator!=(const DTShortIntArray &A,const DTShortIntArray &B);

extern void MemoryCopy(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from);
extern void MemoryCopy(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from,const DTRange &range);
extern void MemoryCopyColumns(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from,const DTRange &range);
extern void MemoryMove(DTMutableShortIntArray &into,ssize_t intoLocation,const DTRange &range);
extern void MemoryMoveColumns(DTMutableShortIntArray &into,ssize_t intoLocation,const DTRange &range);

extern void Swap(DTMutableShortIntArray &,DTMutableShortIntArray &);
extern void Swap(DTShortIntArray &,DTShortIntArray &);
extern void CopyValues(DTMutableShortIntArray &into,const DTShortIntArray &from);
extern DTMutableShortIntArray ExtractIndices(const DTShortIntArray &,const DTRange &); // A single list, range in offsets

// Changing the size of an array
extern DTMutableShortIntArray TruncateSize(const DTShortIntArray &A,ssize_t length);
extern DTMutableShortIntArray IncreaseSize(const DTShortIntArray &A,ssize_t addLength);
extern DTMutableShortIntArray FlipJ(const DTShortIntArray &A);

extern short int Minimum(const DTShortIntArray &);
extern short int Maximum(const DTShortIntArray &);

#endif
