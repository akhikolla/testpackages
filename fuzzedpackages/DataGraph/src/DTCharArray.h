// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTCharArray_Header
#define DTCharArray_Header

#include <iostream>
#include <unistd.h>
#include "DTLock.h"

// By default, range check is turned on.
#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

// An array of char numbers.  See comments inside DTDoubleArray for more information.
class DTCharArrayStorage {
public:
    DTCharArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov);
    ~DTCharArrayStorage();

    DTLock accessLock;
    ssize_t m,n,o,mn,length;
    int referenceCount;
    char *Data;
    
private:
    DTCharArrayStorage(const DTCharArrayStorage &);
    DTCharArrayStorage &operator=(const DTCharArrayStorage &);
};

class DTMutableCharArray;
class DTIndex;
class DTCharArrayRegion;
struct DTRange;

class DTCharArray {

public:
    DTCharArray() : Storage(new DTCharArrayStorage(0,0,0)), invalidEntry(0) {}
    virtual ~DTCharArray() {Storage->accessLock.Lock(); int refCnt = (--Storage->referenceCount); Storage->accessLock.Unlock(); if (refCnt==0) delete Storage;}
    DTCharArray(const DTCharArray &A) : Storage(A.Storage), invalidEntry(0) {Storage->accessLock.Lock(); Storage->referenceCount++;Storage->accessLock.Unlock(); }
    DTCharArray &operator=(const DTCharArray &A);

protected:
    // If you get a notice that this is protected, change DTCharArray to DTMutableCharArray
    explicit DTCharArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : Storage(new DTCharArrayStorage(mv,nv,ov)), invalidEntry(0) {}

public:
    DTMutableCharArray Copy() const;

    // Size information.
    ssize_t m() const {return Storage->m;}
    ssize_t n() const {return Storage->n;}
    ssize_t o() const {return Storage->o;}
    ssize_t Length() const {return Storage->length;}
    bool IsEmpty() const {return (Storage->length==0);}
    bool NotEmpty() const {return (Storage->length!=0);}
    ssize_t MemoryUsed(void) const {return Length();}

    // Low level access
    int ReferenceCount() const {Storage->accessLock.Lock(); int refCnt = Storage->referenceCount; Storage->accessLock.Unlock(); return refCnt;}
    const char *Pointer() const {return Storage->Data;}

    // Allow A(i) and A(i,j), but check each access.
    char operator()(ssize_t i) const {if (i<0 || i>=Storage->length) {PrintErrorMessage(i); return invalidEntry;} return Storage->Data[i];}
    char operator()(ssize_t i,ssize_t j) const {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n) {PrintErrorMessage(i,j); return invalidEntry;} return Storage->Data[i+j*Storage->m];}
    char operator()(ssize_t i,ssize_t j,ssize_t k) const {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o) {PrintErrorMessage(i,j,k); return invalidEntry;} return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Debug functions, since gdb can't call the () operator.
    char e(int i) const;
    char e(int i,int j) const;
    char e(int i,int j,int k) const;
    void pinfo(void) const;
    void pi(int i) const; // (i,:)
    void pj(int j) const; // (:,j)
    void pall(void) const;  // Uses the same layout as DataTank in the variable monitor.
	void psigns(void) const; // -,0,+
    
    // Support for subregions
    const DTCharArrayRegion operator()(DTIndex) const;
    const DTCharArrayRegion operator()(DTIndex,DTIndex) const;
    const DTCharArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
protected:
    DTCharArrayStorage *Storage;
    char invalidEntry;

    // Error messages for index access.
    void PrintErrorMessage(ssize_t i) const;
    void PrintErrorMessage(ssize_t i,ssize_t j) const;
    void PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const;
};

class DTMutableCharArray : public DTCharArray
{
public:
    DTMutableCharArray() : DTCharArray() {}
    explicit DTMutableCharArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : DTCharArray(mv,nv,ov) {}
    DTMutableCharArray(const DTMutableCharArray &A) : DTCharArray(A) {}

    DTMutableCharArray &operator=(const DTMutableCharArray &A) {DTCharArray::operator=(A); return *this;}

    // Assignment
    DTMutableCharArray &operator=(char a);

    // Raw access
    char *Pointer() {return Storage->Data;}
    const char *Pointer() const {return Storage->Data;}

    // High level access
    char operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    char &operator()(ssize_t i)
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
        return Storage->Data[i];}
    char operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    char &operator()(ssize_t i,ssize_t j)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
        return Storage->Data[i+j*Storage->m];}
    char operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    char &operator()(ssize_t i,ssize_t j,ssize_t k)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
        return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    
    // Support for subregions
    const DTCharArrayRegion operator()(DTIndex) const;
    const DTCharArrayRegion operator()(DTIndex,DTIndex) const;
    const DTCharArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
    DTCharArrayRegion operator()(DTIndex);
    DTCharArrayRegion operator()(DTIndex,DTIndex);
    DTCharArrayRegion operator()(DTIndex,DTIndex,DTIndex);
};

extern bool operator==(const DTCharArray &A,const DTCharArray &B);
extern bool operator!=(const DTCharArray &A,const DTCharArray &B);

// Changing the size of an array
extern DTMutableCharArray TruncateSize(const DTCharArray &A,ssize_t length);
extern DTMutableCharArray IncreaseSize(const DTCharArray &A,ssize_t addLength);

// Misc utilities
extern DTMutableCharArray Transpose(const DTCharArray &A);
extern DTMutableCharArray FlipJ(const DTCharArray &A);
extern void Swap(DTMutableCharArray &,DTMutableCharArray &);
extern void Swap(DTCharArray &,DTCharArray &);
extern void CopyValues(DTMutableCharArray &into,const DTCharArray &from);

extern DTMutableCharArray ExtractIndices(const DTCharArray &,const DTRange &); // A single list, range in offsets

#endif
