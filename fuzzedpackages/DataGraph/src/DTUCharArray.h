// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTUCharArray_Header
#define DTUCharArray_Header

// By default, range check is turned on.
#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

#include <iostream>
#include <unistd.h>
#include "DTLock.h"

// An array of unsigned char numbers.  See comments inside DTDoubleArray for more information.
class DTUCharArrayStorage {
public:
    DTUCharArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov);
    ~DTUCharArrayStorage();

    DTLock accessLock;
    ssize_t m,n,o,mn,length;
    int referenceCount;
    unsigned char *Data;
    
private:
    DTUCharArrayStorage(const DTUCharArrayStorage &);
    DTUCharArrayStorage &operator=(const DTUCharArrayStorage &);
};

class DTMutableUCharArray;
class DTIndex;
class DTUCharArrayRegion;

class DTUCharArray {

public:
    DTUCharArray() : Storage(new DTUCharArrayStorage(0,0,0)), invalidEntry(0) {}
    virtual ~DTUCharArray() {Storage->accessLock.Lock(); int refCnt = (--Storage->referenceCount); Storage->accessLock.Unlock(); if (refCnt==0) delete Storage;}
    DTUCharArray(const DTUCharArray &A) : Storage(A.Storage), invalidEntry(0) {Storage->accessLock.Lock(); Storage->referenceCount++;Storage->accessLock.Unlock(); }
    DTUCharArray &operator=(const DTUCharArray &A);

protected:
    // If you get a notice that this is protected, change DTUCharArray to DTMutableUCharArray
    explicit DTUCharArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : Storage(new DTUCharArrayStorage(mv,nv,ov)), invalidEntry(0) {}

public:
    DTMutableUCharArray Copy() const;

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
    const unsigned char *Pointer() const {return Storage->Data;}

    // Allow A(i) and A(i,j), but check each access.
    unsigned char operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    unsigned char operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    unsigned char operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}

    // Debug functions, since gdb can't call the () operator.
    unsigned char e(int i) const;
    unsigned char e(int i,int j) const;
    unsigned char e(int i,int j,int k) const;
    void pinfo(void) const;
    void pi(int i) const; // (i,:)
    void pj(int j) const; // (:,j)
    void pall(void) const;  // Uses the same layout as DataTank in the variable monitor.
    
    // Support for subregions
    const DTUCharArrayRegion operator()(DTIndex) const;
    const DTUCharArrayRegion operator()(DTIndex,DTIndex) const;
    const DTUCharArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
protected:
    DTUCharArrayStorage *Storage;
    unsigned char invalidEntry;

    // Error messages for index access.
    void PrintErrorMessage(ssize_t i) const;
    void PrintErrorMessage(ssize_t i,ssize_t j) const;
    void PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const;
};

class DTMutableUCharArray : public DTUCharArray
{
public:
    DTMutableUCharArray() : DTUCharArray() {}
    explicit DTMutableUCharArray(ssize_t mv,ssize_t nv=1,ssize_t ov=1) : DTUCharArray(mv,nv,ov) {}
    DTMutableUCharArray(const DTMutableUCharArray &A) : DTUCharArray(A) {}

    DTMutableUCharArray &operator=(const DTMutableUCharArray &A) {DTUCharArray::operator=(A); return *this;}

    // Assignment
    DTMutableUCharArray &operator=(unsigned char a);

    // Raw access
    unsigned char *Pointer() {return Storage->Data;}
    const unsigned char *Pointer() const {return Storage->Data;}

    // High level access
    unsigned char operator()(ssize_t i) const
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
            return Storage->Data[i];}
    unsigned char &operator()(ssize_t i)
        {if (i<0 || i>=Storage->length)
            {PrintErrorMessage(i); return invalidEntry;}
        return Storage->Data[i];}
    unsigned char operator()(ssize_t i,ssize_t j) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
            return Storage->Data[i+j*Storage->m];}
    unsigned char &operator()(ssize_t i,ssize_t j)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n)
            {PrintErrorMessage(i,j); return invalidEntry;}
        return Storage->Data[i+j*Storage->m];}
    unsigned char operator()(ssize_t i,ssize_t j,ssize_t k) const
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
            return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    unsigned char &operator()(ssize_t i,ssize_t j,ssize_t k)
        {if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o)
            {PrintErrorMessage(i,j,k); return invalidEntry;}
        return Storage->Data[i+j*Storage->m+k*Storage->mn];}
    
    // Support for subregions
    const DTUCharArrayRegion operator()(DTIndex) const;
    const DTUCharArrayRegion operator()(DTIndex,DTIndex) const;
    const DTUCharArrayRegion operator()(DTIndex,DTIndex,DTIndex) const;
    
    DTUCharArrayRegion operator()(DTIndex);
    DTUCharArrayRegion operator()(DTIndex,DTIndex);
    DTUCharArrayRegion operator()(DTIndex,DTIndex,DTIndex);
};

// Misc
extern DTMutableUCharArray Transpose(const DTUCharArray &A);
extern bool operator==(const DTUCharArray &A,const DTUCharArray &B);
extern bool operator!=(const DTUCharArray &A,const DTUCharArray &B); 

extern unsigned char Minimum(const DTUCharArray &);
extern unsigned char Maximum(const DTUCharArray &);

extern void Swap(DTMutableUCharArray &,DTMutableUCharArray &);
extern void Swap(DTUCharArray &,DTUCharArray &);
extern void CopyValues(DTMutableUCharArray &into,const DTUCharArray &from);

// Changing the size of an array
extern DTMutableUCharArray TruncateSize(const DTUCharArray &A,ssize_t length);
extern DTMutableUCharArray IncreaseSize(const DTUCharArray &A,ssize_t addLength);
extern DTMutableUCharArray FlipJ(const DTUCharArray &A);

#endif
