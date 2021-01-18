// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTUShortIntArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"
#include "DTUtilities.h"

#include <cstring>
#include <iostream>
#include <string>

DTUShortIntArrayStorage::DTUShortIntArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
{
    // Check if it's called correctly.
    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    referenceCount = 1;
    mn = m*n;

    Data = length==0 ? NULL : new unsigned short int[(size_t)length];
}

DTUShortIntArrayStorage::~DTUShortIntArrayStorage()
{
    delete [] Data;
}

DTUShortIntArray &DTUShortIntArray::operator=(const DTUShortIntArray &A)
{
    // Allow A = A
    if (Storage==A.Storage) return *this;
    
    Storage->accessLock.Lock();
    A.Storage->accessLock.Lock();
    Storage->referenceCount--;
    int refCnt = Storage->referenceCount;
    Storage->accessLock.Unlock();
    if (refCnt==0) delete Storage;
    Storage = A.Storage;
    Storage->referenceCount++;
    Storage->accessLock.Unlock();
    
    return *this;
}

DTMutableUShortIntArray DTUShortIntArray::Copy() const
{
    DTMutableUShortIntArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(unsigned short int));
    return CopyInto;
}

unsigned short int DTUShortIntArray::e(int i) const
{
    if (i<0 || i>=Storage->length) {
#ifndef DG_NOSTDErrOut
        std::cerr << "Out of bounds" << std::endl;
#endif
        return invalidEntry;
    }
    else
        return Storage->Data[i];
}

unsigned short int DTUShortIntArray::e(int i,int j) const
{
    if (i<0 || i>=Storage->m || j<0 || j>=Storage->n) {
#ifndef DG_NOSTDErrOut
        std::cerr << "Out of bounds" << std::endl;
#endif
        return invalidEntry;
    }
    else
        return Storage->Data[i+j*Storage->m];
}

unsigned short int DTUShortIntArray::e(int i,int j,int k) const
{
    if (i<0 || i>=Storage->m || j<0 || j>=Storage->n || k<0 || k>=Storage->o) {
#ifndef DG_NOSTDErrOut
        std::cerr << "Out of bounds" << std::endl;
#endif
        return invalidEntry;
    }
    else
        return Storage->Data[i+j*Storage->m+k*Storage->mn];
}

void DTUShortIntArray::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    if (o()==0)
        std::cerr << "Empty" << std::endl;
    else if (o()==1) {
        if (n()==1)
            std::cerr << m() << " entries" << std::endl;
        else
            std::cerr << m() << " x " << n() << " array" << std::endl;
    }
    else
        std::cerr << m() << " x " << n() << " x " << o() << " array" << std::endl;
    std::cerr << std::flush;
#endif
}

void DTUShortIntArray::pi(int i) const
{
#ifndef DG_NOSTDErrOut
    if (i<0 || i>=m()) {
        std::cerr << "Out of bounds." << std::endl;
    }
    else {
        ssize_t howMany = n();
        ssize_t j;
        for (j=0;j<howMany-1;j++) std::cerr << operator()(i,j) << ", ";
        if (howMany>0) std::cerr << operator()(i,howMany-1);
        std::cerr << std::endl;
    }
#endif
}

void DTUShortIntArray::pj(int j) const
{
#ifndef DG_NOSTDErrOut
    if (j<0 || j>=n()) {
        std::cerr << "Out of bounds." << std::endl;
    }
    else {
        ssize_t howMany = m();
        ssize_t i;
        for (i=0;i<howMany-1;i++) std::cerr << operator()(i,j) << ", ";
        if (howMany>0) std::cerr << operator()(howMany-1,j);
        std::cerr << std::endl;
    }
#endif
}

void DTUShortIntArray::pall(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t mv = m();
    ssize_t nv = n();
    ssize_t i,j;
    if (mv==0) {
        std::cerr << "Empty" << std::endl;
    }
    else {
        for (j=0;j<nv;j++) {
            for (i=0;i<mv-1;i++) std::cerr << operator()(i,j) << ", ";
            std::cerr << operator()(mv-1,j);
            std::cerr << std::endl;
        }
    }
#endif
}

DTMutableUShortIntArray TruncateSize(const DTUShortIntArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableUShortIntArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableUShortIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableUShortIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableUShortIntArray();
        }
        newM = A.m();
        newN = length/A.m();
        newO = 1;
    }
    else {
        newM = length;
        newN = 1;
        newO = 1;
    }

    DTMutableUShortIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(unsigned short int));
    return toReturn;
}

DTMutableUShortIntArray IncreaseSize(const DTUShortIntArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableUShortIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableUShortIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableUShortIntArray();
        }
        newM = A.m();
        newN = A.n() + addLength/A.m();
        newO = 1;
    }
    else {
        newM = A.m() + addLength;
        newN = 1;
        newO = 1;
    }

    DTMutableUShortIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(unsigned short int));
    return toReturn;
}

void DTUShortIntArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTUShortIntArray",i,Storage->length);
}

void DTUShortIntArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTUShortIntArray",i,j,Storage->m,Storage->n);
}

void DTUShortIntArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTUShortIntArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableUShortIntArray &DTMutableUShortIntArray::operator=(unsigned short int a)
{
    const ssize_t howManyNumbers = Length();
    ssize_t i;
    unsigned short int *Data = Pointer();
    for (i=0;i<howManyNumbers;i++)
        Data[i] = a;
    
    return *this;
}

bool operator==(const DTUShortIntArray &A,const DTUShortIntArray &B)
{
    return DTOperatorArrayEqualsArray<DTUShortIntArray,unsigned short int>(A,B);
}

bool operator!=(const DTUShortIntArray &A,const DTUShortIntArray &B)
{
    return !(A==B);
}

DTMutableUShortIntArray Transpose(const DTUShortIntArray &A)
{
    return DTTransposeArray<DTUShortIntArray,DTMutableUShortIntArray,unsigned short int>(A);
}

DTMutableUShortIntArray FlipJ(const DTUShortIntArray &A)
{
    return DTArrayFlipJ<DTUShortIntArray,DTMutableUShortIntArray,unsigned short int>(A);
}

void Swap(DTMutableUShortIntArray &A,DTMutableUShortIntArray &B)
{
	DTMutableUShortIntArray C = A;
	A = B;
	B = C;
}

void Swap(DTUShortIntArray &A,DTUShortIntArray &B)
{
	DTUShortIntArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableUShortIntArray &into,const DTUShortIntArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableUShortIntArray,UShortIntArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length()*sizeof(unsigned short int));
	}
}

DTMutableUShortIntArray ExtractIndices(const DTUShortIntArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(DTUShortIntArray,Range)","Range is out of bounds");
        return DTMutableUShortIntArray();
    }
    
    DTMutableUShortIntArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start, (size_t)r.length*sizeof(unsigned short int));
    return toReturn;
}

unsigned short int Minimum(const DTUShortIntArray &A)
{
    ssize_t len = A.Length();
    unsigned short int minV = 32767;
    
    unsigned short int v;
    ssize_t i;
    
    const unsigned short int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }
    
    return minV;
}

unsigned short int Maximum(const DTUShortIntArray &A)
{
    ssize_t len = A.Length();
    unsigned short int maxV = 0;
    
    unsigned short int v;
    ssize_t i;
    
    const unsigned short int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV);
    }
    
    return maxV;
}

DTMutableUShortIntArray operator+(const DTUShortIntArray &A,const DTUShortIntArray &B)
{
    return DTAddArrays<DTUShortIntArray,DTMutableUShortIntArray,unsigned short>("UShortIntArray+UShortIntArray",A,B);
}

