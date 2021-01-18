// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTUCharArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"

#include <cstring>

DTUCharArrayStorage::DTUCharArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
{
    // Check if it's called correctly.
    if (mv<0 || nv<0 || ov<0) DTErrorMessage("DTMutableUCharArray", "Negative index in constructor");
    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    referenceCount = 1;
    mn = m*n;

    Data = length==0 ? NULL : new unsigned char[(size_t)length];
}

DTUCharArrayStorage::~DTUCharArrayStorage()
{
    delete [] Data;
}

DTUCharArray &DTUCharArray::operator=(const DTUCharArray &A)
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

DTMutableUCharArray DTUCharArray::Copy() const
{
    DTMutableUCharArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(unsigned char));
    return CopyInto;
}

unsigned char DTUCharArray::e(int i) const
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

unsigned char DTUCharArray::e(int i,int j) const
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

unsigned char DTUCharArray::e(int i,int j,int k) const
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

void DTUCharArray::pinfo(void) const
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

void DTUCharArray::pi(int i) const
{
#ifndef DG_NOSTDErrOut
    if (i<0 || i>=m()) {
        std::cerr << "Out of bounds." << std::endl;
    }
    else {
        ssize_t howMany = n();
        ssize_t j;
        for (j=0;j<howMany-1;j++) std::cerr << (int)operator()(i,j) << ", ";
        if (howMany>0) std::cerr << (int)operator()(i,howMany-1);
        std::cerr << std::endl;
    }
#endif
}

void DTUCharArray::pj(int j) const
{
#ifndef DG_NOSTDErrOut
    if (j<0 || j>=n()) {
        std::cerr << "Out of bounds." << std::endl;
    }
    else {
        ssize_t howMany = m();
        ssize_t i;
        for (i=0;i<howMany-1;i++) std::cerr << (int)operator()(i,j) << ", ";
        if (howMany>0) std::cerr << (int)operator()(howMany-1,j);
        std::cerr << std::endl;
    }
#endif
}

void DTUCharArray::pall(void) const
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
            for (i=0;i<mv-1;i++) std::cerr << (int)operator()(i,j) << ", ";
            std::cerr << (int)operator()(mv-1,j);
            std::cerr << std::endl;
        }
    }
#endif
}

DTMutableUCharArray TruncateSize(const DTUCharArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableUCharArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableUCharArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableUCharArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableUCharArray();
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

    DTMutableUCharArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(unsigned char));
    return toReturn;
}

DTMutableUCharArray IncreaseSize(const DTUCharArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableUCharArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableUCharArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableUCharArray();
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

    DTMutableUCharArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(unsigned char));
    return toReturn;
}

void DTUCharArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTUCharArray",i,Storage->length);
}

void DTUCharArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTUCharArray",i,j,Storage->m,Storage->n);
}

void DTUCharArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTUCharArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableUCharArray &DTMutableUCharArray::operator=(unsigned char a)
{
    const ssize_t howManyNumbers = Length();
    ssize_t i;
    unsigned char *Data = Pointer();
    if (a==0) {
        memset(Data,0,(size_t)howManyNumbers);
    }
    else {
        for (i=0;i<howManyNumbers;i++)
            Data[i] = a;
    }
    return *this;
}

bool operator==(const DTUCharArray &A,const DTUCharArray &B)
{
    return DTOperatorArrayEqualsArray<DTUCharArray,unsigned char>(A,B);
}

bool operator!=(const DTUCharArray &A,const DTUCharArray &B)
{
    return !(A==B);
}

DTMutableUCharArray Transpose(const DTUCharArray &A)
{
    return DTTransposeArray<DTUCharArray,DTMutableUCharArray,unsigned char>(A);
}

DTMutableUCharArray FlipJ(const DTUCharArray &A)
{
    return DTArrayFlipJ<DTUCharArray,DTMutableUCharArray,unsigned char>(A);
}

void Swap(DTMutableUCharArray &A,DTMutableUCharArray &B)
{
	DTMutableUCharArray C = A;
	A = B;
	B = C;
}

void Swap(DTUCharArray &A,DTUCharArray &B)
{
	DTUCharArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableUCharArray &into,const DTUCharArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableUCharArray,UCharArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length());
	}
}

unsigned char Minimum(const DTUCharArray &A)
{
    ssize_t len = A.Length();
    unsigned char minV = 255;
    
    unsigned char v;
    ssize_t i;
    
    const unsigned char *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }
    
    return minV;
}

unsigned char Maximum(const DTUCharArray &A)
{
    ssize_t len = A.Length();
    unsigned char maxV = 0;
    
    unsigned char v;
    ssize_t i;
    
    const unsigned char *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV);
    }
    
    return maxV;
}


