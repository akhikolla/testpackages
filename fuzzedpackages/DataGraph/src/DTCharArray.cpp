// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTCharArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"
#include "DTUtilities.h"

#include <cstring>

DTCharArrayStorage::DTCharArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
: accessLock(), m(0), n(0), o(0), mn(0), length(0), referenceCount(1), Data(NULL) {
    // Check if it's called correctly.
    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    mn = m*n;

    Data = length==0 ? NULL : new char[(size_t)length];
}

DTCharArrayStorage::~DTCharArrayStorage()
{
    delete [] Data;
}

DTCharArray &DTCharArray::operator=(const DTCharArray &A)
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

DTMutableCharArray DTCharArray::Copy() const
{
    DTMutableCharArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(char));
    return CopyInto;
}

char DTCharArray::e(int i) const
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

char DTCharArray::e(int i,int j) const
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

char DTCharArray::e(int i,int j,int k) const
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

void DTCharArray::pinfo(void) const
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

void DTCharArray::pi(int i) const
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

void DTCharArray::pj(int j) const
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

void DTCharArray::pall(void) const
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

void DTCharArray::psigns(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t mv = m();
    ssize_t nv = n();
    ssize_t i,j;
    if (mv==0) {
        std::cerr << "Empty" << std::endl;
    }
    else if (o()>1) {
        std::cerr << "only works for 2D arrays" << std::endl;
    }
    else {
		char val;
        for (j=0;j<nv;j++) {
            for (i=0;i<mv;i++) {
				val = operator()(i,j);
				if (val<0)
					std::cerr << "-";
				else if (val==0)
					std::cerr << "0";
				else
					std::cerr << "+";
			}
            std::cerr << std::endl;
        }
    }
#endif
}

DTMutableCharArray TruncateSize(const DTCharArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableCharArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableCharArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableCharArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableCharArray();
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

    DTMutableCharArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(char));
    return toReturn;
}

DTMutableCharArray IncreaseSize(const DTCharArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableCharArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableCharArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableCharArray();
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

    DTMutableCharArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(char));
    return toReturn;
}

void DTCharArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTCharArray",i,Storage->length);
}

void DTCharArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTCharArray",i,j,Storage->m,Storage->n);
}

void DTCharArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTCharArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableCharArray &DTMutableCharArray::operator=(char a)
{
    const ssize_t howManyNumbers = Length();
    memset(Pointer(),a,(size_t)howManyNumbers);
    return *this;
}

bool operator==(const DTCharArray &A,const DTCharArray &B)
{
    return DTOperatorArrayEqualsArray<DTCharArray,char>(A,B);
}

bool operator!=(const DTCharArray &A,const DTCharArray &B)
{
    return !(A==B);
}

DTMutableCharArray Transpose(const DTCharArray &A)
{
    return DTTransposeArray<DTCharArray,DTMutableCharArray,char>(A);
}

DTMutableCharArray FlipJ(const DTCharArray &A)
{
    return DTArrayFlipJ<DTCharArray,DTMutableCharArray,char>(A);
}

void Swap(DTMutableCharArray &A,DTMutableCharArray &B)
{
	DTMutableCharArray C = A;
	A = B;
	B = C;
}

void Swap(DTCharArray &A,DTCharArray &B)
{
	DTCharArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableCharArray &into,const DTCharArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableCharArray,CharArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length());
	}
}

DTMutableCharArray ExtractIndices(const DTCharArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(DTCharArray,Range)","Range is out of bounds");
        return DTMutableCharArray();
    }
    
    DTMutableCharArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start,(size_t)r.length);
    return toReturn;
}

