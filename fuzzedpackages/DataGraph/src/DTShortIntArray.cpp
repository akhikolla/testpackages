// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTShortIntArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"
#include "DTUtilities.h"

#include <cstring>

DTShortIntArrayStorage::DTShortIntArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
{
    // Check if it's called correctly.
    if (mv<0 || nv<0 || ov<0) DTErrorMessage("DTMutableShortIntArray", "Negative index in constructor");

    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    referenceCount = 1;
    mn = m*n;

    Data = length==0 ? NULL : new short int[(size_t)length];
}

DTShortIntArrayStorage::~DTShortIntArrayStorage()
{
    delete [] Data;
}

DTShortIntArray &DTShortIntArray::operator=(const DTShortIntArray &A)
{
    // Allow A = A
    if (Storage!=A.Storage) {
        Storage->accessLock.Lock();
        A.Storage->accessLock.Lock();
        Storage->referenceCount--;
        int refCnt = Storage->referenceCount;
        Storage->accessLock.Unlock();
        if (refCnt==0) delete Storage;
        Storage = A.Storage;
        Storage->referenceCount++;
        Storage->accessLock.Unlock();
    }
    
    return *this;
}

DTMutableShortIntArray DTShortIntArray::Copy() const
{
    DTMutableShortIntArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(short int));
    return CopyInto;
}

short int DTShortIntArray::e(int i) const
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

short int DTShortIntArray::e(int i,int j) const
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

short int DTShortIntArray::e(int i,int j,int k) const
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

void DTShortIntArray::pinfo(void) const
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

void DTShortIntArray::pi(int i) const
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

void DTShortIntArray::pj(int j) const
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

void DTShortIntArray::pall(void) const
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

DTMutableShortIntArray TruncateSize(const DTShortIntArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableShortIntArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableShortIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableShortIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableShortIntArray();
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

    DTMutableShortIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(short int));
    return toReturn;
}

DTMutableShortIntArray IncreaseSize(const DTShortIntArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableShortIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableShortIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableShortIntArray();
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

    DTMutableShortIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(short int));
    return toReturn;
}

void DTShortIntArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTShortIntArray",i,Storage->length);
}

void DTShortIntArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTShortIntArray",i,j,Storage->m,Storage->n);
}

void DTShortIntArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTShortIntArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableShortIntArray &DTMutableShortIntArray::operator=(short int a)
{
    const ssize_t howManyNumbers = Length();
    ssize_t i;
    short int *Data = Pointer();
    for (i=0;i<howManyNumbers;i++)
        Data[i] = a;
    
    return *this;
}

bool operator==(const DTShortIntArray &A,const DTShortIntArray &B)
{
    return DTOperatorArrayEqualsArray<DTShortIntArray,short int>(A,B);
}

bool operator!=(const DTShortIntArray &A,const DTShortIntArray &B)
{
    return !(A==B);
}

DTMutableShortIntArray Transpose(const DTShortIntArray &A)
{
    return DTTransposeArray<DTShortIntArray,DTMutableShortIntArray,short int>(A);
}

DTMutableShortIntArray FlipJ(const DTShortIntArray &A)
{
    return DTArrayFlipJ<DTShortIntArray,DTMutableShortIntArray,short int>(A);
}

void Swap(DTMutableShortIntArray &A,DTMutableShortIntArray &B)
{
	DTMutableShortIntArray C = A;
	A = B;
	B = C;
}

void Swap(DTShortIntArray &A,DTShortIntArray &B)
{
	DTShortIntArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableShortIntArray &into,const DTShortIntArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableShortIntArray,ShortIntArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length()*sizeof(short int));
	}
}

void MemoryCopy(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from)
{
    MemoryCopy(into,intoLocation,from,DTRange(0,from.Length()));
}

void MemoryCopy(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from,const DTRange &range)
{
    // Check validty
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.Length()) {
        DTErrorMessage("MemoryCopy","Invalid destination range");
        return;
    }
    if (range.end()>from.Length()) {
        DTErrorMessage("MemoryCopy","Invalid source range");
        return;
    }
    if (range.length) std::memcpy(into.Pointer()+intoLocation,from.Pointer()+range.start,(size_t)range.length*sizeof(short int));
}

void MemoryCopyColumns(DTMutableShortIntArray &into,ssize_t intoLocation,const DTShortIntArray &from,const DTRange &range)
{
    // Check validty
    if (into.m()!=from.m()) {
        DTErrorMessage("MemoryCopyColumns","Need to be the same number of columns");
        return;
    }
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.n()) {
        DTErrorMessage("MemoryCopyColumns","Invalid destination range");
        return;
    }
    if (range.end()>from.n()) {
        DTErrorMessage("MemoryCopyColumns","Invalid source range");
        return;
    }
    if (range.length) std::memcpy(into.Pointer()+into.m()*intoLocation,from.Pointer()+into.m()*range.start,(size_t)(range.length*into.m())*sizeof(short int));
}

void MemoryMove(DTMutableShortIntArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid destination range");
        return;
    }
    if (range.end()>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+intoLocation,into.Pointer()+range.start,(size_t)range.length*sizeof(short int));
}

void MemoryMoveColumns(DTMutableShortIntArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid destination range");
        return;
    }
    if (range.end()>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+into.m()*intoLocation,into.Pointer()+into.m()*range.start,(size_t)(range.length*into.m())*sizeof(short int));
}

DTMutableShortIntArray ExtractIndices(const DTShortIntArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(DTShortIntArray,Range)","Range is out of bounds");
        return DTMutableShortIntArray();
    }
    
    DTMutableShortIntArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start, (size_t)r.length*sizeof(short int));
    return toReturn;
}

short int Minimum(const DTShortIntArray &A)
{
    ssize_t len = A.Length();
    short int minV = 32767;
    
    short int v;
    ssize_t i;
    
    const short int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }
    
    return minV;
}

short int Maximum(const DTShortIntArray &A)
{
    ssize_t len = A.Length();
    short int maxV = 0;
    
    short int v;
    ssize_t i;
    
    const short int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV);
    }
    
    return maxV;
}


