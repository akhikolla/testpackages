// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

// Needed to implement INT32_MAX in a more portable way
#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "DTIntArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"
#include "DTUtilities.h"

#include <cstring>
#include <algorithm>


DTIntArrayStorage::DTIntArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
{
    // Check if it's called correctly.
    if (mv<0 || nv<0 || ov<0) DTErrorMessage("DTMutableDoubleArray", "Negative index in constructor");

    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    referenceCount = 1;
    mn = m*n;

    Data = length==0 ? NULL : new int[(size_t)length];
}

DTIntArrayStorage::~DTIntArrayStorage()
{
    delete [] Data;
}

DTIntArray &DTIntArray::operator=(const DTIntArray &A)
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

DTMutableIntArray DTIntArray::Copy() const
{
    DTMutableIntArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(int));
    return CopyInto;
}

int DTIntArray::e(int i) const
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

int DTIntArray::e(int i,int j) const
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

int DTIntArray::e(int i,int j,int k) const
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

void DTIntArray::pinfo(void) const
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

void DTIntArray::pi(int i) const
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

void DTIntArray::pj(int j) const
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

void DTIntArray::pall(void) const
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

void DTIntArray::prange(int startIndex,int endIndex) const
{
#ifndef DG_NOSTDErrOut
    // Offsets [s,e], both included.
    ssize_t i;
    if (startIndex<0) {
        std::cerr << "start out of bounds" << std::endl;
        return;
    }
    if (endIndex>=Length()) {
        std::cerr << "end out of bounds" << std::endl;
        return;
    }
    for (i=startIndex;i<endIndex;i++) {
        std::cerr << operator()(i) << ", ";
    }
    if (startIndex<=endIndex) {
        std::cerr << operator()(i) << std::endl;
    }
#endif
}

void DTIntArray::pcrange(int startIndex,int endIndex) const
{
#ifndef DG_NOSTDErrOut
    // Print (:,[s:e]), otherwise like pall().
    ssize_t mv = m();
    ssize_t nv = n();
    ssize_t i,j;
    if (startIndex<0) {
        std::cerr << "start out of bounds" << std::endl;
        return;
    }
    if (endIndex>=nv) {
        std::cerr << "end out of bounds" << std::endl;
        return;
    }

    for (j=startIndex;j<=endIndex;j++) {
        for (i=0;i<mv-1;i++) std::cerr << operator()(i,j) << ", ";
        std::cerr << operator()(mv-1,j);
        std::cerr << std::endl;
    }
#endif
}

void DTIntArray::psigns(void) const
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
		int val;
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

ssize_t DTIntArray::Find(int v) const
{
    const int *D = Pointer();
    ssize_t len = Length();
    ssize_t i;
    for (i=0;i<len;i++) {
        if (D[i]==v) break;
    }
    
    return (i<len ? i : -1);
}

DTMutableIntArray CombineRows(const DTIntArray &First,const DTIntArray &Second)
{
    if (First.IsEmpty()) return Second.Copy();
    if (Second.IsEmpty()) return First.Copy();
    if (First.n()!=Second.n()) {
        DTErrorMessage("CombineRows(A,B)","A and B have to have the same number of columns.");
        return DTMutableIntArray();
    }
    if (First.o()!=1 || Second.o()!=1) {
        DTErrorMessage("CombineRows(A,B)","A and B have to be one or two dimensional.");
        return DTMutableIntArray();
    }
    
    DTMutableIntArray toReturn(First.m()+Second.m(),First.n());
    if (First.n()==1) {
        std::memcpy(toReturn.Pointer(),First.Pointer(),(size_t)First.Length()*sizeof(int));
        std::memcpy(toReturn.Pointer()+First.Length(),Second.Pointer(),(size_t)Second.Length()*sizeof(int));
    }
    else {
        ssize_t j,n = First.n();
        const int *FirstD = First.Pointer();
        const int *SecondD = Second.Pointer();
        int *toReturnD = toReturn.Pointer();
        ssize_t mFirst = First.m();
        ssize_t mSecond = Second.m();
        ssize_t stride = toReturn.m();
        for (j=0;j<n;j++) {
            std::memcpy(toReturnD+j*stride,FirstD+j*mFirst,(size_t)mFirst*4);
            std::memcpy(toReturnD+j*stride+mFirst,SecondD+j*mSecond,(size_t)mSecond*4);
        }
    }
    
    return toReturn;
}

DTMutableIntArray CombineColumns(const DTIntArray &First,const DTIntArray &Second)
{
    if (First.IsEmpty()) return Second.Copy();
    if (Second.IsEmpty()) return First.Copy();
    if (First.m()!=Second.m()) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to have the same number of rows.");
        return DTMutableIntArray();
    }
    if (First.o()!=1 || Second.o()!=1) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to be two dimensional.");
        return DTMutableIntArray();
    }
    
    DTMutableIntArray toReturn(First.m(),First.n()+Second.n());
    std::memcpy(toReturn.Pointer(),First.Pointer(),(size_t)First.Length()*sizeof(int));
    std::memcpy(toReturn.Pointer()+First.Length(),Second.Pointer(),(size_t)Second.Length()*sizeof(int));
    
    return toReturn;
}

DTMutableIntArray TruncateSize(const DTIntArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableIntArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableIntArray();
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

    DTMutableIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(int));
    return toReturn;
}

DTMutableIntArray IncreaseSize(const DTIntArray &A)
{
    return IncreaseSize(A,A.Length());
}

DTMutableIntArray IncreaseSize(const DTIntArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableIntArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableIntArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableIntArray();
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

    DTMutableIntArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(int));
    return toReturn;
}

void DTIntArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTIntArray",i,Storage->length);
}

void DTIntArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTIntArray",i,j,Storage->m,Storage->n);
}

void DTIntArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTIntArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableIntArray &DTMutableIntArray::operator=(int a)
{
    const ssize_t howManyNumbers = Length();
    ssize_t i;
    int *Data = Pointer();
    if (a==0) {
        memset(Data,0,(size_t)howManyNumbers*sizeof(int));
    }
    else {
        for (i=0;i<howManyNumbers;i++)
        Data[i] = a;
    }
    
    return *this;
}

void DTMutableIntArray::operator+=(int v)
{
    DTPlusEqualsScalar<DTMutableIntArray,int>(*this,v);
}

void DTMutableIntArray::operator-=(int v)
{
    DTMinusEqualsScalar<DTMutableIntArray,int>(*this,v);
}

bool operator==(const DTIntArray &A,const DTIntArray &B)
{
    return DTOperatorArrayEqualsArray<DTIntArray,int>(A,B);
}

bool operator!=(const DTIntArray &A,const DTIntArray &B)
{
    return !(A==B);
}

DTMutableIntArray Transpose(const DTIntArray &A)
{
    return DTTransposeArray<DTIntArray,DTMutableIntArray,int>(A);
}

DTMutableIntArray Reshape(const DTIntArray &A,ssize_t m,ssize_t n,ssize_t o)
{
    /*
    if (m<0 || n<0 || o<0) {
        DTErrorMessage("Reshape(DTIntArray,...)","One of the new dimensions is negative.");
        return DTMutableIntArray();
    }
     */
    if (m*n*o!=A.Length()) {
        DTErrorMessage("Reshape(DTIntArray,...)","Size before and after need to be the same.");
        return DTMutableIntArray();
    }

    DTMutableIntArray toReturn(m,n,o);
    if (toReturn.Length()) {
        std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(int));
    }

    return toReturn;
}

DTMutableIntArray Sort(const DTIntArray &A)
{
    DTMutableIntArray toReturn = Reshape(A,A.Length());;
    std::sort(toReturn.Pointer(),toReturn.Pointer()+toReturn.Length());
    return toReturn;
}

DTMutableIntArray FlipJ(const DTIntArray &A)
{
    return DTArrayFlipJ<DTIntArray,DTMutableIntArray,int>(A);
}

ssize_t FindEntry(const DTIntArray &A,int val)
{
    ssize_t Len = A.Length();
    ssize_t i;
    for (i=0;i<Len;i++) {
        if (A(i)==val) break;
    }
    return (i==Len ? -1 : i);
}

ssize_t FindEntryInSorted(const DTIntArray &A,int val)
{
    if (A.Length()==0 || A(0)>val || A(A.Length()-1)<val)
        return -1;
    
    ssize_t StrictlyBefore = A.Length();
    ssize_t AfterOrEqual = 0;
    ssize_t LookAt;
    while (StrictlyBefore-AfterOrEqual>1) {
        LookAt = (AfterOrEqual+StrictlyBefore)/2;
        if (val<A(LookAt))
            StrictlyBefore = LookAt;
        else
            AfterOrEqual = LookAt;
    }
    
    if (A(AfterOrEqual)==val)
        return AfterOrEqual;
    else
        return -1;
}

void Swap(DTMutableIntArray &A,DTMutableIntArray &B)
{
	DTMutableIntArray C = A;
	A = B;
	B = C;
}

void Swap(DTIntArray &A,DTIntArray &B)
{
	DTIntArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableIntArray &into,const DTIntArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableIntArray,IntArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length()*sizeof(int));
	}
}

void CopyValuesIntoAndAdd(DTMutableIntArray &into,ssize_t offset,const DTIntArray &from,ssize_t add)
{
	if (offset<0 || offset+from.Length()>into.Length()) {
		DTErrorMessage("CopyValuesInto(array,offset,array)","Copying outside the valid range.");
		return;
	}
	int *D = into.Pointer()+offset;
	std::memcpy(D,from.Pointer(),sizeof(int)*(size_t)from.Length());
	ssize_t i,howMany = from.Length();
	for (i=0;i<howMany;i++)
		D[i] += add;
}

void CopyIntoColumn(DTMutableIntArray &into,const DTIntArray &list,ssize_t j)
{
    if (into.m()!=list.Length()) {
		DTErrorMessage("CopyIntoColumns(into,list,j)","Need into.m()==list.Length()");
    }
    else if (into.o()!=1) {
		DTErrorMessage("CopyIntoColumns(into,list,j)","into is a 3D array (into.o()>1)");
    }
    else if (j<0 || j>into.n()) {
        DTErrorMessage("CopyIntoColumns(into,list,j)","j out of bounds");
    }
    else {
        std::memcpy(into.Pointer()+j*into.m(),list.Pointer(),(size_t)into.m()*sizeof(int));
    }
}

void CopyIntoColumns(DTMutableIntArray &into,const DTRange &intoRange,const DTIntArray &from,const DTRange &fromRange)
{
    if (into.o()!=1 || from.o()!=1) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","into is a 3D array (into.o()>1 or from.o()>1)");
    }
    else if (into.n()<intoRange.end() || from.n()<fromRange.end()) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","Out of bounds");
    }
    else if (intoRange.length!=fromRange.length) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","fromRange.length!=toRange.length");
    }
    else if (into.m()!=from.m()) {
        DTErrorMessage("CopyIntoColumns(into,range,from,range)","from.m()!=to.m()");
    }
    else {
        std::memcpy(into.Pointer()+intoRange.start*into.m(),from.Pointer()+fromRange.start*from.m(),(size_t)(from.m()*intoRange.length)*sizeof(int));
    }
}

void MemoryCopy(DTMutableIntArray &into,ssize_t intoLocation,const DTIntArray &from,ssize_t fromLocation,ssize_t numberOfEntries)
{
    // Check validty
    if (numberOfEntries==0) return;
    if (numberOfEntries<0) {
        DTErrorMessage("MemoryCopy","Invalid number of entries to copy");
        return;
    }
    if (intoLocation<0 || intoLocation+numberOfEntries>into.Length()) {
        DTErrorMessage("MemoryCopy","Invalid destination range");
        return;
    }
    if (fromLocation<0 || fromLocation+numberOfEntries>from.Length()) {
        DTErrorMessage("MemoryCopy","Invalid source range");
        return;
    }
    std::memcpy(into.Pointer()+intoLocation,from.Pointer()+fromLocation,(size_t)numberOfEntries*sizeof(int));
}

void MemoryCopyColumns(DTMutableIntArray &into,ssize_t intoLocation,const DTIntArray &from,ssize_t fromLocation,ssize_t numberOfColumns)
{
    // Check validty
    if (numberOfColumns==0) return;
    if (into.m()!=from.m()) {
        DTErrorMessage("MemoryCopyColumns","Need to be the same number of columns");
        return;
    }
    if (numberOfColumns<0) {
        DTErrorMessage("MemoryCopyColumns","Invalid number of entries to copy");
        return;
    }
    if (intoLocation<0 || intoLocation+numberOfColumns>into.n()) {
        DTErrorMessage("MemoryCopyColumns","Invalid destination range");
        return;
    }
    if (fromLocation<0 || fromLocation+numberOfColumns>from.n()) {
        DTErrorMessage("MemoryCopyColumns","Invalid source range");
        return;
    }
    std::memcpy(into.Pointer()+into.m()*intoLocation,from.Pointer()+into.m()*fromLocation,(size_t)(numberOfColumns*into.m())*sizeof(int));
}

void MemoryMove(DTMutableIntArray &into,ssize_t intoLocation,ssize_t fromLocation,ssize_t numberOfEntries)
{
    if (numberOfEntries==0) return;
    if (numberOfEntries<0) {
        DTErrorMessage("MemoryMove","Invalid number of entries to copy");
        return;
    }
    if (intoLocation<0 || intoLocation+numberOfEntries>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid destination range");
        return;
    }
    if (fromLocation<0 || fromLocation+numberOfEntries>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid source range");
        return;
    }
    memmove(into.Pointer()+intoLocation,into.Pointer()+fromLocation,(size_t)numberOfEntries*sizeof(int));
}

void MemoryMoveColumns(DTMutableIntArray &into,ssize_t intoLocation,ssize_t fromLocation,ssize_t numberOfColumns)
{
    if (numberOfColumns==0) return;
    if (numberOfColumns<0) {
        DTErrorMessage("MemoryMoveColumns","Invalid number of entries to copy");
        return;
    }
    if (intoLocation<0 || intoLocation+numberOfColumns>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid destination range");
        return;
    }
    if (fromLocation<0 || fromLocation+numberOfColumns>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid source range");
        return;
    }
    memmove(into.Pointer()+into.m()*intoLocation,into.Pointer()+into.m()*fromLocation,(size_t)(numberOfColumns*into.m())*sizeof(int));
}

DTMutableIntArray ExtractColumns(const DTIntArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableIntArray();
        }
        else {
            DTErrorMessage("ExtractColumns(IntArray,IntArray)","Int array is empty");
            return DTMutableIntArray();
        }
    }
    
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(IntArray,IntArray)","Does not work for 3D arrays");
        return DTMutableIntArray();
    }
    
    ssize_t i,j,howMany = indices.Length();
    ssize_t index;
    ssize_t m = A.m();
    ssize_t n = A.n();
    
    DTMutableIntArray toReturn(A.m(),howMany);
    bool outOfBounds = false;
    for (j=0;j<howMany;j++) {
        index = indices(j);
        if (index<0 || index>=n) {
            outOfBounds = true; // Set break point here
            for (i=0;i<m;i++) toReturn(i,j) = INT32_MAX;
        }
        else {
            for (i=0;i<m;i++) toReturn(i,j) = A(i,index);
        }
    }
    if (outOfBounds) {
        DTErrorMessage("ExtractColumns(IntArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

DTMutableIntArray ExtractColumns(const DTIntArray &A,const DTRange &r)
{
    if (r.end()>A.n()) {
        DTErrorMessage("ExtractColumns(IntArray,Range)","Range is out of bounds");
        return DTMutableIntArray();
    }
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(IntArray,Range)","Does not work for 3D arrays");
        return DTMutableIntArray();
    }
    
    DTMutableIntArray toReturn(A.m(),r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start*A.m(),(size_t)(r.length*A.m())*sizeof(int));
    return toReturn;
}

DTMutableIntArray ExtractIndices(const DTIntArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableIntArray();
        }
        else {
            DTErrorMessage("ExtractIndices(IntArray,IntArray)","Int array is empty");
            return DTMutableIntArray();
        }
    }
    
    ssize_t len = A.Length();
    ssize_t i,howMany = indices.Length();
    ssize_t index;
    DTMutableIntArray toReturn(howMany);
    
    bool outOfBounds = false;
    for (i=0;i<howMany;i++) {
        index = indices(i);
        
        if (index<0 || index>=len) {
            outOfBounds = true; // Set break point here
            toReturn(i) = INT32_MAX;
        }
        else {
            toReturn(i) = A(index);
        }
    }
    if (outOfBounds) {
        DTErrorMessage("ExtractIndices(IntArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

DTMutableIntArray ExtractIndices(const DTIntArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(IntArray,Range)","Range is out of bounds");
        return DTMutableIntArray();
    }
    
    DTMutableIntArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start, (size_t)r.length*sizeof(int));
    return toReturn;
}

int Minimum(const DTIntArray &A)
{
    ssize_t len = A.Length();
    int minV = INT32_MAX;
    
    int v;
    ssize_t i;
    
    const int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }
    
    return minV;
}

int Maximum(const DTIntArray &A)
{
    ssize_t len = A.Length();
    int maxV = INT32_MIN;
    
    int v;
    ssize_t i;
    
    const int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV);
    }
    
    return maxV;
}

void Range(const DTIntArray &A,int &minVRet,int &maxVRet)
{
    ssize_t len = A.Length();
    int minV = INT32_MAX;
    int maxV = INT32_MIN;
    
    int v;
    ssize_t i;
    
    const int *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
        maxV = (maxV < v ? v : maxV);
    }
    
    minVRet = minV;
    maxVRet = maxV;
}

