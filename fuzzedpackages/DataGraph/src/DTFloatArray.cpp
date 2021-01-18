// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTFloatArray.h"
#include "DTError.h"
#include "DTArrayTemplates.h"
#include "DTIntArray.h"
#include "DTUtilities.h"

#include <math.h>
#include <cstring>
#include <algorithm>
#include <cmath>

DTFloatArrayStorage::DTFloatArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
{
    // Check if it's called correctly.
    if (mv<0 || nv<0 || ov<0) DTErrorMessage("DTMutableFloatArray", "Negative index in constructor");

    m = mv>0 ? mv : 0;
    n = nv>0 ? nv : 0;
    o = ov>0 ? ov : 0;
    length = m*n*o;
    if (length==0) m = n = o = 0;
    referenceCount = 1;
    mn = m*n;

    Data = length==0 ? NULL : new float[(size_t)length];
}

DTFloatArrayStorage::~DTFloatArrayStorage()
{
    delete [] Data;
}

DTFloatArray &DTFloatArray::operator=(const DTFloatArray &A)
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

DTMutableFloatArray DTFloatArray::Copy() const
{
    DTMutableFloatArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)Length()*sizeof(float));
    return CopyInto;
}

float DTFloatArray::e(int i) const
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

float DTFloatArray::e(int i,int j) const
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

float DTFloatArray::e(int i,int j,int k) const
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

void DTFloatArray::pinfo(void) const
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

void DTFloatArray::pi(int i) const
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

void DTFloatArray::pj(int j) const
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

void DTFloatArray::pall(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t mv = m();
    ssize_t nv = n();
    ssize_t ov = o();
    if (mv*nv*ov>1000) {
        std::cerr << "More than 1000 numbers, save it to a file instead." << std::endl;
    }
    else if (mv==0) {
        std::cerr << "Empty" << std::endl;
    }
    else {
        ssize_t i,j,k;
        if (ov==1) {
            for (j=0;j<nv;j++) {
                for (i=0;i<mv-1;i++) std::cerr << operator()(i,j) << ", ";
                std::cerr << operator()(mv-1,j);
                std::cerr << std::endl;
            }
        }
        else {
            for (k=0;k<ov;k++) {
                std::cerr << "k = " << k << ":" << std::endl;
                for (j=0;j<nv;j++) {
                    std::cerr << "  ";
                    for (i=0;i<mv-1;i++) std::cerr << operator()(i,j,k) << ", ";
                    std::cerr << operator()(mv-1,j,k);
                    std::cerr << std::endl;
                }
            }
        }
    }
#endif
}

void DTFloatArray::prange(int startIndex,ssize_t endIndex) const
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

void DTFloatArray::pcrange(int startIndex,int endIndex) const
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

void DTFloatArray::psigns(void) const
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
		float val;
        for (j=0;j<nv;j++) {
            for (i=0;i<mv;i++) {
				val = operator()(i,j);
				if (isfinite(val)==0) {
                    if (std::isnan(val))
						std::cerr << "*";
					else if (val<0)
						std::cerr << "P";
					else if (val>0)
						std::cerr << "M";
				}
				else if (val<0)
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

ssize_t DTFloatArray::Find(float v) const
{
    const float *D = Pointer();
    ssize_t len = Length();
    ssize_t i;
    for (i=0;i<len;i++) {
        if (D[i]==v) break;
    }
    
    return (i<len ? i : -1);
}

DTMutableFloatArray TruncateSize(const DTFloatArray &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return DTMutableFloatArray();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return DTMutableFloatArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableFloatArray();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return DTMutableFloatArray();
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

    DTMutableFloatArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(float));
    return toReturn;
}

DTMutableFloatArray IncreaseSize(const DTFloatArray &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return DTMutableFloatArray();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return DTMutableFloatArray();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return DTMutableFloatArray();
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

    DTMutableFloatArray toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(float));
    return toReturn;
}

void DTFloatArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTFloatArray",i,Storage->length);
}

void DTFloatArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTFloatArray",i,j,Storage->m,Storage->n);
}

void DTFloatArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTFloatArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableFloatArray &DTMutableFloatArray::operator=(float a)
{
    const ssize_t howManyNumbers = Length();
    ssize_t i;
    float *Data = Pointer();
    for (i=0;i<howManyNumbers;i++)
        Data[i] = a;
    
    return *this;
}

void DTMutableFloatArray::operator*=(float v)
{
    DTTimesEqualsScalar<DTMutableFloatArray,float>(*this,v);
}

void DTMutableFloatArray::operator/=(float v)
{
    DTDivideEqualsScalar<DTMutableFloatArray,float>(*this,v);
}

void DTMutableFloatArray::operator+=(float v)
{
    DTPlusEqualsScalar<DTMutableFloatArray,float>(*this,v);
}

void DTMutableFloatArray::operator-=(float v)
{
    DTMinusEqualsScalar<DTMutableFloatArray,float>(*this,v);
}

void DTMutableFloatArray::operator+=(const DTFloatArray &B)
{
    DTPlusEqualsArray<DTFloatArray,DTMutableFloatArray,float>(*this,B);
}

void DTMutableFloatArray::operator-=(const DTFloatArray &B)
{
    DTMinusEqualsArray<DTFloatArray,DTMutableFloatArray,float>(*this,B);
}

void DTMutableFloatArray::operator*=(const DTFloatArray &B)
{
    DTTimesEqualsArray<DTFloatArray,DTMutableFloatArray,float>(*this,B);
}

void DTMutableFloatArray::operator/=(const DTFloatArray &B)
{
    DTDivideEqualsArray<DTFloatArray,DTMutableFloatArray,float>(*this,B);
}

bool operator==(const DTFloatArray &A,const DTFloatArray &B)
{
    return DTOperatorArrayEqualsArray<DTFloatArray,float>(A,B);
}

bool operator!=(const DTFloatArray &A,const DTFloatArray &B)
{
    return !(A==B);
}

DTMutableFloatArray Transpose(const DTFloatArray &A)
{
    if (A.IsEmpty()) return DTMutableFloatArray();

    const ssize_t m = A.m();
    const ssize_t n = A.n();
    const ssize_t o = A.o();
    ssize_t i,j,k;

    DTMutableFloatArray toReturn;
    float *toReturnD;
    const float *AD = A.Pointer();

    if (A.o()!=1) {
        toReturn = DTMutableFloatArray(o,n,m);
        toReturnD = toReturn.Pointer();
        ssize_t ijkNew,ijkOld;
        ssize_t no = n*o;
        for (k=0;k<o;k++) {
            for (j=0;j<n;j++) {
                ijkNew = k + j*o;
                ijkOld = j*m + k*m*n;
                for (i=0;i<m;i++) {
                    toReturnD[ijkNew] = AD[ijkOld]; // toReturn(k,j,i) = A(i,j,k)
                    ijkNew += no;
                    ijkOld++;
                }
            }
        }
    }
    else {
        toReturn = DTMutableFloatArray(n,m);
        toReturnD = toReturn.Pointer();
        ssize_t ijNew, ijOld;
        if (m==1 || n==1) {
            std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)(m*n)*sizeof(float));
        }
        else {
            for (j=0;j<n;j++) {
                ijNew = j;
                ijOld = j*m;
                for (i=0;i<m;i++) {
                    toReturnD[ijNew] = AD[ijOld]; // toReturn(j,i) = A(i,j)
                    ijNew += n;
                    ijOld++;
                }
            }
        }
    }

    return toReturn;
}

DTMutableFloatArray Reshape(const DTFloatArray &A,ssize_t m,ssize_t n,ssize_t o)
{
    if (m<0 || n<0 || o<0) {
        DTErrorMessage("Reshape(DTFloatArray,...)","One of the new dimensions is negative.");
        return DTMutableFloatArray();
    }
    if (m*n*o!=A.Length()) {
        DTErrorMessage("Reshape(DTFloatArray,...)","Size before and after need to be the same.");
        return DTMutableFloatArray();
    }

    DTMutableFloatArray toReturn(m,n,o);
    if (toReturn.Length()) {
        std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(float));
    }

    return toReturn;
}

DTMutableFloatArray Sort(const DTFloatArray &A)
{
    DTMutableFloatArray toReturn = Reshape(A,A.Length());;
    std::sort(toReturn.Pointer(),toReturn.Pointer()+toReturn.Length());
    return toReturn;
}

DTMutableFloatArray FlipJ(const DTFloatArray &A)
{
    return DTArrayFlipJ<DTFloatArray,DTMutableFloatArray,float>(A);
}

void Swap(DTMutableFloatArray &A,DTMutableFloatArray &B)
{
	DTMutableFloatArray C = A;
	A = B;
	B = C;
}

void Swap(DTFloatArray &A,DTFloatArray &B)
{
	DTFloatArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableFloatArray &into,const DTFloatArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableFloatArray,FloatArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),(size_t)into.Length()*sizeof(float));
	}
}

void CopyIntoColumns(DTMutableFloatArray &into,const DTRange &intoRange,const DTFloatArray &from,const DTRange &fromRange)
{
    if (into.o()!=1 || from.o()!=1) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","into is a 3D array (into.o()>1 or from.o()>1)");
    }
    else if (into.n()<ssize_t(intoRange.start+intoRange.length) || from.n()<ssize_t(fromRange.start+fromRange.length)) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","Out of bounds");
    }
    else if (intoRange.length!=fromRange.length) {
		DTErrorMessage("CopyIntoColumns(into,range,from,range)","fromRange.length!=toRange.length");
    }
    else if (into.m()!=from.m()) {
        DTErrorMessage("CopyIntoColumns(into,range,from,range)","from.m()!=to.m()");
    }
    else {
        std::memcpy(into.Pointer()+intoRange.start*into.m(),from.Pointer()+fromRange.start*from.m(),(size_t)(from.m()*intoRange.length)*sizeof(float));
    }
}

void MemoryCopy(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from)
{
    MemoryCopy(into,intoLocation,from,DTRange(0,from.Length()));
}

void MemoryCopy(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from,const DTRange &range)
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
    if (range.length) std::memcpy(into.Pointer()+intoLocation,from.Pointer()+range.start,(size_t)range.length*sizeof(float));
}

void MemoryCopyColumns(DTMutableFloatArray &into,ssize_t intoLocation,const DTFloatArray &from,const DTRange &range)
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
    if (range.length) std::memcpy(into.Pointer()+into.m()*intoLocation,from.Pointer()+into.m()*range.start,(size_t)(range.length*into.m())*sizeof(float));
}

void MemoryMove(DTMutableFloatArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid destination range");
        return;
    }
    if (range.end()>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+intoLocation,into.Pointer()+range.start,(size_t)range.length*sizeof(float));
}

void MemoryMoveColumns(DTMutableFloatArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid destination range");
        return;
    }
    if (range.end()>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+into.m()*intoLocation,into.Pointer()+into.m()*range.start,(size_t)(range.length*into.m())*sizeof(float));
}

DTMutableFloatArray Region(const DTFloatArray &A,const DTRange &iRange,const DTRange &jRange,const DTRange &kRange)
{
	if (iRange.end()>A.m() || jRange.end()>A.n() || kRange.end()>A.o()) {
		DTErrorMessage("SubArray(A,....)","Index out of range");
	}
	if (iRange.length<=0 || jRange.length<=0 || kRange.length<=0) return DTMutableFloatArray();
    ssize_t imin = iRange.start;
    ssize_t icount = iRange.length;
    ssize_t jmin = jRange.start;
    ssize_t jcount = jRange.length;
    ssize_t kmin = kRange.start;
    ssize_t kcount = kRange.length;
    DTMutableFloatArray toReturn(icount,jcount,kcount);
	int j,k;
	float *into = toReturn.Pointer();
	const float *from = A.Pointer();
	const ssize_t m = A.m();
	const ssize_t n = A.n();
	for (k=0;k<kcount;k++) {
		for (j=0;j<jcount;j++) {
			if (icount==1)
				toReturn(0,j,k) = A(imin,jmin+j,kmin+k);
			else
				std::memcpy(into+j*icount+k*icount*jcount,from + imin + (jmin+j)*m+(kmin+k)*m*n, (size_t)icount*sizeof(float));
		}
	}
	
	return toReturn;
}

DTMutableFloatArray Region(const DTFloatArray &A,const DTRange &iRange,const DTRange &jRange)
{
    return Region(A,iRange,jRange,DTRange(0,1));
}

DTMutableFloatArray ExtractIndices(const DTFloatArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableFloatArray();
        }
        else {
            DTErrorMessage("ExtractIndices(FloatArray,IntArray)","Float array is empty");
            return DTMutableFloatArray();
        }
    }
    
    ssize_t len = A.Length();
    ssize_t i,howMany = indices.Length();
    int index;
    DTMutableFloatArray toReturn(howMany);
    
    bool outOfBounds = false;
    for (i=0;i<howMany;i++) {
        index = indices(i);
        
        if (index<0 || index>=len) {
            outOfBounds = true; // Set break point here
            toReturn(i) = NAN;
        }
        else {
            toReturn(i) = A(index);
        }
    }
    if (outOfBounds) {
        DTErrorMessage("ExtractIndices(FloatArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

DTMutableFloatArray ExtractIndices(const DTFloatArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(FloatArray,Range)","Range is out of bounds");
        return DTMutableFloatArray();
    }
    
    DTMutableFloatArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start,(size_t)r.length*sizeof(float));
    return toReturn;
}

DTMutableFloatArray ExtractColumns(const DTFloatArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableFloatArray();
        }
        else {
            DTErrorMessage("ExtractColumns(FloatArray,IntArray)","Float array is empty");
            return DTMutableFloatArray();
        }
    }
    
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(FloatArray,IntArray)","Does not work for 3D arrays");
        return DTMutableFloatArray();
    }
    
    ssize_t i,j,howMany = indices.Length();
    int index;
    ssize_t m = A.m();
    ssize_t n = A.n();
    
    DTMutableFloatArray toReturn(A.m(),howMany);
    bool outOfBounds = false;
    for (j=0;j<howMany;j++) {
        index = indices(j);
        if (index<0 || index>=n) {
            outOfBounds = true; // Set break point here
            for (i=0;i<m;i++) toReturn(i,j) = NAN;
        }
        else {
            for (i=0;i<m;i++) toReturn(i,j) = A(i,index);
        }
    }
    if (outOfBounds) {
        DTErrorMessage("ExtractColumns(FloatArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

DTMutableFloatArray ExtractColumns(const DTFloatArray &A,const DTRange &r)
{
    if (r.end()>A.n()) {
        DTErrorMessage("ExtractColumns(FloatArray,Range)","Range is out of bounds");
        return DTMutableFloatArray();
    }
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(FloatArray,Range)","Does not work for 3D arrays");
        return DTMutableFloatArray();
    }
    
    DTMutableFloatArray toReturn(A.m(),r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start*A.m(), (size_t)(r.length*A.m())*sizeof(float));
    return toReturn;
}

float Minimum(const DTFloatArray &A)
{
    ssize_t len = A.Length();
    float minV = INFINITY;
    
    float v;
    ssize_t i;
    
    const float *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }
	
    return minV;
}

float Maximum(const DTFloatArray &A)
{
    ssize_t len = A.Length();
    float maxV = -INFINITY;
    
    float v;
    ssize_t i;
    
    const float *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV);
    }
    
    return maxV;
}

float Mean(const DTFloatArray &A)
{
    ssize_t len = A.Length();
    
    float v;
    ssize_t i;
    
    const float *D = A.Pointer();
    float sum = 0;
    for (i=0;i<len;i++) {
        v = D[i];
        sum += v;
    }
    
    return sum/len;
}

DTMutableFloatArray Minimum(const DTFloatArray &A,const DTFloatArray &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("Minimum(DTFloatArray,DTFloatArray)","Need to have the same size");
        return DTMutableFloatArray();
    }

    float va,vb;
    ssize_t i, len = A.Length();

    const float *AD = A.Pointer();
    const float *BD = B.Pointer();

    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    float *RD = toReturn.Pointer();

    for (i=0;i<len;i++) {
        va = AD[i];
        vb = BD[i];
        RD[i] = (va<vb ? va : vb);
    }

    return toReturn;
}

DTMutableFloatArray Maximum(const DTFloatArray &A,const DTFloatArray &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("Maximum(DTFloatArray,DTFloatArray)","Need to have the same size");
        return DTMutableFloatArray();
    }

    float va,vb;
    ssize_t i, len = A.Length();

    const float *AD = A.Pointer();
    const float *BD = B.Pointer();

    DTMutableFloatArray toReturn(A.m(),A.n(),A.o());
    float *RD = toReturn.Pointer();

    for (i=0;i<len;i++) {
        va = AD[i];
        vb = BD[i];
        RD[i] = (va>vb ? va : vb);
    }
    
    return toReturn;
}

DTMutableFloatArray CombineColumns(const DTFloatArray &First,const DTFloatArray &Second)
{
    return CombineColumns(First,Second,Second.n());
}

DTMutableFloatArray CombineColumns(const DTFloatArray &First,const DTFloatArray &Second,ssize_t fromSecond)
{
    if (First.m()!=Second.m()) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to have the same number of rows.");
        return DTMutableFloatArray();
    }
    if (First.IsEmpty())
        return DTMutableFloatArray();
    if (First.o()!=1 || Second.o()!=1) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to be two dimensional.");
        return DTMutableFloatArray();
    }
    if (fromSecond>Second.n()) {
        DTErrorMessage("CombineColumns(A,B,fromSecond)","Too many columns specified.");
        return DTMutableFloatArray();
    }
    
    DTMutableFloatArray toReturn(First.m(),First.n()+fromSecond);
    std::memcpy(toReturn.Pointer(),First.Pointer(),(size_t)First.Length()*sizeof(float));
    std::memcpy(toReturn.Pointer()+First.Length(),Second.Pointer(),(size_t)(Second.m()*fromSecond)*sizeof(float));
    
    return toReturn;
}
