// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTDoubleArray.h"
#include "DTError.h"
#include "DTIntArray.h"
#include "DTList.h"
#include "DTUtilities.h"

#include "DTArrayTemplates.h"

#include <math.h>
#include <cstring>
#include <algorithm>
#include <limits>
#include <cmath>

#if !defined(INFINITY)
#if defined(WIN32)
#define INFINITY std::numeric_limits<double>::infinity();
#else
#define INFINITY 1e500
#endif
#endif

DTDoubleArrayStorage::DTDoubleArrayStorage(ssize_t mv,ssize_t nv,ssize_t ov)
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
    mutableReferences = 0;
    
    Data = length==0 ? NULL : new double[(size_t)length];
}

DTDoubleArrayStorage::~DTDoubleArrayStorage()
{
    delete [] Data;
}

DTDoubleArray::~DTDoubleArray()
{
    if (Storage) {
        Storage->accessLock.Lock();
        int refCnt = (--Storage->referenceCount);
        Storage->accessLock.Unlock();
        if (refCnt==0) delete Storage;
    }
}

DTDoubleArray::DTDoubleArray(const DTDoubleArray &A) 
{
    Storage = A.Storage;
    Storage->accessLock.Lock();
    Storage->referenceCount++;
    Storage->accessLock.Unlock();
}

DTDoubleArray &DTDoubleArray::operator=(const DTDoubleArray &A)
{
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

int DTDoubleArray::ReferenceCount() const
{
    Storage->accessLock.Lock(); 
    int refCnt = Storage->referenceCount;
    Storage->accessLock.Unlock();
    return refCnt;
}

int DTDoubleArray::MutableReferences() const
{
    Storage->accessLock.Lock(); 
    int toReturn = Storage->mutableReferences;
    Storage->accessLock.Unlock();
    return toReturn;
}

ssize_t DTDoubleArray::m() const
{
    return Storage->m;
}

ssize_t DTDoubleArray::n() const
{
    return Storage->n;
}

ssize_t DTDoubleArray::o() const
{
    return Storage->o;
}

ssize_t DTDoubleArray::Length() const
{
    return Storage->length;
}

bool DTDoubleArray::IsEmpty() const
{
    return (Storage->length==0);
}

bool DTDoubleArray::NotEmpty() const
{
    return (Storage->length!=0);
}

double DTDoubleArray::e(int i) const
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

double DTDoubleArray::e(int i,int j) const
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

double DTDoubleArray::e(int i,int j,int k) const
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

void DTDoubleArray::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    if (o()==0)
        std::cerr << "Empty" << std::endl;
    else if (o()==1) {
        if (n()==1) {
            if (m()<10) {
                std::cerr << "(";
                for (ssize_t i=0;i<m();i++) {
                    if (i>0) std::cerr << ",";
                    std::cerr << operator()(i);
                }
                std::cerr << ")" << std::endl;
            }
            else {
                std::cerr << "double vector with " << m() << " entries" << std::endl;
            }
        }
        else
            std::cerr << m() << " x " << n() << " double array" << std::endl;
    }
    else
        std::cerr << m() << " x " << n() << " x " << o() << " double array" << std::endl;
    
    std::cerr.flush();
#endif
}

void DTDoubleArray::pi(int i) const
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

void DTDoubleArray::pj(int j) const
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

void DTDoubleArray::pall(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t mv = m();
    ssize_t nv = n();
    size_t i,j;
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

void DTDoubleArray::pslice(int k) const
{
#ifndef DG_NOSTDErrOut
    ssize_t mv = m();
    ssize_t nv = n();
    size_t i,j;
    if (mv==0) {
        std::cerr << "Empty" << std::endl;
    }
    else if (k<0 || k>=o()) {
        std::cerr << "Out of bounds, o = " << o() << std::endl;
    }
    else {
        for (j=0;j<nv;j++) {
            for (i=0;i<mv-1;i++) std::cerr << operator()(i,j) << ", ";
            std::cerr << operator()(mv-1,j,k);
            std::cerr << std::endl;
        }
    }
#endif
}

void DTDoubleArray::prange(int startIndex,int endIndex) const
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

void DTDoubleArray::pcrange(int startIndex,int endIndex) const
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

void DTDoubleArray::psigns(void) const
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
		double val;
        for (j=0;j<nv;j++) {
            for (i=0;i<mv;i++) {
				val = operator()(i,j);
                if (std::isfinite(val)==0) {
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

ssize_t DTDoubleArray::Find(double v) const
{
    const double *D = Pointer();
    ssize_t len = Length();
    ssize_t i;
    for (i=0;i<len;i++) {
        if (D[i]==v) break;
    }

    return (i<len ? i : -1);
}

DTMutableDoubleArray DTDoubleArray::Copy() const
{
    DTMutableDoubleArray CopyInto(m(),n(),o());
    // Check that the allocation worked.
    if (CopyInto.Length()!=Length()) return CopyInto; // Failed.  Already printed an error message.
    std::memcpy(CopyInto.Pointer(),Pointer(),(size_t)(Length()*sizeof(double)));
    return CopyInto;
}

void DTDoubleArray::PrintErrorMessage(ssize_t i) const
{
    DTErrorOutOfRange("DTDoubleArray",i,Storage->length);
}

void DTDoubleArray::PrintErrorMessage(ssize_t i,ssize_t j) const
{
    DTErrorOutOfRange("DTDoubleArray",i,j,Storage->m,Storage->n);
}

void DTDoubleArray::PrintErrorMessage(ssize_t i,ssize_t j,ssize_t k) const
{
    DTErrorOutOfRange("DTDoubleArray",i,j,k,Storage->m,Storage->n,Storage->o);
}

DTMutableDoubleArray::DTMutableDoubleArray(const DTMutableDoubleArray &A)
: DTDoubleArray(A)
{
    Storage->accessLock.Lock();
    Storage->mutableReferences++;
    Storage->accessLock.Unlock();
}

DTMutableDoubleArray::~DTMutableDoubleArray()
{
    Storage->accessLock.Lock();
    int refCnt = (--Storage->referenceCount);
    Storage->mutableReferences--;
    Storage->accessLock.Unlock();
    if (refCnt==0) delete Storage;
    Storage = NULL;
}

DTMutableDoubleArray &DTMutableDoubleArray::operator=(const DTMutableDoubleArray &A)
{
    if (Storage!=A.Storage) {
        Storage->accessLock.Lock();
        A.Storage->accessLock.Lock();
        Storage->referenceCount--;
        int refCnt = Storage->referenceCount;
        Storage->mutableReferences--;
        Storage->accessLock.Unlock();
        if (refCnt==0) delete Storage;
        Storage = A.Storage;
        Storage->referenceCount++;
        Storage->mutableReferences++;
        Storage->accessLock.Unlock();
    }

    return *this;
}

DTMutableDoubleArray &DTMutableDoubleArray::operator=(double a)
{
    const size_t howManyNumbers = Length();
    if (a==0.0) {
        memset(Pointer(),0,sizeof(double)*howManyNumbers);
    }
    else {
//#ifndef _ANSI_SOURCE
//        memset_pattern8(Pointer(),&a,sizeof(double)*howManyNumbers);
//#else
        size_t i;
        double *Data = Pointer();
        for (i=0;i<howManyNumbers;i++)
            Data[i] = a;
//#endif
    }
    
    return *this;
}

void DTMutableDoubleArray::operator*=(double v)
{
    DTTimesEqualsScalar<DTMutableDoubleArray,double>(*this,v);
}

void DTMutableDoubleArray::operator/=(double v)
{
    DTDivideEqualsScalar<DTMutableDoubleArray,double>(*this,v);
}

void DTMutableDoubleArray::operator+=(double v)
{
    DTPlusEqualsScalar<DTMutableDoubleArray,double>(*this,v);
}

void DTMutableDoubleArray::operator-=(double v)
{
    DTMinusEqualsScalar<DTMutableDoubleArray,double>(*this,v);
}

void DTMutableDoubleArray::operator+=(const DTDoubleArray &B)
{
    DTPlusEqualsArray<DTDoubleArray,DTMutableDoubleArray,double>(*this,B);
}

void DTMutableDoubleArray::operator-=(const DTDoubleArray &B)
{
    DTMinusEqualsArray<DTDoubleArray,DTMutableDoubleArray,double>(*this,B);
}

void DTMutableDoubleArray::operator*=(const DTDoubleArray &B)
{
    DTTimesEqualsArray<DTDoubleArray,DTMutableDoubleArray,double>(*this,B);
}

void DTMutableDoubleArray::operator/=(const DTDoubleArray &B)
{
    DTDivideEqualsArray<DTDoubleArray,DTMutableDoubleArray,double>(*this,B);
}

bool operator==(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return DTOperatorArrayEqualsArray<DTDoubleArray,double>(A,B);
}

bool operator!=(const DTDoubleArray &A,const DTDoubleArray &B)
{
    return !(A==B);
}

DTMutableDoubleArray TruncateSize(const DTDoubleArray &A,ssize_t length)
{
    return DTTruncateArraySize<DTDoubleArray,DTMutableDoubleArray,double>(A,length);
}

DTMutableDoubleArray IncreaseSize(const DTDoubleArray &A)
{
    return IncreaseSize(A,A.Length());
}

DTMutableDoubleArray IncreaseSize(const DTDoubleArray &A,ssize_t addLength)
{
    return DTIncreaseArraySize<DTDoubleArray,DTMutableDoubleArray,double>(A,addLength);
}

DTMutableDoubleArray Transpose(const DTDoubleArray &A)
{
    return DTTransposeArray<DTDoubleArray,DTMutableDoubleArray,double>(A);
}

DTMutableDoubleArray Reshape(const DTDoubleArray &A,ssize_t m,ssize_t n,ssize_t o)
{
    if (m<0 || n<0 || o<0) {
        DTErrorMessage("Reshape(DTDoubleArray,...)","One of the new dimensions is negative.");
        return DTMutableDoubleArray();
    }
    if (m*n*o!=A.Length()) {
        DTErrorMessage("Reshape(DTDoubleArray,...)","Size before and after need to be the same.");
        return DTMutableDoubleArray();
    }
    
    DTMutableDoubleArray toReturn(m,n,o);
    if (toReturn.Length()) {
        std::memcpy(toReturn.Pointer(),A.Pointer(),A.Length()*sizeof(double));
    }
    
    return toReturn;
}

DTMutableDoubleArray Sort(const DTDoubleArray &A)
{
    DTMutableDoubleArray toReturn = Reshape(A,A.Length());
    std::sort(toReturn.Pointer(),toReturn.Pointer()+toReturn.Length());
    return toReturn;
}

struct DTArraySortingPair {
    bool operator<(const DTArraySortingPair &a) const {return (v<a.v);}
    double v;
	size_t i;
};

DTMutableIntArray SortedOrder(const DTDoubleArray &A)
{
	// A(returned(i)) = i'th entry in a sorted version of A.  So A(SortedOrder(A)) = Sort(A)
	size_t i,howMany = A.Length();
	DTMutableList<DTArraySortingPair> list(howMany);
	for (i=0;i<howMany;i++) {
		list(i).v = A(i);
		list(i).i = i;
	}
    std::sort(list.Pointer(),list.Pointer()+howMany);
	DTMutableIntArray toReturn(howMany);
	for (i=0;i<howMany;i++) {
		toReturn(i) = int(list(i).i);
	}
	return toReturn;
}

DTMutableDoubleArray FlipJ(const DTDoubleArray &A)
{
    return DTArrayFlipJ<DTDoubleArray,DTMutableDoubleArray,double>(A);
}

void Swap(DTMutableDoubleArray &A,DTMutableDoubleArray &B)
{
	DTMutableDoubleArray C = A;
	A = B;
	B = C;
}

void Swap(DTDoubleArray &A,DTDoubleArray &B)
{
	DTDoubleArray C = A;
	A = B;
	B = C;
}

void CopyValues(DTMutableDoubleArray &into,const DTDoubleArray &from)
{
	if (into.m()!=from.m() || into.n()!=from.n() || into.o()!=from.o()) {
		DTErrorMessage("CopyValues(MutableDoubleArray,DoubleArray)","Incompatible sizes");
	}
	else if (into.NotEmpty()) {
		std::memcpy(into.Pointer(),from.Pointer(),into.Length()*sizeof(double));
	}
}

void CopyIntoColumn(DTMutableDoubleArray &into,const DTDoubleArray &list,ssize_t j)
{
    if (into.m()!=list.Length()) {
		DTErrorMessage("CopyIntoColumn(into,list,j)","Need into.m()==list.Length()");
    }
    else if (into.o()!=1) {
		DTErrorMessage("CopyIntoColumn(into,list,j)","into is a 3D array (into.o()>1)");
    }
    else if (j<0 || j>into.n()) {
        DTErrorMessage("CopyIntoColumn(into,list,j)","j out of bounds");
    }
    else {
        std::memcpy(into.Pointer()+j*into.m(),list.Pointer(),into.m()*sizeof(double));
    }
}

void CopyIntoColumns(DTMutableDoubleArray &into,const DTRange &intoRange,const DTDoubleArray &from,const DTRange &fromRange)
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
        std::memcpy(into.Pointer()+intoRange.start*into.m(),from.Pointer()+fromRange.start*from.m(),(size_t)(from.m()*intoRange.length*sizeof(double)));
    }
}

void MemoryCopy(DTMutableDoubleArray &into,ssize_t intoLocation,const DTDoubleArray &from)
{
    MemoryCopy(into,intoLocation,from,DTRange(0,from.Length()));
}

void MemoryCopy(DTMutableDoubleArray &into,ssize_t intoLocation,const DTDoubleArray &from,const DTRange &range)
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
    if (range.length) std::memcpy(into.Pointer()+intoLocation,from.Pointer()+range.start,(size_t)(range.length*sizeof(double)));
}

void MemoryCopyColumns(DTMutableDoubleArray &into,ssize_t intoLocation,const DTDoubleArray &from,const DTRange &range)
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
    if (range.length) std::memcpy(into.Pointer()+into.m()*intoLocation,from.Pointer()+into.m()*range.start,(size_t)(range.length*into.m()*sizeof(double)));
}

void MemoryMove(DTMutableDoubleArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid destination range");
        return;
    }
    if (range.end()>into.Length()) {
        DTErrorMessage("MemoryMove","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+intoLocation,into.Pointer()+range.start,(size_t)(range.length*sizeof(double)));
}

void MemoryMoveColumns(DTMutableDoubleArray &into,ssize_t intoLocation,const DTRange &range)
{
    if (intoLocation<0 || ssize_t(intoLocation+range.length)>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid destination range");
        return;
    }
    if (range.end()>into.n()) {
        DTErrorMessage("MemoryMoveColumns","Invalid source range");
        return;
    }
    if (range.length) memmove(into.Pointer()+into.m()*intoLocation,into.Pointer()+into.m()*range.start,(size_t)(range.length*into.m()*sizeof(double)));
}

DTMutableDoubleArray SubArray(const DTDoubleArray &A,ssize_t imin,ssize_t icount,ssize_t jmin,ssize_t jcount,ssize_t kmin,ssize_t kcount)
{
	if (imin<0 || imin+icount>A.m() || jmin<0 || jmin+jcount>A.n() || kmin<0 || kmin+kcount>A.o()) {
		DTErrorMessage("SubArray(A,....)","Index out of range");
	}
	if (icount<=0 || jcount<=0 || kcount<=0) return DTMutableDoubleArray();
	DTMutableDoubleArray toReturn(icount,jcount,kcount);
	ssize_t j,k;
	double *into = toReturn.Pointer();
	const double *from = A.Pointer();
	const ssize_t m = A.m();
	const ssize_t n = A.n();
	for (k=0;k<kcount;k++) {
		for (j=0;j<jcount;j++) {
			if (icount==1)
				toReturn(0,j,k) = A(imin,jmin+j,kmin+k);
			else
				std::memcpy(into+j*icount+k*icount*jcount,from + imin + (jmin+j)*m+(kmin+k)*m*n,(size_t)(icount*sizeof(double)));
		}
	}
	
	return toReturn;
}

DTMutableDoubleArray Region(const DTDoubleArray &A,const DTRange &iRange,const DTRange &jRange,const DTRange &kRange)
{
    if (iRange.end()>A.m() || jRange.end()>A.n() || kRange.end()>A.o()) {
        DTErrorMessage("SubArray(A,....)","Index out of range");
    }
    if (iRange.length<=0 || jRange.length<=0 || kRange.length<=0) return DTMutableDoubleArray();
    ssize_t imin = iRange.start;
    ssize_t icount = iRange.length;
    ssize_t jmin = jRange.start;
    ssize_t jcount = jRange.length;
    ssize_t kmin = kRange.start;
    ssize_t kcount = kRange.length;
    DTMutableDoubleArray toReturn(icount,jcount,kcount);
    int j,k;
    double *into = toReturn.Pointer();
    const double *from = A.Pointer();
    const ssize_t m = A.m();
    const ssize_t n = A.n();
    for (k=0;k<kcount;k++) {
        for (j=0;j<jcount;j++) {
            if (icount==1)
                toReturn(0,j,k) = A(imin,jmin+j,kmin+k);
            else
                std::memcpy(into+j*icount+k*icount*jcount,from + imin + (jmin+j)*m+(kmin+k)*m*n,(size_t)(icount*sizeof(double)));
        }
    }
    
    return toReturn;
}

DTMutableDoubleArray Region(const DTDoubleArray &A,const DTRange &iRange,const DTRange &jRange)
{
    return Region(A,iRange,jRange,DTRange(0,1));
}

DTMutableDoubleArray ExtractColumns(const DTDoubleArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableDoubleArray();
        }
        else {
            DTErrorMessage("ExtractColumns(DoubleArray,IntArray)","Double array is empty");
            return DTMutableDoubleArray();
        }
    }
    
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(DoubleArray,IntArray)","Does not work for 3D arrays");
        return DTMutableDoubleArray();
    }
    
    ssize_t i,j,howMany = indices.Length();
    ssize_t index;
    ssize_t m = A.m();
    ssize_t n = A.n();
    
    DTMutableDoubleArray toReturn(A.m(),howMany);
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
        DTErrorMessage("ExtractColumns(DoubleArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

DTMutableDoubleArray ExtractRows(const DTDoubleArray &A,const DTRange &r)
{
    if (r.end()>A.m()) {
        DTErrorMessage("ExtractRows(DoubleArray,Range)","Range is out of bounds");
        return DTMutableDoubleArray();
    }
    if (A.o()>1) {
        DTErrorMessage("ExtractRows(DoubleArray,Range)","Does not work for 3D arrays");
        return DTMutableDoubleArray();
    }
    
    DTMutableDoubleArray toReturn(r.length,A.n());
    ssize_t m = A.m();
    ssize_t n = A.n();
    const double *fromD = A.Pointer();
    double *toD = toReturn.Pointer();
    ssize_t j;
    for (j=0;j<n;j++) {
        std::memcpy(toD+j*r.length,fromD+r.start+j*m,r.length*sizeof(double));
    }

    return toReturn;
}

DTMutableDoubleArray ExtractColumns(const DTDoubleArray &A,const DTRange &r)
{
    if (r.end()>A.n()) {
        DTErrorMessage("ExtractColumns(DoubleArray,Range)","Range is out of bounds");
        return DTMutableDoubleArray();
    }
    if (A.o()>1) {
        DTErrorMessage("ExtractColumns(DoubleArray,Range)","Does not work for 3D arrays");
        return DTMutableDoubleArray();
    }

    DTMutableDoubleArray toReturn(A.m(),r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start*A.m(), r.length*A.m()*sizeof(double));
    return toReturn;
}

DTMutableDoubleArray ExtractIndices(const DTDoubleArray &A,const DTIntArray &indices)
{
    if (A.IsEmpty()) {
        if (indices.IsEmpty()) {
            return DTMutableDoubleArray();
        }
        else {
            DTErrorMessage("ExtractIndices(DoubleArray,IntArray)","Double array is empty");
            return DTMutableDoubleArray();
        }
    }
    
    ssize_t len = A.Length();
    ssize_t i,howMany = indices.Length();
    ssize_t index;
    DTMutableDoubleArray toReturn(howMany);
    
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
        DTErrorMessage("ExtractIndices(DoubleArray,IntArray)","Index out of bounds");
    }
    
    return toReturn;
}

void AddToColumnRange(DTMutableDoubleArray &A,const DTRange &r,const DTDoubleArray &B)
{
    // A(range) += B
    if (A.m()!=B.m() || A.o()!=1 || B.o()!=1 || A.n()<r.end() || B.n()!=ssize_t(r.length)) {
        DTErrorMessage("AddToColumnRange(MutableDoubleArray,Range,DoubleArray)","Incompatible sizes");
        return;
    }

    size_t howMany = B.Length();
    size_t i;
    double *AD = A.Pointer() + A.m()*r.start;
    const double *BD = B.Pointer();
    for (i=0;i<howMany;i++) {
        AD[i] += BD[i];
    }
}

void AddToColumnRange(DTMutableDoubleArray &A,const DTRange &ar,const DTDoubleArray &B,const DTRange br,double scale)
{
    // A(range) += B(range)*scale
    if (A.m()!=B.m() || A.o()!=1 || B.o()!=1 || A.n()<ar.end() || B.n()<br.end() || ar.length!=br.length) {
        DTErrorMessage("AddToColumnRange(MutableDoubleArray,Range,DoubleArray,Range,double)","Incompatible sizes");
        return;
    }

    size_t howMany = A.m()*ar.length;
    size_t i;
    double *AD = A.Pointer() + A.m()*ar.start;
    const double *BD = B.Pointer() + B.m()*br.start;
    for (i=0;i<howMany;i++) {
        AD[i] += BD[i]*scale;
    }
}

DTMutableDoubleArray ExtractIndices(const DTDoubleArray &A,const DTRange &r)
{
    if (r.end()>A.Length()) {
        DTErrorMessage("ExtractIndices(DoubleArray,Range)","Range is out of bounds");
        return DTMutableDoubleArray();
    }
    
    DTMutableDoubleArray toReturn(r.length);
    std::memcpy(toReturn.Pointer(), A.Pointer()+r.start, r.length*sizeof(double));
    return toReturn;
}

double InfinityNorm(const DTDoubleArray &A)
{
    size_t len = A.Length();
    double maxV = 0;
    
    double v;
    size_t i;
    
    const double *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = fabs(D[i]);
        if (v>maxV) maxV = v;
    }
    
    return maxV;
}

double Minimum(const DTDoubleArray &A)
{
    size_t len = A.Length();
    double minV = INFINITY;
    
    double v;
    size_t i;
    
    const double *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        minV = (v < minV ? v : minV);
    }

    return minV;
}

double Minimum(const DTDoubleArray &A,ssize_t &index)
{
    size_t len = A.Length();
    double minV = INFINITY;
    index = -1;

    double v;
    size_t i;

    const double *D = A.Pointer();

    for (i=0;i<len;i++) {
        v = D[i];
        if (v < minV) {
            minV = v;
            index = i;
        }
    }

    return minV;
}

double Maximum(const DTDoubleArray &A)
{
    size_t len = A.Length();
    double maxV = -INFINITY;
    
    double v;
    size_t i;
    
    const double *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        maxV = (maxV < v ? v : maxV); // Will ignore NAN values
    }
    
    return maxV;
}

DTMutableDoubleArray Minimum(const DTDoubleArray &A,const DTDoubleArray &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("Minimum(DTDoubleArray,DTDoubleArray)","Need to have the same size");
        return DTMutableDoubleArray();
    }

    double va,vb;
    size_t i, len = A.Length();

    const double *AD = A.Pointer();
    const double *BD = B.Pointer();

    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    double *RD = toReturn.Pointer();

    for (i=0;i<len;i++) {
        va = AD[i];
        vb = BD[i];
        RD[i] = (va<vb ? va : vb);
    }

    return toReturn;
}

DTMutableDoubleArray Minimum(const DTDoubleArray &A,double vb)
{
    double va;
    size_t i, len = A.Length();
    
    const double *AD = A.Pointer();
    
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    double *RD = toReturn.Pointer();
    
    for (i=0;i<len;i++) {
        va = AD[i];
        RD[i] = (va<vb ? va : vb);
    }
    
    return toReturn;
}

DTMutableDoubleArray Maximum(const DTDoubleArray &A,const DTDoubleArray &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("Maximum(DTDoubleArray,DTDoubleArray)","Need to have the same size");
        return DTMutableDoubleArray();
    }

    double va,vb;
    size_t i, len = A.Length();

    const double *AD = A.Pointer();
    const double *BD = B.Pointer();

    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    double *RD = toReturn.Pointer();

    for (i=0;i<len;i++) {
        va = AD[i];
        vb = BD[i];
        RD[i] = (va>vb ? va : vb);
    }

    return toReturn;
}

bool ContainsNonFinite(const DTDoubleArray &A)
{
    const double *AD = A.Pointer();
    ssize_t i,howMany = A.Length();
    for (i=0;i<howMany;i++) {
        if (isfinite(AD[i])==0) break;
    }
    return (i<howMany);
}

DTMutableDoubleArray Maximum(const DTDoubleArray &A,double vb)
{
    double va;
    size_t i, len = A.Length();
    
    const double *AD = A.Pointer();
    
    DTMutableDoubleArray toReturn(A.m(),A.n(),A.o());
    double *RD = toReturn.Pointer();
    
    for (i=0;i<len;i++) {
        va = AD[i];
        RD[i] = (va>vb ? va : vb);
    }
    
    return toReturn;
}

ssize_t FindIndexOfMaximum(const DTDoubleArray &A)
{
    size_t len = A.Length();
    double maxV = -INFINITY;
    
    double v;
    size_t i;
    ssize_t toReturn = -1;
    
    const double *D = A.Pointer();
    
    for (i=0;i<len;i++) {
        v = D[i];
        if (maxV<v) {
            maxV = v;
            toReturn = i;
        }
    }
    
    return toReturn;
}

double Mean(const DTDoubleArray &A)
{
    size_t len = A.Length();
    
    double v;
    size_t i;
    
    const double *D = A.Pointer();
    double sum = 0;
    for (i=0;i<len;i++) {
        v = D[i];
        sum += v;
    }
    
    return sum/len;
}

DTMutableDoubleArray CombineColumns(const DTDoubleArray &First,const DTDoubleArray &Second)
{
    if (First.m()!=Second.m()) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to have the same number of rows.");
        return DTMutableDoubleArray();
    }
    if (First.IsEmpty())
        return DTMutableDoubleArray();
    if (First.o()!=1 || Second.o()!=1) {
        DTErrorMessage("CombineColumns(A,B)","A and B have to be two dimensional.");
        return DTMutableDoubleArray();
    }
    
    DTMutableDoubleArray toReturn(First.m(),First.n()+Second.n());
    std::memcpy(toReturn.Pointer(),First.Pointer(),First.Length()*sizeof(double));
    std::memcpy(toReturn.Pointer()+First.Length(),Second.Pointer(),Second.Length()*sizeof(double));
    
    return toReturn;
}

