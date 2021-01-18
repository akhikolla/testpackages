// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

// Template functions to implement array operators.

#include "DTError.h"
#include <unistd.h>
#include <cstring>

template <class T,class TM,class Td>
TM DTAddArrays(const char *name,const T &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage(name,"Incompatible sizes.");
        return TM();
    }

    TM toReturn(A.m(),A.n(),A.o());
    ssize_t len = A.Length();
    const Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    ssize_t i;
    for (i=0;i<len;i++)
        retP[i] = AP[i]+BP[i];

    return toReturn;
}

template <class T,class TM,class Td>
TM DTSubtractArrays(const char *name,const T &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage(name,"Incompatible sizes.");
        return TM();
    }
    
    TM toReturn(A.m(),A.n(),A.o());
    ssize_t len = A.Length();
    const Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    ssize_t i;
    for (i=0;i<len;i++)
        retP[i] = AP[i]-BP[i];

    return toReturn;
}

template <class T,class TM,class Td>
TM DTMultiplyArrays(const char *name,const T &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage(name,"Incompatible sizes.");
        return TM();
    }

    TM toReturn(A.m(),A.n(),A.o());
    ssize_t len = A.Length();
    const Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    ssize_t i;
    for (i=0;i<len;i++)
        retP[i] = AP[i]*BP[i];

    return toReturn;
}

template <class T,class TM,class Td>
TM DTDivideArrays(const char *name,const T &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage(name,"Incompatible sizes.");
        return TM();
    }

    TM toReturn(A.m(),A.n(),A.o());
    ssize_t len = A.Length();
    const Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    ssize_t i;
    for (i=0;i<len;i++)
        retP[i] = AP[i]/BP[i];

    return toReturn;
}

template <class T,class TM,class Td>
TM DTArrayPlusNumber(const T &A,Td b)
{
    TM toReturn(A.m(),A.n(),A.o());
    ssize_t i,len = A.Length();
    const Td *AP = A.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = AP[i]+b;
    return toReturn;
}

template <class T,class TM,class Td>
TM DTArrayTimesNumber(const T &A,Td b)
{
    TM toReturn(A.m(),A.n(),A.o());
    ssize_t i,len = A.Length();
    const Td *AP = A.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = AP[i]*b;
    return toReturn;
}

template <class T,class TM,class Td>
TM DTArrayDivideByNumber(const T &A,Td b)
{
    TM toReturn(A.m(),A.n(),A.o());
    ssize_t i,len = A.Length();
    const Td *AP = A.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = AP[i]/b;
    return toReturn;
}

template <class T,class TM,class Td>
TM DTNumberMinusArray(Td a,const T &B)
{
    TM toReturn(B.m(),B.n(),B.o());
    ssize_t i,len = B.Length();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = a-BP[i];
    return toReturn;
}

template <class T,class TM,class Td>
TM DTNumberDividedByArray(Td a,const T &B)
{
    TM toReturn(B.m(),B.n(),B.o());
    ssize_t i,len = B.Length();
    const Td *BP = B.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = a/BP[i];
    return toReturn;
}

template <class T,class TM,class Td>
TM DTNegateArray(const T &A)
{
    TM toReturn(A.m(),A.n(),A.o());
    ssize_t i,len = A.Length();
    const Td *AP = A.Pointer();
    Td *retP = toReturn.Pointer();
    for (i=0;i<len;i++) retP[i] = -AP[i];
    return toReturn;
}

template <class T,class TM,class Td>
void DTPlusEqualsArray(TM &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("A+=B","Incompatible sizes.");
        return;
    }
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    for (i=0;i<len;i++) AP[i] += BP[i];
}

template <class T,class TM,class Td>
void DTMinusEqualsArray(TM &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("A-=B","Incompatible sizes.");
        return;
    }
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    for (i=0;i<len;i++) AP[i] -= BP[i];
}

template <class T,class TM,class Td>
void DTTimesEqualsArray(TM &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("A*=B","Incompatible sizes.");
        return;
    }
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    for (i=0;i<len;i++) AP[i] *= BP[i];
}

template <class T,class TM,class Td>
void DTDivideEqualsArray(TM &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o()) {
        DTErrorMessage("A/=B","Incompatible sizes.");
        return;
    }
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    const Td *BP = B.Pointer();
    for (i=0;i<len;i++) AP[i] /= BP[i];
}

template <class TM,class Td>
void DTPlusEqualsScalar(TM &A,Td b)
{
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    for (i=0;i<len;i++) AP[i] += b;
}

template <class TM,class Td>
void DTMinusEqualsScalar(TM &A,Td b)
{
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    for (i=0;i<len;i++) AP[i] -= b;
}

template <class TM,class Td>
void DTTimesEqualsScalar(TM &A,Td b)
{
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    for (i=0;i<len;i++) AP[i] *= b;
}

template <class TM,class Td>
void DTDivideEqualsScalar(TM &A,Td b)
{
    ssize_t i,len = A.Length();
    Td *AP = A.Pointer();
    for (i=0;i<len;i++) AP[i] /= b;
}

template <class T,class TM,class Td>
TM DTTruncateArraySize(const T &A,ssize_t length)
{
    // New length needs to fit as a MxNxO array
    // where MNO = length and
    // if o>1, length = m*n*k
    // if o=1 and n>1 length = m*k
    // if o=1 and n=1 everything is ok.

    if (length==0) return TM();
    if (A.IsEmpty()) {
        DTErrorMessage("TruncateSize(Array,Length)","Array is empty.");
        return TM();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (length%(A.m()*A.n())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return TM();
        }
        newM = A.m();
        newN = A.n();
        newO = length/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (length%(A.m())!=0) {
            DTErrorMessage("TruncateSize(Array,Length)","Invalid new dimension");
            return TM();
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

    TM toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)length*sizeof(Td));
    return toReturn;
}

template <class T,class TM,class Td>
TM DTIncreaseArraySize(const T &A,ssize_t addLength)
{
    if (addLength<0) {
        DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be >0.");
        return TM();
    }

    ssize_t newM,newN,newO;
    if (A.o()>1) {
        if (addLength%(A.m()*A.n())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m*n");
            return TM();
        }
        newM = A.m();
        newN = A.n();
        newO = A.o() + addLength/(A.m()*A.n());
    }
    else if (A.n()>1) {
        if (addLength%(A.m())!=0) {
            DTErrorMessage("IncreaseSize(Array,Length)","Length needs to be a multiple of m");
            return TM();
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

    TM toReturn(newM,newN,newO);
    std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)A.Length()*sizeof(Td));
    return toReturn;
}

template <class T,class td>
bool DTOperatorArrayEqualsArray(const T &A,const T &B)
{
    if (A.m()!=B.m() || A.n()!=B.n() || A.o()!=B.o())
        return false;
    if (A.Pointer()==B.Pointer())
        return true;
    
    ssize_t len = A.Length();
    const td *AP = A.Pointer();
    const td *BP = B.Pointer();
    
    return (memcmp(AP,BP,(size_t)len*sizeof(td))==0);
}

template <class T,class TM,class Td>
TM DTTransposeArray(const T &A)
{
    if (A.IsEmpty()) return TM();

    const ssize_t m = A.m();
    const ssize_t n = A.n();
    const ssize_t o = A.o();
    ssize_t i,j,k;

    TM toReturn;
    Td *toReturnD;
    const Td *AD = A.Pointer();

    if (A.o()!=1) {
        toReturn = TM(o,n,m);
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
        toReturn = TM(n,m);
        toReturnD = toReturn.Pointer();
        ssize_t ijNew, ijOld;
        if (m==1 || n==1) {
            std::memcpy(toReturn.Pointer(),A.Pointer(),(size_t)(m*n)*sizeof(Td));
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

template <class T,class TM,class Td>
TM DTArrayFlipJ(const T &A)
{
    TM toReturn(A.m(),A.n(),A.o());
    
    ssize_t m = A.m();
    ssize_t n = A.n();
    ssize_t o = A.o();
    ssize_t mn = m*n;
    
    const Td *fromP = A.Pointer();
    Td *toP = toReturn.Pointer();
    
    ssize_t j,k;
    for (k=0;k<o;k++) {
        for (j=0;j<n;j++) {
            std::memcpy(toP+j*m+k*mn,fromP+(n-1-j)*m+k*mn,(size_t)m*sizeof(Td));
        }
    }
    
    return toReturn;
}

