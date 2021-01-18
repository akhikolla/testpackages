// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTListT_Header
#define DTListT_Header

#include "DTError.h"

// Simple list construct.  If you have a big data structure, consider using the STL library
// instead.  This array will not resize, will call the constructor for every element, and
// contains one additional entry - the outOfRange entry.
//
// The benefit is that memory is contiguous, random list access is order 1 and the list
// can be returned from a function, since the storage is reference counted.  It also has the
// same const/mutable convention that the Array class does.

// To allocate N entries
//   DTMutableList<TypeName> list(N);
//
// List access
//   list(3) = list(2);

#ifndef DTRangeCheck
#define DTRangeCheck 1
#endif

template <class T>
class DTList {
public:
    DTList() : Data(0), length(0), refCount(new int(1)), outOfRange() {}
    DTList(const DTList<T> &A) : Data(A.Data), length(A.length), refCount(A.refCount), outOfRange() {++(*refCount);}
protected:
    DTList(ssize_t len) : Data(len<=0 ? 0 : new T[(size_t)len]), length(len<=0 ? 0 : len), refCount(new int(1)), outOfRange() {}
public:
    
    virtual ~DTList() {
        --(*refCount);
        if (*refCount==0) {delete [] Data; delete refCount;}
        Data = 0; refCount = 0; length=0;
    }

    DTList<T> &operator=(const DTList<T> &A) {
        if (A.refCount!=refCount) { // Otherwise doing A=A.
            --(*refCount);
            if (*refCount==0) {delete [] Data; delete refCount;}
            refCount = A.refCount;
            ++(*refCount);
            length = A.length;
            Data = A.Data;
        }
        return *this;
    }

    ssize_t MemoryUsed(void) const {return length*sizeof(T);}

    const T *Pointer(void) const {return Data;}

    bool IsEmpty(void) const {return (Data==0);}
    ssize_t Length(void) const {return length;}

#if DTRangeCheck
    const T operator()(ssize_t i) const  {if (i<0 || i>=length) {DTErrorOutOfRange("DTList<T>",i,length); return outOfRange;} return Data[i];}
#else
    const T operator()(ssize_t i) const  {return Data[i];}
#endif

protected:
    T *Data;
    ssize_t length;
    int *refCount;

    // Should be static.
    T outOfRange;
};

template <class T>
class DTMutableList : public DTList<T> {
public:
    DTMutableList() : DTList<T>() {}
    DTMutableList(ssize_t len) : DTList<T>(len) {}
    DTMutableList(const DTMutableList<T> &A) : DTList<T>(A) {}

    DTMutableList<T> &operator=(const DTMutableList<T> &A) {DTList<T>::operator=(A); return *this;}

    T *Pointer(void) {return DTList<T>::Data;}
    const T *Pointer(void) const {return DTList<T>::Data;}

#if DTRangeCheck
    T &operator()(ssize_t i) {if (i<0 || i>=DTList<T>::length) {DTErrorOutOfRange("DTList<T>",i,DTList<T>::length); return DTList<T>::outOfRange;} return DTList<T>::Data[i];}
    T operator()(ssize_t i) const  {if (i<0 || i>= DTList<T>::length) {DTErrorOutOfRange("DTList<T>",i,DTList<T>::length); return DTList<T>::outOfRange;} return DTList<T>::Data[i];}
#else
    T &operator()(ssize_t i) {return DTList<T>::Data[i];}
    T operator()(ssize_t i) const  {return DTList<T>::Data[i];}
#endif

    DTMutableList<T> &operator=(T v) {for (ssize_t i=0;i<DTList<T>::length;i++) DTList<T>::Data[i] = v; return *this;}
};

template <class T> DTMutableList<T> TruncateSize(const DTList<T> &A,ssize_t length)
{
    if (length>A.Length()) length = A.Length();
    DTMutableList<T> toReturn(length);
    const T *fromP = A.Pointer();
    T *toP = toReturn.Pointer();
    for (ssize_t i=0;i<length;i++) toP[i] = fromP[i];
    return toReturn;
}

template <class T> DTMutableList<T> IncreaseSize(const DTList<T> &A,ssize_t addLength)
{
    DTMutableList<T> toReturn(A.Length()+(addLength>=0 ? addLength : 0));
    ssize_t len = A.Length();
    const T *fromP = A.Pointer();
    T *toP = toReturn.Pointer();
    for (ssize_t i=0;i<len;i++) toP[i] = fromP[i];
    return toReturn;
}

template <class T> DTMutableList<T> Copy(const DTList<T> &A)
{
    DTMutableList<T> toReturn(A.Length());
    ssize_t len = A.Length();
    const T *fromP = A.Pointer();
    T *toP = toReturn.Pointer();
    for (ssize_t i=0;i<len;i++) toP[i] = fromP[i];
    return toReturn;
}

#endif
