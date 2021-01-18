// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef _DTPointer_Header
#define _DTPointer_Header

// A very simple class that handles a reference counted pointer.
// The benefit of this is that this class will "own" the pointer and free
// the memory when the reference count hits 0.

// Used very much like a pointer.
// For example
//    DTPointer<DTAudio> theAudioFile(new DTAudio(....));
// allows you to pass an audio variable into and out of functions, even though
// the DTAudio file does not allow assignment.
//
// To call a member function, use
//    theAudioFile->NumberOfChannels();
// just as if it was a pointer.

#ifndef DTUseThreads
#define DTUseThreads 0
#endif

#if DTUseThreads
#include <pthread.h>
#else
#include <stdio.h>
#endif

template <class T>
class DTPointer {
public:
    // Functions
#if DTUseThreads
    DTPointer() : ref(new int(1)), mutexLock(new pthread_mutex_t()), Value(NULL) {pthread_mutex_init(mutexLock,NULL);} // Dangerous.
    DTPointer(T *Va) : ref(new int(1)), mutexLock(new pthread_mutex_t()), Value(Va) {pthread_mutex_init(mutexLock,NULL);}
    explicit DTPointer(T Va) : ref(new int(1)), mutexLock(new pthread_mutex_t()), Value(new T(Va)) {pthread_mutex_init(mutexLock,NULL);}
    DTPointer(const DTPointer<T> &ToC) : ref(NULL), mutexLock(NULL), Value(NULL) {
        pthread_mutex_lock(ToC.mutexLock);
        ref = ToC.ref;
        mutexLock = ToC.mutexLock; 
        Value = ToC.Value; 
        (*ref)++; 
        pthread_mutex_unlock(mutexLock);
    }
    virtual ~DTPointer() {
        pthread_mutex_lock(mutexLock);
        if (--(*ref)==0) {
            pthread_mutex_unlock(mutexLock);
            if (Value) delete Value;
            pthread_mutex_destroy(mutexLock);
            delete mutexLock;
            delete ref;
        }
        else {
            pthread_mutex_unlock(mutexLock);
        }
    }
#else
    DTPointer() : ref(new int(1)), Value(NULL) {} // Dangerous.
    DTPointer(T *Va) : ref(new int(1)), Value(Va) {}
    // explicit DTPointer(T Va) : ref(new int(1)), Value(new T(Va)) {}
    DTPointer(const DTPointer<T> &ToC) : ref(NULL), Value(NULL) {
        ref = ToC.ref;
        Value = ToC.Value; 
        (*ref)++; 
    }
    virtual ~DTPointer() {
        if (--(*ref)==0) {
            if (Value) delete Value;
            delete ref;
        }
    }
#endif
    
    operator bool() const {return (Value!=NULL);}

#if DTUseThreads
    DTPointer<T> &operator=(const DTPointer<T> &ToC) {
        if (ref!=ToC.ref) {
            pthread_mutex_lock(ToC.mutexLock);
            pthread_mutex_lock(mutexLock);
            if (--(*ref)==0) {
                pthread_mutex_unlock(mutexLock);
                if (Value) delete Value;
                pthread_mutex_destroy(mutexLock);
                delete ref;
            }
            else {
                pthread_mutex_unlock(mutexLock);
            }
            ref = ToC.ref;
            mutexLock = ToC.mutexLock;
            Value = ToC.Value;
            (*ref)++;
            pthread_mutex_unlock(mutexLock);
        }
        return *this;
    }
#else
    DTPointer<T> &operator=(const DTPointer<T> &ToC) {
        if (ref!=ToC.ref) {
            if (--(*ref)==0) {
                if (Value) delete Value;
                delete ref;
            }
            ref = ToC.ref;
            Value = ToC.Value;
            (*ref)++;
        }
        return *this;
    }
#endif
    
    T *operator->() const {return Value;}

    const T &operator*() const {return *Value;}

#if DTUseThreads
    int ReferenceCount(void) const {
        pthread_mutex_lock(mutexLock);
        int toReturn = *ref;
        pthread_mutex_unlock(mutexLock);
        return toReturn;
    }
#else
    int ReferenceCount(void) const {
        int toReturn = *ref;
        return toReturn;
    }
#endif

    const T *Data() const {return Value;}

protected:

    // Data
    int *ref;    // count how many use this structure.
#if DTUseThreads
    pthread_mutex_t *mutexLock;
#endif
    T *Value;
};

template <class T>
class DTMutablePointer : public DTPointer<T> {
public:

    DTMutablePointer() : DTPointer<T>() {}
    DTMutablePointer(T *Va) : DTPointer<T>(Va) {}
    // explicit DTMutablePointer(T Va) : DTPointer<T>(Va) {}
    DTMutablePointer(const DTMutablePointer<T> &A) : DTPointer<T>(A) {}

    DTMutablePointer<T> &operator=(const DTMutablePointer<T> &A) {DTPointer<T>::operator=(A); return *this;}
    
    T *operator->() {return DTPointer<T>::Value;}
    T *operator->() const {return DTPointer<T>::Value;}

    T &operator*() {return *DTPointer<T>::Value;}
    const T &operator*() const {return *DTPointer<T>::Value;}

    T *Data() {return DTPointer<T>::Value;}
    const T *Data() const {return DTPointer<T>::Value;}

};
    
#endif
