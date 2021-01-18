// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef _DTLock_Header
#define _DTLock_Header

// Assignment and copy constructors are not thread safe for this class.

#ifndef DTUseThreads
#define DTUseThreads 0
#endif

#if DTUseThreads
#include <pthread.h>
#endif

class DTLock {
public:
#if DTUseThreads
    DTLock() : mutexLock(NULL) {mutexLock = new pthread_mutex_t(); pthread_mutex_init(mutexLock,NULL);}
    ~DTLock() {pthread_mutex_destroy(mutexLock); delete mutexLock;}
    
    bool operator==(const DTLock &L) const {return (mutexLock==L.mutexLock);}

    void Lock(void) const {pthread_mutex_lock(mutexLock);}
    bool TryLock(void) const {return (pthread_mutex_trylock(mutexLock)==0);}
    void Unlock(void) const {pthread_mutex_unlock(mutexLock);}
#else
    DTLock() {}
    ~DTLock() {}
    
    bool operator==(const DTLock &) const {return false;}

    void Lock(void) const {}
    bool TryLock(void) const {return true;}
    void Unlock(void) const {}
#endif

private:
    // Can not copy this object.
    DTLock(const DTLock &ToC);
    DTLock &operator=(const DTLock &ToC);
    
#if DTUseThreads
    pthread_mutex_t *mutexLock;
#endif
};

#endif

