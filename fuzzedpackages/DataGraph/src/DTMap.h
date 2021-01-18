// Part of DTSource. Copyright 2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTMap_h
#define DTMap_h

#include "DTError.h"
#include "DTPointer.h"
#include "DTList.h"

// A simple key value container.  Key is a string, the value is a template class.  Use

// DTMutableMap<DTDataFile> dataFileByName;
// dataFileByName(pathName) = DTDataFile(pathName,DTFile::NewReadWrite);
// dataFileByName(pathName).Save(array,"name");

template <class T>
struct DTMapStorage {
    std::map<std::string,T> map;
    
    T notFound;
};

template <class T>
class DTMap {
public:
    DTMap() : storage(DTPointer<DTMapStorage<T> >(new DTMapStorage<T>())) {}
    
    T operator()(const std::string &s) const {
        typename std::map<std::string,T>::const_iterator it = storage->map.find(s);
        if (it==storage->map.end()) {
            DTErrorMessage("DTMutableMap(string)","key not found");
            return storage->notFound;
        }
        else {
            return it->second;
        }
    }
    
    size_t NumberOfKeys(void) const {return storage->map.size();}
    DTMutableList<std::string> Keys(void) const {
        typename std::map<std::string,T>::const_iterator it = storage->map.begin();
        typename std::map<std::string,T>::const_iterator theEnd = storage->map.end();
        int pos = 0;
        DTMutableList<std::string> toReturn(storage->map.size());
        for (it=storage->map.begin(); it!=theEnd; ++it)
            toReturn(pos++) = it->first;
        std::sort(toReturn.Pointer(),toReturn.Pointer()+toReturn.Length());
        return toReturn;
    }
    
    bool Contains(const std::string &s) const {
        typename std::map<std::string,T>::const_iterator it = storage->map.find(s);
        return (it!=storage->map.end());
    }
    
protected:
    DTPointer<DTMapStorage<T> > storage;
};

template <class T>
class DTMutableMap : public DTMap<T> {
public:
    DTMutableMap() : mutableStorage(DTMutablePointer<DTMapStorage<T> >(new DTMapStorage<T>())) {DTMap<T>::storage = mutableStorage;}
    
    T operator()(const std::string &s) const {
        typename std::map<std::string,T>::const_iterator it = mutableStorage->map.find(s);
        if (it==mutableStorage->map.end()) {
            DTErrorMessage("DTMutableMap(string)","key not found");
            return mutableStorage->notFound;
        }
        else {
            return it->second;
        }
    }
    T &operator()(const std::string &s) {
        return mutableStorage->map[s];
    }
    
    void Erase(const std::string &s) {mutableStorage->map.erase(s);}
    
private:
    DTMutablePointer<DTMapStorage<T> > mutableStorage;
};

#endif /* DTMap_h */
