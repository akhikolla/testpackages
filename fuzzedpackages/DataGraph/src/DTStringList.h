// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTStringList_H
#define DTStringList_H

#include "DTCharArray.h"
#include "DTIntArray.h"
#include "DTDataStorage.h"
#include "DTList.h"

struct DTRange;

class DTStringList {
public:
    explicit DTStringList() : characters(), offsets() {}
    DTStringList(const DTCharArray &); // Compute the offsets
    DTStringList(const DTCharArray &,const DTIntArray &);
    explicit DTStringList(const DTList<std::string> &entries);
    
    size_t NumberOfStrings(void) const {return offsets.Length();}
    std::string operator()(ssize_t) const;
    
    DTCharArray Characters(void) const {return characters;}
    DTIntArray Offsets(void) const {return offsets;}
    
    void pinfo(void) const;
    void pall(void) const;
    
private:
    DTCharArray characters;
    DTIntArray offsets;
};

extern DTStringList ExtractIndices(const DTStringList &A,const DTRange &);

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTStringList &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTStringList &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTStringList &toWrite); // One time value, self documenting.

#endif
