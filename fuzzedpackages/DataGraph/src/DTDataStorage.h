// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTDataStorage_Header
#define DTDataStorage_Header

/*
 Base class for DTDataContainer (memory based storage), DTDataFile (file based) and DTMatlabDataFile (file based).
 */

#include "DTPointer.h"
#include "DTList.h"

class DTDoubleArray;
class DTFloatArray;
class DTIntArray;
class DTCharArray;
class DTUCharArray;
class DTShortIntArray;
class DTUShortIntArray;

#include <string>

class DTDataStorage {
public:
    virtual ~DTDataStorage() {}

    virtual DTMutablePointer<DTDataStorage> AsPointer() const = 0;
    
    DTList<std::string> AllVariableNamesWithPrefix(const std::string &) const;
    virtual DTList<std::string> AllVariableNames(void) const = 0;
    virtual bool Contains(const std::string &name) const = 0;
    virtual bool IsReadOnly(void) const = 0;

    // Saving data.
    virtual void Save(int,const std::string &name) = 0;
    virtual void Save(double,const std::string &name) = 0;
    virtual void Save(const DTDoubleArray &A,const std::string &name) = 0;
    virtual void Save(const DTFloatArray &A,const std::string &name) = 0;
    virtual void Save(const DTIntArray &A,const std::string &name) = 0;
    virtual void Save(const DTCharArray &A,const std::string &name) = 0;
    virtual void Save(const DTUCharArray &A,const std::string &name) = 0;
    virtual void Save(const DTShortIntArray &A,const std::string &name) = 0;
    virtual void Save(const DTUShortIntArray &A,const std::string &name) = 0;
    virtual void Save(const std::string &theString,const std::string &name) = 0;

    virtual bool SavedAsCharacter(const std::string &name) const = 0; // Signed or unsigned
    virtual bool SavedAsShort(const std::string &name) const = 0; // Signed or unsigned
    virtual bool SavedAsInt(const std::string &name) const = 0;
    virtual bool SavedAsFloat(const std::string &name) const = 0;
    virtual bool SavedAsDouble(const std::string &name) const = 0;
    virtual bool SavedAsString(const std::string &name) const = 0;

    virtual void Flush(void) const;

    // Reading data.
    virtual DTDoubleArray ReadDoubleArray(const std::string &name) const = 0;
    virtual DTFloatArray ReadFloatArray(const std::string &name) const = 0;
    virtual DTIntArray ReadIntArray(const std::string &name) const = 0;
    virtual DTCharArray ReadCharArray(const std::string &name) const = 0;
    virtual DTUCharArray ReadUCharArray(const std::string &name) const = 0;
    virtual DTShortIntArray ReadShortIntArray(const std::string &name) const = 0;
    virtual DTUShortIntArray ReadUShortIntArray(const std::string &name) const = 0;
    virtual double ReadNumber(const std::string &name) const;
    virtual int ReadInt(const std::string &name) const;
    virtual std::string ReadString(const std::string &name) const = 0;

    std::string ResolveName(const std::string &name) const;

    // Debug
    void pinfo() const;

protected:
    virtual void printInfo(void) const;
};

extern void Read(const DTDataStorage &input,const std::string &name,double &toReturn);
extern void Read(const DTDataStorage &input,const std::string &name,int &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,double theVar);
extern void Read(const DTDataStorage &input,const std::string &name,std::string &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const std::string &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTDoubleArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTDoubleArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTFloatArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTFloatArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTIntArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTIntArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTCharArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTCharArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTUCharArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTUCharArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTShortIntArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTShortIntArray &theVar);
extern void Read(const DTDataStorage &input,const std::string &name,DTUShortIntArray &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTUShortIntArray &theVar);

// Each variable type has a WriteOne(...) function, which will save the variable to a data file 
// along with a type description.  These functions are similar.  However, this will save
// the variable as a "List Of Numbers" if you hand it a vector, but an "Array" if the second dimension is >1.
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTDoubleArray &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTFloatArray &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTIntArray &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTCharArray &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const std::string &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,double theVar);

#endif
