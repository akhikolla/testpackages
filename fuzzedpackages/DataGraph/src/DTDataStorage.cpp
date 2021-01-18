// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTDataStorage.h"

#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTCharArray.h"
#include "DTUCharArray.h"
#include "DTFloatArray.h"
#include "DTShortIntArray.h"
#include "DTUShortIntArray.h"
#include "DTError.h"

#include <set>

std::string DTDataStorage::ResolveName(const std::string &name) const
{
    if (SavedAsString(name)==false)
        return name;
    
    std::string theName = ReadString(name);
    
    // First check if this is a one step redirect
    if (SavedAsString(theName)==false)
        return theName;
    
    // Deeper redirect, need to avoid circular references
    std::set<std::string> soFar;
    soFar.insert(name);
    while (SavedAsString(theName) && soFar.count(theName)==0) {
        soFar.insert(theName);
        theName = ReadString(theName);
    }
    if (soFar.count(theName)) {
        DTErrorMessage("DTDataStorage::ResolveName","Circular reference for "+name);
        return name;
    }
    else if (Contains(theName)==false)
        return name;
    else
        return theName;
}

void DTDataStorage::Flush(void) const
{
    
}

DTList<std::string> DTDataStorage::AllVariableNamesWithPrefix(const std::string &prefix) const
{
    DTList<std::string> allEntries = AllVariableNames();
    
    // Now go through the entries that start with the given prefix
    // This can clearly be made more efficient, and will be if people request it.
    // Right now, call the virtual function to get all of the entries, and pick from that.

    size_t howMany = allEntries.Length();
    size_t i;
    DTMutableIntArray whichAreIncluded(howMany);
    size_t pos = 0;
    size_t lenOfPrefix = prefix.length();

    for (i=0;i<howMany;i++) {
        if ((allEntries(i).compare(0,lenOfPrefix,prefix)==0))
            whichAreIncluded(pos++) = int(i);
    }
    
    DTMutableList<std::string> toReturn(pos);
    for (i=0;i<pos;i++) {
        toReturn(i) = allEntries(whichAreIncluded(i));
    }
    
    return toReturn;
}

void DTDataStorage::pinfo(void) const
{
    printInfo();
}

void DTDataStorage::printInfo(void) const
{
    // Overwrite to print content.
}

double DTDataStorage::ReadNumber(const std::string &name) const
{
    DTDoubleArray theArr = ReadDoubleArray(name);
    if (theArr.IsEmpty() || theArr.Length()!=1)
        return 0.0;

    return theArr(0);
}

int DTDataStorage::ReadInt(const std::string &name) const
{
    DTIntArray theArr = ReadIntArray(name);
    if (theArr.IsEmpty() || theArr.Length()!=1)
        return 0;

    return theArr(0);
}

void Read(const DTDataStorage &input,const std::string &name,double &toReturn)
{
    toReturn = input.ReadNumber(name);
}

void Read(const DTDataStorage &input,const std::string &name,int &toReturn)
{
    double temp;
    Read(input,name,temp);
    toReturn = int(temp);
}

void Write(DTDataStorage &output,const std::string &name,double theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,std::string &toReturn)
{
    toReturn = input.ReadString(name);
}

void Write(DTDataStorage &output,const std::string &name,const std::string &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTDoubleArray &toReturn)
{
    toReturn = input.ReadDoubleArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTDoubleArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTFloatArray &toReturn)
{
    toReturn = input.ReadFloatArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTFloatArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTIntArray &toReturn)
{
    toReturn = input.ReadIntArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTIntArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTCharArray &toReturn)
{
    toReturn = input.ReadCharArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTCharArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTUCharArray &toReturn)
{
    toReturn = input.ReadUCharArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTUCharArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTShortIntArray &toReturn)
{
    toReturn = input.ReadShortIntArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTShortIntArray &theVar)
{
    output.Save(theVar,name);
}

void Read(const DTDataStorage &input,const std::string &name,DTUShortIntArray &toReturn)
{
    toReturn = input.ReadUShortIntArray(name);
}

void Write(DTDataStorage &output,const std::string &name,const DTUShortIntArray &theVar)
{
    output.Save(theVar,name);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTDoubleArray &theVar)
{
    output.Save(theVar,name);
    if (theVar.n()>1)
        output.Save("Array","Seq_"+name);
    else
        output.Save("NumberList","Seq_"+name);
    output.Flush();
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTFloatArray &theVar)
{
    output.Save(theVar,name);
    if (theVar.n()>1)
        output.Save("Array","Seq_"+name);
    else
        output.Save("NumberList","Seq_"+name);
    output.Flush();
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTIntArray &theVar)
{
    output.Save(theVar,name);
    if (theVar.n()>1)
        output.Save("Array","Seq_"+name);
    else
        output.Save("NumberList","Seq_"+name);
    output.Flush();
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTCharArray &theVar)
{
    output.Save(theVar,name);
    if (theVar.n()>1)
        output.Save("Array","Seq_"+name);
    else
        output.Save("NumberList","Seq_"+name);
    output.Flush();
}

void WriteOne(DTDataStorage &output,const std::string &name,const std::string &theVar)
{
    output.Save(theVar,name);
    output.Save("String","Seq_"+name);
    output.Flush();
}

void WriteOne(DTDataStorage &output,const std::string &name,double theVar)
{
    output.Save(theVar,name);
    output.Save("Real Number","Seq_"+name);
    output.Flush();
}
