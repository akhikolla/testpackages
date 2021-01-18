// Part of DTSource. Copyright 2012. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTDictionary.h"
#include "DTError.h"
#include <math.h>
#include "DTDataStorage.h"
#include "DTIntArray.h"
#include "DTUtilities.h"

 
DTDictionary::DTDictionary(void)
{
	content = DTPointer<DTDictionaryStorage>(new DTDictionaryStorage());
}

DTDictionaryAccess DTDictionary::operator()(const std::string &name) const
{
    return DTDictionaryAccess(*this,name);
}

double DTDictionary::GetNumber(const std::string &s) const
{
	const DTDictionaryStorage *p = content.Data();
    std::map<std::string,double>::const_iterator it = p->numberDictionary.find(s);
	if (it==p->numberDictionary.end()) {
		DTErrorMessage("dictionary.Number(string)","key \"" + s + "\" not found, returning NAN");
		return NAN;
	}
	else {
		return it->second;
	}
}

double DTDictionary::GetNumber(const std::string &s,double ifNotDefined) const // If not in the dictionary, quietly return the second argument.
{
	const DTDictionaryStorage *p = content.Data();
	std::map<std::string,double>::const_iterator it = p->numberDictionary.find(s);
	if (it==p->numberDictionary.end()) {
        return ifNotDefined;
	}
	else {
		return it->second;
	}
}

DTDoubleArray DTDictionary::GetArray(const std::string &s) const
{
	const DTDictionaryStorage *p = content.Data();
	std::map<std::string,DTDoubleArray>::const_iterator it = p->arrayDictionary.find(s);
	if (it==p->arrayDictionary.end()) {
		DTErrorMessage("dictionary.GetArray(string)","key not found, returning an empty array.");
		return DTDoubleArray();
	}
	else {
		return it->second;
	}
}

std::string DTDictionary::GetString(const std::string &s) const
{
	const DTDictionaryStorage *p = content.Data();
	std::map<std::string,std::string>::const_iterator it = p->stringDictionary.find(s);
	if (it==p->stringDictionary.end()) {
		DTErrorMessage("dictionary.GetString(string)","key not found, returning an empty string.");
		return "";
	}
	else {
		return it->second;
	}
}

DTDictionary DTDictionary::GetDictionary(const std::string &s) const
{
	const DTDictionaryStorage *p = content.Data();
	std::map<std::string,DTDictionary>::const_iterator it = p->dictionaryDictionary.find(s);
	if (it==p->dictionaryDictionary.end()) {
		DTErrorMessage("dictionary.GetDictionary(string)","key not found, returning an empty dictionary.");
		return DTDictionary();
	}
	else {
		return it->second;
	}
}

DTDictionary::ValueType DTDictionary::TypeOf(const std::string &s) const
{
    const DTDictionaryStorage *p = content.Data();
    if (p->numberDictionary.count(s)) return DTDictionary::Number;
    if (p->arrayDictionary.count(s)) return DTDictionary::Array;
    if (p->stringDictionary.count(s)) return DTDictionary::String;
    if (p->dictionaryDictionary.count(s)) return DTDictionary::Dictionary;
    return DTDictionary::NotFound;
}

bool DTDictionary::Contains(const std::string &s) const
{
    return (TypeOf(s)!=DTDictionary::NotFound);
}

size_t DTDictionary::NumberOfKeys(void) const
{
    const DTDictionaryStorage *p = content.Data();
    return (p->numberDictionary.size() + p->arrayDictionary.size() + p->stringDictionary.size() + p->dictionaryDictionary.size());
}

DTList<std::string> DTDictionary::AllNumberKeys(void) const
{
    const DTDictionaryStorage *p = content.Data();
    size_t howMany = p->numberDictionary.size();
    DTMutableList<std::string> toReturn(howMany);
    std::map<std::string,double>::const_iterator dit = p->numberDictionary.begin();
    size_t pos;
    for (pos=0;pos<howMany;pos++) {
        toReturn(pos) = dit->first;
        dit++;
    }
    return toReturn;
}

DTMutableDoubleArray DTDictionary::NumbersForKeys(const DTList<std::string> &keys) const
{
    std::map<std::string,double> &numbers = content->numberDictionary;

    size_t howMany = keys.Length();
    DTMutableDoubleArray toReturn(howMany);

    std::map<std::string,double>::const_iterator it;
    std::map<std::string,double>::const_iterator theEnd = numbers.end();
    size_t i;
    for (i=0;i<howMany;i++) {
        it = numbers.find(keys(i));
        if (it==theEnd) {
            DTErrorMessage("dictionary.NumbersForKeys(list)","key not found, returning NAN");
            toReturn(i) = NAN;
        }
        else {
            toReturn(i) = it->second;
        }
    }
    return toReturn;
}

DTMutableDictionary DTDictionary::Copy(void) const
{
    DTMutableDictionary toReturn;
    toReturn+=*this;
    return toReturn;
}

void DTDictionary::pinfoWithPrefix(std::string prefix) const
{
#ifndef DG_NOSTDErrOut
    int howMany = 0;
    const DTDictionaryStorage *from = content.Data();
    std::map<std::string,double>::const_iterator nit = from->numberDictionary.begin();
    std::map<std::string,double>::const_iterator nitend = from->numberDictionary.end();
    while (nit!=nitend) {
        std::cerr << prefix << nit->first << " = " << nit->second << std::endl;
        nit++;
        howMany++;
    }

    std::map<std::string,std::string>::const_iterator sit = from->stringDictionary.begin();
    std::map<std::string,std::string>::const_iterator sitend = from->stringDictionary.end();
    while (sit!=sitend) {
        std::cerr << prefix << sit->first << " = " << sit->second << std::endl;
        sit++;
        howMany++;
    }

    std::map<std::string,DTDoubleArray>::const_iterator ait = from->arrayDictionary.begin();
    std::map<std::string,DTDoubleArray>::const_iterator aitend = from->arrayDictionary.end();
    while (ait!=aitend) {
        std::cerr << prefix << ait->first << " = ";
        ait->second.pinfo();
        ait++;
        howMany++;
    }

    std::map<std::string,DTDictionary>::const_iterator dit = from->dictionaryDictionary.begin();
    std::map<std::string,DTDictionary>::const_iterator ditend = from->dictionaryDictionary.end();
    while (dit!=ditend) {
        std::cerr << prefix << dit->first << " = \n";
        dit->second.pinfoWithPrefix(prefix+"   ");
        dit++;
        howMany++;
    }
    
    if (howMany==0)
        std::cerr << prefix << "empty\n";
#endif
}

void DTDictionary::pinfo(void) const
{
    pinfoWithPrefix("");
}

DTMutableDictionary::DTMutableDictionary()
{
	mutableContent = DTMutablePointer<DTDictionaryStorage>(new DTDictionaryStorage());
    content = mutableContent;
}

DTMutableDictionary::DTMutableDictionary(DTList<std::string> &keys,DTDoubleArray &numbers)
{
	mutableContent = DTMutablePointer<DTDictionaryStorage>(new DTDictionaryStorage());
    content = mutableContent;

    if (keys.Length()!=numbers.Length()) {
        DTErrorMessage("DTMutableDoubleArray(keys,DTDoubleArray)","Lengths have to match");
        return;
    }
    std::map<std::string,double> &numberDict = mutableContent->numberDictionary;
    size_t i,howMany = keys.Length();
    for (i=0;i<howMany;i++) {
        numberDict[keys(i)] = numbers(i);
    }
}

void DTMutableDictionary::Remove(const std::string &s)
{
	DTDictionaryStorage *p = mutableContent.Data();
	p->numberDictionary.erase(p->numberDictionary.find(s));
	p->arrayDictionary.erase(p->arrayDictionary.find(s));
}

void DTMutableDictionary::operator+=(const DTDictionary &from)
{
    // Add both the numbers and arrays
    std::map<std::string,double> &addToNumbers = mutableContent->numberDictionary;
    std::map<std::string,double>::const_iterator nit = from.content->numberDictionary.begin();
    std::map<std::string,double>::const_iterator nitend = from.content->numberDictionary.end();
    while (nit!=nitend) {
        addToNumbers[nit->first] = nit->second;
        nit++;
    }
    
    // Arrays
    std::map<std::string,DTDoubleArray> &addToArrays = mutableContent->arrayDictionary;
    std::map<std::string,DTDoubleArray>::const_iterator ait = from.content->arrayDictionary.begin();
    std::map<std::string,DTDoubleArray>::const_iterator aitend = from.content->arrayDictionary.end();
    while (ait!=aitend) {
        addToArrays[ait->first] = ait->second;
        ait++;
    }
    
    // strings
    std::map<std::string,std::string> &addToStrings = mutableContent->stringDictionary;
    std::map<std::string,std::string>::const_iterator sit = from.content->stringDictionary.begin();
    std::map<std::string,std::string>::const_iterator sitend = from.content->stringDictionary.end();
    while (sit!=sitend) {
        addToStrings[sit->first] = sit->second;
        sit++;
    }
    
    // dictionaries
    std::map<std::string,DTDictionary> &addToDictionaries = mutableContent->dictionaryDictionary;
    std::map<std::string,DTDictionary>::const_iterator dit = from.content->dictionaryDictionary.begin();
    std::map<std::string,DTDictionary>::const_iterator ditend = from.content->dictionaryDictionary.end();
    while (dit!=ditend) {
        addToDictionaries[dit->first] = dit->second;
        dit++;
    }
}

void DTMutableDictionary::AddNumbers(const DTDoubleArray &numbers,const DTList<std::string> &keys)
{
    if (keys.Length()!=numbers.Length()) {
        DTErrorMessage("DTMutableDictionary::AddNumbers(DTDoubleArray,keys)","Lengths have to match");
        return;
    }
    std::map<std::string,double> &numberDict = mutableContent->numberDictionary;
    size_t i,howMany = keys.Length();
    for (i=0;i<howMany;i++) {
        numberDict[keys(i)] = numbers(i);
    }
}

DTMutableDictionaryAssignment DTMutableDictionary::operator()(std::string name)
{
    return DTMutableDictionaryAssignment(*this,name);
}

DTDictionaryAccess DTMutableDictionary::operator()(const std::string &name) const
{
    return DTDictionaryAccess(*this,name);
}

void DTMutableDictionaryAssignment::operator=(const DTDictionaryAccess &v)
{
    switch (v.dict.TypeOf(v.name)) {
        case DTDictionary::Number:
            dict.SetNumber(v.dict.GetNumber(v.name),name);
            break;
        case DTDictionary::Array:
            dict.SetArray(v.dict.GetArray(v.name),name);
            break;
        case DTDictionary::String:
            dict.SetString(v.dict.GetString(v.name),name);
            break;
        case DTDictionary::Dictionary:
            dict.SetDictionary(v.dict.GetDictionary(v.name),name);
            break;
        case DTDictionary::NotFound:
            DTErrorMessage("Dictionary(name) = dictionary(name)","Name not found");
            break;
    }
}

DTMutableDictionary operator+(const DTDictionary &A,const DTDictionary &B)
{
    DTMutableDictionary toReturn = A.Copy();
    toReturn+=B;
    return toReturn;
}
                              
void Read(const DTDataStorage &input,const std::string &name,DTDictionary &toReturn)
{
    DTIntArray lengths = input.ReadIntArray(name);
    if (lengths.Length()!=4) {
        DTErrorMessage("Read(DTDataStorage,name,DTDictionary","Invalid format.");
        toReturn = DTDictionary();
        return;
    }
    
    int i,howMany;
    DTMutableDictionary readInto;
    
    howMany = lengths(0);
    for (i=0;i<howMany;i++) readInto(input.ReadString(name+"_n"+DTInt2String(i)+"_name")) = input.ReadNumber(name+"_n"+DTInt2String(i)+"_value");
    
    howMany = lengths(1);
    for (i=0;i<howMany;i++) readInto(input.ReadString(name+"_s"+DTInt2String(i)+"_name")) = input.ReadString(name+"_s"+DTInt2String(i)+"_value");
    
    howMany = lengths(2);
    for (i=0;i<howMany;i++) readInto(input.ReadString(name+"_a"+DTInt2String(i)+"_name")) = input.ReadDoubleArray(name+"_a"+DTInt2String(i)+"_value");
    
    howMany = lengths(3);
    DTDictionary temp;
    for (i=0;i<howMany;i++) {
        Read(input,name+"_d"+DTInt2String(i)+"_value",temp);
        readInto(input.ReadString(name+"_d"+DTInt2String(i)+"_name")) = temp;
    }
    
    toReturn = readInto;
}

void Write(DTDataStorage &output,const std::string &name,const DTDictionary &theDictionary)
{
    // Structure is
    // name - An array.  Four numbers, how many numbers, strings, arrays and dictionaries there are
    // name_n#_name - name of scalar
    // name_n#_value - value for scalar
    // name_s#_name - name of string
    // name_s#_value
    // name_a#_name - name of array
    // name_a#_value
    // name_d#_name - name of dictionary
    // name_d#_value - The dictionary, another DTDictionary variable.
    
    const DTDictionaryStorage &storage = *theDictionary.content;
    
    DTMutableIntArray lengths(4);
    
    map<string,double>::const_iterator nit = storage.numberDictionary.begin();
    map<string,double>::const_iterator nitend = storage.numberDictionary.end();
    int count = 0;
    while (nit!=nitend) {
        Write(output,name+"_n"+DTInt2String(count)+"_name",nit->first);
        Write(output,name+"_n"+DTInt2String(count)+"_value",nit->second);
        nit++;
        count++;
    }
    lengths(0) = count;
    
    map<string,std::string>::const_iterator sit = storage.stringDictionary.begin();
    map<string,std::string>::const_iterator sitend = storage.stringDictionary.end();
    count = 0;
    while (sit!=sitend) {
        Write(output,name+"_s"+DTInt2String(count)+"_name",sit->first);
        Write(output,name+"_s"+DTInt2String(count)+"_value",sit->second);
        sit++;
        count++;
    }
    lengths(1) = count;
    
    map<string,DTDoubleArray>::const_iterator ait = storage.arrayDictionary.begin();
    map<string,DTDoubleArray>::const_iterator aitend = storage.arrayDictionary.end();
    count = 0;
    while (ait!=aitend) {
        Write(output,name+"_a"+DTInt2String(count)+"_name",ait->first);
        Write(output,name+"_a"+DTInt2String(count)+"_value",ait->second);
        ait++;
        count++;
    }
    lengths(2) = count;
    
    map<string,DTDictionary>::const_iterator dit = storage.dictionaryDictionary.begin();
    map<string,DTDictionary>::const_iterator ditend = storage.dictionaryDictionary.end();
    count = 0;
    while (dit!=ditend) {
        Write(output,name+"_d"+DTInt2String(count)+"_name",dit->first);
        Write(output,name+"_d"+DTInt2String(count)+"_value",dit->second);
        dit++;
        count++;
    }
    lengths(3) = count;
    
    Write(output,name,lengths);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTDictionary &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"Dictionary");
    output.Flush();
}


