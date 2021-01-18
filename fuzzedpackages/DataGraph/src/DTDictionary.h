// Part of DTSource. Copyright 2012. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTDictionary_H
#define DTDictionary_H

// This is a method to group multiple variables into a single object.
// C and C++ don't allow named arguments like some other programming languages.
// This means that when you have a long list of arguments of the same type it can
// get fairly brittle.  For example if you swap arguments accidentally.  And variable number of
// arguments and default values is also tricky.  The dictionary class is intended to simplify
// this.  It is based on the STL map class, but makes it easy to pass in and out, save and load
// from file etc.

// Some examples of how to use it:
//   DTMutableDictionary input;
//   //..populate input
//   input("a") = 5.03;
//   input("b") = someArray; // of type DTDoubleArray
//   input("c") = "string"
//   input("d") = anotherDictionary("e"); // Keeps the same type
//   input("more") = dict; // dict is a DTDictionary
//   DTDictionary output = foo(input); // All inputs handled with a dictionary, and outputs are put into a dictionary.
// where foo extracts the input as
//      double a = input("a");
//      DTDoubleArray b = input("b");
//      double optional = input.GetNumber("optional",2.3); // 2.3 is the value if "optional" is not defined.
//      DTDictionary dict = input("more");
// If the key doesn't exist or doesn't have the same type, this is a run-time error and you should set breakpoints in DTError.cpp
// or DTDictionary.cpp to catch this case.

// DTDictionary is an immutable (read only) dictionary.  Since DTMutableDictionary is derived
// from DTDictionary you can hand in a DTMutableDictionary to a function that expects a DTDictionary but
// not vice versa.  The dictionary is a smart pointer so that copies don't copy the members.

#include "DTPointer.h"
#include "DTDoubleArray.h"
#include "DTList.h"
#include <string>
#include <map>

class DTMutableDictionary;
class DTDataStorage;
class DTDictionary;

struct DTDictionaryStorage {
	std::map<std::string,double> numberDictionary;
	std::map<std::string,DTDoubleArray> arrayDictionary;
	std::map<std::string,std::string> stringDictionary;
	std::map<std::string,DTDictionary> dictionaryDictionary;
};

class DTMutableDictionaryAssignment;
class DTDictionaryAccess;

class DTDictionary {
public:
	DTDictionary();


    // Low level access, get one number or array by name.  Will display an error message if not found.
    enum ValueType {NotFound=0,Number,Array,String,Dictionary};
    DTDictionary::ValueType TypeOf(const std::string &) const;
    
    bool Contains(const std::string &) const;

    // value = dict("a") works for both types.  Run time error if the variable isn't defined or has the wrong type.
    DTDictionaryAccess operator()(const std::string &) const;
    // Dictionary can be used to pass multiple variables into a function.  If the argument has to exist use value = dict(name)
    // but if it doesn't have to exist you can supply a default value to use if the entry is not defined.
    // For example if eps is an optional threshold you can use
    // double epsilon = inputArguments.GetNumber("epsilon",1e-9);
    
    double GetNumber(const std::string &,double ifNotDefined) const; // If not in the dictionary, quietly return the second argument.
    
    // List based access
    bool NotEmpty(void) const {return (NumberOfKeys()!=0);}
    size_t NumberOfKeys(void) const; // both for numbers and arrays.
    DTList<std::string> AllNumberKeys(void) const;
    DTMutableDoubleArray NumbersForKeys(const DTList<std::string> &) const;
    
    void pinfo(void) const;
    
    DTMutableDictionary Copy(void) const; // If you want a copy.

    // Functions that are used by the accessor class.  double v = dict(name) is cleaner than double v = dict.GetNumber(name);
    // Instead of calling these, use the handy syntax
    // value = dictionary("key") and let the compiler automatically figure out which type you want based on the value type
    double GetNumber(const std::string &) const;
    DTDoubleArray GetArray(const std::string &) const;
    std::string GetString(const std::string &) const;
    DTDictionary GetDictionary(const std::string &) const;
    
protected:
    void pinfoWithPrefix(std::string) const;
    
	DTPointer<DTDictionaryStorage> content;
    friend class DTMutableDictionary;
    friend void Write(DTDataStorage &output,const std::string &name,const DTDictionary &);
};


class DTMutableDictionary : public DTDictionary {
public:
	DTMutableDictionary();
    DTMutableDictionary(DTList<std::string> &keys,DTDoubleArray &numbers);

	// Set values by using the syntax dict.Number("a") = 3.24;
    DTMutableDictionaryAssignment operator()(std::string name);
    DTDictionaryAccess operator()(const std::string &) const;
    
    // Remove entries from a dictionary
	void Remove(const std::string &);
    
    // Add multiple entries in one call
    void operator+=(const DTDictionary &);
    void AddNumbers(const DTDoubleArray &numbers,const DTList<std::string> &keys);

    // Accessor entries, used by the assignment class below.
    void SetNumber(double v,const std::string &s) {mutableContent->numberDictionary[s] = v;}
    void SetString(std::string v,const std::string &s) {mutableContent->stringDictionary[s] = v;}
    void SetArray(const DTDoubleArray &v,const std::string &s) {mutableContent->arrayDictionary[s] = v;}
    void SetDictionary(const DTDictionary &v,const std::string &s) {mutableContent->dictionaryDictionary[s] = v;}
    
private:
	DTMutablePointer<DTDictionaryStorage> mutableContent;
};

// Combine two dictionaries. The second dictionary overwrites any entries that are also defined in the first.
extern DTMutableDictionary operator+(const DTDictionary &,const DTDictionary &);

extern void Read(const DTDataStorage &input,const std::string &name,DTDictionary &);
extern void Write(DTDataStorage &output,const std::string &name,const DTDictionary &);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTDictionary &);

////////////////////////////////////////////////////////////////////////////////////
// Below this are classes that are used to make the assignment syntax work.
// Don't need to know them since you don't create them directly.
// Since C++ only overloads based on the argument need to create these classes.
// they are "invisible" in that you don't assign them to local variables,
// the compiler will just create them as temporary variables.
////////////////////////////////////////////////////////////////////////////////////
class DTMutableDictionaryAssignment {
public:
    DTMutableDictionaryAssignment(const DTMutableDictionary &d,const std::string &s) : dict(d), name(s) {}
    void operator=(double v) {dict.SetNumber(v,name);}
    void operator=(std::string v) {dict.SetString(v,name);}
    void operator=(const DTDoubleArray &v) {dict.SetArray(v,name);}
    void operator=(const DTDictionary &v) {dict.SetDictionary(v.Copy(),name);}
    void operator=(const DTDictionaryAccess &);
    operator double(void) {return dict.GetNumber(name);}
    operator DTDoubleArray(void) {return dict.GetArray(name);}
    operator std::string(void) {return dict.GetString(name);}
    operator DTDictionary(void) {return dict.GetDictionary(name);}
    
private:
    DTMutableDictionary dict;
    std::string name;
};

class DTDictionaryAccess {
public:
    DTDictionaryAccess(const DTDictionary &d,const std::string &s) : dict(d), name(s) {}
    
    operator double(void) {return dict.GetNumber(name);}
    operator DTDoubleArray(void) {return dict.GetArray(name);}
    operator std::string(void) {return dict.GetString(name);}
    operator DTDictionary(void) {return dict.GetDictionary(name);}
    
protected:
    friend class DTMutableDictionaryAssignment;

    DTDictionary dict;
    std::string name;
};


#endif
