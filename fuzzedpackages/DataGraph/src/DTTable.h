// Part of DTSource. Copyright 2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTTable_H
#define DTTable_H

#include "DTIntArray.h"
#include "DTDoubleArray.h"
#include "DTFloatArray.h"
#include "DTIntArray.h"
#include "DTShortIntArray.h"
#include "DTCharArray.h"
#include "DTStringList.h"
#include "DTList.h"
#include "DTPointer.h"

class DTDataStorage;
class DTTable;
class DTTableColumn;
struct DTTableStructure;
struct DTRange;
class DTTableColumnBase;


// A column can be a number, text or date, and you can also have sub-tables.
struct DTTableColumnStructure
{
    std::string name;
    std::string type; // Number, Text, Date, Point2D, Table
    DTPointer<DTTableStructure> tableStructure;
    
    bool operator!=(const DTTableColumnStructure &v) const;
    bool operator==(const DTTableColumnStructure &v) const;
};

class DTTableColumn
{
public:
    enum StorageType {Empty,Double, Float, Int, Short, Character, Date, Point2D, Point3D, IndexedText, Text, Table, SingleNumber, SingleText, SingleDate};
            
    DTTableColumn();
    explicit DTTableColumn(const std::string &nm);
    DTTableColumn(const DTPointer<DTTableColumnBase> &,const std::string &);
    
    void SetMask(const DTCharArray &);
    
    // StorageType Type(void) const {return type;}
    DTPointer<DTTableColumnBase> ContentPointer(void) const {return contentPointer;}
    
    static DTTableColumn TextColumn(const std::string &,const DTIntArray &,const DTCharArray &);
    static DTTableColumn TextColumn(const std::string &,const DTIntArray &,const DTStringList &);
    static DTTableColumn TextColumn(const std::string &,const DTStringList &);

    static DTTableColumn NumberColumn(const std::string &,const DTIntArray &,const DTCharArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTIntArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTDoubleArray &,const DTCharArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTDoubleArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTShortIntArray &,const DTCharArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTShortIntArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTCharArray &,const DTCharArray &);
    static DTTableColumn NumberColumn(const std::string &,const DTCharArray &);
    
    static DTTableColumn DateColumn(const std::string &,const DTDoubleArray &,const DTCharArray &);
    static DTTableColumn DateColumn(const std::string &,const DTDoubleArray &);
    
    static DTTableColumn SingleNumberColumn(const std::string &,double);
    static DTTableColumn SingleDateColumn(const std::string &,double);
    static DTTableColumn SingleTextColumn(const std::string &,const std::string &);
    
    static DTTableColumn TableColumn(const std::string &,const DTTable &);
    
    bool IsEmpty(void) const;
    bool NotEmpty(void) const;
    bool IsNumberColumn(void) const;
    bool IsTextColumn(void) const;
    bool IsDateColumn(void) const;
    bool IsTable(void) const;
    
    // Get data
    ssize_t NumberOfRows(void) const;
    DTTableColumnStructure Structure(void) const;
    // DTDoubleArray NumberArray(void) const;
    
    DTTableColumn ExtractRows(const DTRange &) const;
    DTTableColumn ChangeName(const std::string &) const;
    std::string Name(void) const {return name;}

    void pinfoWithIndent(const std::string &) const;
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    
    static DTTableColumn Read(const DTDataStorage &output,const std::string &name);
    void Write(DTDataStorage &output,const std::string &name) const;
    void WriteSingle(DTDataStorage &output,const std::string &name) const;
    void WriteStructure(DTDataStorage &dataFile,const std::string &name) const;
    void ReadFrom(const DTDataStorage &input,const std::string &name);

private:
    std::string name;
    DTPointer<DTTableColumnBase> contentPointer;
    DTCharArray mask;

    friend class DTTable;
};

// Check the structure of a table, so that it can be checked.
struct DTTableStructure
{
    DTList<DTTableColumnStructure> columns;

    bool operator!=(const DTTableStructure &v) const;
    bool operator==(const DTTableStructure &v) const;
};

class DTTable
{
public:
    DTTable() : _numberOfRows(0) {}
    DTTable(const DTList<DTTableColumn> &);
    
    enum ColumnType {OutOfBounds,Empty,Numerical, Date, Text, Point2D, Table,SingleNumber,SingleDate,SingleText};
    
    bool IsEmpty(void) const {return (_columns.Length()==0);}
    ssize_t NumberOfColumns(void) const {return _columns.Length();}
    ssize_t NumberOfRows(void) const {return _numberOfRows;}
    
    ColumnType Type(int col) const;
    DTTableColumn Column(int col) const;
    DTTableColumn operator()(const std::string &) const;
    DTTableStructure Structure(void) const;
    
    DTTable ExtractRows(const DTRange &) const;
    
    void pinfo(void) const;
    void pall(void) const;
    
    void WriteStructure(DTDataStorage &,const std::string &) const;
    
private:
    void pinfoWithIndent(const std::string &) const;
    void WriteStructureInternal(DTDataStorage &dataFile,const std::string &name) const;
    
    DTList<DTTableColumn> _columns;
    ssize_t _numberOfRows;

    friend class DTTableColumn;
};

// Reading and writing
extern void Read(const DTDataStorage &input,const std::string &name,DTTable &toReturn);
extern void Write(DTDataStorage &output,const std::string &name,const DTTable &theVar);
extern void WriteOne(DTDataStorage &output,const std::string &name,const DTTable &toWrite); // One time value, self documenting.


// *************************************************************************************
// Content for each table.
// *************************************************************************************

class DTTableColumnBase
{
public:
    DTTableColumnBase();
    virtual ~DTTableColumnBase() {}
    ssize_t NumberOfRows() const {return numberOfRows;}
    virtual std::string Type(void) const = 0;
    virtual DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const = 0;
    virtual void WriteToFile(DTDataStorage &output,const std::string &name) const = 0;
    virtual void ReadFrom(const DTDataStorage &output,const std::string &name) = 0;
    
    virtual DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const = 0;

protected:
    ssize_t numberOfRows;
};

class DTTableColumnTable : public DTTableColumnBase
{
public:
    DTTableColumnTable() {}
    explicit DTTableColumnTable(const DTTable &t) : table(t) {numberOfRows = t.NumberOfRows();}
    ~DTTableColumnTable() {};

    DTTable Table(void) const {return table;}

    std::string Type(void) const {return "Table";}
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);

private:
    DTTable table;
};

class DTTableColumnNumber : public DTTableColumnBase
{
public:
    DTTableColumnNumber() {}
    explicit DTTableColumnNumber(const DTDoubleArray &);
    explicit DTTableColumnNumber(const DTFloatArray &);
    explicit DTTableColumnNumber(const DTIntArray &);
    explicit DTTableColumnNumber(const DTShortIntArray &);
    explicit DTTableColumnNumber(const DTCharArray &);
    ~DTTableColumnNumber() {};

    DTDoubleArray Values(void) const;

    std::string Type(void) const {return "Number";}
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);

private:
    DTDoubleArray doubleVersion;
    DTFloatArray floatVersion;
    DTIntArray intVersion;
    DTShortIntArray shortVersion;
    DTCharArray charVersion;
};

class DTTableColumnDate : public DTTableColumnBase
{
public:
    DTTableColumnDate() {}
    explicit DTTableColumnDate(const DTDoubleArray &);
    ~DTTableColumnDate() {};

    DTDoubleArray Values(void) const;

    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    std::string Type(void) const {return "Date";}
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);

private:
    DTDoubleArray doubleVersion;
};

class DTTableColumnText : public DTTableColumnBase
{
public:
    DTTableColumnText() : isIndexed(false) {}
    explicit DTTableColumnText(const DTStringList &);
    DTTableColumnText(const DTStringList &,const DTIntArray &);
    ~DTTableColumnText() {};

    DTStringList StringList(void) const;

    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    std::string Type(void) const {return "Text";}
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);

private:
    DTStringList stringList;
    
    // Two ways to store this.  Either as a list of string objects, or as
    // an indexed list.  In the second case the indexed array is an index into the stringList, and the
    // length of the column is based on the length of the indexec list.
    // An offset of -1 means that the string is empty.
    bool isIndexed;
    DTIntArray indexed;
};

#include "DTPointCollection3D.h"

class DTTableColumnPoint3D : public DTTableColumnBase
{
public:
    DTTableColumnPoint3D() {}
    explicit DTTableColumnPoint3D(const DTPointCollection3D &);
    ~DTTableColumnPoint3D() {};

    DTPointCollection3D Points(void) const {return pointList;}

    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    std::string Type(void) const {return "Point3D";}
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);

private:
    DTPointCollection3D pointList;
};

#include "DTSurface3D.h"

class DTTableColumnSurface : public DTTableColumnBase
{
public:
    DTTableColumnSurface() {}
    DTTableColumnSurface(const DTSurface3D &,const DTIntArray &);
    ~DTTableColumnSurface() {};

    std::string Type(void) const {return "Surface";}
    DTPointer<DTTableColumnBase> ExtractRows(const DTRange &) const;
    DTList<std::string> StringVersionForPinfo(ssize_t &maxLen) const;
    void WriteToFile(DTDataStorage &output,const std::string &name) const;
    void ReadFrom(const DTDataStorage &output,const std::string &name);
    
    DTSurface3D Surface(void) const {return surface;}
    DTIntArray StartsOfIntervals(void) const {return startsOfIntervals;}

private:
    DTSurface3D surface;
    DTIntArray startsOfIntervals;
};


#endif
