//
//  DTTable.cpp
//  Chromosomes
//
//  Created by David Adalsteinsson on 4/2/17.
//
//

#include "DTTable.h"
#include "DTUtilities.h"
#include "DTArrayConversion.h"

#include <set>

// http://docs.rexamine.com/R-devel/Rinternals_8h_source.html

bool DTTableColumnStructure::operator!=(const DTTableColumnStructure &v) const
{
    return (!this->operator==(v));
}

bool DTTableColumnStructure::operator==(const DTTableColumnStructure &v) const
{
    if (name!=v.name) return false;
    if (type!=v.type) return false;
    if (type=="Table")
        return (*tableStructure==*v.tableStructure);
    
    return true;
}

#pragma mark Table Column

DTTableColumn::DTTableColumn()
{
    contentPointer = DTPointer<DTTableColumnBase>(new DTTableColumnNumber());
}

DTTableColumn::DTTableColumn(const std::string &nm)
{
    name = nm;
    contentPointer = DTPointer<DTTableColumnBase>(new DTTableColumnNumber());
}

DTTableColumn::DTTableColumn(const DTPointer<DTTableColumnBase> &ptr,const std::string &nm)
{
    // type = DTTableColumn::Empty;
    contentPointer = ptr;
    if (!contentPointer) {
        DTErrorMessage("DTTableColumn(pointer,name)", "Pointer is NULL");
        contentPointer = DTPointer<DTTableColumnBase>(new DTTableColumnNumber());
    }
    name = nm;
}

void DTTableColumn::SetMask(const DTCharArray &m)
{
    if (m.NotEmpty()) {
        if (m.Length()!=contentPointer->NumberOfRows()) {
            DTErrorMessage("DTTableColumn::SetMask", "Mask length is not valid");
            return;
        }
    }
    mask = m;
}

DTTableColumn DTTableColumn::TextColumn(const std::string &nm,const DTIntArray &off,const DTCharArray &strings)
{
    // This is indexed
    return DTTableColumn::TextColumn(nm,off,DTStringList(strings));
}

DTTableColumn DTTableColumn::TextColumn(const std::string &nm,const DTIntArray &off,const DTStringList &strings)
{
    DTTableColumn col(DTPointer<DTTableColumnBase>(new DTTableColumnText(strings,off)),nm);
    return col;
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTIntArray &a,const DTCharArray &m)
{
    DTTableColumn col(DTPointer<DTTableColumnBase>(new DTTableColumnNumber(a)),nm);
    col.SetMask(m);
    return col;
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTIntArray &a)
{
    return DTTableColumn(new DTTableColumnNumber(a),nm);
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTDoubleArray &a,const DTCharArray &m)
{
    DTTableColumn col(DTPointer<DTTableColumnBase>(new DTTableColumnNumber(a)),nm);
    col.SetMask(m);
    return col;
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTDoubleArray &a)
{
    return DTTableColumn(new DTTableColumnNumber(a),nm);
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTShortIntArray &a,const DTCharArray &m)
{
    return DTTableColumn(new DTTableColumnNumber(a),nm);
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTShortIntArray &a)
{
    return DTTableColumn(new DTTableColumnNumber(a),nm);
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTCharArray &a,const DTCharArray &m)
{
    DTTableColumn col(DTPointer<DTTableColumnBase>(new DTTableColumnNumber(a)),nm);
    col.SetMask(m);
    return col;
}

DTTableColumn DTTableColumn::NumberColumn(const std::string &nm,const DTCharArray &a)
{
    return DTTableColumn(new DTTableColumnNumber(a),nm);
}

DTTableColumn DTTableColumn::DateColumn(const std::string &nm,const DTDoubleArray &a,const DTCharArray &m)
{
    DTTableColumn col(DTPointer<DTTableColumnBase>(new DTTableColumnDate(a)),nm);
    col.SetMask(m);
    return col;
}

DTTableColumn DTTableColumn::DateColumn(const std::string &nm,const DTDoubleArray &a)
{
    return DTTableColumn(new DTTableColumnDate(a),nm);
}

DTTableColumn DTTableColumn::TextColumn(const std::string &nm,const DTStringList &a)
{
    return DTTableColumn(new DTTableColumnText(a),nm);
}

DTTableColumn DTTableColumn::TableColumn(const std::string &nm,const DTTable &t)
{
    return DTTableColumn(new DTTableColumnTable(t),nm);
    //DTTableColumn toReturn(nm);
    // toReturn.type = DTTableColumn::Table;
    //toReturn.table = DTPointer<DTTable>(new DTTable(t));
    //return toReturn;
}

bool DTTableColumn::IsEmpty(void) const
{
    return (contentPointer->NumberOfRows()==0);
}

bool DTTableColumn::NotEmpty(void) const
{
    return (contentPointer->NumberOfRows()!=0);
}

bool DTTableColumn::IsNumberColumn(void) const
{
    return (contentPointer->Type()=="Number");
}

bool DTTableColumn::IsTextColumn(void) const
{
    return (contentPointer->Type()=="Text");
}

bool DTTableColumn::IsDateColumn(void) const
{
    return (contentPointer->Type()=="Date");
}

bool DTTableColumn::IsTable(void) const
{
    return (contentPointer->Type()=="Table");
}

/*
DTTableColumn DTTableColumn::SingleNumberColumn(const std::string &nm,double v)
{
    DTTableColumn toReturn(nm);
    toReturn.type = DTTableColumn::SingleNumber;
    DTMutableDoubleArray single(1);
    single(0) = v;
    toReturn.doubleVersion = single;
    return toReturn;
}

DTTableColumn DTTableColumn::SingleDateColumn(const std::string &nm,double v)
{
    DTTableColumn toReturn(nm);
    toReturn.type = DTTableColumn::SingleDate;
    DTMutableDoubleArray single(1);
    single(0) = v;
    toReturn.doubleVersion = single;
    return toReturn;
}

DTTableColumn DTTableColumn::SingleTextColumn(const std::string &nm,const std::string &v)
{
    DTTableColumn toReturn(nm);
    toReturn.type = DTTableColumn::SingleText;
    DTMutableList<std::string> single(1);
    single(0) = v;
    toReturn.stringVersion = DTStringList(single);
    return toReturn;
}
 
 */

DTTableColumn DTTableColumn::ExtractRows(const DTRange &r) const
{
    DTRange finalRange = Intersection(r,DTRange(0,NumberOfRows()));
    DTTableColumn toReturn(name);
    // toReturn.type = type;
    if (finalRange.length==0) return toReturn;
    return DTTableColumn(contentPointer->ExtractRows(r),Name());
    
    /*
    switch (type) {
        case Empty:
            break;
        case Double:
        case Date:
            toReturn.doubleVersion = ExtractIndices(doubleVersion,final);
            break;
        case Float:
            toReturn.floatVersion = ExtractIndices(floatVersion,final);
            break;
        case Int:
            toReturn.intVersion = ExtractIndices(intVersion,final);
            break;
        case Short:
            toReturn.shortVersion = ExtractIndices(shortVersion,final);
            break;
        case Point2D:
            toReturn.doubleVersion = ExtractColumns(doubleVersion,final);
            break;
        case Character:
            toReturn.charVersion = ExtractIndices(charVersion,final);
            break;
        case IndexedText:
            DTErrorMessage("DTTableColumn::ExtractRows","Not done yet - IndexedText");
            break;
        case Text:
            DTErrorMessage("DTTableColumn::ExtractRows","Not done yet - Text");
            break;
        case Table:
            DTErrorMessage("DTTableColumn::ExtractRows","Not done yet - Table");
            break;
        case SingleDate:
        case SingleNumber:
            toReturn.doubleVersion = doubleVersion;
            break;
        case SingleText:
            toReturn.stringVersion = stringVersion;
            break;
    }
*/

    return toReturn;
}

ssize_t DTTableColumn::NumberOfRows(void) const
{
    return contentPointer->NumberOfRows();

    /*
    ssize_t toReturn = 0;
    switch (type) {
        case Empty:
            toReturn = contentPointer->NumberOfRows();
            break;
        case Double:
        case Date:
            toReturn = doubleVersion.Length();
            break;
        case Float:
            toReturn = floatVersion.Length();
            break;
        case Int:
            toReturn = intVersion.Length();
            break;
        case Short:
            toReturn = shortVersion.Length();
            break;
        case Character:
            toReturn = charVersion.Length();
            break;
        case IndexedText:
            toReturn = intVersion.Length();
            break;
        case Point2D:
            toReturn = doubleVersion.n();
            break;
        case Point3D:
            toReturn = doubleVersion.n();
            break;
        case Text:
            toReturn = stringVersion.NumberOfStrings();
            break;
        case Table:
            return (table ? table->NumberOfRows() : 0);
            break;
        case SingleDate:
        case SingleNumber:
        case SingleText:
            return 0;
            break;
    }
    return toReturn;
     */
}

DTTableColumn DTTableColumn::ChangeName(const std::string &nm) const
{
    DTTableColumn toReturn = *this;
    toReturn.name = nm;
    return toReturn;
}

/*
DTDoubleArray DTTableColumn::NumberArray(void) const
{
    DTDoubleArray toReturn;
    switch (type) {
        case Empty:
            break;
        case Double:
            toReturn = doubleVersion;
            break;
        case Date:
            DTErrorMessage("TableColumn::NumberArray","This is a date column");
            break;
        case Float:
            toReturn = ConvertToDouble(floatVersion);
            break;
        case Int:
            toReturn = ConvertToDouble(intVersion);
            break;
        case Short:
            toReturn = ConvertToDouble(shortVersion);
            break;
        case Character:
            toReturn = ConvertToDouble(charVersion);
            break;
        case IndexedText:
        case Text:
            DTErrorMessage("TableColumn::NumberArray","This is a text column");
            break;
        case Table:
            DTErrorMessage("TableColumn::NumberArray","This is a table");
            break;
        case Point2D:
            DTErrorMessage("TableColumn::NumberArray","This is a list of points");
            break;
        case SingleDate:
        case SingleNumber:
        case SingleText:
            DTErrorMessage("TableColumn::NumberArray","This is a single value");
            break;
    }
    return toReturn;
}
*/

DTTableColumnStructure DTTableColumn::Structure(void) const
{
    DTTableColumnStructure toReturn;
    toReturn.name = name;
    toReturn.type = contentPointer->Type();
    if (toReturn.type=="Table") {
        
        // toReturn.tableStructure = DTPointer<DTTableStructure>(new DTTableStructure(table->Structure()))
    }
    
    /*
    switch (type) {
        case Empty:
            toReturn.type = "Empty";
            break;
        case Double:
        case Float:
        case Int:
        case Short:
        case Character:
            toReturn.type = "Number";
            break;
        case Date:
            toReturn.type = "Date";
            break;
        case IndexedText:
        case Text:
            toReturn.type = "Text";
            break;
        case Table:
            toReturn.type = "Table";
            toReturn.tableStructure = DTPointer<DTTableStructure>(new DTTableStructure(table->Structure()));
            break;
        case SingleDate:
            toReturn.type = "Single Date";
            break;
        case SingleNumber:
            toReturn.type = "Single Number";
            break;
        case SingleText:
            toReturn.type = "Single Text";
            break;
        case Point2D:
            toReturn.type = "Point2D";
            break;
        case Point3D:
            toReturn.type = "Point3D";
            break;
    }
     */
    
    return toReturn;
}

void DTTableColumn::pinfoWithIndent(const std::string &indent) const
{
#ifndef DG_NOSTDErrOut
    std::string padding = ".................................";

    std::string desc = name + " ";
    if (desc.length()<30) desc = desc + string(padding,0,30-desc.length());
    if (desc.length()>30) desc = string(desc,0,30);
    
    std::cerr << indent << desc << " - ";
    
    std::string t = contentPointer->Type();
    ssize_t len = t.length();
    std::string e;
    if (len<10) {
        e.append(10-len,' ');
        e.append(t);
    }
    e.append(" - ");
    std::cerr << e << NumberOfRows() << " rows";
    std::cerr << std::endl;
    
    if (contentPointer->Type()=="Table") {
        (((DTTableColumnTable *)contentPointer.Data())->Table()).pinfoWithIndent(indent+"    ");
    }
    /*
    switch (type) {
        case Empty:
        {
            std::string t = contentPointer->Type();
            ssize_t len = t.length();
            std::string e;
            if (len<10) {
                e.append(10-len,' ');
                e.append(t);
            }
            e.append(" - ");
            std::cerr << e << NumberOfRows() << " rows";
        }
            break;
        case Double:
            std::cerr << "double - " << doubleVersion.Length() << " rows";
            break;
        case Date:
            std::cerr << "  date - " << doubleVersion.Length() << " rows";
            break;
        case Float:
            std::cerr << " float - " << floatVersion.Length() << " rows";
            break;
        case Int:
            std::cerr << "   int - " << intVersion.Length() << " rows";
            break;
        case Short:
            std::cerr << " short - " << shortVersion.Length() << " rows";
            break;
        case Character:
            std::cerr << "  char - " << charVersion.Length() << " rows";
            break;
        case IndexedText:
            std::cerr << "Text (ind) - " << intVersion.Length() << " rows";
            break;
        case Text:
            std::cerr << "  Text - ";
            break;
        case Point2D:
            std::cerr << "Point2D - " << doubleVersion.n() << " points";
            break;
        case Point3D:
            std::cerr << "Point3D - " << doubleVersion.n() << " points";
            break;
        case Table:
            std::cerr << " Table - \n";
            table->pinfoWithIndent(indent+"   ");
            break;
        case SingleDate:
            std::cerr << " date - " << doubleVersion(0) << endl;
            break;
        case SingleNumber:
            std::cerr << " number - " << doubleVersion(0) << endl;
            break;
        case SingleText:
            std::cerr << " text - " << stringVersion(0) << endl;
            break;
    }
    
     */
#endif
}

DTList<std::string> DTTableColumn::StringVersionForPinfo(ssize_t &maxLen) const
{
    return contentPointer->StringVersionForPinfo(maxLen);

    /*
    DTMutableList<std::string> toReturn(NumberOfRows());
    maxLen = name.length();
    ssize_t i,howMany,thisLen;
    std::string str;
    
    switch (type) {
        case Empty:
            if (contentPointer) {
                return contentPointer->StringVersionForPinfo(maxLen);
            }
            else {
                DTErrorMessage("Not implemented");
            }
            break;
        case Double:
            howMany = doubleVersion.Length();
            for (i=0;i<howMany;i++) {
                str = DTFloat2StringShort(doubleVersion(i));
                thisLen = str.length();
                if (thisLen>maxLen) maxLen = thisLen;
                toReturn(i) = str;
            }
            break;
        default:
            DTErrorMessage("DTTableColumn::StringVersionForPinfo","Missing for type");
    }
    
    return toReturn;
     */
}

DTTableColumn DTTableColumn::Read(const DTDataStorage &input,const std::string &name)
{
    std::string nm = input.ReadString(name+"N");
    std::string varName = name+"V";
    std::string type = input.ReadString(varName+"_T");
    
    DTTableColumn col;
    if (type=="Number") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnNumber()),nm);
    }
    else if (type=="UTF8" || type=="Text") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnText()),nm);
    }
    else if (type=="Date") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnDate()),nm);
    }
    else if (type=="Table") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnTable()),nm);
    }
    else if (type=="Point3D") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnPoint3D()),nm);
    }
    else if (type=="Surface") {
        col = DTTableColumn(DTPointer<DTTableColumnBase>(new DTTableColumnSurface()),nm);
    }
    else {
        DTErrorMessage("DTTableColumn::Read", "Unexpected type");
    }
    col.ReadFrom(input,varName);
    if (input.Contains(varName+"_mask")) {
        DTCharArray mask;
        mask = input.ReadCharArray(varName+"_mask");
        col.SetMask(mask);
    }
    
    return col;
}

void DTTableColumn::WriteStructure(DTDataStorage &output,const std::string &saveAs) const
{
    output.Save(name,saveAs+"N");
    output.Save(contentPointer->Type(),saveAs+"T");
    if (contentPointer->Type()=="Table") {
        (((DTTableColumnTable *)contentPointer.Data())->Table()).WriteStructureInternal(output,saveAs+"T");
    }
    /*
    switch (type) {
        case Empty:
            output.Save(contentPointer->Type(),saveAs+"T");
            break;
        case Double:
        case Float:
        case Int:
        case Short:
        case Character:
            output.Save("Number",saveAs+"T");
            break;
        case Date:
            output.Save("Date",saveAs+"T");
            break;
        case Point2D:
            output.Save("Point2D",saveAs+"T");
            break;
        case IndexedText:
        case Text:
            output.Save("Text",saveAs+"T");
            break;
        case Table:
            output.Save("Table",saveAs+"T");
            table->WriteStructureInternal(output,saveAs+"T");
            break;
        case SingleDate:
            output.Save("Single Date",saveAs+"T");
            break;
        case SingleNumber:
            output.Save("Single Number",saveAs+"T");
            break;
        case SingleText:
            output.Save("Single Text",saveAs+"T");
            break;
    }
     */
}

void DTTableColumn::Write(DTDataStorage &output,const std::string &saveAs) const
{
    output.Save(name,saveAs+"N");
    if (mask.NotEmpty())
        output.Save(mask,saveAs+"V_mask");
    
    output.Save(contentPointer->Type(),saveAs+"V_T");
    contentPointer->WriteToFile(output,saveAs+"V");

    /*
    switch (type) {
        case Empty:
            output.Save(contentPointer->Type(),saveAs+"V_T");
            contentPointer->WriteToFile(output,saveAs+"V");
            break;
        case Double:
            output.Save("Number",saveAs+"V_T");
            output.Save(doubleVersion,saveAs+"V");
            break;
        case Float:
            output.Save("Number",saveAs+"V_T");
            output.Save(floatVersion,saveAs+"V");
            break;
        case Int:
            output.Save("Number",saveAs+"V_T");
            output.Save(intVersion,saveAs+"V");
            break;
        case Short:
            output.Save("Number",saveAs+"V_T");
            output.Save(shortVersion,saveAs+"V");
            break;
        case Character:
            output.Save("Number",saveAs+"V_T");
            output.Save(charVersion,saveAs+"V");
            break;
        case Date:
            output.Save("Date",saveAs+"V_T");
            output.Save(doubleVersion,saveAs+"V");
            break;
        case Point2D:
            output.Save("Point2D",saveAs+"V_T");
            output.Save(doubleVersion,saveAs+"V");
            break;
        case IndexedText:
            output.Save("Indexed",saveAs+"V_T");
            output.Save(charVersion,saveAs+"V_S");
            output.Save(intVersion,saveAs+"V");
            break;
        case Text:
            output.Save("UTF8",saveAs+"V_T");
            output.Save(stringVersion.Characters(),saveAs+"V");
            break;
        case Table:
            output.Save("Table",saveAs+"V_T");
            ::Write(output,saveAs+"V",*table);
            break;
        case SingleDate:
            output.Save("Single Date",saveAs+"V_T");
            output.Save(doubleVersion,saveAs+"V");
            break;
        case SingleNumber:
            output.Save("Single Number",saveAs+"V_T");
            output.Save(doubleVersion,saveAs+"V");
            break;
        case SingleText:
            output.Save("Single Text",saveAs+"V_T");
            output.Save(stringVersion(0),saveAs+"V");
            break;
    }
     */
}

void DTTableColumn::WriteSingle(DTDataStorage &output,const std::string &saveAs) const
{
    if (mask.NotEmpty())
        output.Save(mask,saveAs+"_mask");
    
    output.Save(contentPointer->Type(),saveAs+"_T");
    contentPointer->WriteToFile(output,saveAs);
    
    /*
    switch (type) {
        case Empty:
            output.Save(contentPointer->Type(),saveAs+"_T");
            contentPointer->WriteToFile(output,saveAs);
            break;
        case Double:
            output.Save("Number",saveAs+"_T");
            output.Save(doubleVersion,saveAs);
            break;
        case Float:
            output.Save("Number",saveAs+"_T");
            output.Save(floatVersion,saveAs);
            break;
        case Int:
            output.Save("Number",saveAs+"_T");
            output.Save(intVersion,saveAs);
            break;
        case Short:
            output.Save("Number",saveAs+"_T");
            output.Save(shortVersion,saveAs);
            break;
        case Character:
            output.Save("Number",saveAs+"_T");
            output.Save(charVersion,saveAs);
            break;
        case Date:
            output.Save("Date",saveAs+"_T");
            output.Save(doubleVersion,saveAs);
            break;
        case Point2D:
            output.Save("Point2D",saveAs+"_T");
            output.Save(doubleVersion,saveAs);
            break;
        case IndexedText:
            output.Save("Indexed",saveAs+"_T");
            output.Save(charVersion,saveAs+"_S");
            output.Save(intVersion,saveAs);
            break;
        case Text:
            output.Save("UTF8",saveAs+"_T");
            output.Save(stringVersion.Characters(),saveAs);
            break;
        case Table:
            ::Write(output,saveAs,*table);
            break;
        case SingleDate:
            output.Save("Single Date",saveAs+"_T");
            output.Save(doubleVersion,saveAs);
            break;
        case SingleNumber:
            output.Save("Single Number",saveAs+"_T");
            output.Save(doubleVersion,saveAs);
            break;
        case SingleText:
            output.Save("Single Text",saveAs+"_T");
            output.Save(stringVersion(0),saveAs);
            break;
    }
    */
}

void DTTableColumn::ReadFrom(const DTDataStorage &input,const std::string &varName)
{
    contentPointer->ReadFrom(input,varName);
}

#pragma mark DTTableStructure

bool DTTableStructure::operator!=(const DTTableStructure &v) const
{
    ssize_t i,howMany = columns.Length();
    if (howMany!=v.columns.Length()) return true;
    
    for (i=0;i<howMany;i++) {
        if (columns(i)!=v.columns(i)) return true;
    }
    return false;
}

bool DTTableStructure::operator==(const DTTableStructure &v) const
{
    return !(this->operator!=(v));
}

#pragma mark Table

DTTable::DTTable(const DTList<DTTableColumn> &list)
{
    _numberOfRows = 0;
    
    ssize_t i,howMany = list.Length();
    
    // Require the names to be unique.
    std::set<std::string> theNames;
    std::string colName;
    for (i=0;i<howMany;i++) {
        colName = list(i).Name();
        if (theNames.count(colName)) {
            DTErrorMessage("DTTable(Columns)","Column names have to be unique");
            break;
        }
        theNames.insert(colName);
    }
    if (i<howMany) {
        return;
    }
    
    _columns = list;
    // Find the maximum row length
    ssize_t maxLen = 0;
    ssize_t thisLen;
    for (i=0;i<howMany;i++) {
        thisLen = list(i).NumberOfRows();
        if (thisLen>maxLen) maxLen = thisLen;
    }
    _numberOfRows = maxLen;
}

/*
DTTable::ColumnType DTTable::Type(int col) const
{
    if (col<0 || col>=_columns.Length()) {
        DTErrorOutOfRange("DTTable", col,_columns.Length());
        return DTTable::OutOfBounds;
    }
    
    const DTTableColumn &column = _columns(col);
    DTTable::ColumnType toReturn = DTTable::Empty;
    
    switch (column.type) {
        case DTTableColumn::Empty:
            break;
        case DTTableColumn::Double:
        case DTTableColumn::Float:
        case DTTableColumn::Int:
        case DTTableColumn::Short:
        case DTTableColumn::Character:
            toReturn = DTTable::Numerical;
            break;
        case DTTableColumn::Date:
            toReturn = DTTable::Date;
            break;
        case DTTableColumn::IndexedText:
        case DTTableColumn::Text:
            toReturn = DTTable::Text;
            break;
        case DTTableColumn::Table:
            toReturn = DTTable::Table;
            break;
        case DTTableColumn::SingleDate:
            toReturn = DTTable::SingleDate;
            break;
        case DTTableColumn::SingleNumber:
            toReturn = DTTable::SingleNumber;
            break;
        case DTTableColumn::SingleText:
            toReturn = DTTable::SingleText;
            break;
        case DTTableColumn::Point2D:
            toReturn = DTTable::Point2D;
            break;
    }
    
    return toReturn;
}
 */

DTTableColumn DTTable::Column(int col) const
{
    if (col<0 || col>=_columns.Length()) {
        DTErrorOutOfRange("DTTable", col,_columns.Length());
        return DTTableColumn();
    }

    return _columns(col);
}

DTTableColumn DTTable::operator()(const std::string &nm) const
{
    ssize_t col,howMany = _columns.Length();
    for (col=0;col<howMany;col++) {
        if (_columns(col).Name()==nm) {
            return _columns(col);
        }
    }
    DTErrorMessage("Table(name)","Column not found");
    return DTTableColumn(nm);
}

DTTableStructure DTTable::Structure(void) const
{
    ssize_t i,howMany = _columns.Length();
    DTMutableList<DTTableColumnStructure> columnStructure(howMany);
    for (i=0;i<howMany;i++) {
        columnStructure(i) = _columns(i).Structure();
    }
    DTTableStructure toReturn;
    toReturn.columns = columnStructure;
    return toReturn;
}

DTTable DTTable::ExtractRows(const DTRange &range) const
{
    ssize_t i,howMany = _columns.Length();
    DTMutableList<DTTableColumn> newColumns(howMany);
    for (i=0;i<howMany;i++) {
        newColumns(i) = _columns(i).ExtractRows(range);
    }
    return DTTable(newColumns);
}

void DTTable::WriteStructureInternal(DTDataStorage &dataFile,const std::string &writeName) const
{
    int howManyColumns = int(NumberOfColumns());
    for (int i=0;i<howManyColumns;i++) {
        Column(i).WriteStructure(dataFile,writeName+"_"+DTInt2String(i+1));
    }
    
    dataFile.Save(howManyColumns,writeName+"_N");
}

void DTTable::WriteStructure(DTDataStorage &dataFile,const std::string &writeName) const
{
    std::string infoName = "SeqInfo_" + writeName;
    WriteStructureInternal(dataFile,infoName);
}

void DTTable::pinfoWithIndent(const std::string &indent) const
{
    ssize_t howMany = NumberOfColumns();
    for (int i=0;i<howMany;i++) {
        Column(i).pinfoWithIndent(indent);
    }
}

void DTTable::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t howMany = NumberOfColumns();
    std::cerr << howMany << " columns and " << NumberOfRows() << " rows" << std::endl;
    pinfoWithIndent("  ");
#endif
}

void DTTable::pall(void) const
{
#ifndef DG_NOSTDErrOut
    ssize_t howManyColumns = NumberOfColumns();
    ssize_t howManyRows = NumberOfRows();
    
    DTMutableList<DTList<std::string> > content(howManyColumns);
    DTMutableIntArray maximumLengths(howManyColumns);
    int colNumber;
    ssize_t entryLength;
    ssize_t maxStrLen = 0;
    ssize_t totalWidth = 0;
    for (colNumber=0;colNumber<howManyColumns;colNumber++) {
        content(colNumber) = Column(colNumber).StringVersionForPinfo(maxStrLen);
        std::string headerName = Column(colNumber).Name();
        entryLength = headerName.length();
        if (entryLength>maxStrLen) {
            maxStrLen = entryLength;
        }
        maximumLengths(colNumber) = int(maxStrLen);
        totalWidth += maxStrLen+1;
    }
    
    std::string line;
    ssize_t rowNumber;

    line.append("  # |");
    for (colNumber=0;colNumber<howManyColumns;colNumber++) {
        line.append(" ");
        std::string headerName = Column(colNumber).Name();
        entryLength = headerName.length();
        if (entryLength+1<maximumLengths(colNumber)) {
            headerName.append(maximumLengths(colNumber)-entryLength,' ');
        }
        line.append(headerName);
        line.append(" |");
    }
    ssize_t padLength;
    
    std::cerr << line << std::endl;
    line.clear();
    line.append("----|");
    for (colNumber=0;colNumber<howManyColumns;colNumber++) {
        entryLength = maximumLengths(colNumber);
        line.append(entryLength+2,'-');
        line.append("|");
    }
    std::cerr << line << std::endl;
    line.clear();

    // Now the content
    std::string entry;
    for (rowNumber=0;rowNumber<howManyRows;rowNumber++) {
        entry = DTInt2String(int(rowNumber)+1);
        padLength = 3-entry.length();
        if (padLength>0) {
            line.append(padLength,' ');
        }
        line.append(entry);
        line.append(" |");
        for (colNumber=0;colNumber<howManyColumns;colNumber++) {
            line.append(" ");
            entryLength = maximumLengths(colNumber);
            DTList<std::string> &columnContent = content(colNumber);
            if (rowNumber<columnContent.Length()) {
                entry = columnContent(rowNumber);
                line.append(entry);
                padLength = entryLength-entry.length();
                if (padLength>0) line.append(padLength,' ');
            }
            else {
                line.append(entryLength,' ');
            }
            line.append(" |");
        }
        std::cerr << line << std::endl;
        line.clear();
    }
    std::cerr << std::endl;
#endif
}

void Read(const DTDataStorage &input,const std::string &readName,DTTable &toReturn)
{
    if (input.Contains(readName)==false) {
        toReturn = DTTable();
        return;
    }
    int howManyColumns = input.ReadInt(readName);
    if (howManyColumns<=0) {
        toReturn = DTTable();
        return;
    }
    
    DTMutableList<DTTableColumn> columns(howManyColumns);
    
    int i;
    for (i=0;i<howManyColumns;i++) {
        columns(i) = DTTableColumn::Read(input,readName+"_"+DTInt2String(i));
    }
 
    toReturn = DTTable(columns);
}

void Write(DTDataStorage &output,const std::string &name,const DTTable &theVar)
{
    ssize_t howMany = theVar.NumberOfColumns();
    int i;
    for (i=0;i<howMany;i++) {
        theVar.Column(i).Write(output,name+"_"+DTInt2String(i));
    }
    
    output.Save(int(theVar.NumberOfRows()),name+"_R");
    output.Save(int(theVar.NumberOfColumns()),name);
}

void WriteOne(DTDataStorage &output,const std::string &name,const DTTable &toWrite)
{
    Write(output,name,toWrite);
    Write(output,"Seq_"+name,"Table");
    toWrite.WriteStructure(output,name);
    output.Flush();
}

#pragma mark Columns

DTTableColumnBase::DTTableColumnBase()
{
    numberOfRows = 0;
}

#pragma mark Numbers

DTTableColumnNumber::DTTableColumnNumber(const DTDoubleArray &list)
{
    doubleVersion = list;
    numberOfRows = doubleVersion.Length();
}

DTTableColumnNumber::DTTableColumnNumber(const DTFloatArray &list)
{
    floatVersion = list;
    numberOfRows = floatVersion.Length();
}

DTTableColumnNumber::DTTableColumnNumber(const DTIntArray &list)
{
    intVersion = list;
    numberOfRows = intVersion.Length();
}

DTTableColumnNumber::DTTableColumnNumber(const DTShortIntArray &list)
{
    shortVersion = list;
    numberOfRows = shortVersion.Length();
}

DTTableColumnNumber::DTTableColumnNumber(const DTCharArray &list)
{
    charVersion = list;
    numberOfRows = charVersion.Length();
}

DTDoubleArray DTTableColumnNumber::Values(void) const
{
    if (doubleVersion.NotEmpty()) {
        return doubleVersion;
    }
    else if (floatVersion.NotEmpty()) {
        return ConvertToDouble(floatVersion);
    }
    else if (intVersion.NotEmpty()) {
        return ConvertToDouble(intVersion);
    }
    else if (shortVersion.NotEmpty()) {
        return ConvertToDouble(shortVersion);
    }
    else if (charVersion.NotEmpty()) {
        return ConvertToDouble(charVersion);
    }
    else {
        return DTDoubleArray();
    }
}

DTPointer<DTTableColumnBase> DTTableColumnNumber::ExtractRows(const DTRange &r) const
{
    DTRange finalRange = Intersection(r,DTRange(0,NumberOfRows()));

    if (doubleVersion.NotEmpty()) {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(doubleVersion, finalRange)));
    }
    else if (floatVersion.NotEmpty()) {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(floatVersion, finalRange)));
    }
    else if (intVersion.NotEmpty()) {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(intVersion, finalRange)));
    }
    else if (shortVersion.NotEmpty()) {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(shortVersion, finalRange)));
    }
    else if (charVersion.NotEmpty()) {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(charVersion, finalRange)));
    }
    else {
        return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(DTDoubleArray()));
    }
}

DTList<std::string> DTTableColumnNumber::StringVersionForPinfo(ssize_t &maxLen) const
{
    ssize_t i,thisLen;
    DTMutableList<std::string> toReturn(numberOfRows);
    std::string str;
    maxLen = 0;
    if (doubleVersion.Length()==numberOfRows) {
        for (i=0;i<numberOfRows;i++) {
            str = DTFloat2StringShort(doubleVersion(i));
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }
    else if (floatVersion.Length()==numberOfRows) {
        for (i=0;i<numberOfRows;i++) {
            str = DTFloat2StringShort(floatVersion(i));
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }
    else if (intVersion.Length()==numberOfRows) {
        for (i=0;i<numberOfRows;i++) {
            str = DTFloat2StringShort(intVersion(i));
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }
    else if (shortVersion.Length()==numberOfRows) {
        for (i=0;i<numberOfRows;i++) {
            str = DTFloat2StringShort(shortVersion(i));
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }
    else if (charVersion.Length()==numberOfRows) {
        for (i=0;i<numberOfRows;i++) {
            str = DTFloat2StringShort(charVersion(i));
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }

    return toReturn;
}

void DTTableColumnNumber::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    if (doubleVersion.NotEmpty()) {
        Write(output, name, doubleVersion);
    }
    else if (floatVersion.NotEmpty()) {
        Write(output, name, floatVersion);
    }
    else if (intVersion.NotEmpty()) {
        Write(output, name, intVersion);
    }
    else if (shortVersion.NotEmpty()) {
        Write(output, name, shortVersion);
    }
    else if (charVersion.NotEmpty()) {
        Write(output, name, charVersion);
    }
    else {
        Write(output, name, DTDoubleArray());
    }
}

void DTTableColumnNumber::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    if (input.SavedAsDouble(name)) {
        doubleVersion = input.ReadDoubleArray(name);
        numberOfRows = doubleVersion.Length();
    }
    else if (input.SavedAsFloat(name)) {
        floatVersion = input.ReadFloatArray(name);
        numberOfRows = floatVersion.Length();
    }
    else if (input.SavedAsInt(name)) {
        intVersion = input.ReadIntArray(name);
        numberOfRows = intVersion.Length();
    }
    else if (input.SavedAsShort(name)) {
        shortVersion = input.ReadShortIntArray(name);
        numberOfRows = shortVersion.Length();
    }
    else if (input.SavedAsCharacter(name)) {
        charVersion = input.ReadCharArray(name);
        numberOfRows = charVersion.Length();
    }
    else {
        doubleVersion = input.ReadDoubleArray(name);
        numberOfRows = doubleVersion.Length();
    }
}

#pragma mark Date

DTTableColumnDate::DTTableColumnDate(const DTDoubleArray &list)
{
    doubleVersion = list;
    numberOfRows = doubleVersion.Length();
}

DTDoubleArray DTTableColumnDate::Values(void) const
{
    return doubleVersion;
}

DTPointer<DTTableColumnBase> DTTableColumnDate::ExtractRows(const DTRange &r) const
{
    DTRange finalRange = Intersection(r,DTRange(0,NumberOfRows()));
    return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(ExtractIndices(doubleVersion, finalRange)));
}

DTList<std::string> DTTableColumnDate::StringVersionForPinfo(ssize_t &maxLen) const
{
    ssize_t i,thisLen;
    DTMutableList<std::string> toReturn(numberOfRows);
    std::string str;
    maxLen = 0;
    for (i=0;i<numberOfRows;i++) {
        str = DTFloat2StringShort(doubleVersion(i));
        toReturn(i) = str;
        thisLen = str.length();
        if (thisLen>maxLen) maxLen = thisLen;
    }

    return toReturn;
}

void DTTableColumnDate::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    Write(output, name, doubleVersion);
}

void DTTableColumnDate::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    doubleVersion = input.ReadDoubleArray(name);
    numberOfRows = doubleVersion.Length();
}

#pragma mark Text

DTTableColumnText::DTTableColumnText(const DTStringList &list)
{
    stringList = list;
    numberOfRows = stringList.NumberOfStrings();
    isIndexed = false;
}

DTTableColumnText::DTTableColumnText(const DTStringList &entries,const DTIntArray &index)
{
    // Validity check
    ssize_t i,howMany = index.Length();
    int howManyEntries = entries.NumberOfStrings();
    int val;
    for (i=0;i<howMany;i++) {
        val = index(i);
        if (val<-1) {
            DTErrorMessage("DTTableColumnText(StringList,IntArray)","Negative index");
            break;
        }
        else if (val>=howManyEntries) {
            DTErrorMessage("DTTableColumnText(StringList,IntArray)","Index out of bounds");
            break;
        }
    }
    if (i==howMany) {
        stringList = entries;
        indexed = index;
        numberOfRows = indexed.Length();
        isIndexed = true;
    }
    else {
        numberOfRows = 0;
        isIndexed = false;
    }
}

DTStringList DTTableColumnText::StringList(void) const
{
    return stringList;
}

DTPointer<DTTableColumnBase> DTTableColumnText::ExtractRows(const DTRange &r) const
{
    DTRange finalRange = Intersection(r,DTRange(0,NumberOfRows()));
    if (isIndexed) {
        DTMutableIntArray subIndices = ExtractIndices(indexed,finalRange);
        ssize_t howManyLabels = stringList.NumberOfStrings();
        DTMutableCharArray found(howManyLabels);
        found = 0;
        int ind;
        ssize_t i,newLength = subIndices.Length();
        for (i=0;i<newLength;i++) {
            ind = subIndices(i);
            if (ind>=0) found(ind) = 1;
        }
        // Remapper
        ssize_t pos = 0;
        DTMutableIntArray newNumber(howManyLabels);
        newNumber = -1;
        for (i=0;i<howManyLabels;i++) {
            if (found(i)) {
                newNumber(i) = pos++;
            }
        }
        // Extract these entries from the string list
        DTMutableList<std::string> newEntries(pos);
        int posInNew = 0;
        for (i=0;i<howManyLabels;i++) {
            if (found(i)) {
                newEntries(posInNew++) = stringList(i);
            }
        }
        
        DTMutableIntArray newIndex(newLength);
        for (i=0;i<newLength;i++) {
            ind = subIndices(i);
            if (ind>=0) {
                subIndices(i) = newNumber(ind);
            }
        }
        return DTPointer<DTTableColumnBase>(new DTTableColumnText(DTStringList(newEntries),subIndices));
    }
    else {
        DTStringList subset = ExtractIndices(stringList,finalRange);
        return DTPointer<DTTableColumnBase>(new DTTableColumnText(subset));
    }
}

DTList<std::string> DTTableColumnText::StringVersionForPinfo(ssize_t &maxLen) const
{
    ssize_t i,thisLen;
    DTMutableList<std::string> toReturn(numberOfRows);
    std::string str;
    maxLen = 0;
    if (isIndexed) {
        int ind;
        for (i=0;i<numberOfRows;i++) {
            ind = indexed(i);
            if (ind>=0) {
                str = stringList(ind);
                toReturn(i) = str;
                thisLen = str.length();
                if (thisLen>maxLen) maxLen = thisLen;
            }
        }
    }
    else {
        for (i=0;i<numberOfRows;i++) {
            str = stringList(i);
            toReturn(i) = str;
            thisLen = str.length();
            if (thisLen>maxLen) maxLen = thisLen;
        }
    }
    
    return toReturn;
}

void DTTableColumnText::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    if (isIndexed) {
        Write(output,name+"_S",stringList);
        Write(output,name,indexed);
    }
    else {
        Write(output,name,stringList);
    }
}

void DTTableColumnText::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    if (input.Contains(name+"_S")) {
        Read(input,name+"_S",stringList);
        Read(input,name,indexed);
        numberOfRows = indexed.Length();
        isIndexed = true;
    }
    else {
        Read(input,name,stringList);
        numberOfRows = stringList.NumberOfStrings();
        isIndexed = false;
    }
}

#pragma mark Points in 3D

DTTableColumnPoint3D::DTTableColumnPoint3D(const DTPointCollection3D &list)
{
    pointList = list;
    numberOfRows = pointList.NumberOfPoints();
}

DTPointer<DTTableColumnBase> DTTableColumnPoint3D::ExtractRows(const DTRange &) const
{
    DTErrorMessage("Not defined yet");
    return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(DTDoubleArray()));
}

DTList<std::string> DTTableColumnPoint3D::StringVersionForPinfo(ssize_t &maxLen) const
{
    ssize_t i,thisLen,howMany = pointList.NumberOfPoints();
    DTMutableList<std::string> toReturn(howMany);
    std::string str;
    DTPoint3D p;
    maxLen = 0;
    for (i=0;i<howMany;i++) {
        p = pointList(i);
        str = DTFloat2StringShort(p.x);
        str.append(", ");
        str.append(DTFloat2StringShort(p.y));
        str.append(", ");
        str.append(DTFloat2StringShort(p.z));
        thisLen = str.length();
        if (thisLen>maxLen) maxLen = thisLen;
        toReturn(i) = str;
    }

    return toReturn;
}

void DTTableColumnPoint3D::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    Write(output, name, pointList);
}

void DTTableColumnPoint3D::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    Read(input, name,pointList);
    numberOfRows = pointList.NumberOfPoints();
}

#pragma mark Surfaces

DTTableColumnSurface::DTTableColumnSurface(const DTSurface3D &surf,const DTIntArray &startEnd)
{
    surface = surf;
    startsOfIntervals = startEnd;
    numberOfRows = startsOfIntervals.Length()-1;
}

DTPointer<DTTableColumnBase> DTTableColumnSurface::ExtractRows(const DTRange &) const
{
    DTErrorMessage("Not defined yet");
    return DTPointer<DTTableColumnBase>(new DTTableColumnNumber(DTDoubleArray()));
}

DTList<std::string> DTTableColumnSurface::StringVersionForPinfo(ssize_t &maxLen) const
{
    ssize_t i,thisLen,howMany = numberOfRows;
    DTMutableList<std::string> toReturn(howMany);
    std::string str;
    DTPoint3D p;
    maxLen = 0;
    for (i=0;i<howMany;i++) {
        int numberOfTriangles = startsOfIntervals(i+1)-startsOfIntervals(i);
        str = DTInt2String(int(numberOfTriangles)) + " tri";
        thisLen = str.length();
        if (thisLen>maxLen) maxLen = thisLen;
        toReturn(i) = str;
    }

    return toReturn;
}

void DTTableColumnSurface::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    Write(output,name+"_St",startsOfIntervals);
    Write(output, name,surface);
}

void DTTableColumnSurface::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    Read(input, name,surface);
    Read(input,name+"_St",startsOfIntervals);
    numberOfRows = startsOfIntervals.Length()-1;
}

#pragma mark Table

DTList<std::string> DTTableColumnTable::StringVersionForPinfo(ssize_t &maxLen) const
{
    maxLen = 0;
    return DTList<std::string>();
}

DTPointer<DTTableColumnBase> DTTableColumnTable::ExtractRows(const DTRange &) const
{
    DTErrorMessage("Not defined yet");
    return DTPointer<DTTableColumnBase>(new DTTableColumnTable(DTTable()));
}

void DTTableColumnTable::WriteToFile(DTDataStorage &output,const std::string &name) const
{
    ::Write(output,name,table);
}

void DTTableColumnTable::ReadFrom(const DTDataStorage &input,const std::string &name)
{
    Read(input, name,table);
}
