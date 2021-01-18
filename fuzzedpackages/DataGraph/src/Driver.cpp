#include <Rcpp.h>
using namespace Rcpp;

#include "DTDataFile.h"
#include "DTUtilities.h"
#include "DTDoubleArray.h"
#include "DTCharArray.h"
#include "DTIntArray.h"
#include "DTTable.h"
#include "DTArrayConversion.h"
#include "DTDoubleArrayOperators.h"
#include "DTDictionary.h"
#include "DTMap.h"

// http://kbroman.org/pkg_primer/pages/cran.html
// Submission policies
// https://cran.r-project.org/web/packages/policies.html

DTTable ConvertToTable(DataFrame df);
DTTable ConvertFromTimeSeries(const std::string &name,SEXP x);
DTTable ConvertFromMatrix(const std::string &name,SEXP x);
bool ConvertToTableIfPossible(const std::string &name,SEXP x,DTTable &returnTable);


DTTableColumn ConvertFromLogicalColumn(const std::string &name,SEXP x)
{
    int *data = INTEGER(x);
    int i,n = Rf_length(x);
    DTMutableCharArray values(n);
    DTMutableCharArray mask;
    
    bool saveMask = false;
    for (i=0;i<n;i++) {
        if (data[i]==NA_LOGICAL) {
            if (saveMask==false) {
                mask = DTMutableCharArray(n);
                mask = 1;
                saveMask = true;
            }
            mask(i) = 0;
            values(i) = 0;
        }
        else {
            values(i) = (data[i]!=0);
        }
    }
    
    if (saveMask) {
        return DTTableColumn::NumberColumn(name,values,mask);
    }
    else {
        return DTTableColumn::NumberColumn(name,values);
    }
}

DTCharArray UTF8BufferFrom(SEXP x)
{
    if (TYPEOF(x) != STRSXP) {
        Rcout << "Internal error, needs to be a string column, the type is " << Rf_type2char(TYPEOF(x)) << " (" << TYPEOF(x) << ")" << endl;
        return DTCharArray();
    }

    int lengthOfBuffer = 1000;
    DTMutableCharArray utf8Buffer(lengthOfBuffer);
    int posInBuffer = 0;
    
    int n = Rf_length(x);
    for (int i = 0; i < n; ++i) {
        SEXP xi = STRING_ELT(x, i);
        if (xi==NA_STRING) {
            if (posInBuffer+1>lengthOfBuffer) {
                utf8Buffer = IncreaseSize(utf8Buffer,lengthOfBuffer);
                lengthOfBuffer = utf8Buffer.Length();
            }
            utf8Buffer(posInBuffer++) = 0;
        }
        else {
            const char* utf8 = Rf_translateCharUTF8(xi);
            int length = strlen(utf8);
            if (posInBuffer+length+1>lengthOfBuffer) {
                utf8Buffer = IncreaseSize(utf8Buffer,posInBuffer+length+1);
                lengthOfBuffer = utf8Buffer.Length();
            }
            memcpy(utf8Buffer.Pointer()+posInBuffer,utf8,length+1);
            posInBuffer+=length+1;
        }
    }
    utf8Buffer = TruncateSize(utf8Buffer,posInBuffer);
    
    return utf8Buffer;
}

DTTableColumn ConvertFromIndexedStrings(const std::string &name,SEXP x)
{
    if (TYPEOF(x) != INTSXP) {
        Rcout << "The column " << name << " is corrupt (type)" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
    
    SEXP x_levels = Rf_getAttrib(x, Rf_install("levels"));
    if (TYPEOF(x_levels) != STRSXP) {
        Rcout << "The column " << name << " is corrupt (levels)" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
    
    int n = Rf_length(x);
    DTMutableIntArray offsets(n);
    // int *intValues = INTEGER(x);
    memcpy(offsets.Pointer(),INTEGER(x),n*sizeof(int));
    int *intValues = offsets.Pointer();
    for (int i=0;i<n;i++) {
        if (intValues[i]==NA_INTEGER) {
            intValues[i] = -1;
        }
        else {
            intValues[i]--;
        }
        // Rcout << intValues[i] << endl;
    }
    
    // auto values = factorCodesToPrimitiveArray(x);
    // auto levels = chrToPrimitiveArray(x_levels);
    
    // bool ordered = Rf_inherits(x, "ordered");
    
    // Now extract the string objects
    DTCharArray utf8Buffer = UTF8BufferFrom(x_levels);
    
    return DTTableColumn::TextColumn(name,offsets,utf8Buffer);
}

void ConvertToDoubleArray(SEXP x,DTMutableDoubleArray &da,DTCharArray &m)
{
    double *data = REAL(x);
    int i,n = Rf_length(x);
    DTMutableDoubleArray values(n);
    DTMutableCharArray mask;
    
    bool saveMask = false;
    memcpy(values.Pointer(),data,sizeof(double)*n);
    for (i=0;i<n;i++) {
        if (R_IsNA(data[i])) {
            if (saveMask==false) {
                mask = DTMutableCharArray(n);
                mask = 1;
                saveMask = true;
            }
            mask(i) = 0;
        }
    }
    
    da = values;
    if (saveMask) m = mask;
}

DTTableColumn ConvertFromRealColumn(const std::string &name,SEXP x)
{
    DTMutableDoubleArray values;
    DTCharArray mask;
    ConvertToDoubleArray(x,values,mask);
    
    if (mask.Length()) {
        return DTTableColumn::NumberColumn(name,values,mask);
    }
    else {
        return DTTableColumn::NumberColumn(name,values);
    }
}

void ConvertToIntArray(SEXP x,DTMutableIntArray &ia,DTCharArray &m)
{
    int *data = INTEGER(x);
    int i,n = Rf_length(x);
    DTMutableIntArray values(n);
    DTMutableCharArray mask;
    
    bool saveMask = false;
    memcpy(values.Pointer(),data,sizeof(int)*n);
    for (i=0;i<n;i++) {
        if (data[i]==NA_INTEGER) {
            if (saveMask==false) {
                mask = DTMutableCharArray(n);
                mask = 1;
                saveMask = true;
            }
            mask(i) = 0;
        }
    }
    
    ia = values;
    if (saveMask) m = mask;
}

DTTableColumn ConvertFromIntegerColumn(const std::string &name,SEXP x)
{
    DTMutableIntArray values;
    DTCharArray mask;
    ConvertToIntArray(x,values,mask);
    
    if (mask.Length()) {
        return DTTableColumn::NumberColumn(name,values,mask);
    }
    else {
        return DTTableColumn::NumberColumn(name,values);
    }
}

DTTableColumn ConvertFromStringColumn(const std::string &name,SEXP x)
{
    DTCharArray utf8Buffer = UTF8BufferFrom(x);
    return DTTableColumn::TextColumn(name,DTStringList(utf8Buffer));
}

DTTableColumn ConvertFromDateColumn(const std::string &name,SEXP x)
{
    DTCharArray mask;
    
    if (TYPEOF(x)==REALSXP) {
        DTMutableDoubleArray da;
        ConvertToDoubleArray(x,da,mask);
        
        da *= 3600*24; // Since it is days from 1970
        if (mask.Length()) {
            return DTTableColumn::DateColumn(name,da,mask);
        }
        else {
            return DTTableColumn::DateColumn(name,da);
        }
    }
    else if (TYPEOF(x)==INTSXP) {
        DTMutableIntArray ia;
        ConvertToIntArray(x,ia,mask);
        DTMutableDoubleArray da = ConvertToDouble(ia);
        da *= 3600*24; // Since it is days from 1970
        
        if (mask.Length()) {
            return DTTableColumn::DateColumn(name,da,mask);
        }
        else {
            return DTTableColumn::DateColumn(name,da);
        }
    }
    else {
        Rcout << "The column " << name << " is not a properly saved date column" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
}

DTTableColumn ConvertFromTimeStampColumn(const std::string &name,SEXP x)
{
    DTCharArray mask;

    if (TYPEOF(x)==REALSXP) {
        DTMutableDoubleArray da;
        ConvertToDoubleArray(x,da,mask);
        
        if (mask.Length()) {
            return DTTableColumn::DateColumn(name,da,mask);
        }
        else {
            return DTTableColumn::DateColumn(name,da);
        }
    }
    else if (TYPEOF(x)==INTSXP) {
        DTMutableIntArray ia;
        ConvertToIntArray(x,ia,mask);
        DTMutableDoubleArray da = ConvertToDouble(ia);
        da *= 3600*24; // Since it is days from 1970
        
        if (mask.Length()) {
            return DTTableColumn::DateColumn(name,da,mask);
        }
        else {
            return DTTableColumn::DateColumn(name,da);
        }
    }
    else {
        Rcout << "The column " << name << " is a time stamp, but using an unexpected number format" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
}

DTTableColumn ConvertSingleColumn(const std::string &name,SEXP x)
{
    int typeOfX = TYPEOF(x);
    
    if (typeOfX==REALSXP) {
        return ConvertFromRealColumn(name,x);
    }
    else if (typeOfX==LGLSXP) {
        return ConvertFromLogicalColumn(name,x);
    }
    else if (typeOfX==INTSXP) {
        return ConvertFromIntegerColumn(name,x);
    }
    else if (typeOfX==STRSXP) {
        return ConvertFromStringColumn(name,x);
    }
    else if (typeOfX==RAWSXP) {
        Rcout << name << " : is a raw byte object that can't be saved.  Left blank." << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
    else {
        Rcout << name << " : can not be converted, please report (" << TYPEOF(x) << ") " << Rf_type2char(TYPEOF(x)) << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    }
}

bool ConvertToTableIfPossible(const std::string &name,SEXP x,DTTable &returnTable)
{
    // If this is a table, convert it and return true, otherwise return false
    // int typeOfX = TYPEOF(x);
    
    if (Rf_inherits(x,"ts")) {
        returnTable = ConvertFromTimeSeries(name,x);
        return true;
    }
    if (Rf_isMatrix(x)) {
        returnTable = ConvertFromMatrix(name,x);
        return true;
    }
    if (Rf_inherits(x,"table")) {
        Rcout << "Column " << name << " is a table.  Not supported at this time. Please report" << endl;
        returnTable = DTTable();
        return true;
    }
    
    if (Rf_isFrame(x)) {
        returnTable = ConvertToTable(DataFrame(x));
        return true;
    }
    
    SEXP nameList = Rf_getAttrib(x,Rf_install("names"));

    if (Rf_isNumeric(x) && nameList && TYPEOF(nameList)==STRSXP) {
        // Since names are defined, need to extract two columns and put it into a table
        
        // Rcout << name << " has string names" << endl;
        // See if I have a names column specified
        DTMutableList<DTTableColumn> twoColumns(2);
        twoColumns(0) = ConvertSingleColumn("name",nameList);
        twoColumns(1) = ConvertSingleColumn("value",x);
        returnTable = DTTable(twoColumns);
        return true;
    }
    
    // Not a table
    return false;
}

DTTableColumn ConvertToColumn(const std::string &name,SEXP x)
{
    int typeOfX = TYPEOF(x);
    
    // Might be a table
    DTTable tableVersion;
    if (ConvertToTableIfPossible(name,x,tableVersion)) {
        return DTTableColumn::TableColumn(name,tableVersion);
    }
    
    if (Rf_inherits(x, "factor")) {
        return ConvertFromIndexedStrings(name,x);
    } else if (Rf_inherits(x, "Date")) {
        return ConvertFromDateColumn(name,x);
        // Rcout << "Date" << endl;
    } else if (Rf_inherits(x, "time") || Rf_inherits(x, "hms")) {
        Rcout << "Haven't implemented time or hms format yet, please report" << endl;
        return DTTableColumn(name);
        // Rcout << "time" << endl;
    } else if (Rf_inherits(x, "POSIXct")) {
        return ConvertFromTimeStampColumn(name,x);
    } else if (Rf_inherits(x, "POSIXlt")) {
        Rcout << name << " : Is a POSIXlt array that needs to be converted to POSIXct, saving a blank column at this time" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    } else if (Rf_inherits(x, "dist")) {
        Rcout << "Can't save a dist class yet.  Not clear what it should map to in DataGraph" << endl;
        return DTTableColumn::NumberColumn(name,DTDoubleArray());
    } else {
        if (typeOfX==19) {
            if (Rf_isFrame(x)) {
                DTTable subTable = ConvertToTable(DataFrame(x));
                return DTTableColumn::TableColumn(name,subTable);
                // Rcout << name << " is a frame" << endl;
            }
            else if (Rf_isList(x)) {
                Rcout << name << " : is a list, not supported yet, saving an empty column" << endl;
                return DTTableColumn::NumberColumn(name,DTDoubleArray());
            }
            else {
                Rcout << name << " : has an unknown type saving an empty column" << endl;
                return DTTableColumn::NumberColumn(name,DTDoubleArray());
            }
            /*
            List list(x);
            std::string st = list.attr() ;
            Rcout << "attr = " << st << endl;
             */
            /*
            DataFrame df(x);
            CharacterVector names = df.names();
            int howManyColumns = df.size();
            Rcout << name << " is a list : ";
            for (int i=0;i<howManyColumns;i++) {
                Rcout << std::string(names[i]) << " ";
            }
            Rcout << endl;
             */
        }
        else {
            return ConvertSingleColumn(name,x);
        }
    }
    
    return DTTableColumn::NumberColumn(name,DTDoubleArray());
}


/*
 #define NILSXP       0    // nil = NULL
 #define SYMSXP       1    // symbols
 #define LISTSXP      2    // lists of dotted pairs
 #define CLOSXP       3    // closures
 #define ENVSXP       4    // environments
 #define PROMSXP      5    // promises: [un]evaluated closure arguments
 #define LANGSXP      6    // language constructs (special lists)
 #define SPECIALSXP   7    // special forms
 #define BUILTINSXP   8    // builtin non-special forms
 #define CHARSXP      9    // "scalar" string type (internal only)
 #define LGLSXP      10    // logical vectors
 // 11 and 12 were factors and ordered factors in the 1990s
 #define INTSXP      13    // integer vectors
 #define REALSXP     14    // real variables
 #define CPLXSXP     15    // complex variables
 #define STRSXP      16    // string vectors
 #define DOTSXP      17    // dot-dot-dot object
 #define ANYSXP      18    // make "any" args work.
 Used in specifying types for symbol
 registration to mean anything is okay
 #define VECSXP      19    // generic vectors
 #define EXPRSXP     20    // expressions vectors
 #define BCODESXP    21    // byte code
 #define EXTPTRSXP   22    // external pointer
 #define WEAKREFSXP  23    // weak reference
 #define RAWSXP      24    // raw bytes
 #define S4SXP       25    // S4, non-vector
 
 // used for detecting PROTECT issues in memory.c
 #define NEWSXP      30    // fresh node creaed in new page
 #define FREESXP     31    // node released by GC
 
 #define FUNSXP      99    // Closure or Builtin or Special
 */

/*
 
 library(devtools)
 install()
 library(DataGraph)

 
 openDTBin("/tmp/t")
 addDTBin("/tmp/t","iris",iris)
 closeDTBin("/tmp/t")

 */

DTTable ConvertToTable(DataFrame df)
{
    DTTableColumn rowNameColumn;
    SEXP rowNames = Rf_getAttrib(df, Rf_install("row.names"));
    if (TYPEOF(rowNames)==STRSXP) {
        rowNameColumn = ConvertToColumn("row.names",rowNames);
    }
    // Rcout << "type = " << Rf_type2char(TYPEOF(rown)) << endl;

    CharacterVector names = df.names();
    std::string name;
    int howManyColumnsInFrame = df.size();
    int howManyColumns = howManyColumnsInFrame;
    if (rowNameColumn.NotEmpty()) {
        howManyColumns++;
    }
    DTMutableList<DTTableColumn> columns(howManyColumns);
    int posInColumns = 0;
    if (rowNameColumn.NotEmpty()) {
        columns(posInColumns++) = rowNameColumn;
    }

    
    for (int i=0;i<howManyColumnsInFrame;i++) {
        name = std::string(names[i]);
        columns(posInColumns++) = ConvertToColumn(name,df[i]);
    }
    
    DTTable theTable(columns);
    
    return theTable;
}

int ComputeDaysSinceJan1st1970(int year,int month,int day)
{
    int m = (month + 9)%12;                /* mar=0, feb=11 */
    int y = year - m/10;                     /* if Jan/Feb, year-- */
    return (y*365 + y/4 - y/100 + y/400 + (m*306 + 5)/10 + (day - 1)) - 719468;
}

DTTable ConvertFromMatrix(const std::string &name,SEXP x)
{
    DTMutableCharArray mask;
    
    // Might be a multivariate time series
    int rowsInMatrix = 0;
    int columnsInMatrix = 0;
    
    SEXP dims = Rf_getAttrib(x, Rf_install("dim"));
    DTMutableIntArray ia;
    ConvertToIntArray(dims,ia,mask);
    if (ia.Length()!=2) {
        // I don't think this will happen, but just in case
        Rcout << "Only support a two dimensional matrices (" << name <<")" << endl;
        return DTTable();
    }
    rowsInMatrix = ia(0);
    columnsInMatrix = ia(1);
    
    // SEXP dimnames = Rf_getAttrib(x,Rf_install("dimnames"));
    // Require the rows and column names to be specified
    SEXP dimN = Rf_getAttrib(x,Rf_install("dimnames"));
    if (TYPEOF(dimN)!=VECSXP) {
        // Don't want to save a matrix that doesn't have any labels.
        Rcout << "Can only save matrices that have dimension names defined.  The entry " << name << " will be saved as an empty table (" << TYPEOF(dimN) << ")" << endl;
        return DTTable();
    }
    List dimNames(dimN);
    
    SEXP rowLabels = dimNames[0];
    DTCharArray utf8Buffer = UTF8BufferFrom(rowLabels);
    DTStringList rowNameList = DTStringList(utf8Buffer);
    
    SEXP labels = dimNames[1];
    utf8Buffer = UTF8BufferFrom(labels);
    DTStringList columnNameList = DTStringList(utf8Buffer);
    
    int typeOfX = TYPEOF(x);
    
    DTTableColumn col;
    if (typeOfX==REALSXP) {
        col = ConvertFromRealColumn("value",x);
    }
    else if (typeOfX==LGLSXP) {
        col = ConvertFromLogicalColumn("value",x);
    }
    else if (typeOfX==INTSXP) {
        col = ConvertFromIntegerColumn("value",x);
    }
    else if (typeOfX==STRSXP) {
        col = ConvertFromStringColumn("value",x);
    }
    else {
        Rcout << "Could not convert the variable " << name << " into a table.  Please report" << endl;
        return DTTable();
    }

    // Chop the matrix up
    DTMutableList<DTTableColumn> columns;
    if (columnsInMatrix==0) {
        columns = DTMutableList<DTTableColumn>(1);
        columns(0) = col;
    }
    else {
        columns = DTMutableList<DTTableColumn>(columnsInMatrix+1);
        int i;
        columns(0) = DTTableColumn::TextColumn("row names",rowNameList);
        for (i=0;i<columnsInMatrix;i++) {
            columns(i+1) = col.ExtractRows(DTRange(i*rowsInMatrix,rowsInMatrix)).ChangeName(columnNameList(i));
        }
    }
    
    DTTable tableToReturn(columns);
    
    return tableToReturn;
}

DTTable ConvertFromTimeSeries(const std::string &name,SEXP x)
{
    DTMutableCharArray mask;

    // Might be a multivariate time series
    int rowsInMatrix = 0;
    int columnsInMatrix = 0;
    
    DTStringList columnNameList;
    
    if (Rf_inherits(x,"matrix")) {
        SEXP dims = Rf_getAttrib(x, Rf_install("dim"));
        DTMutableIntArray ia;
        ConvertToIntArray(dims,ia,mask);
        if (ia.Length()!=2) {
            Rcout << "Only support a two dimensional time series matrix (" << name <<")" << endl;
            return DTTable();
        }
        rowsInMatrix = ia(0);
        columnsInMatrix = ia(1);
        
        // SEXP dimnames = Rf_getAttrib(x,Rf_install("dimnames"));
        List dimNames(Rf_getAttrib(x,Rf_install("dimnames")));
        SEXP labels = dimNames[1];
        
        DTCharArray utf8Buffer = UTF8BufferFrom(labels);
        columnNameList = DTStringList(utf8Buffer);
    }

    int typeOfX = TYPEOF(x);

    DTTableColumn col;
    if (typeOfX==REALSXP) {
        col = ConvertFromRealColumn("value",x);
    }
    else if (typeOfX==LGLSXP) {
        col = ConvertFromLogicalColumn("value",x);
    }
    else if (typeOfX==INTSXP) {
        col = ConvertFromIntegerColumn("value",x);
    }
    else if (typeOfX==STRSXP) {
        col = ConvertFromStringColumn("value",x);
    }
    else {
        Rcout << "The column is a time column, but can't be converted" << endl;
        return DTTable();
    }

    // The time values are saved
    SEXP info = Rf_getAttrib(x, Rf_install("tsp"));
    DTMutableDoubleArray da;
    ConvertToDoubleArray(info,da,mask);
    if (da.Length()!=3) {
        Rcout << "The column is a time column, but can't be converted" << endl;
        return DTTable();
    }
    
    double start = da(0);
    // double end = da(1);
    double howMany = da(2);
    double stride = 1.0/howMany;
    
    int howManyRows = col.NumberOfRows();
    if (rowsInMatrix==0)
        rowsInMatrix = howManyRows;

    DTMutableDoubleArray timeValues(rowsInMatrix);
    timeValues = NAN;
    double t;
    int i,year;
    if (howMany==12 || howMany==4 || howMany==2 || howMany==6 || howMany==3) {
        // This is divisible by 12
        int month;
        for (i=0;i<rowsInMatrix;i++) {
            t = start + stride*i+1e-5;
            year = int(floor(t));
            month = int(round((t-floor(t))*12)+1);
            timeValues(i) = ComputeDaysSinceJan1st1970(year,month,1)*86400;
        }
    }
    else if (howMany==1) {
        for (i=0;i<rowsInMatrix;i++) {
            t = start + i+1e-5;
            year = int(round(t));
            timeValues(i) = ComputeDaysSinceJan1st1970(year,1,1)*86400;
        }
    }
    else {
        double fraction;
        double startDate = NAN, endDate = NAN;
        int yearForStart = -1122352;
        for (i=0;i<rowsInMatrix;i++) {
            t = start + stride*i+1e-5;
            year = int(floor(t));
            fraction = t-floor(t);
            // start + fraction*(end-start) = start*(1-fraction) + end*fraction
            if (year!=yearForStart) {
                startDate = ComputeDaysSinceJan1st1970(year,1,1)*86400;
                endDate =  ComputeDaysSinceJan1st1970(year+1,1,1)*86400;
                yearForStart = year;
            }
            timeValues(i) = startDate*(1-fraction) + endDate*fraction;
        }
    }
    DTTableColumn dateColumn = DTTableColumn::DateColumn("date",timeValues,DTCharArray());
    
    DTMutableList<DTTableColumn> columns;
    if (columnsInMatrix==0) {
        columns = DTMutableList<DTTableColumn>(2);
        columns(0) = dateColumn;
        columns(1) = col;
    }
    else {
        columns = DTMutableList<DTTableColumn>(columnsInMatrix+1);
        columns(0) = dateColumn;
        for (i=0;i<columnsInMatrix;i++) {
            columns(i+1) = col.ExtractRows(DTRange(i*rowsInMatrix,rowsInMatrix)).ChangeName(columnNameList(i));
        }
    }
    
    DTTable tableToReturn(columns);
    
    return tableToReturn;
}

#pragma mark Drivering routines

// [[Rcpp::export]]
void writeDTable(const std::string& path,SEXP data)
{
    // Write a single table. If this is not a table complain to the user and don't write anything.
    DTTable theTable;
    if (ConvertToTableIfPossible("Input",data,theTable)==false) {
        Rcout << "The input argument is not a table.";
        return;
    }
    
    // make sure that the path ends with a dtable
    std::string pathStd = path;
    std::string::size_type location = path.find_last_of(".");
    if (location==std::string::npos || path.substr(location+1)!="dtable") {
        pathStd = path+".dtable";
    }
    
    DTDataFile outputFile(pathStd,DTFile::NewReadWrite);
    
    
    WriteOne(outputFile,"Var",theTable);
    
    outputFile.SaveIndex();
}

// **************************************************************************************************************
#pragma mark DataTable

struct DGGlobalTableStorage
{
    // DGGlobalTableStorage(const DTDataFile &df) : howManySaved(0), dataFile(df) {}
    
    DTMutableDictionary information;
    DTTableStructure tableStructure;
    DTDataFile dataFile;
};

struct DGGlobalStorageClass
{
    // std::map<std::string,DTDataFile> binFiles;
    DTMutableMap<DGGlobalTableStorage> tableFiles;
    // std::map<std::string,DTPointer<DGGlobalTableStorage> > tableFiles;
};

static DGGlobalStorageClass *DGGlobalStorage = NULL;

// [[Rcpp::export]]
void openDTable(const std::string& path)
{
    if (DGGlobalStorage==NULL) {
        DGGlobalStorage = new DGGlobalStorageClass();
    }
    
    std::string pathStd = path;
    std::string::size_type location = pathStd.find_last_of(".");
    if (location==std::string::npos || pathStd.substr(location+1)!="dtable") {
        pathStd = pathStd+".dtable";
    }
    
    DTDataFile dataFile;
    if (DGGlobalStorage->tableFiles.Contains(pathStd)) {
        Rcout << "The file " << pathStd << " is already open" << endl;
        return;
    }

    DGGlobalTableStorage storage;
    storage.dataFile = DTDataFile(pathStd,DTFile::NewReadWrite);
    storage.dataFile.SaveIndex();
    storage.information("Count") = 0;
    DGGlobalStorage->tableFiles(pathStd) = storage;
}

// [[Rcpp::export]]
void syncDTable(const std::string& path)
{
    std::string pathStd = path;
    std::string::size_type location = pathStd.find_last_of(".");
    if (location==std::string::npos || pathStd.substr(location+1)!="dtable") {
        pathStd = pathStd+".dtable";
    }
    
    if (DGGlobalStorage==NULL) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    if (DGGlobalStorage->tableFiles.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    DGGlobalTableStorage &storage = DGGlobalStorage->tableFiles(pathStd);
    storage.dataFile.Sync();
}

// [[Rcpp::export]]
void closeDTable(const std::string& path)
{
    std::string pathStd = path;
    std::string::size_type location = pathStd.find_last_of(".");
    if (location==std::string::npos || pathStd.substr(location+1)!="dtable") {
        pathStd = pathStd+".dtable";
    }
    
    if (DGGlobalStorage==NULL) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    if (DGGlobalStorage->tableFiles.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    DGGlobalStorage->tableFiles.Erase(pathStd);
}

// [[Rcpp::export]]
void addDTable(const std::string& path,SEXP data)
{
    // make sure that the path ends with a dtable
    std::string pathStd = path;
    std::string::size_type location = path.find_last_of(".");
    if (location==std::string::npos || path.substr(location+1)!="dtable") {
        pathStd = pathStd+".dtable";
    }
    
    if (DGGlobalStorage==NULL) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    if (DGGlobalStorage->tableFiles.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    DGGlobalTableStorage &storage = DGGlobalStorage->tableFiles(pathStd);
    
    DTTable theTable;
    if (ConvertToTableIfPossible("Input",data,theTable)==false) {
        Rcout << "The input argument is not a table." << endl;
        return;
    }
    
    int howManySaved = storage.information("Count");
    if (howManySaved==0) {
        theTable.WriteStructure(storage.dataFile,"Var");
        storage.dataFile.Save("Table","Seq_Var");
        storage.tableStructure = theTable.Structure();
    }
    else {
        if (theTable.Structure()!=storage.tableStructure) {
            Rcout << "All tables have to have the same structure as the first table." << endl;
            return;
        }
    }
    
    Write(storage.dataFile,"Var_" + DTInt2String(howManySaved),theTable);
    howManySaved++;
    storage.information("Count") = howManySaved;
}

// **************************************************************************************************************
#pragma mark DTBin support

struct DGVariableInfo
{
    DGVariableInfo() : howManySaved(0), lastTimeValue(-1) {}
    
    std::string type;
    DTTableStructure tableStructure;
    int howManySaved;
    double lastTimeValue;
};

struct DGGlobalDTBinStorage
{
    DTMutableMap<DGVariableInfo> variableInfo;

    DTDataFile dataFile;
};

struct DGGlobalStorageClassDTBin
{
    // std::map<std::string,DTDataFile> binFiles;
    DTMutableMap<DGGlobalDTBinStorage> map;
    // std::map<std::string,DTPointer<DGGlobalDTBinStorage> > dtbinFiles;
};

static DGGlobalStorageClassDTBin *DGGlobalStorageDTBin = NULL;

extern std::string StandardizeDTBinName(const std::string &path);

extern std::string StandardizeDTBinName(const std::string &pathIn)
{
    std::string path = pathIn;
    std::string::size_type location = path.find_last_of(".");
    if (location==std::string::npos || path.substr(location+1)!="dtbin") {
        path = path+".dtbin";
    }
    return path;
}

// [[Rcpp::export]]
void openDTBin(const std::string& path)
{
    if (DGGlobalStorageDTBin==NULL) {
        DGGlobalStorageDTBin = new DGGlobalStorageClassDTBin();
    }
    
    std::string pathStd = StandardizeDTBinName(path);
    
    DTDataFile dataFile;
    if (DGGlobalStorageDTBin->map.Contains(pathStd)) {
        Rcout << "The file " << pathStd << " is already open" << endl;
        return;
    }
    
    DGGlobalDTBinStorage addThis;
    addThis.dataFile = DTDataFile(pathStd,DTFile::NewReadWrite);
    addThis.dataFile.SaveIndex();
    
    DGGlobalStorageDTBin->map(pathStd) = addThis;
}

// [[Rcpp::export]]
void closeDTBin(const std::string &path)
{
    std::string pathStd = StandardizeDTBinName(path);
    
    if (DGGlobalStorageDTBin==NULL || DGGlobalStorageDTBin->map.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    DGGlobalStorageDTBin->map.Erase(pathStd);
}

// [[Rcpp::export]]
void syncDTBin(const std::string &path)
{
    std::string pathStd = StandardizeDTBinName(path);
    
    if (DGGlobalStorageDTBin==NULL || DGGlobalStorageDTBin->map.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened" << endl;
        return;
    }
    
    DGGlobalDTBinStorage &information = DGGlobalStorageDTBin->map(pathStd);
    information.dataFile.Sync();
}

// [[Rcpp::export]]
void infoDTBin(const std::string &path="")
{
    int i, howMany;
    std::string padding = ". . . . . . . . . . . . . . . . . . . . . . . . . . .";
    std::string spaces = "                ";
    

    if (path.length()==0) {
        if (DGGlobalStorageDTBin==NULL) {
            Rcout << "No file is open" << endl;
            return;
        }
        // Display a list of files that are open
        DTList<std::string> keys = DGGlobalStorageDTBin->map.Keys();
        howMany = keys.Length();
        std::string pathName;
        
        Rcout << "***********************************************************" << endl;
        Rcout << "File name                              # of entries        " << endl;
        Rcout << "-----------------------------------------------------------" << endl;

        for (i=0;i<howMany;i++) {
            pathName = keys(i);
            std::string desc = pathName + " ";
            int len = desc.length();
            // If it is too long, trim from the start
            if (len>40) {
                desc = string(desc,len-37,len);
                desc = "..."+desc;
                len = desc.length();
            }
            if (len%2) {
                desc = desc + " ";
                len++;
            }
            if (len<46) desc = desc + string(padding,0,46-len);
            if (len>46) desc = string(desc,0,46);
            Rcout << desc;

            // Display how many entries are stored
            DGGlobalDTBinStorage &information = DGGlobalStorageDTBin->map(pathName);
            int howManyEntries = information.variableInfo.NumberOfKeys();
            if (howManyEntries==0) {
                Rcout << "empty" << endl;
            }
            else if (howManyEntries==1) {
                Rcout << "1 variable" << endl;
            }
            else {
                Rcout << howManyEntries << " variables" << endl;
            }
        }
        Rcout << "***********************************************************" << endl;
        return;
    }
    
    std::string pathStd = StandardizeDTBinName(path);
    
    if (DGGlobalStorageDTBin==NULL || DGGlobalStorageDTBin->map.Contains(pathStd)==false) {
        Rcout << "The file " << pathStd << " has not been opened";
        return;
    }
    
    DGGlobalDTBinStorage &information = DGGlobalStorageDTBin->map(pathStd);
    DTList<std::string> keys = information.variableInfo.Keys();
    howMany = keys.Length();
    
    Rcout << "***********************************************************" << endl;
    Rcout << "File : " << pathStd << endl;
    if (howMany==0)
        Rcout << "No variables" << endl;
    else if (howMany==1)
        Rcout << "One variable" << endl;
    else
        Rcout << howMany << " variables" << endl;
    Rcout << "***********************************************************" << endl;
    Rcout << "Name                           Type              # times   " << endl;
    Rcout << "-----------------------------------------------------------" << endl;
    for (i=0;i<howMany;i++) {
        DGVariableInfo &info = information.variableInfo(keys(i));
        std::string desc = keys(i) + " ";
        int len = desc.length();
        if (len%2) {
            desc = desc + " ";
            len++;
        }
        if (len<30) desc = desc + string(padding,0,30-len);
        if (len>30) desc = string(desc,0,30);
        Rcout << desc << " ";
        desc = info.type + " ";
        if (desc.length()<17) desc = desc + string(spaces,0,17-desc.length());
        if (desc.length()>17) desc = string(desc,0,17);
        Rcout << desc;
        Rcout << " " << info.howManySaved;
        Rcout << endl;
    }
    Rcout << "***********************************************************" << endl;
}

// Add a table
// [[Rcpp::export]]
void addDTBin(const std::string& path,const std::string &name,SEXP data,double time=NA_REAL)
{
    std::string pathStd = StandardizeDTBinName(path);
    if (DGGlobalStorageDTBin==NULL || DGGlobalStorageDTBin->map.Contains(pathStd)==false) {
        Rcout << "addDBin: The file " << pathStd << " has not been opened";
        return;
    }
    
    if (isnan(time)==false && time<0) {
        Rcout << "addDBin: The time value needs to be >= 0";
    }
    if (name.size()==0) {
        Rcout << "addDBin: name needs to be non-empty";
    }
    DGGlobalDTBinStorage &information = DGGlobalStorageDTBin->map(pathStd);
    
    // See if this should be converted to a table.  Various reasons for that,
    // data frames, time series and matrices that have named dimensions should be saved into a table.
    DTTable theTable;
    if (ConvertToTableIfPossible(name,data,theTable)) {
        if (information.variableInfo.Contains(name)) {
            // Argument validation
            if (information.variableInfo(name).type!="Table") {
                Rcout << "addDBin: The entry " << name << " inside " << pathStd << " needs to be a table";
                return;
            }
            if (information.variableInfo(name).tableStructure!=theTable.Structure()) {
                Rcout << "addDBin: The table " << name << " inside " << pathStd << " has a different structure";
                return;
            }
        }
        else {
            // The variable didn't exist, create it
            DGVariableInfo newInfo;
            newInfo.type = "Table";
            newInfo.tableStructure = theTable.Structure();
            
            theTable.WriteStructure(information.dataFile,name);
            information.dataFile.Save("Table","Seq_"+name);
            
            information.variableInfo(name) = newInfo;
        }
        
        DGVariableInfo &info = information.variableInfo(name);
        
        if (isnan(time)) {
            time = info.lastTimeValue+1.0;
        }
        else if (time<=info.lastTimeValue) {
            Rcout << "addDBin: The table " << name << "inside " << pathStd << "has a time value that is >= the one you specified";
            return;
        }
        
        Write(information.dataFile,name+"_" + DTInt2String(info.howManySaved),theTable);
        Write(information.dataFile,name+"_" + DTInt2String(info.howManySaved) + "_time",time);
        info.lastTimeValue = time;
        info.howManySaved++;
    }
    else {
        // Interpret it as a column
        DTTableColumn column = ConvertToColumn(name,data);
        std::string type;
        if (column.IsNumberColumn()) {
            type = "List of Numbers";
        }
        else if (column.IsTextColumn()) {
            type = "List of Strings";
        }
        else if (column.IsDateColumn()) {
            type = "List of Dates";
        }
        else {
            Rcout << "addDBin: The variable could not be converted : " << TYPEOF(data) << " type : " << Rf_type2char(TYPEOF(data)) << endl;
            return;
        }
        
        // Rcout << "type = " << type << endl;
        
        if (information.variableInfo.Contains(name)) {
            // Argument validation
            if (information.variableInfo(name).type!=type) {
                Rcout << "addDBin: The entry " << name << " inside " << pathStd << " has a different type";
                return;
            }
        }
        else {
            // The variable didn't exist, create it
            DGVariableInfo newInfo;
            newInfo.type = type;
            information.dataFile.Save(type,"Seq_"+name);
            information.variableInfo(name) = newInfo;
        }

        DGVariableInfo &info = information.variableInfo(name);
        
        if (isnan(time)) {
            time = info.lastTimeValue+1.0;
        }
        else if (time<=info.lastTimeValue) {
            Rcout << "addDBin: The table " << name << "inside " << pathStd << "has a time value that is >= the one you specified";
            return;
        }
        
        // Save name_0 for the values then name_0_T etc.
        column.WriteSingle(information.dataFile,name+"_" + DTInt2String(info.howManySaved));
        Write(information.dataFile,name+"_" + DTInt2String(info.howManySaved) + "_time",time);
        info.lastTimeValue = time;
        info.howManySaved++;
    }
}

