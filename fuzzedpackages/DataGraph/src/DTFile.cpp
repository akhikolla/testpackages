// Part of DTSource. Copyright 2004-2015. David A. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTFile.h"

#include "DTError.h"
#include <iostream>

#include <string>
#include <math.h>

#include "DTEndianSwap.h"
#include "DTDoubleArray.h"
#include "DTFloatArray.h"
#include "DTIntArray.h"
#include "DTCharArray.h"
#include "DTUCharArray.h"
#include "DTShortIntArray.h"
#include "DTUShortIntArray.h"
#include "DTUtilities.h"
#include <limits>

struct DTFileStorage
{
    DTFileStorage();
    ~DTFileStorage();
    
    std::string name;
    FILE *file;
    bool readOnly;
    int retainCount;
    DTFilePosition lengthOfFile;
    
    DTFile::Endian endian;
    bool swapEndian;
};

DTFileStorage::DTFileStorage()
{
    retainCount = 1;
    readOnly = true;
    lengthOfFile = -1;
    file = NULL;
    endian = DTFile::Native;
    swapEndian = false;
}

DTFileStorage::~DTFileStorage()
{
    if (file) fclose(file);
}

DTFile::DTFile()
{
    storage = new DTFileStorage();
}

DTFile::DTFile(const std::string &nm,OpenType openT)
{
    storage = new DTFileStorage();
    storage->name = nm;

    if (openT==DTFile::ReadOnly) {
        storage->file = fopen(nm.c_str(),"rb");
        storage->readOnly = true;
    }
    else if (openT==DTFile::ExistingReadWrite) {
        storage->file = fopen(nm.c_str(),"r+b");
        storage->readOnly = false;
    }
    else {
        // Should delete the file first, to preserve hard links.
        remove(nm.c_str());
        storage->file = fopen(nm.c_str(),"w+b");
        storage->readOnly = false;
    }

    if (storage->file==NULL) {
        std::string msg = "Could not open the file \"";
        msg = msg + nm + "\"";
        DTErrorMessage("DTFile(name,type)",msg);
    }
}

DTFile::DTFile(const std::string &nm,Endian endian,OpenType openT)
{
    storage = new DTFileStorage();
    storage->name = nm;
    storage->endian = endian;
    
    // figure out what the endian ness of this machine is
    if (RunningOnBigEndianMachine()) {
        storage->swapEndian = (endian==DTFile::LittleEndian);
    }
    else {
        storage->swapEndian = (endian==DTFile::BigEndian);
    }

    if (openT==DTFile::ReadOnly) {
        storage->file = fopen(nm.c_str(),"rb");
        storage->readOnly = true;
    }
    else if (openT==DTFile::ExistingReadWrite) {
        storage->file = fopen(nm.c_str(),"r+b");
        storage->readOnly = false;
    }
    else {
        // Should delete the file first, to preserve hard links.
        remove(storage->name.c_str());
        storage->file = fopen(nm.c_str(),"w+b");
        storage->readOnly = false;
    }

    if (storage->file==NULL) {
        std::string msg = "Could not open the file \"";
        msg = msg + nm + "\"";
        DTErrorMessage("DTFile(name,type)",msg);
    }
}

DTFile::~DTFile()
{
    if (--(storage->retainCount)==0) {
        delete storage;
    }
}

DTFile::DTFile(const DTFile &C)
{
    storage = C.storage;
    storage->retainCount++;
}

bool DTFile::RunningOnBigEndianMachine(void)
{
    int fourBytes[1];
    fourBytes[0] = 128912422;
    short int *asTwoShorts = (short int *)fourBytes;
    return (asTwoShorts[0]==1967);
}

DTFile::Endian DTFile::EndianForMachine(void)
{
    int fourBytes[1];
    fourBytes[0] = 128912422;
    short int *asTwoShorts = (short int *)fourBytes;
    return (asTwoShorts[0]==1967 ? DTFile::BigEndian : DTFile::LittleEndian);
}

DTFile &DTFile::operator=(const DTFile &C)
{
    if (storage!=C.storage) {
        if (--(storage->retainCount)==0) {
            delete storage;
        }
        storage = C.storage;
        storage->retainCount++;
    }
    
    return *this;
}

bool DTFile::CanOpen(const std::string &name,OpenType)
{
    bool toReturn = false;
    FILE *tempFile = NULL;
    
    tempFile = fopen(name.c_str(),"rb");
    toReturn = (tempFile!=NULL);
    if (tempFile) fclose(tempFile);
    
    return toReturn;
}

DTFilePosition DTFile::Position(void) const
{
    if (!storage->file) return 0;
#ifdef WIN32
    DTFilePosition toReturn = ftell(storage->file);
#else
    DTFilePosition toReturn = ftello(storage->file);
#endif
    return toReturn;
}

void DTFile::SetPosition(DTFilePosition pos) const
{
    if (!storage->file) return;
#ifdef WIN32
    fseek(storage->file,pos,SEEK_SET);
#else
    fseeko(storage->file,pos,SEEK_SET);
#endif
}

void DTFile::MovePosition(DTFilePosition pos) const
{
    if (!storage->file || pos==0) return;
#ifdef WIN32
    fseek(storage->file,pos,SEEK_CUR);
#else
    fseeko(storage->file,pos,SEEK_CUR);
#endif
}

void DTFile::MoveToEnd(void) const
{
    if (!storage->file) return;
#ifdef WIN32
    fseek(storage->file,0,SEEK_END);
#else
    fseeko(storage->file,0,SEEK_END);
#endif
}

bool DTFile::AtEndOfFile(void) const
{
    if (!storage->file) return true;
    return (feof(storage->file)!=0);
}

void DTFile::Flush(void) const
{
    if (storage->file) fflush(storage->file);
}

FILE *DTFile::FILEForReading(void) const
{
    return storage->file;
}

FILE *DTFile::FILEForWriting(void)
{
    storage->lengthOfFile = -1;
    return storage->file;
}

bool DTFile::IsOpen(void) const
{
    return (storage->file!=NULL);
}

bool DTFile::IsReadOnly(void) const
{
    return storage->readOnly;
}

string DTFile::Name(void) const
{
    return storage->name;
}

DTFile::Endian DTFile::EndianType(void) const
{
    return storage->endian;
}

DTFilePosition DTFile::Length(void) const
{
    if (storage->lengthOfFile>=0)
        return storage->lengthOfFile;
    
    DTFilePosition nowAt = Position();
    MoveToEnd();
    DTFilePosition toReturn = Position();
    SetPosition(nowAt);
    
    storage->lengthOfFile = toReturn;
    
    return toReturn;
}

bool DTFile::Find(char c) const
{
    // Find a specific character in a file.  If the character is found, place
    // the current read position at that character and return true.
    // Otherwise go back to the position that the file had at the start of the read and return false.
    DTFilePosition howLong = Length();
    DTFilePosition nowAt = Position();
    DTFilePosition startsAt = nowAt;
    
    DTMutableCharArray buffer(1024);
    DTFilePosition i,howMuchToRead;
    bool foundIt = false;
    while (1) {
        howMuchToRead = 1024;
        if (nowAt+howMuchToRead>howLong)
            howMuchToRead = howLong-nowAt;
        if (howMuchToRead==0)
            break;
        if (!ReadBinary(buffer,howMuchToRead))
            break;
        nowAt += howMuchToRead;
        for (i=0;i<howMuchToRead;i++) {
            if (buffer(i)==c) {
                foundIt = true;
                MovePosition(i-howMuchToRead);
                break;
            }
        }
        if (i<howMuchToRead)
            break;
    }
    
    if (foundIt==false)
        SetPosition(startsAt);
    
    return foundIt;
}

string DTFile::ReadLine(ssize_t maxLen) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadLine()","No file");
        return std::string();
    }

    FILE *theFile = FILEForReading();

    // Read until we hit \n or \r.  Swallow that character.
    DTMutableCharArray buffer(80);
    int temp;
    ssize_t locInBuffer = 0;
    while ((temp = fgetc(theFile))!=EOF && (maxLen<0 || locInBuffer < maxLen)) {
        if (temp=='\n' || temp=='\r' || temp=='\0')
            break;
        if (locInBuffer==buffer.Length()-1)
            buffer = IncreaseSize(buffer,buffer.Length());
        buffer(locInBuffer++) = (char)temp;
    }
    
    // Check if this is a windows text file, where every line ends with a return+new line.
    if (temp=='\r') {
        temp = fgetc(theFile);
        if (temp!='\n') {
            // Step one step back.
            fseek(storage->file,-1,SEEK_CUR);
        }
    }

    buffer(locInBuffer) = '\0';

    std::string toReturn(buffer.Pointer());

    return toReturn;
}

string DTFile::ReadString(size_t length) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadString(length)","No file");
        return std::string();
    }

    if (length<=0)
        return std::string();

    DTMutableCharArray A(length);
    if (fread(A.Pointer(),1,A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadString(length)",
                       "Could not read the required number of characters from the file");
        return std::string();
    }

    std::string toReturn(A.Pointer(),length);

    return toReturn;
}

string DTFile::NextWord(void) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::NextWord()","No file");
        return std::string();
    }
    
    FILE *theFile = FILEForReading();
    
    // Read until we hit \n or \r.  Swallow that character.
    DTMutableCharArray buffer(80);
    int temp;
    size_t locInBuffer = 0;
    size_t lenBuffer = buffer.Length();
    bool atStart = true;
    while ((temp = fgetc(theFile))!=EOF) {
        // Skip over any leading spaces
        if (temp==32 && atStart) continue;
        atStart = false;
        if (temp<33 || temp>126)
            break;
        if (locInBuffer==lenBuffer-1) {
            buffer = IncreaseSize(buffer,buffer.Length());
            lenBuffer = buffer.Length();
        }
        buffer(locInBuffer++) = char(temp);
    }
    
    // Check if this is a windows text file, where every line ends with a return+new line.
    if (temp=='\r') {
        temp = fgetc(theFile);
        if (temp!='\n') {
            // Step one step back.
            fseek(storage->file,-1,SEEK_CUR);
        }
    }
    
    buffer(locInBuffer) = '\0';
    
    std::string toReturn(buffer.Pointer());
    
    return toReturn;
}

double DTFile::ReadNumberWithSeparator(char c,bool &empty,bool &endOfLine) const
{
    endOfLine = false;
    empty = true;
    if (!IsOpen()) {
        DTErrorMessage("DTFile::NextWord()","No file");
        return NAN;
    }
    
    FILE *theFile = FILEForReading();
    
    int temp = EOF;
    // Skip over spaces in the beginning
    while ((temp = fgetc(theFile))!=EOF) {
        if (temp!=' ') break;
    }
    
    if (temp==EOF) {
        return NAN;
    }
    else if (temp==c) {
        return NAN;
    }
    else if (temp=='\r' || temp=='\n') {
        endOfLine = true;
        return NAN;
    }
    
    if (!(temp=='-' || temp=='+' || (temp>='0' && temp<='9') || temp=='.')) {
        // Not a valid number.
        return NAN;
    }
    
    // This is likely a valid number, read from the file until
    // we hit a value that is not a part of a number.
    char buffer[41];
    int posInBuffer = 0;
    
    buffer[posInBuffer++] = (char)temp;
    
    // Read into the buffer until I get an ending
    bool periodAllowed = (temp!='.');
    bool validNumber = true;
    bool inExponent = false;
    bool signInExponentSeen = false;
    while ((temp = fgetc(theFile))!=EOF && posInBuffer<40) {
        if (temp==c || temp==' ' || temp=='\n' || temp=='\r') {
            // Terminate the number
            break;
        }
        if (temp=='.') {
            if (periodAllowed==false) {
                // Not a valid number
                validNumber = false;
                break;
            }
            else {
                buffer[posInBuffer++] = char(temp);
                periodAllowed = false;
            }
        }
        else if (temp>='0' && temp<='9') {
            buffer[posInBuffer++] = char(temp);
        }
        else if (temp=='e' || temp=='E') {
            if (inExponent) {
                validNumber = false; // Can't start two exponents
                break;
            }
            // Start of the exponent
            buffer[posInBuffer++] = char(temp);
            inExponent = true;
            periodAllowed = false;
        }
        else if (temp=='-' || temp=='+') {
            if (inExponent && signInExponentSeen==false) {
                buffer[posInBuffer++] = char(temp);
                signInExponentSeen = true;
            }
            else {
                validNumber = false;
                break;
            }
        }
        else {
            validNumber = false;
            break;
        }
    }
    
    if (validNumber==false) {
        return NAN;
    }
    
    if (temp==' ' && c!=' ') {
        // Absorb spaces until the separator or end of line
        while ((temp = fgetc(theFile))!=EOF && temp==' ') {}
        if (temp!=EOF && temp!='\n' && temp!='\r' && temp!=c) {
            // Not a valid number, likely something like "5A" which should not be accepted
            return NAN;
        }
    }
    
    if (temp=='\r' || temp=='\n') {
        // Absorb the new line
        endOfLine = true;
    }
    
    buffer[posInBuffer] = 0;
    
    char *end;
    double toReturn = strtod(buffer, &end);
    if (end!=buffer+posInBuffer) {
        // Neeed to be able to understand this number all the way to the end
        return NAN;
    }
    
    empty = false;
    
    return toReturn;
}

std::string DTFile::ReadStringWithSeparator(char c,bool &empty,bool &endOfLine) const
{
    endOfLine = false;
    empty = true;
    if (!IsOpen()) {
        DTErrorMessage("DTFile::NextWord()","No file");
        return std::string();
    }
    
    FILE *theFile = FILEForReading();
    
    int temp = EOF;
    // Skip over spaces in the beginning
    while ((temp = fgetc(theFile))!=EOF) {
        if (temp!=' ') break;
    }
    
    if (temp==EOF) {
        return std::string();
    }
    else if (temp==c) {
        empty = false;
        return std::string();
    }
    else if (temp=='\r' || temp=='\n') {
        endOfLine = true;
        return std::string();
    }
    
    // This is likely a valid number, read from the file until
    // we hit a value that is not a part of a number.
    int posInBuffer = 0;
    ssize_t lengthOfBuffer = 40;
    
    DTMutableCharArray buffer(lengthOfBuffer);
    bool insideQuote = (temp=='"');
    if (insideQuote==false) {
        // Otherwise should ignore it
        buffer(posInBuffer++) = char(temp);
    }
    
    // Read into the buffer until I get an ending
    while ((temp = fgetc(theFile))!=EOF) {
        if (insideQuote) {
            if (temp=='"') {
                // End of the string, otherwise just continue
                break;
            }
        }
        else if (temp==c || temp==' ' || temp=='\n' || temp=='\r') {
            // Terminate
            break;
        }
        if (posInBuffer==lengthOfBuffer) {
            buffer = IncreaseSize(buffer,lengthOfBuffer);
            lengthOfBuffer = buffer.Length();
        }
        buffer(posInBuffer++) = char(temp);
    }
    
    if (temp==' ' && c!=' ') {
        // Absorb spaces until the separator or end of line
        while ((temp = fgetc(theFile))!=EOF && temp==' ') {}
        if (temp!=EOF && temp!='\n' && temp!='\r' && temp!=c) {
            return std::string();
        }
    }
    
    if (temp=='\r' || temp=='\n') {
        // Absorb the new line
        endOfLine = true;
    }
    
    empty = false;
    
    if (posInBuffer==0) {
        // This is the empty string, which is OK
        return std::string();
    }
    std::string toReturn(buffer.Pointer(),posInBuffer);
    
    return toReturn;
}

bool DTFile::ReadCharacters(char *chars,size_t howMuchToRead) const
{
    if (fread(chars,1,howMuchToRead,FILEForReading())!=howMuchToRead) {
        DTErrorMessage("DTFile::ReadCharacters(chars,howMuchToRead)","Could not read the required number of values from the file");
        return false;
    }
    return true;
}

#pragma mark ---------- Read in a single number

unsigned short int DTFile::ReadUnsignedShort() const
{
    unsigned short int toReturn = 0;
    
    if (!IsOpen())
        DTErrorMessage("DTFile::ReadUnsignedShort()","No file");
    else if (fread(&toReturn,sizeof(unsigned short int),1,FILEForReading())!=1)
        DTErrorMessage("DTFile::ReadUnsignedShort()","Could not read the number");
    
    // Endian swap
    
    return toReturn;
}

int DTFile::Read_int32(Endian endian) const
{
    int toReturn = 0;

    if (!IsOpen())
        DTErrorMessage("DTFile::Read_int32(Endian)","No file");
    else if (fread(&toReturn,4,1,FILEForReading())!=1)
        DTErrorMessage("DTFile::Read_int32(Endian)","Could not read the number");
    else {
        if (endian!=DTFile::EndianForMachine()) DTSwap4Bytes((unsigned char *)&toReturn,4); // Swap
    }

    return toReturn;
}

ssize_t DTFile::Read_int64(Endian endian) const
{
    ssize_t toReturn = 0;
    
    if (!IsOpen())
        DTErrorMessage("DTFile::Read_int64(Endian)","No file");
    else if (fread(&toReturn,8,1,FILEForReading())!=1)
        DTErrorMessage("DTFile::Read_int64(Endian)","Could not read the number");
    else {
        if (endian!=DTFile::EndianForMachine()) DTSwap8Bytes((unsigned char *)&toReturn,8); // Swap
    }
    
    return toReturn;
}

unsigned int DTFile::Read_uint32(Endian endian) const
{
    unsigned int toReturn = 0;

    if (!IsOpen())
        DTErrorMessage("DTFile::Read_uint32(Endian)","No file");
    else if (fread(&toReturn,4,1,FILEForReading())!=1)
        DTErrorMessage("DTFile::Read_uint32(Endian)","Could not read the number");
    else {
        if (endian!=DTFile::EndianForMachine()) DTSwap4Bytes((unsigned char *)&toReturn,4); // Swap
    }

    return toReturn;
}

float DTFile::Read_float(Endian endian) const
{
    float toReturn = 0;

    if (!IsOpen())
        DTErrorMessage("DTFile::Read_float(Endian)","No file");
    else if (fread(&toReturn,4,1,FILEForReading())!=1)
        DTErrorMessage("DTFile::Read_float(Endian)","Could not read the number");
    else {
        if (endian!=DTFile::EndianForMachine()) DTSwap4Bytes((unsigned char *)&toReturn,4); // Swap
    }

    return toReturn;
}

double DTFile::Read_double(Endian endian) const
{
    double toReturn = 0;

    if (!IsOpen())
        DTErrorMessage("DTFile::Read_double(Endian)","No file");
    else if (fread(&toReturn,8,1,FILEForReading())!=1)
        DTErrorMessage("DTFile::Read_double(Endian)","Could not read the number");
    else {
        if (endian!=DTFile::EndianForMachine()) DTSwap8Bytes((unsigned char *)&toReturn,8); // Swap
    }

    return toReturn;
}

#pragma mark ---------- Read in an array

bool DTFile::ReadBinary(DTMutableDoubleArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(DoubleArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (fread(A.Pointer(),sizeof(double),A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(DoubleArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableDoubleArray &A,const DTRange &range) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(DoubleArray,Range)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (range.end()>A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(DoubleArray,Range)","Range out of bounds");
        return false;
    }

    if (ssize_t(fread(A.Pointer()+range.start,sizeof(double),range.length,FILEForReading()))!=range.length) {
        DTErrorMessage("DTFile::ReadBinary(DoubleArray,Range)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableFloatArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(FloatArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (fread(A.Pointer(),sizeof(float),A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(FloatArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableIntArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(IntArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (fread(A.Pointer(),sizeof(int),A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(IntArray)",
                       "Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableIntArray &A,const DTRange &range) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(IntArray,Range)","No file");
        return false;
    }
    
    if (A.IsEmpty())
        return true;

    if (range.end()>A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(IntArray,Range)","Range out of bounds");
        return false;
    }

    if (ssize_t(fread(A.Pointer()+range.start,sizeof(int),range.length,FILEForReading()))!=range.length) {
        DTErrorMessage("DTFile::ReadBinary(IntArray,Range)","Could not read the required number of values from the file");
        return false;
    }
    
    return true;
}

bool DTFile::ReadBinary(DTMutableShortIntArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(ShortIntArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (fread(A.Pointer(),2,A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(ShortIntArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}


bool DTFile::ReadBinary(DTMutableUShortIntArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(UShortIntArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (fread(A.Pointer(),2,A.Length(),FILEForReading())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::ReadBinary(UShortIntArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableUCharArray &A) const
{
    return ReadBinary(A,0,A.Length());
}

bool DTFile::ReadBinary(DTMutableUCharArray &A,ssize_t len) const
{
    if (A.Length()<len) {
        DTErrorMessage("DTFile::ReadBinary(UCharArray,int)","Invalid length");
        return false;
    }

    return ReadBinary(A,0,len);
}

bool DTFile::ReadBinary(DTMutableUCharArray &A,ssize_t startAt,ssize_t howMuchToRead) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(UCharArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (A.Length()<startAt+howMuchToRead) {
        DTErrorMessage("DTFile::ReadBinary(UCharArray,int start,int length)","Invalid range");
        return false;
    }

    if (howMuchToRead==0)
        return true;

    if (fread(A.Pointer()+startAt,1,howMuchToRead,FILEForReading())!=(unsigned int)howMuchToRead) {
        DTErrorMessage("DTFile::ReadBinary(UCharArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadBinary(DTMutableCharArray &A) const
{
    return ReadBinary(A,A.Length());
}

bool DTFile::ReadBinary(DTMutableCharArray &A,ssize_t len) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadBinary(CharArray)","No file");
        return false;
    }

    if (A.IsEmpty())
        return true;

    if (A.Length()<len) {
        DTErrorMessage("DTFile::ReadBinary(CharArray,int)","Invalid length");
        return false;
    }

    if (fread(A.Pointer(),1,len,FILEForReading())!=(unsigned int)len) {
        DTErrorMessage("DTFile::ReadBinary(CharArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

bool DTFile::ReadAscii(DTMutableDoubleArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadAscii(DoubleArray)","No file");
        return false;
    }

    FILE *theFile = FILEForReading();

    if (A.IsEmpty())
        return true;

    size_t pos = 0;
    size_t len = A.Length();
    int howMany = -1;
    char singleChar;

    while (pos<len) {
        howMany = fscanf(theFile,"%lf",&A(pos));
        if (howMany==0) {
            // Try to skip over one character and try again.
            if (fread(&singleChar,1,1,theFile)!=1)
                break; // end of file
            howMany = fscanf(theFile,"%lf",&A(pos));
            if (howMany<=0)
                break;
        }
        else if (howMany==-1) {
            break;
        }
        pos++;
    }

    if (pos<len) {
        if (howMany==-1) {
            DTErrorMessage("DTFile::ReadAscii(DoubleArray)","Could not read the required number of values from the file");
        }
        return false;
    }

    return true;
}

DTMutableDoubleArray DTFile::ReadFortranBinary()
{
    return DTMutableDoubleArray();
}

bool DTFile::ReadAscii(DTMutableFloatArray &A) const
{
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadAscii(FloatArray)","No file");
        return false;
    }

    FILE *theFile = FILEForReading();

    if (A.IsEmpty())
        return true;

    size_t pos = 0;
    size_t len = A.Length();
    size_t howMany;
    char singleChar;

    while (pos<len) {
        howMany = fscanf(theFile,"%f",&A(pos));
        if (howMany==0) {
            // Try to skip over one character and try again.
            if (fread(&singleChar,1,1,theFile)!=1)
                break; // end of file
            continue;
        }
        pos++;
    }

    if (pos<len) {
        DTErrorMessage("DTFile::ReadAscii(FloatArray)","Could not read the required number of values from the file");
        return false;
    }

    return true;
}

double DTFile::ReadAsciiNumber(void) const
{
#if defined(WIN32) && !defined(NAN)
#define NAN std::numeric_limits<float>::quiet_NaN();
#endif
    if (!IsOpen()) {
        DTErrorMessage("DTFile::ReadAsciiNumber(DTFile)","No file");
        return NAN;
    }

    double toReturn = NAN;
    double temp;
    char singleChar;
    FILE *theFile = FILEForReading();

    while (1) {
        if (fscanf(theFile,"%lf",&temp)==0) {
            // Try to skip over one character and try again.
            if (fread(&singleChar,1,1,theFile)!=1)
                break; // end of file
            continue;
        }
        toReturn = temp;
        break;
    }

    return toReturn;
}

char DTFile::CharacterAtCurrentPosition(void) const
{
    char toReturn = char(getc(FILEForReading()));
    MovePosition(-1);
    return toReturn;
}

bool DTFile::CheckWriteErrorState(const char *errStr) const
{
    if (!IsOpen()) {
        DTErrorMessage(errStr,"No file");
        return true;
    }
    if (storage->readOnly) {
        DTErrorMessage(errStr,"Read only");
        return true;
    }

    return false;
}

bool DTFile::WriteString(string theStr)
{
    if (CheckWriteErrorState("DTFile::WriteString(string)"))
        return false;
    
    const char *cStr = theStr.c_str();
    size_t len = theStr.length();
    
    if (fwrite(cStr,1,len,FILEForWriting())!=(unsigned int)len) {
        DTErrorMessage("DTFile::WriteString(string)","Could not write the string to the file.");
        return false;
    }
    
    
    return true;
}

bool DTFile::WriteStringWithZero(string theStr)
{
    if (CheckWriteErrorState("DTFile::WriteStringWithZero(string)"))
        return false;
    
    const char *cStr = theStr.c_str();
    size_t len = theStr.length()+1;
    
    if (fwrite(cStr,1,len,FILEForWriting())!=(unsigned int)len) {
        DTErrorMessage("DTFile::WriteStringWithZero(string)","Could not write the string to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteUnsignedShort(unsigned short int v)
{
    if (CheckWriteErrorState("DTFile::WriteUnsignedShort(value)"))
        return false;
    
    if (fwrite(&v,sizeof(unsigned short int),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::WriteUnsignedShort(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}

#if defined(WIN32)
bool DTFile::Write8ByteInt(__int64 v)
{
    if (CheckWriteErrorState("DTFile::Write8ByteInt(value)"))
        return false;
    
    if (fwrite(&v,sizeof(__int64),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::Write8ByteInt(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}
#else
bool DTFile::Write8ByteInt(int64_t v)
{
    if (CheckWriteErrorState("DTFile::Write8ByteInt(value)"))
        return false;
    
    if (fwrite(&v,sizeof(int64_t),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::Write8ByteInt(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}
#endif

bool DTFile::WriteFloat(float v)
{
    if (CheckWriteErrorState("DTFile::WriteFloat(value)"))
        return false;

    if (fwrite(&v,sizeof(float),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::WriteFloat(value)","Could not write the number to the file.");
        return false;
    }

    return true;
}

bool DTFile::WriteDouble(double v)
{
    if (CheckWriteErrorState("DTFile::WriteDouble(value)"))
        return false;

    if (fwrite(&v,sizeof(double),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::WriteDouble(value)","Could not write the number to the file.");
        return false;
    }

    return true;
}

bool DTFile::Write4ByteInt(int v)
{
    if (CheckWriteErrorState("DTFile::Write8ByteInt(value)"))
        return false;
    
    if (fwrite(&v,sizeof(int),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::Write8ByteInt(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::Write2ByteInt(short int v)
{
    if (CheckWriteErrorState("DTFile::Write8ByteInt(value)"))
        return false;
    
    if (fwrite(&v,sizeof(short int),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::Write8ByteInt(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::Write1ByteInt(char v)
{
    if (CheckWriteErrorState("DTFile::Write8ByteInt(value)"))
        return false;
    
    if (fwrite(&v,sizeof(char),1,FILEForWriting())!=1) {
        DTErrorMessage("DTFile::Write8ByteInt(value)","Could not write the number to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteRaw(const char *ptr,ssize_t howMany)
{
    if (CheckWriteErrorState("DTFile::WriteRaw(value)"))
        return false;
    
    if (int(fwrite(ptr,1,howMany,FILEForWriting()))!=howMany) {
        DTErrorMessage("DTFile::WriteRaw(ptr,length)","Could not write the data to the file.");
        return false;
    }
    
    return true;
}

#pragma mark ---------- Write an array

bool DTFile::WriteBinary(const DTDoubleArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTDoubleArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),sizeof(double),A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTDoubleArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTFloatArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTFloatArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),sizeof(float),A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTFloatArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTIntArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTIntArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),sizeof(int),A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTIntArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTShortIntArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTShortIntArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),sizeof(short int),A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTShortIntArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTUShortIntArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTUShortIntArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),sizeof(unsigned short int),A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTUShortIntArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTCharArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTCharArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;
    
    if (fwrite(A.Pointer(),1,A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTCharArray)",
                       "Could not write the array to the file.");
        return false;
    }
    
    return true;
}

bool DTFile::WriteBinary(const DTUCharArray &A)
{
    if (CheckWriteErrorState("DTFile::WriteBinary(DTUCharArray)"))
        return false;
    
    if (A.IsEmpty())
        return true;

    if (fwrite(A.Pointer(),1,A.Length(),FILEForWriting())!=(unsigned int)A.Length()) {
        DTErrorMessage("DTFile::WriteBinary(DTUCharArray)",
                       "Could not write the array to the file.");
        return false;
    }

    return true;
}

void DTFile::pinfo(void) const
{
#ifndef DG_NOSTDErrOut
    std::cerr << "File = " << storage->name;
    if (IsReadOnly()) std::cerr << " [read only]";
    std::cerr << std::flush;
#endif
}

DTFolder DTFolder::AppendFolderName(const std::string &n)
{
    return DTFolder(name+"/"+n);
}

DTFile DTFolder::AppendFileName(const std::string &n)
{
    return DTFile(name+"/"+n);
}
