// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see http://www.visualdatatools.com/DTSource/license.html for more information.

#include "DTDataFile.h"

#include "DTDoubleArray.h"
#include "DTIntArray.h"
#include "DTCharArray.h"
#include "DTUCharArray.h"
#include "DTFloatArray.h"
#include "DTShortIntArray.h"
#include "DTUShortIntArray.h"
#include "DTArrayConversion.h"
#include "DTCharArray.h"
#include "DTUtilities.h"
#include "DTEndianSwap.h"
#include "DTLock.h"

#include "DTError.h"

#include <string>
#include <map>
#include <algorithm>

#if defined(WIN32)
typedef __int64 int64_t;
#endif
 
struct DTDataEntry {
    DTDataEntry() : m(0), n(0), o(0), type(0), location(-1) {}
    int m,n,o,type;
    off_t location; // -1 if not found.
};

struct DTDataFileStructure {
    DTDataFileStructure() : blockLength(0), type(0), m(0), n(0), o(0), nameLength(0) {}
    DTDataFileStructure(size_t mv,size_t nv,size_t ov,size_t bl,int tp,
                        size_t nl) : blockLength(int(bl)), type(tp), m(int(mv)), n(int(nv)), o(int(ov)), nameLength(int(nl)) {}
    
    int64_t blockLength;  // Total length of this block (including this entry).
    int type;
    int m;
    int n;
    int o;
    int nameLength;             // Includes the ending \0
};

class DTDataFileContent
{
public:
    DTDataFileContent(DTFile file);
    ~DTDataFileContent();
    
    void ReadInContent();
    
    void Lock(void) const {accessLock.Lock();}
    void Unlock(void) const {accessLock.Unlock();}
    
    DTLock accessLock;
    int referenceCount;
    map<string,DTDataEntry> content;
    
    DTFile file;
    
    bool isAtEnd; // Where the read-head is.
    bool swapBytes;
    bool haveSavedIndexBlock;
    
    // Index file is not saved after each addition, typically you save this just before you close the file.
    bool saveIndex;
    bool saveIndexWhenClosing;
    
    void SaveIndexBlock(void);
    
private:
    // Don't allow assignment
    DTDataFileContent(const DTDataFileContent &);
    DTDataFileContent &operator=(const DTDataFileContent &);
};

struct DTDataFileIndexEntry {
    int64_t location;
    char type;
    int m;
    int n;
    int o;
    int nameLength;
};

DTDataFileContent::~DTDataFileContent()
{
    if (saveIndexWhenClosing) {
        SaveIndexBlock();
    }
}

DTDataFileContent::DTDataFileContent(DTFile f)
{
    referenceCount = 1;
    file = f;
    saveIndex = false;
    isAtEnd = false;
    
    if (f.EndianType()==DTFile::Native) {
        swapBytes = false;
    }
    else {
        swapBytes = (DTFile::EndianForMachine()!=f.EndianType());
    }
    
    // Read in the content of the file, so subsequent access is much faster.
    ReadInContent();
}

void DTDataFileContent::ReadInContent(void)
{
    DTDataFileStructure TheHeader;
    DTDataEntry singleEntry;

    if (!file.IsOpen())
        return;

    off_t EndsAt = file.Length();
    size_t howMuchRead;
    
    // Check if there is an index saved at the end of this file.
    if (EndsAt>24) {
        file.SetPosition(EndsAt-23);
        DTMutableCharArray readInto(23);
        file.ReadBinary(readInto);
        
        if (memcmp(readInto.Pointer()," DT Index end ",14)==0) {
            // Read the last two bytes as an offset
            off_t offset;
            memcpy(&offset,readInto.Pointer()+15,sizeof(off_t));
            if (swapBytes) {
                DTSwap8Bytes(((unsigned char *)&offset),sizeof(off_t));
            }
            
            DTMutableCharArray headerBlock;
            
            // Check to make sure that this points to the start of the variable.
            if (offset>0 && offset<EndsAt-25) {
                file.SetPosition(offset);
                DTDataFileStructure header;
                howMuchRead = fread(&header,1,28,file.FILEForReading());
                if (howMuchRead!=28) {
                    return;
                }
                if (swapBytes) {
                    DTSwap8Bytes(((unsigned char *)&header),8);
                    DTSwap4Bytes(((unsigned char *)&header)+8,20);
                }
                if (header.n==1 && header.o==1 && header.nameLength==17) {
                    // Make sure that the name is correct
                    file.ReadBinary(readInto,17);
                    if (memcmp(readInto.Pointer()," DT Index block ",16)==0) {
                        headerBlock = DTMutableCharArray(header.m-23);
                        file.ReadBinary(headerBlock);
                    }
                }
            }
            
            if (headerBlock.Length() && headerBlock(0)==1) {
                ssize_t indexLength = headerBlock.Length();
                const char *bufferP = headerBlock.Pointer();
                
                ssize_t m,n,o,nameLength,entrySize;
                ssize_t position = 1;
                bool foundAProblem = false;
                const char *buffer;
                std::string name;
                
                while (position<indexLength) {
                    buffer = bufferP + position;
                    singleEntry.type = buffer[0];
                    m = singleEntry.m = *((int *)(buffer+1));
                    n = singleEntry.n = *((int *)(buffer+5));
                    o = singleEntry.o = *((int *)(buffer+9));
                    singleEntry.location = *((off_t *)(buffer+13));
                    position += 21;
                    buffer = bufferP + position;
                    nameLength = 0;
                    while (buffer[nameLength] && position+nameLength<indexLength) nameLength++;
                    
                    if (nameLength==0 || m<0 || n<0 || o<0 || nameLength>10000 || singleEntry.location<0 || singleEntry.location>=EndsAt) {
                        foundAProblem = true;
                        break;
                    }
                    
                    name = std::string(buffer);
                    
                    switch (singleEntry.type) {
                        case 1:
                            entrySize = 8;
                            break;
                        case 2:
                        case 8:
                            entrySize = 4;
                            break;
                        case 9:
                        case 10:
                            entrySize = 2;
                            break;
                        case 11:
                        case 12:
                        case 20:
                            entrySize = 1;
                            break;
                        default:
                            foundAProblem = true;
                            entrySize = 0;
                    }
                    if (foundAProblem) {
                        break;
                    }
                    
                    position += nameLength+1;
                    
                    if (singleEntry.location+entrySize*m*n*o>EndsAt) {
                        // This entry isn't completely stored in the file, ignore it.
                        DTErrorMessage("DataFile::ReadContent","Index has a problem");
                        foundAProblem = true;
                        break;
                    }
                    
                    content[name] = singleEntry;
                }
                
                if (foundAProblem==false) {
                    return;
                }
                else {
                    // Fallback to reading the index from the file.  Clear what is here already
                    content.clear();
                }
            }
        }
        file.SetPosition(0);
    }
    
    // See if there is an index file saved.
    std::string indexName;
    std::string thisFile = file.Name();
    if (thisFile.length()>6 && thisFile.substr(thisFile.length()-6,6)==".dtbin") {
        indexName = thisFile.substr(0,thisFile.length()-6)+".index";
    }
    else {
        indexName = thisFile+".index";
    }
    if (DTFile::CanOpen(indexName,DTFile::ReadOnly)) {
        // Read in the entire index file
        
        // Read in the portion of the file that is valid.  The index file might
        // include entries that haven't been completely saved in the data file, and
        // those entries should be quietly skipped.
        DTFile tempIndex(indexName,DTFile::ReadOnly);
        DTFilePosition indexLength = tempIndex.Length();
        
        const char *shouldBe = "dtbin index file";
        if (indexLength<(ssize_t)strlen(shouldBe)+1) {
            // Too small
            indexLength = 0;
        }

        FILE *indexFilePointer = tempIndex.FILEForReading();
        if (indexFilePointer==NULL) indexLength = 0;
        
        DTMutableCharArray bufferArray((ssize_t)(indexLength ? indexLength+1 : 0));
        
        char *bufferP = bufferArray.Pointer();
        if (indexLength && bufferP) {
            bufferP[indexLength] = 0;
            if ((ssize_t)fread(bufferP,1,(ssize_t)indexLength,indexFilePointer)!=indexLength) {
                indexLength = 0;
                bufferArray = DTMutableCharArray();
                bufferP = 0;
            }
        }
        
        // Check if this is the right type
        if (bufferP!=NULL && memcmp(bufferP, shouldBe,strlen(shouldBe)+1)!=0) {
            bufferArray = DTMutableCharArray();
            bufferP = 0;
        }
        
        if (bufferP) {
            ssize_t m,n,o,nameLength,entrySize;
            ssize_t position = strlen(shouldBe)+1;
            bool foundAProblem = false;
            char *buffer;
            std::string name;
            
            while (position<indexLength) {
                buffer = bufferP + position;
                singleEntry.location = *((int64_t *)buffer);
                singleEntry.type = buffer[8]%100;
                m = singleEntry.m = ((int *)(buffer+9))[0];
                n = singleEntry.n = ((int *)(buffer+9))[1];
                o = singleEntry.o = ((int *)(buffer+9))[2];
                position += 8+1+12;
                buffer = bufferP + position;
                nameLength = 0;
                while (buffer[nameLength]) nameLength++;
                
                if (nameLength==0 || m<0 || n<0 || o<0) {
                    foundAProblem = true;
                    break;
                }
                
                name = std::string(buffer);
                
                switch (singleEntry.type) {
                    case 1:
                        entrySize = 8;
                        break;
                    case 2:
                    case 8:
                        entrySize = 4;
                        break;
                    case 9:
                    case 10:
                        entrySize = 2;
                        break;
                    case 11:
                    case 12:
                    case 20:
                        entrySize = 1;
                        break;
                    default:
                        foundAProblem = true;
                        entrySize = 0;
                }
                if (foundAProblem) {
                    break;
                }
                
                position += nameLength+1;

                if (singleEntry.location+entrySize*m*n*o>EndsAt) {
                    // This entry isn't completely stored in the file, ignore it.
                    DTErrorMessage("DataFile::ReadContent","Index has a problem");
                    foundAProblem = true;
                    break;
                }
                
                content[name] = singleEntry;
            }

            if (foundAProblem==false) {
                return;
            }
            else {
                // Fallback to reading the index from the file.  Clear what is here already
                content.clear();
            }
        }
    }
    
    // The index file wasn't found or not valid, so read the content from the file.
    FILE *theFile = file.FILEForReading();
    
    size_t howMuchWasRead;
    char tempString[255];
    off_t StartsAt;
    
    // The file always start with an identification string.
    // DataTank Binary File v1\0  - Always a big endian file (before Apple's Intel plans)
    // DataTank Binary File LE\0  - Little endian
    // DataTank Binary File BE\0  - Big endian
    
    // that is 24 bytes.
    const char *identifierOld = "DataTank Binary File v1\0";
    const char *identifierLE = "DataTank Binary File LE\0";
    const char *identifierBE = "DataTank Binary File BE\0";
    
    size_t howManyRead = fread(tempString,1,24,theFile);
    if (howManyRead==0) {
        return; // Empty is ok.
    }
    
    if (howManyRead!=24) {
        DTErrorMessage("DTDataFile::ReadInContent","Not a valid DataTank binary format.");
        return;
    }
    
    // Overwrite the swapBytes
    if (strncmp(identifierOld,tempString,24)==0) {
        swapBytes = (DTFile::EndianForMachine()==DTFile::LittleEndian);
    }
    else if (strncmp(identifierLE,tempString,24)==0) {
        swapBytes = (DTFile::EndianForMachine()==DTFile::BigEndian);
    }
    else if (strncmp(identifierBE,tempString,24)==0) {
        swapBytes = (DTFile::EndianForMachine()==DTFile::LittleEndian);
    }
    else {
        DTErrorMessage("DTDataFile::ReadInContent","Not a valid DataTank binary format.");
        return;
    }

    
    bool IsOK = true;
    while (IsOK) {
        StartsAt = file.Position();

        howMuchWasRead = fread(&TheHeader,28,1,theFile);
        
        if (swapBytes) {
            DTSwap8Bytes(((unsigned char *)&TheHeader),8);
            DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
        }
        if (howMuchWasRead==0) {
            break;
        }
        if (TheHeader.blockLength==0 ||
            howMuchWasRead<1 || TheHeader.type>2000 || TheHeader.nameLength>255) {
            DTErrorMessage("Reading In File Content","Invalid file format.");
            break;
        }
        if (int(TheHeader.blockLength)>EndsAt-StartsAt)
            break; // Incomplete data.  Most likely an incomplete write.

        howMuchWasRead = fread(tempString,1,TheHeader.nameLength,theFile);

        // Add this entry to the content list.
        singleEntry.location = StartsAt+28+TheHeader.nameLength;
        singleEntry.m = TheHeader.m;
        singleEntry.n = TheHeader.n;
        singleEntry.o = TheHeader.o;
        singleEntry.type = TheHeader.type;
        content[string(tempString)] = singleEntry;
        
        // Get ready for the next entry.
        file.SetPosition(StartsAt+TheHeader.blockLength);
    }
    isAtEnd = false;
}

void DTDataFileContent::SaveIndexBlock(void)
{
    // Makes it possible to read in the table of contents really quickly.
    // Saves a block as a variable, called " DT Index block ", which means that you should not use this name for a variable.
    // The structure of the content is that it contains a byte stream with one block for each entry, and ends with
    // the text " DT Index end \0" followed by the byte offset where the variable is stored.
    if (content.size()==0) return; // Nothing to save
        
    FILE *theFile = file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        return;
    }
    
    if (!isAtEnd) {
        file.MoveToEnd();
        isAtEnd = true;
    }
    long long int startOfVariable = file.Position();
    
    ssize_t posInBuffer = 0, lengthOfBuffer = 10000;
    DTMutableCharArray buffer(lengthOfBuffer);
    char *bufferD = buffer.Pointer();
    
    map<string,DTDataEntry>::const_iterator mapIterator;
    DTDataEntry fileEntry;
    string name;
    ssize_t strLen;
    
    buffer(0) = 1; // Version number
    posInBuffer = 1;
    
    for (mapIterator=content.begin();mapIterator!=content.end();++mapIterator) {
        name = mapIterator->first;
        fileEntry = mapIterator->second;
        
        strLen = strlen(name.c_str());
        
        // Logarithmic growth
        if (posInBuffer + strLen + 100 > lengthOfBuffer) {
            buffer = IncreaseSize(buffer,buffer.Length()+strLen);
            bufferD = buffer.Pointer();
            lengthOfBuffer = buffer.Length();
        }
        // Order is type - 1 byte
        // m,n,o - 4 bytes
        // location - 8 bytes
        // name - 0 terminated
        
        bufferD[posInBuffer] = (char)fileEntry.type;
        *((int *)(bufferD+posInBuffer+1)) = fileEntry.m;
        *((int *)(bufferD+posInBuffer+5)) = fileEntry.n;
        *((int *)(bufferD+posInBuffer+9)) = fileEntry.o;
        *((off_t *)(bufferD+posInBuffer+13)) = fileEntry.location;
        posInBuffer+=21;
        memcpy(bufferD+posInBuffer,name.c_str(),strLen+1);
        posInBuffer += strLen+1;
    }
    
    std::string endString = " DT Index end ";
    strLen = strlen(endString.c_str());
    memcpy(bufferD+posInBuffer,endString.c_str(),strLen+1);
    posInBuffer += strLen+1;
    
    memcpy(bufferD+posInBuffer,(char *)&startOfVariable,sizeof(long long int));
    posInBuffer += 8;
    
    buffer = TruncateSize(buffer,posInBuffer);
    
    std::string VarName = " DT Index block ";
    DTDataFileStructure TheHeader(buffer.Length(),1,1,
                                  29+VarName.length()+buffer.Length()*sizeof(char),
                                  DTDataFile_Signed8Char,
                                  1+VarName.length());
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = file.Position()+28+TheHeader.nameLength;
    content[VarName] = entry;
    
    if (swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);
    fwrite(buffer.Pointer(),sizeof(char),buffer.Length(),theFile);
}

#pragma mark DTDataFile

DTDataFile::DTDataFile(DTFile file)
: DTDataStorage(), content(NULL)
{
    if (file.Length()==0) {
        // Need to remove the index file if it exists
        std::string indexName;
        std::string name = file.Name();
        if (name.length()>4 && name.substr(name.length()-4,4)==".mat") {
            indexName = name.substr(0,name.length()-4)+".index";
        }
        else if (name.length()>6 && name.substr(name.length()-6,6)==".dtbin") {
            indexName = name.substr(0,name.length()-6)+".index";
        }
        else {
            indexName = name+".index";
        }
        unlink(indexName.c_str());
    }
    content = new DTDataFileContent(file);
}

DTDataFile::DTDataFile()
: DTDataStorage(), content(NULL)
{
    DTFile emptyFile;
    content = new DTDataFileContent(emptyFile);
}

DTDataFile::DTDataFile(const std::string &name,DTFile::OpenType oType)
: DTDataStorage(), content(NULL)
{
    if (oType==DTFile::NewReadWrite) {
        // Need to remove the index file if it exists
        std::string indexName;
        if (name.length()>4 && name.substr(name.length()-4,4)==".mat") {
            indexName = name.substr(0,name.length()-4)+".index";
        }
        else if (name.length()>6 && name.substr(name.length()-6,6)==".dtbin") {
            indexName = name.substr(0,name.length()-6)+".index";
        }
        else {
            indexName = name+".index";
        }
        unlink(indexName.c_str());
    }
    content = new DTDataFileContent(DTFile(name,oType));
}

DTDataFile::DTDataFile(const DTDataFile &C)
: DTDataStorage(C), content(C.content)
{
    C.content->Lock();
    content = C.content;
    content->referenceCount++;
    C.content->Unlock();
}

DTDataFile &DTDataFile::operator=(const DTDataFile &C)
{
    if (content==C.content) return *this; // Slight thread safety issue might come up, if someone is reassigning the content.
    content->Lock();
    C.content->Lock();
    content->referenceCount--;
    if (content->referenceCount==0) {
        content->Unlock();
        delete content;
    }
    else {
        content->Unlock();
    }
    content = C.content;
    content->referenceCount++;
    content->Unlock();
    
    return *this;
}

DTDataFile::~DTDataFile()
{
    content->Lock();
    content->referenceCount--;
    if (content->referenceCount==0) {
        content->Unlock();
        delete content;
    }
    else {
        content->Unlock();
    }
}

DTMutablePointer<DTDataStorage> DTDataFile::AsPointer() const
{
    return DTMutablePointer<DTDataStorage>(new DTDataFile(*this));
}

void DTDataFile::Save(int v,const std::string &name)
{
    DTMutableIntArray temp(1);
    temp(0) = v;
    Save(temp,name);
}

void DTDataFile::WriteHeaderIfNecessary(void)
{
    // Called inside a lock
    if (content->content.size()>0)
        return;
    
    const char *identifierLE = "DataTank Binary File LE\0";
    const char *identifierBE = "DataTank Binary File BE\0";
    const char *identifier;
    
    if (content->swapBytes) {
        if (DTFile::RunningOnBigEndianMachine())
            identifier = identifierLE;
        else
            identifier = identifierBE;
    }
    else {
        if (DTFile::RunningOnBigEndianMachine())
            identifier = identifierBE;
        else
            identifier = identifierLE;
    }
    
    fwrite(identifier,1,24, content->file.FILEForWriting());
    
}

void DTDataFile::Save(double v,const std::string &name)
{
    DTMutableDoubleArray temp(1);
    temp(0) = v;
    Save(temp,name);
}

void DTDataFile::Save(const DTDoubleArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(double),
                                  DTDataFile_Double,
                                  1+VarName.length());
    
    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;

    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);

    if (A.Length()) {
        if (content->swapBytes) {
            DTMutableDoubleArray Temp = A.Copy();
            SwapEndian(Temp);
            fwrite(Temp.Pointer(),sizeof(double),Temp.Length(),theFile);
        }
        else {
            fwrite(A.Pointer(),sizeof(double),A.Length(),theFile);
        }
    }
    
    content->Unlock();
}

void DTDataFile::Save(const DTFloatArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(float),
                                  DTDataFile_Single,
                                  1+VarName.length());

    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();

    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;
    
    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);

    if (A.Length()) {
        if (content->swapBytes) {
            DTMutableFloatArray Temp = A.Copy();
            SwapEndian(Temp);
            fwrite(Temp.Pointer(),sizeof(float),Temp.Length(),theFile);
        }
        else {
            fwrite(A.Pointer(),sizeof(float),A.Length(),theFile);
        }
    }

    content->Unlock();
}

void DTDataFile::Save(const DTIntArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(int),
                                  DTDataFile_Signed32Int,
                                  1+VarName.length());
    
    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;
    
    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);

    if (A.Length()) {
        if (content->swapBytes) {
            DTMutableIntArray Temp = A.Copy();
            SwapEndian(Temp);
            fwrite(Temp.Pointer(),sizeof(int),Temp.Length(),theFile);
        }
        else {
            fwrite(A.Pointer(),sizeof(int),A.Length(),theFile);
        }
    }

    content->Unlock();
}

void DTDataFile::Save(const DTUCharArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(char),
                                  DTDataFile_Unsigned8Char,
                                  1+VarName.length());

    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;
    
    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);
    if (A.Length()) fwrite(A.Pointer(),sizeof(char),A.Length(),theFile);

    content->Unlock();
}

void DTDataFile::Save(const DTCharArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(char),
                                  DTDataFile_Signed8Char,
                                  1+VarName.length());

    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;

    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);
    if (A.Length()) fwrite(A.Pointer(),sizeof(char),A.Length(),theFile);

    content->Unlock();
}

void DTDataFile::Save(const DTShortIntArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(short),
                                  DTDataFile_Short,
                                  1+VarName.length());
    
    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;
    
    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);
    
    if (A.Length()) {
        if (content->swapBytes) {
            DTMutableShortIntArray Temp = A.Copy();
            SwapEndian(Temp);
            fwrite(Temp.Pointer(),sizeof(short int),Temp.Length(),theFile);
        }
        else {
            fwrite(A.Pointer(),sizeof(short int),A.Length(),theFile);
        }
    }

    content->Unlock();
}

void DTDataFile::Save(const DTUShortIntArray &A,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(A.m(),A.n(),A.o(),
                                  29+VarName.length()+A.Length()*sizeof(short),
                                  DTDataFile_UnsignedShort,
                                  1+VarName.length());

    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;

    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);

    if (A.Length()) {
        if (content->swapBytes) {
            DTMutableUShortIntArray Temp = A.Copy();
            SwapEndian(Temp);
            fwrite(Temp.Pointer(),sizeof(unsigned short int),Temp.Length(),theFile);
        }
        else {
            fwrite(A.Pointer(),sizeof(unsigned short int),A.Length(),theFile);
        }
    }

    content->Unlock();
}

void DTDataFile::Save(const std::string &theString,const std::string &VarName)
{
    content->Lock();
    if (IsReadOnly()) {
        DTErrorMessage("DTDataFile::Save","File is read only.");
        content->Unlock();
        return;
    }
    FILE *theFile = content->file.FILEForWriting();
    if (theFile==NULL) {
        DTErrorMessage("DTDataFile::Save","Empty File.");
        content->Unlock();
        return;
    }
    
    DTDataFileStructure TheHeader(theString.length()+1,1,1,
                                  29+VarName.length()+theString.length()+1,
                                  DTDataFile_String,
                                  1+VarName.length());

    if (!content->isAtEnd) {
        content->file.MoveToEnd();
        content->isAtEnd = true;
    }
    
    WriteHeaderIfNecessary();
    
    DTDataEntry entry;
    entry.m = TheHeader.m;
    entry.n = TheHeader.n;
    entry.o = TheHeader.o;
    entry.type = TheHeader.type;
    entry.location = content->file.Position()+28+TheHeader.nameLength;
    content->content[VarName] = entry;
    
    if (content->swapBytes) {
        DTSwap8Bytes(((unsigned char *)&TheHeader),8);
        DTSwap4Bytes(((unsigned char *)&TheHeader)+8,20);
    }
    fwrite(&TheHeader,28, 1, theFile);
    fwrite(VarName.c_str(),sizeof(char),1+VarName.length(),theFile);
    fwrite(theString.c_str(),sizeof(char),theString.length()+1,theFile);

    content->Unlock();
}

void DTDataFile::Sync(void) const
{
    if (IsReadOnly()) return;
    
    content->Lock();
    content->file.Flush();
    
    std::string fileN = content->file.Name();
    std::string::size_type location = fileN.find_last_of(".");
    if (location==std::string::npos) {
        // Can't save a sync if you don't specify the name
        content->Unlock();
        return;
    }
    
    std::string ending = fileN.substr(location+1);
    if (ending.find_last_of("/")!=std::string::npos) {
        // The period is not in the last path component
        content->Unlock();
        return;
    }
    if (location==0 || fileN[location-1]=='/') {
        // The . is at the beginning of the file name, i.e. this is a hidden file
        content->Unlock();
        return;
    }

    std::string syncName = fileN.substr(0,location)+".sync";
    DTFilePosition currentLength = content->file.Length();
    
    remove(syncName.c_str());
    FILE *syncFile = fopen(syncName.c_str(),"w+b");
    fwrite((void *)&currentLength, sizeof(DTFilePosition), 1, syncFile);
    fclose(syncFile);

    content->Unlock();
}

void DTDataFile::Flush(void) const
{
    if (IsReadOnly()) return;
    
    content->Lock();
    content->file.Flush();
    content->Unlock();
}

DTDataEntry DTDataFile::FindVariable(const std::string &name) const
{
    // Inside a lock
    // should be in the content->content list
    map<string,DTDataEntry>::const_iterator searchResult = content->content.find(name);
    
    if (searchResult==content->content.end()) {
        return DTDataEntry();
    }
    else {
        return searchResult->second;
    }
}

DTList<std::string> DTDataFile::AllVariableNames(void) const
{
    content->Lock();
    DTMutableList<std::string> toReturn(content->content.size());
    
    map<string,DTDataEntry>::const_iterator mapIterator;
    int pos = 0;
    DTDataEntry fileEntry;
    
    for (mapIterator=content->content.begin();mapIterator!=content->content.end();++mapIterator) {
        toReturn(pos++) = mapIterator->first;
    }
    
    sort(toReturn.Pointer(),toReturn.Pointer()+toReturn.Length());
    
    content->Unlock();
    return toReturn;
}

struct DTDataFilePosString {
    DTFilePosition pos;
    std::string description;

    bool operator<(const DTDataFilePosString &A) const {return (pos<A.pos);}
};

void DTDataFile::printInfo(void) const
{
#ifndef DG_NOSTDErrOut
    vector<DTDataFilePosString> list;
    DTDataFilePosString entry;
    std::string desc,stringValue;
    DTDataEntry fileEntry;

    std::cerr << "------------------------------------------------------------------------" << std::endl;
    std::cerr << "Content of \"" << content->file.Name() << "\" - ";
    size_t howMany = content->content.size();
    if (howMany==0)
        std::cerr << "empty" << std::endl;
    else if (howMany==1)
        std::cerr << "1 entry" << std::endl;
    else
        std::cerr << howMany << " entries" << std::endl;
    std::cerr << "------------------------------------------------------------------------" << std::endl;
    
    std::string padding = ".................................";
    
    map<string,DTDataEntry>::const_iterator mapIterator;
    for (mapIterator=content->content.begin();mapIterator!=content->content.end();++mapIterator) {
        fileEntry = mapIterator->second;
        entry.pos = fileEntry.location;
        desc = mapIterator->first + " ";
        // Pad to make it 30 characters
        if (desc.length()<30)
            desc = desc + string(padding,0,30-desc.length());
        switch (fileEntry.type) {
            case DTDataFile_Double:
                desc += " - double - ";
                break;
            case DTDataFile_Single:
                desc += " -  float - ";
                break;
            case DTDataFile_Signed32Int:
                desc += " -    int - ";
                break;
            case DTDataFile_UnsignedShort:
                desc += " - UShort - ";
                break;
            case DTDataFile_Short:
                desc += " -  short - ";
                break;
            case DTDataFile_Unsigned8Char:
                desc += " -  UChar - ";
                break;
            case DTDataFile_Signed8Char:
                desc += " -   char - ";
                break;
            case DTDataFile_String:
                desc += " - string - ";
                break;
            default:
                desc += " - ?????? - ";
                break;
        }
        // Dimension.
        if (fileEntry.type==DTDataFile_String) {
            stringValue = ReadString(mapIterator->first);
            if (stringValue.length()>25) {
                desc += "\""+string(stringValue,0,15) + "...\" - " + DTInt2String(fileEntry.m*fileEntry.n*fileEntry.o) + " characters";
            }
            else {
                desc += "\""+stringValue+"\"";
            }
        }
        else {
            if (fileEntry.m==0)
                desc += "Empty";
            else if (fileEntry.m==1 && fileEntry.n==1 && fileEntry.o==1)
                desc += DTFloat2StringShort(ReadNumber(mapIterator->first));
            else if (fileEntry.n==1 && fileEntry.o==1)
                desc += DTInt2String(fileEntry.m) + " numbers";
            else if (fileEntry.o==1)
                desc += DTInt2String(fileEntry.m) + " x " + DTInt2String(fileEntry.n) + " array";
            else
                desc += DTInt2String(fileEntry.m) + " x " + DTInt2String(fileEntry.n) + " x " + DTInt2String(fileEntry.o) + " array";
        }
        entry.description = desc;
        list.push_back(entry);
    }
    
    sort(list.begin(),list.end());

    // Print the content
    size_t howLong = list.size();
    size_t pos = 0;
    vector<DTDataFilePosString>::iterator iter;
	if (howLong<400) {
		for (iter=list.begin();iter!=list.end();++iter) {
			cerr << iter->description << std::endl;
		}
	}
	else {
		// Skip from the middle
		for (iter=list.begin();pos<350;++iter) {
            std::cerr << iter->description << std::endl;
			pos++;
		}
		cerr << "Skipping " << howLong-30-pos <<  " entries" << std::endl;
		while (pos<howLong-30) {
			++iter;
			pos++;
		}
		while (iter!=list.end()) {
            std::cerr << iter->description << std::endl;
			++iter;
		}
	}
    std::cerr << flush;
#endif
}

bool DTDataFile::IsOpen(void) const
{
    return content->file.IsOpen();
}

bool DTDataFile::Contains(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    return (entry.location>=0);
}

bool DTDataFile::SizeOf(const std::string &name,int &m,int &n,int &o) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) {
        m = n = o = 0;
        return false;
    }
    m = entry.m;
    n = entry.n;
    o = entry.o;
    return true;
}

bool DTDataFile::IsReadOnly() const
{
    return content->file.IsReadOnly();
}

bool DTDataFile::SavedAsCharacter(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_Signed8Char || entry.type==DTDataFile_Unsigned8Char);
}

bool DTDataFile::SavedAsShort(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_Short || entry.type==DTDataFile_UnsignedShort);
}

bool DTDataFile::SavedAsInt(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_Signed32Int);
}

bool DTDataFile::SavedAsFloat(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_Single);
}

bool DTDataFile::SavedAsDouble(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_Double);
}

bool DTDataFile::SavedAsString(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    content->Unlock();
    if (entry.location<0) return false;
    return (entry.type==DTDataFile_String);
}

void DTDataFile::SaveIndex(void)
{
    content->saveIndexWhenClosing = true; // Not saved until the file is closed.
}

DTDoubleArray DTDataFile::ReadDoubleArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadDoubleArray(name)",msg);
        content->Unlock();
        return DTDoubleArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;

    // Now read the array.
    DTMutableDoubleArray toReturn(m,n,o);

    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        content->file.ReadBinary(toReturn);
        if (content->swapBytes) SwapEndian(toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is a float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        // This is an short arrray.
        DTMutableCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadDoubleArray(name)","Trying to read in a string.");
        toReturn = DTMutableDoubleArray();
    }

    content->Unlock();
    return toReturn;
}

DTFloatArray DTDataFile::ReadFloatArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadFloatArray(name)",msg);
        content->Unlock();
        return DTFloatArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;
    
    // Now read the array.
    DTMutableFloatArray toReturn(m,n,o);

    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        content->file.ReadBinary(toReturn);
        if (content->swapBytes) SwapEndian(toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        // This is an short arrray.
        DTMutableCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadFloatArray(name)","Trying to read in a string.");
        toReturn = DTMutableFloatArray();
    }
    content->Unlock();
    
    return toReturn;
}

DTIntArray DTDataFile::ReadIntArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadIntArray(name)",msg);
        content->Unlock();
        return DTIntArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;
    
    // Now read the array.
    DTMutableIntArray toReturn(m,n,o);
    
    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is an float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        content->file.ReadBinary(toReturn);
        if (content->swapBytes) SwapEndian(toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        // This is an short arrray.
        DTMutableCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadIntArray(name)","Trying to read in a string.");
        toReturn = DTMutableIntArray();
    }
    content->Unlock();

    return toReturn;
}

DTCharArray DTDataFile::ReadCharArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadCharArray(name)",msg);
        content->Unlock();
        return DTCharArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;

    // Now read the array.
    DTMutableCharArray toReturn(m,n,o);

    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is an float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        content->file.ReadBinary(toReturn);
    }
    else if (entry.type==DTDataFile_String || entry.type==DTDataFile_Unsigned8Char) {
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadCharArray(name)",
                       "Haven't taken care of this case.");
    }
    content->Unlock();
    
    return toReturn;
}

DTUCharArray DTDataFile::ReadUCharArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadUCharArray(name)",msg);
        content->Unlock();
        return DTUCharArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;
    
    // Now read the array.
    DTMutableUCharArray toReturn(m,n,o);
    
    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is an float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_String || entry.type==DTDataFile_Unsigned8Char || entry.type==DTDataFile_Signed8Char) {
        content->file.ReadBinary(toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadCharArray(name)",
                       "Haven't taken care of this case.");
        toReturn = DTMutableUCharArray();
    }
    content->Unlock();
    
    return toReturn;
}

DTShortIntArray DTDataFile::ReadShortIntArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadShortIntArray(name)",msg);
        content->Unlock();
        return DTShortIntArray();
    }

    DTFilePosition StartsAt = entry.location;

    int m = entry.m;
    int n = entry.n;
    int o = entry.o;

    // Now read the array.
    DTMutableShortIntArray toReturn(m,n,o);

    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is an float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        // This is an unsigned short arrray.
        DTMutableUShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        content->file.ReadBinary(toReturn);
        if (content->swapBytes) SwapEndian(toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        // This is an short arrray.
        DTMutableCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadUShortIntArray(name)","Trying to read in a string.");
        toReturn = DTMutableShortIntArray();
    }
    content->Unlock();
    
    return toReturn;
}

DTUShortIntArray DTDataFile::ReadUShortIntArray(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the variable \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadUShortIntArray(name)",msg);
        content->Unlock();
        return DTUShortIntArray();
    }
    
    DTFilePosition StartsAt = entry.location;
    
    int m = entry.m;
    int n = entry.n;
    int o = entry.o;
    
    // Now read the array.
    DTMutableUShortIntArray toReturn(m,n,o);
    
    content->file.SetPosition(StartsAt);
    content->isAtEnd = false;
    if (entry.type==DTDataFile_Double) {
        // This is a double arrray.
        DTMutableDoubleArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Single) {
        // This is an float arrray.
        DTMutableFloatArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed32Int) {
        // This is an int arrray.
        DTMutableIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_UnsignedShort) {
        content->file.ReadBinary(toReturn);
        if (content->swapBytes) SwapEndian(toReturn);
    }
    else if (entry.type==DTDataFile_Short) {
        // This is an short arrray.
        DTMutableShortIntArray temp(m,n,o);
        content->file.ReadBinary(temp);
        if (content->swapBytes) SwapEndian(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Unsigned8Char) {
        // This is an unsigned short arrray.
        DTMutableUCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else if (entry.type==DTDataFile_Signed8Char) {
        // This is an short arrray.
        DTMutableCharArray temp(m,n,o);
        content->file.ReadBinary(temp);
        ConvertArray(temp,toReturn);
    }
    else {
        DTErrorMessage("dataFile.ReadUShortIntArray(name)","Trying to read in a string.");
        toReturn = DTMutableUShortIntArray();
    }
    content->Unlock();
    
    return toReturn;
}

std::string DTDataFile::ReadString(const std::string &name) const
{
    content->Lock();
    DTDataEntry entry = FindVariable(name);
    if (entry.location<0) {
        std::string msg = string("Did not find the string \"") + name + "\" inside the datafile.";
        DTErrorMessage("dataFile.ReadString(name)",msg);
        content->Unlock();
        return std::string();
    }
    if (entry.type!=DTDataFile_String) {
        std::string msg = string("The variable \"") + name + "\" is not a string.";
        DTErrorMessage("dataFile.ReadString(name)",msg);
        content->Unlock();
        return std::string();
    }

    if (entry.m==0) {
        content->Unlock();
        return std::string();
    }
    
    int m = entry.m;
    int n = entry.n;
    int o = entry.o;

    // Now read the array.
    content->file.SetPosition(entry.location);
    content->isAtEnd = false;

    DTMutableCharArray temp(m,n,o);
    content->file.ReadBinary(temp);

    std::string toReturn;
    if (temp(temp.Length()-1)=='\0')
        toReturn = temp.Pointer();
    else
        toReturn = string(temp.Pointer(),temp.Length());
    content->Unlock();
        
    return toReturn;
}

