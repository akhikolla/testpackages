
// Part of DTSource. Copyright 2004-2017. David Adalsteinsson.
// see https://www.visualdatatools.com/DTSource/license.html for more information.

#ifndef DTFile_Header
#define DTFile_Header

// Wrapper around a FILE pointer.  You can extract the pointer, or use the member functions
// to read the data, set file locations etc.

// This class will "own" the pointer and will call fclose() when the last reference is deleted.

// This is used in the data files, but handles the lowest level access.
// Also contains a number of ascii input routines.

#include <stdio.h>
#include <string>
#include <sys/types.h>
#include <stdlib.h>
#include <stdint.h>

using namespace std;

class DTDoubleArray;
class DTFloatArray;
class DTIntArray;
class DTShortIntArray;
class DTUShortIntArray;
class DTCharArray;
class DTUCharArray;

class DTMutableDoubleArray;
class DTMutableFloatArray;
class DTMutableIntArray;
class DTMutableCharArray;
class DTMutableUCharArray;
class DTMutableShortIntArray;
class DTMutableUShortIntArray;

struct DTFileStorage;
struct DTRange;

#if defined(WIN32)
typedef __int64 DTFilePosition;
#else
typedef off_t DTFilePosition;
#endif

class DTFile {
public:
    enum OpenType {ReadOnly, ExistingReadWrite, NewReadWrite};
    enum Endian {Native,LittleEndian,BigEndian};

    DTFile();
    DTFile(const std::string &,OpenType=ExistingReadWrite);
    DTFile(const std::string &,Endian,OpenType=ExistingReadWrite);
    
    DTFile(const DTFile &);
    ~DTFile();
    DTFile &operator=(const DTFile &);
    
    static bool RunningOnBigEndianMachine(void);
    static Endian EndianForMachine(void);
    
    static bool CanOpen(const std::string &,OpenType=ExistingReadWrite);

    bool IsOpen(void) const;
    bool IsReadOnly(void) const;
    std::string Name(void) const;
    DTFilePosition Length(void) const;
    Endian EndianType(void) const;
    
    // Position in the file.
    DTFilePosition Position(void) const;
    void SetPosition(DTFilePosition) const;
    void MovePosition(DTFilePosition) const;
    void MoveToEnd(void) const;
    bool AtEndOfFile(void) const;

    void Flush(void) const;
        
    // In order to do anything, need to use the underlying FILE pointer.
    FILE *FILEForReading(void) const;
    FILE *FILEForWriting(void); // Clears the length cache
    
    // Reading a string
    std::string ReadLine(ssize_t maxLen=-1) const; // To the next newline or \0 character
    std::string ReadString(size_t length) const; // Exact length.
    std::string NextWord(void) const; // Read until new line, space or non-printable character.  Swallows the next character.
    
    double ReadNumberWithSeparator(char c,bool &empty,bool &endOfLine) const; // Skips over spaces before and after ends after the separator
    std::string ReadStringWithSeparator(char c,bool &empty,bool &endOfLine) const; // Skips over spaces before and after ends after the separator.  Checks for quotes.
    
    // Searching
    bool Find(char) const; // If found, position=location of character.
    
    // Binary input.  The size is determined by the array size.
    // If you need to change between big and little endian, use SwapEndian()
    // on the array.
    
    bool ReadCharacters(char *chars,size_t howMuchToRead) const;
    
    unsigned short int ReadUnsignedShort() const;

    // Read in a single number (binary)
    int Read_int32(Endian endian=DTFile::LittleEndian) const;
    unsigned int Read_uint32(Endian endian=DTFile::LittleEndian) const;
    ssize_t Read_int64(Endian endian=DTFile::LittleEndian) const;
    float Read_float(Endian endian=DTFile::LittleEndian) const;
    double Read_double(Endian endian=DTFile::LittleEndian) const;

    // Read in an array of entries.
    bool ReadBinary(DTMutableDoubleArray &A) const;
    bool ReadBinary(DTMutableDoubleArray &A,const DTRange &) const;
    bool ReadBinary(DTMutableFloatArray &A) const;
    bool ReadBinary(DTMutableIntArray &A) const;
    bool ReadBinary(DTMutableIntArray &A,const DTRange &) const;
    bool ReadBinary(DTMutableShortIntArray &A) const;
    bool ReadBinary(DTMutableUShortIntArray &A) const;
    bool ReadBinary(DTMutableUCharArray &A) const;
    bool ReadBinary(DTMutableUCharArray &A,ssize_t howMuchToRead) const;
    bool ReadBinary(DTMutableUCharArray &A,ssize_t startAt,ssize_t howMuchToRead) const;
    bool ReadBinary(DTMutableCharArray &) const;
    bool ReadBinary(DTMutableCharArray &,ssize_t howMuchToRead) const;
    
    // Read in Fortran based arrays (binary fortran)
    DTMutableDoubleArray ReadFortranBinary(); // Tries to deduce endian order automatically
    
    char CharacterAtCurrentPosition(void) const;
    
    // Ascii input.
    bool ReadAscii(DTMutableDoubleArray &A) const;
    bool ReadAscii(DTMutableFloatArray &A) const;

    double ReadAsciiNumber() const;

    // Binary output
    
    bool WriteUnsignedShort(unsigned short int);
#if defined(WIN32)
    bool Write8ByteInt(__int64);
#else
    bool Write8ByteInt(int64_t);
#endif
    bool WriteFloat(float v);
    bool WriteDouble(double v);
    bool Write4ByteInt(int);
    bool Write2ByteInt(short int);
    bool Write1ByteInt(char);
    bool WriteRaw(const char *,ssize_t howMany);
    bool WriteString(string); // Will not save \0
    bool WriteStringWithZero(string);

    bool WriteBinary(const DTDoubleArray &);
    bool WriteBinary(const DTFloatArray &);
    bool WriteBinary(const DTIntArray &);
    bool WriteBinary(const DTShortIntArray &);
    bool WriteBinary(const DTUShortIntArray &);
    bool WriteBinary(const DTCharArray &);
    bool WriteBinary(const DTUCharArray &);
    
    // Debugging info.
    void pinfo(void) const;
    
private:
    bool CheckWriteErrorState(const char *) const;    
    
    DTFileStorage *storage;
};

class DTFolder {
public:
    explicit DTFolder(const std::string &nm) : name(nm) {}
    
    std::string Name(void) const {return name;}
    DTFolder AppendFolderName(const std::string &);
    DTFile AppendFileName(const std::string &);
    
private:
    std::string name;
};

#endif
