/*
 * InputFile.cpp
 *
 *  Created on: Mar 2, 2017
 *      Author: mikolajr
 */

#include "ArrayFile.h"

#include "Exception.h"

#include <fstream>

/* CONSTANT DEFINES */
#define PLINK_RAW_NCOL_META 6 // First 6 columns in PLINK .raw file are meta data (FID IID PAT MAT SEX PHENOTYPE)


using namespace std;

namespace SignificantPattern
{
    // Constructor
    ArrayFile::ArrayFile ()
    {
    }
    // Destructor
    ArrayFile::~ArrayFile ()
    {
    }

//    // Copy constructor
//    InputFile::InputFile (const InputFile& other)
//    {
//        // call copy assignment operator
//        *this = other;
//    }
    // Copy assignment operator.
    ArrayFile& ArrayFile::operator=(const ArrayFile& other)
    {
       if (this != &other) {
           if (!other.isInitialised()) {
               cleanupMemory();
           } else {
               reallocArray(other.getArrayDimensions());
               copyArray(other.getArrayPtr(), other.getArrayDimensions());
           }
       }
       return *this;
    }

    longint ArrayFile::getArraySize() const
    {
        const std::vector<longint>& dimensions = getArrayDimensions();
        longint msize = 1;
        for(size_t i = 0; i < dimensions.size(); i++) msize *= dimensions[i];
        return msize;
    }
    void ArrayFile::copyArray(const unsigned char *m, const std::vector<longint>& dimensions)
    {
        if (!isSameSizeArray(dimensions))
            throw Exception("Can't copy memory with inconsistent sizes");
        std::copy(m, m + getArraySize(), getArrayPtr());
    }
    void ArrayFile::initArray()
    {
        #ifdef DEBUG
        fprintf(stderr,"ArrayFile::initArray()\n");
        #endif
        markInit();
    }
    void ArrayFile::allocArray(const std::vector<longint>& dimensions)
    {
        #ifdef DEBUG
        fprintf(stderr,"ArrayFile::allocArray(), dim=%lu\n", dimensions.size());
        #endif
        initArray();
    }
    void ArrayFile::reallocArray(const std::vector<longint>& dimensions)
    {
        // re-allocate memory only if necessary
        if (!isSameSizeArray(dimensions)) {
            cleanupMemory();
            allocArray(dimensions);
        } else {
            initArray();
        }
    }



    bool ArrayFile::canReadFile(const std::string& filename)
    {
        std::ifstream file(filename.c_str()); // opens file for read w/o error
        return file.is_open(); // destructor closes the file
    }

    void ArrayFile::tryOpenFile(const std::string& filename, /*out*/ std::ifstream& file)
    {
        //Try to open file, throwing an error if it fails
        file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
        try
        {
            file.open(filename.c_str());
        }
        catch (const std::ios_base::failure& e)
        {
            throw Exception("Failed opening " + filename);
        }
        file.exceptions ( std::ifstream::goodbit );
    }

//    void InputFile::readFile(const std::string& filename)
//    {
//        checkFile(filename);
//        initialiseMemory();
//        parseFile(filename);
//    }

    void ArrayFile::writeFile(const std::string& filename) const
    {
        if (!isInitialised()) {
            throw Exception("Nothing to write.");
        }
        std::ofstream file(filename.c_str());
        writeFileStream(file);
    }


    // std::set<char> ArrayFile::getArrayItems() {
    //     unsigned char *ptr = getArrayPtr();
    //     return std::set<char>(ptr, ptr+getArraySize());
    // }

} /* namespace SignificantPattern */
