/*
 * InputFile.h
 *
 *  Created on: Mar 2, 2017
 *      Author: mikolajr
 */

#ifndef SRC_INPUTFILE_H_
#define SRC_INPUTFILE_H_

#include <chrono>
// #include <set>
#include <string>
#include <vector>

#include "types.h"

/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of buffer used to read files

namespace SignificantPattern
{

    /// Abstract class for a vector or matrix of input data read from a file
    class ArrayFile
    {
    private:
        // last time the array was initialised (in ms since epoch)
        long long lastInitTime = 0;
        inline void markInit() {
            lastInitTime = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        }

    protected:

        // cleanup memory
        virtual void cleanupMemory() = 0;

        // static array & size getters
        virtual unsigned char *getArrayPtr() const = 0;
        virtual std::vector<longint> getArrayDimensions() const = 0;
        virtual longint getArraySize() const;
        inline bool isSameSizeArray(const std::vector<longint>& dimensions) const {
            return getArrayDimensions() == dimensions;
        }
        inline bool isSameSizeArray(const ArrayFile& other) const {
            return isSameSizeArray(other.getArrayDimensions());
        }

        // allocate array using size of the array from the other instance
        // ATTENTION: call the base class impl _at the end_ of the derived class impl
        virtual void allocArray(const std::vector<longint>& dimensions) = 0;
        void allocArray(const ArrayFile& other) {
            allocArray(other.getArrayDimensions());
        }
        void copyArray(const unsigned char *m, const std::vector<longint>& dimensions);
        void reallocArray(const std::vector<longint>& dimensions);
        virtual void initArray() = 0;

//        // initialise already allocated memory
//        inline void initialiseAllocatedMemory() { cleanupNonreusableMemory(); initArray(); };

        /// read from file
        void tryOpenFile(const std::string& filename, /*out*/ std::ifstream& file);
//        virtual void checkFile(const std::string& filename) = 0;
//        virtual void parseFile(const std::string& filename) = 0;

        /// read from file
        bool canReadFile(const std::string& filename);
//        void readFile(const std::string& filename);
        /// write to file
        virtual void writeFileStream(std::ofstream& file) const = 0;
        void writeFile(const std::string& filename) const;

    public:
        ArrayFile();
        virtual ~ArrayFile();

//        ArrayFile (const ArrayFile& other);
        virtual ArrayFile& operator=(const ArrayFile& other);

        // check if memory is initialised
        virtual int isInitialised() const = 0;
        inline long long getLastInitTime() const { return lastInitTime; }

    //    // get set of all items in the array
    //    std::set<char> getArrayItems();

    };

} /* namespace SignificantPattern */

#endif /* SRC_INPUTFILE_H_ */
