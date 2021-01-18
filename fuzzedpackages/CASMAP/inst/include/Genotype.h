/*
 * Genotype.h
 *
 *  Created on: Aug 25, 2016
 *      Author: mabaker
 */

#ifndef SRC_GENOTYPE_H_
#define SRC_GENOTYPE_H_

#include <string>
#include <vector>

#include "ArrayFile.h"
#include "Phenotype.h"
#include "types.h"

namespace SignificantPattern
{

    class Genotype : public ArrayFile
    {
    private:
        // super class pattern for code independence of changes in inheritance
        typedef ArrayFile super;

        /// Number of features
        longint L;
        /// Number of observations
        longint N;
        // The original dataset matrix stored as a LxN matrix, where L is the sequence length
        // and N the number of observed sequences. This exploits the row-major storage of C matrices.
        unsigned char **X_tr;

        void checkEthDataFile(const std::string& filename, longint N_expected, /*out*/ longint& L);
        void parseEthDataFile(const std::string& filename, /*out*/ unsigned char *data_buf, const std::string& encoding);
        void checkPlinkRawFile(const std::string& filename, const Phenotype& phenotype, /*out*/ longint& L, bool& hasHeader);
        void parsePlinkRawFile(const std::string& filename, bool hasHeader, const Phenotype& phenotype, /*out*/ unsigned char *data_buf);
        void splitRawLine(const std::string& line,
            std::string& fid, std::string& iid,
            std::string& pat, std::string& mat, std::string& sex,
            short& phenotype, std::vector<short>& variants, long lineno);

    protected:
        inline void setNumFeatures(longint l) { L = l; }
        inline void setNumObservations(longint n) { N = n; }
        inline void setMatrixPtr(unsigned char **tr) { X_tr = tr; }

        inline unsigned char *getArrayPtr() const override {
            return isInitialised() ? getMatrixPtr()[0] : NULL;
        }
        inline std::vector<longint> getArrayDimensions() const override {
            return {getNumFeatures(), getNumObservations()};
        }
        void cleanupMemory() override;
        virtual void allocArray(const std::vector<longint>& dimensions) override;
        virtual void initArray() override;


//        // read from file
//        void checkFile(const std::string& filename) override;
//        void parseFile(const std::string& filename) override;
        // write to file
        void writeFileStream(std::ofstream& file) const override;

    public:
        Genotype ();
        Genotype (const Genotype& other);
        virtual Genotype& operator=(const Genotype& other);
        virtual ~Genotype ();

        int isInitialised() const override;

        /// matrix init (mem allocation)
        void initialiseMatrix(longint L, longint N);
        /// read from file
        void readETHFile(const std::string& filename, longint N, const std::string& encoding);
        void readPlinkRawFile(const std::string& filename, const Phenotype& phenotype);
        /// write to file
        void writeETHFile(const std::string& filename);

        /// returns the number of columns of the matrix
        inline longint getNumObservations() const { return N; }
        /// returns the number of rows of the matrix
        inline longint getNumFeatures() const { return L; }
        /// returns the 2-D matrix pointer (to row pointers)
        inline unsigned char **getMatrixPtr() const { return X_tr; }
};

} /* namespace SignificantPattern */

#endif /* SRC_GENOTYPE_H_ */
