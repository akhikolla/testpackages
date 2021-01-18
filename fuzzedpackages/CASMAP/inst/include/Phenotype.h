/*
 * Phenotype.h
 *
 *  Created on: Aug 25, 2016
 *      Author: mabaker
 */

#ifndef SRC_PHENOTYPE_H_
#define SRC_PHENOTYPE_H_

#include <limits> // numeric_limits
#include <map>
#include <string>
#include <vector>

#include "ArrayFile.h"
#include "types.h"

namespace SignificantPattern
{

    /// Vector of phenotype data read from an ETH or PLINK file
    /**
     * Phenotype data is stored as a vector of zeroes and ones
     */
    class Phenotype  : public ArrayFile
    {
    private:
        /**
         * super class pattern for code independence of changes in inheritance
         */
        typedef ArrayFile super;

        /**
         * @name CONSTANTS
         */
        /**
         * Maximum number of classes (labels) for observations. Limited by the
         * unsigned char type.
         */
        static constexpr unsigned short MAX_NUMCLASS = std::numeric_limits<unsigned char>::max();
        /**
         * Mininum number of columns in the PLINK format (FID, IID, father,
         * mother, sex, phenotype, value, ...)
         */
        static constexpr unsigned short PLINK_MINNCOL = 6;
        /**
         * Minimum number of columns in the PLINK short format (FID, IID, value,
         * ...)
         */
        static constexpr unsigned short PLINK_MINNCOL_SHORT = 2;
        ///@}

        /**
         * Number of observations (length of the vector)
         */
        longint N;
        /**
         * Vector of class labels.
         */
        unsigned char *Y_tr;
        /**
         * Number of observations in each class.
         */
        std::vector<longint> nv;
        //
        /**
         * Map of FID + IID to line number to match lines with other input files
         * in the PLINK format.
         *
         * Tech: constructor and destructor for the map is auto-inserted into
         *       object's constructor and destructor, respectively.
         */
        std::map<std::string, longint> fidiid2lineno;

        /**
         * Get FID + IID key string for the #fidiid2lineno mapping.
         *
         * @param fid PLINK FID
         * @param iid PLINK IID
         * @return #fidiid2lineno map key
         */
        inline const std::string getFIDAndIIDKeyStr(const std::string& fid, const std::string& iid) const {
            return fid + " " + iid;
        }

        void checkNumObservations(const std::string &filename, longint N,
                                  longint N_expected);
        unsigned short parseTargetValue(const std::string &value_str,
                                        const long lineno, const long colno,
                                        const unsigned short offset_value,
                                        const unsigned short max_value);

        /**
         * Do a scan of the file containing the class labels to check the file
         * format and to compute the total number of observations for further
         * memory allocation.
         *
         * @param filename File name
         * @param N Input parameter to for `setNumObservations()` call
         * @param d Input parameter for `setNumClasses()` call
         */
        void checkEthLabelsFile(const std::string &filename, /*out*/ longint &N,
                                unsigned short &d);
        void parseEthLabelsFile(const std::string &filename, longint N,
                                /*out*/ unsigned char *labels_buf,
                                std::vector<longint> &nv);

        void readPlinkLabelsFile(const std::string &filename,
                                 longint N_expected,
                                 const Phenotype *const otherOrderingPtr,
                                 const bool short_format,
                                 const unsigned short value_ncol,
                                 const unsigned short offset_value,
                                 const unsigned short max_value);
        void checkPlinkLabelsFile(const std::string &filename,
                                  const Phenotype *const otherOrderingPtr,
                                  const bool short_format,
                                  const unsigned short value_ncol,
                                  const unsigned short offset_value,
                                  const unsigned short max_value,
                                  /*out*/ longint &N, unsigned short &d,
                                  std::map<std::string, longint> &id2n);
        void parsePlinkLabelsFile(const std::string &filename, longint N,
                                  const bool short_format,
                                  const unsigned short value_ncol,
                                  const unsigned short offset_value,
                                  const unsigned short max_value,
                                  /*out*/ unsigned char *labels_buf,
                                  std::vector<longint> &nv);
        /**
         * Split a single line of a PLINK labels file (FAM, or COV format), and
         * parse a targe non-negative (short) integer value.
         *
         * @param line Line contents
         * @param lineno Line number (for errors reporting)
         * @param short_format Boolean flag indicating a short format (w/ at
         *                     least #PLINK_MINNCOL_SHORT columns), or a long
         *                     format (w/ at least #PLINK_MINNCOL columns)
         * @param value_ncol Number of a column w/ a target integer value to
         *                   parse (e.g. #PLINK_MINNCOL to read phenotype
         *                   integer in FAM file)
         * @param offset_value Target value offset, which will be checked and
         *                     the deducted from the target value (e.g. 1 for
         *                     phenotype in FAM file)
         * @param max_value Maximal allowed value; if 0 then maximal value of
         *                  the target type (+ offset_value)
         * @param fid set to FID column string
         * @param iid set to IID column string
         * @param pat set to PAT column string or "" if short_format
         * @param mat set to MAT column string or "" if short_format
         * @param sex set to SEX column string or "" if short_format
         * @param phenotype set to PHENOTYPE column string or "" if short_format
         * @param value set to value_ncol column (short) unsigned integer
         */
        void splitLabelsLine(const std::string &line, const long lineno,
                             const bool short_format,
                             const unsigned short value_ncol,
                             const unsigned short offset_value,
                             const unsigned short max_value,
                             /*out*/ std::string &fid, std::string &iid,
                             std::string &pat, std::string &mat,
                             std::string &sex, std::string &phenotype,
                             unsigned short &value);

      protected:
        inline void setNumObservations(longint n) { N = n; }
        inline void setVectorPtr(unsigned char *tr) { Y_tr = tr; }
        /**
         * Setter of #nv, used for memory preallocation.
         *
         * Tech: using C++11 move semantics, i.e. arg w/o reference (&) and an
         * explicit std::move call; cf. http://stackoverflow.com/a/21227582.
         *
         * @param nv_new New #nv value
         */
        inline void setNumObservationsInClasses(std::vector<longint> nv_new) {
            nv = std::move(nv_new);
        }
        inline std::vector<longint> &getNumObservationsInClassesRef() { return nv; }
        void setNumClasses(unsigned short d);
        inline bool isSetNumClasses() { return getNumClasses() != 0; }
        inline void resetNumClasses() { setNumClasses(getNumClasses()); }
        /**
         * Setter of #fidiid2lineno, used for memory cleanup.
         *
         * Tech: using C++11 move semantics, i.e. arg w/o reference (&) and an
         * explicit std::move call; cf. http://stackoverflow.com/a/21227582.
         *
         * @param fidiid2lineno_new New #fidiid2lineno_new value
         */
        inline void setPlinkFIDAndIIDToLineMap(const std::map<std::string, longint> fidiid2lineno_new) {
            fidiid2lineno = std::move(fidiid2lineno_new);
        }
        inline std::map<std::string, longint> const &getPlinkFIDAndIIDToLineMap() const {
            return fidiid2lineno;
        }
        inline void resetPlinkFIDAndIIDToLineMap() { fidiid2lineno.clear(); }

        inline unsigned char *getArrayPtr() const override {
            return isInitialised() ? getVectorPtr() : NULL;
        }
        inline std::vector<longint> getArrayDimensions() const override {
            return {getNumObservations()};
        }
        void cleanupMemory() override;
        void resetNonreusableMemory();
        void copyNonreusableMemory(const Phenotype& other);
        virtual void allocArray(const std::vector<longint>& dimensions) override;
        virtual void initArray() override;

        /**
         * Write to a file stream.
         *
         * @param file Reference to an opened file stream to write to.
         */
        void writeFileStream(std::ofstream& file) const override;

    public:
        Phenotype ();
        Phenotype (const Phenotype& other);
        virtual Phenotype& operator=(const Phenotype& other);
        virtual ~Phenotype ();

        int isInitialised() const override;

        /**
         * @name Memory allocation
         */
        ///@{
        void initialiseMatrix(longint N);
        void initialiseMatrix(longint N, unsigned short value);
        ///@}
        /**
         * @name Read from file
         */
        ///@{
        void readETHFile(const std::string& filename, longint N_expected = 0);
        void readPlinkFamFile(const std::string& filename, longint N_expected = 0);
        void readPlinkFamFile(const std::string& filename, const Phenotype& otherOrdering);
        void readPlinkCovFile(const std::string& filename, longint N_expected = 0);
        void readPlinkCovFile(const std::string& filename, const Phenotype& otherOrdering);
        ///@}
        /**
         * @name Write to file
         */
        ///@{
        void writeETHFile(const std::string& filename) const;
        ///@}
        /**
         * Query for the PLINK line number given FID and IID. Returns `-1` if
         * line with given IDs has not been found, including case when
         * #hasReadPlink() is `false`.
         *
         * @param fid PLINK FID
         * @param iid PLINK IID
         * @return Line number or `-1`
         */
        inline long getPlinkLineForFIDAndIID(const std::string& fid, const std::string& iid) const {
            auto it = fidiid2lineno.find(getFIDAndIIDKeyStr(fid, iid));
            if (it == fidiid2lineno.end()) return -1;
            else return it->second;
        }
        /**
         * Test if PLINK file has been read.
         *
         * @return `false` if no file, or file in a format different than PLINK
         *         has been read
         */
        inline bool hasReadPlink() const {
            return not fidiid2lineno.empty();
        }

       /**
         * Get length of the labels (phenotype) vector
         */
        inline longint getNumObservations() const { return N; }
        /**
         * Get number of positive observations in the labels (phenotype) vector.
         */
        inline unsigned short getNumClasses() const { return nv.size(); }
        /**
         * Read-only #nv getter.
         *
         * @return constant reference to the vector
         */
        inline std::vector<longint> const &getNumObservationsInClasses() const { return nv; }
        /**
         * Get pointer to the labels (phenotype) vector.
         */
        inline unsigned char *getVectorPtr() const { return Y_tr; }

    };

} /* namespace SignificantPattern */

#endif /* SRC_PHENOTYPE_H_ */
