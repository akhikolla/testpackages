/*
 * Genotype.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: mabaker
 */

#include "Genotype.h"
#include "Exception.h"

#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>

/* CONSTANT DEFINES */
#define PLINK_RAW_NCOL_META 6 // First 6 columns in PLINK .raw file are meta data (FID IID PAT MAT SEX PHENOTYPE)

using namespace std;

namespace SignificantPattern
{
    // Constructor
    Genotype::Genotype ()
    {
        setMatrixPtr(0);
        setNumFeatures(0); setNumObservations(0);
    }
    // Copy constructor
    Genotype::Genotype (const Genotype& other)
    {
        // call copy assignment operator
        *this = other;
    }
    //TODO move constructor/operator
    // Copy assignment operator
    Genotype& Genotype::operator=(const Genotype& other)
    {
        super::operator=(other);
        return *this;
    }
    // Destructor
    Genotype::~Genotype ()
    {
        cleanupMemory();
    }

    void Genotype::cleanupMemory()
    {
        int n_dim = isInitialised();
        if (n_dim > 0) {
            unsigned char **mat = getMatrixPtr();
            if (n_dim > 1) delete [] mat[0];
            delete [] mat;
            setMatrixPtr(0);
        }
        setNumFeatures(0); setNumObservations(0);
    }

    int Genotype::isInitialised() const {
        unsigned char **mat = getMatrixPtr();
        return (mat) ? ((mat[0]) ? 2 : 1) : 0;
    }

    void Genotype::allocArray(const std::vector<longint>& dimensions) {
        #ifdef DEBUG
        fprintf(stderr, "Genotype::allocArray(): BEGIN, current data size: LxN=%ldx%ld, mem addr = %p\n", getNumFeatures(), getNumObservations(), (void *) getArrayPtr());
        #endif
        longint L = dimensions[0], N = dimensions[1];
        setMatrixPtr(new unsigned char *[L]);
        unsigned char **mat = getMatrixPtr();
        mat[0] = new unsigned char[L*N];
        setNumFeatures(L); setNumObservations(N);
        super::allocArray(dimensions);
        #ifdef DEBUG
        fprintf(stderr, "Genotype::allocArray(): END, current data size: LxN=%ldx%ld, mem addr = %p\n", getNumFeatures(), getNumObservations(), (void *) getArrayPtr());
        #endif
    }

    void Genotype::initArray() {
        #ifdef DEBUG
        fprintf(stderr, "Genotype::initArray()\n");
        #endif
        const longint L = getNumFeatures(), N = getNumObservations();
        unsigned char **mat = getMatrixPtr();
        std::fill_n(mat[0], L*N, 0);
        for(longint j=1; j<L; j++)
            mat[j] = mat[0] + j*N;
        super::initArray();
    }

    void Genotype::initialiseMatrix(longint L, longint N)
    {
        reallocArray({L, N});
    }

    void Genotype::readETHFile(const std::string& filename, longint N, const std::string& encoding)
    {
        #ifdef DEBUG
        fprintf(stderr, "Genotype::readETHFile(): BEGIN, current data size: LxN=%ldx%ld\n", getNumFeatures(), getNumObservations());
        #endif
        longint L;
        checkEthDataFile(filename, N, L);
        initialiseMatrix(L, N);
        #ifdef DEBUG
        fprintf(stderr, "Genotype::readETHFile(): pre-parse, current data size: LxN=%ldx%ld\n", getNumFeatures(), getNumObservations());
        #endif
        parseEthDataFile(filename, getArrayPtr(), encoding);
        #ifdef DEBUG
        fprintf(stderr, "Genotype::readETHFile(): END, current data size: LxN=%ldx%ld\n", getNumFeatures(), getNumObservations());
        #endif
    }

    void Genotype::readPlinkRawFile(const std::string& filename, const Phenotype& phenotype)
    {
        longint N = phenotype.getNumObservations();
        bool hasHeader; longint L;
        checkPlinkRawFile(filename, phenotype, /*out*/ L, hasHeader);
        initialiseMatrix(L, N);
        parsePlinkRawFile(filename, hasHeader, phenotype, /*out*/ getArrayPtr());
    }

    void Genotype::checkPlinkRawFile(const std::string& filename, const Phenotype& phenotype, /*out*/ longint& L, bool& hasHeader)
    {
        const longint N_expected = phenotype.getNumObservations();
        hasHeader = false;
        L = 0;
        ifstream file; tryOpenFile(filename, file);
        longint lineno = 0; // actual line number, incl. empty lines, for error reporting
        std::string line;
        std::getline(file, line); lineno++;
        if (line.size() == 0) return;

        longint N_read = 0;
        if (line.substr(0,3) == "FID") {
            hasHeader = true;
        } else {
            N_read++;
        }

        // count white space separated columns of the first line
        const std::string delimiters(" \t\r");
        std::string::size_type pos, lastPos = 0, length = line.length();
        std::stringstream errmsgstream;
           while(lastPos < length + 1)
           {
              pos = line.find_first_of(delimiters, lastPos);
              if(pos == std::string::npos) pos = length;

              if (pos == lastPos) {
                  errmsgstream <<  "In data file, line " << lineno << ", unexpected number of columns: " << L;
                  throw Exception(errmsgstream.str());
              }


              lastPos = pos + 1;
              L++;
           }
        L -= PLINK_RAW_NCOL_META;

        // read and split lines for format check; start w/ read line 1
        std::string fid, iid, pat, mat, sex;
        long j; // target line number
        std::vector<longint> rownrLabelsToData(N_expected, -1);
        short phenotypeCode; vector<short> variants(L);
        if (!hasHeader)
            splitRawLine(line, fid, iid, pat, mat, sex, phenotypeCode, variants, N_read+hasHeader);
        bool skippingTrailingEmptyLines = false;
        while (!file.eof()) {
            std::getline(file, line); lineno++;
            if (line.size() > 0) {
                if (!skippingTrailingEmptyLines) {
                    N_read++;
                    // Rem: if read more than N_expected, then either there is
                    // gonna be duplicate FID, IID entry or no matching one
                    splitRawLine(line, fid, iid, pat, mat, sex, phenotypeCode, variants, N_read+hasHeader);
                    j = phenotype.getPlinkLineForFIDAndIID(fid, iid);
                    // fprintf(stderr, "[TEMP] Genotype::checkPlinkRawFile(), key=%s %s => row=%ld\n", fid.c_str(), iid.c_str(), j);
                    if (j < 0) {
                        errmsgstream << "In data file, line " << lineno << ", no matching FID and IID line in labels file.";
                        throw Exception(errmsgstream.str());
                    }
                    if (rownrLabelsToData[j] >= 0) {
                        errmsgstream << "In data file, line " << lineno << ", duplicate FID and IID : \"" << fid << "\", \"" << iid << "\" .";
                        throw Exception(errmsgstream.str());
                    }
                    rownrLabelsToData[j] = N_read-1;
                } else if (line.size() > 0) {
                    errmsgstream << "In data file, line " << lineno <<  ", non-empty trailing line.";
                    throw Exception(errmsgstream.str());
                }
            } else {
                skippingTrailingEmptyLines = true;
            }
        }
        if (N_read < N_expected) {
            errmsgstream << "In data file, line " << lineno <<  ", less than expected " << N_expected << " lines of samples.";
            throw Exception(errmsgstream.str());
        }

        file.close();
    }

    void Genotype::parsePlinkRawFile(const std::string& filename, bool hasHeader, const Phenotype& phenotype, /*out*/ unsigned char *data_buf) {

        std::ifstream file; tryOpenFile(filename, file);
        std::string line;
        std::string fid, iid, pat, mat, sex;
        long j; // target line number
        vector<short> variants(L);
        short phenotypeCode;
        if (hasHeader) std::getline(file, line);
        for (longint lineno=0; lineno<N; ++lineno)
        {
            std::getline(file, line);
            splitRawLine(line, fid, iid, pat, mat, sex, phenotypeCode, variants, lineno);
            j = phenotype.getPlinkLineForFIDAndIID(fid, iid);
            // should be already checked in checkPlinkRawFile method
            assert(j >= 0 && j < N);
            for (long i=0; i<L; ++i)
                data_buf[N*i+j] = variants[i];
        }
        file.close();

    }


    //LP: add function. To check when PLINK function is ready.
    // void Genotype::encodingProcessedPlinkFile(std::vector<short>& variants, const std::string& encoding){
    // 	std::vector<int>::iterator it = variants.begin();
    // 	char curr_char; short curr_val;
    // 	if (encoding == "dominant"){
    // 	for (; it!= variants.end(); ++it){
    // 		curr_char = *it;
    // 		switch (curr_char) {
    // 		case '0': curr_val = 0;
    //				break;
    //	case '1':
    //	case '2':
    //		curr_val = 1;
    //	        break;
    //	}
    //	variants[it-PLINK_RAW_NCOL_META] = curr_val;
    // 	}
    // 	}

    //	if (encoding == "recessive"){
    //	for (; it!= variants.end(); ++it){
    //		curr_char = *it;
    //		switch (curr_char) {
    //		case '0':
    //		case '1':
    //			curr_val = 0;
    //				break;
    //		case '2':
    //			curr_val = 1;
    //		        break;
    //		}
    //		variants[it-PLINK_RAW_NCOL_META] = curr_val;
    //	}
    //	}
    // }

    void Genotype::splitRawLine(const std::string& line,
        std::string& fid, std::string& iid,
        std::string& pat, std::string& mat, std::string& sex,
        short& phenotype, std::vector<short>& variants, long lineno)
    {
        /**
         * Format reference: https://www.cog-genomics.org/plink2/formats#raw
         * Format limitation: assuming only one field per SNP variant
         */
        const std::string delimiters(" \t\r");
        std::string::size_type pos, lastPos = 0, length = line.length();
        std::stringstream errmsgstream;
           longint i = 0; // value (column) number
           longint L_plus = variants.size()+PLINK_RAW_NCOL_META; // expected number of values (columns) w/ metadata
           char curr_char; short curr_val;
           while(lastPos < length + 1)
           {
              pos = line.find_first_of(delimiters, lastPos);
              if(pos == std::string::npos) pos = length;

              if (pos == lastPos) {
                  errmsgstream << "Unexpected number of values in line " << lineno  << ".  Expected " << L_plus;
                  throw Exception(errmsgstream.str());
              }
            switch (i)
            {
                case 0:
                    fid = string(line.data()+lastPos, (int)pos-lastPos);
                    break;
                case 1:
                    iid = string(line.data()+lastPos, (int)pos-lastPos);
                    break;
                case 2:
                    pat = string(line.data()+lastPos, (int)pos-lastPos);
                    break;
                case 3:
                    mat = string(line.data()+lastPos, (int)pos-lastPos);
                    break;
                case 4:
                    sex = string(line.data()+lastPos, (int)pos-lastPos);
                    break;
                case 5:
                    curr_char = *(line.data()+lastPos);
                    switch (curr_char) {
                        case '1':
                            curr_val = 0;
                            break;
                        case '2':
                            curr_val = 1;
                            break;
                        default:
                            errmsgstream << "In data file, line " << lineno << ", col " << (lastPos+1) << ", invalid phenotype character '" << curr_char << "'";
                            throw Exception(errmsgstream.str());
                            break;
                    }
                    phenotype = curr_val;
                    break;
                default:
                    curr_char = *(line.data()+lastPos);
                    switch (curr_char) {
                        case '0':
                            curr_val = 0;
                            break;
                        case '1':
                        case '2': // collapse 1 and 2 genotype values into 1
                            curr_val = 1;
                            break;
                        default: // valid 'NA' will also raise an error here
                            errmsgstream << "In data file, line " << lineno << ", col " << (lastPos+1) << ", invalid variant character '" << curr_char << "'";
                            throw Exception(errmsgstream.str());
                            break;
                    }
                    variants[i-PLINK_RAW_NCOL_META] = curr_val;
                    break;
            }

              lastPos = pos + 1;
              i++;
              if (i > L_plus) {
                  errmsgstream << "In data file, line " << lineno << ", too many values; expected " << L_plus;
                  throw Exception(errmsgstream.str());
              }
           }
           if (i < L_plus) {
               errmsgstream << "In data file, line " << lineno << ", too few values; expected " << L_plus;
               throw Exception(errmsgstream.str());
           }

    }

    void Genotype::writeFileStream(std::ofstream& file) const
    {
        unsigned char *arr = getArrayPtr();
        longint N = getNumObservations();
        for (longint i=0; i<getNumFeatures(); ++i)
        {
            for (longint j=0; j<N; ++j)
            {
                file << (unsigned short) arr[i*N+j];
                if (j < N-1) file << ' ';
            }
            file << endl;
        }
    }

    void Genotype::writeETHFile(const std::string& filename)
    {
        writeFile(filename);
    }

    void Genotype::checkEthDataFile(const std::string& filename, longint N_expected, /*out*/ longint& L){
        ifstream f_dat;
        int n_read;
        unsigned char char_to_int[256];//Array for converting chars to int fast
        char *read_buf_start ,*read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
        unsigned char read_buf_uchar;
        // We only allow the following chars; everything else is invalid
        // new line
        const char new_line_int = 3;
        unsigned char new_line = '\n';
        // values separators (separators are optional)
        const char sep_int = 2;
        unsigned char sep_space = ' '; unsigned char sep_tab = '\t'; unsigned char sep_comma = ',';
        // values
        const char val_int = 1;
        unsigned char val_zero = '0'; unsigned char val_one = '1';
        // everything else
        const char rest_int = 0;

        tryOpenFile(filename, f_dat);

        //Try to allocate memory for the buffer, giving an error message if it fails
        std::string read_buf;
        read_buf.resize(READ_BUF_SIZ);
        read_buf_start = &read_buf[0];

        //Initialize the char to int converter
        std::fill_n(char_to_int, 256, rest_int);
        char_to_int[val_zero] = val_int; char_to_int[val_one] = val_int;
        char_to_int[new_line] = new_line_int;
        char_to_int[sep_space] = sep_int; char_to_int[sep_tab] = sep_int; char_to_int[sep_comma] = sep_int;

        // Read the entire file, counting the number of lines and checking chars
        L = 1; // line number
        longint N_line = 0; // value number
        longint col_line = 0; // column number
        std::stringstream errmsgstream;
        while(1){
            // Try to read READ_BUF_SIZ chars from the file containing the class labels
            f_dat.read(read_buf_start, READ_BUF_SIZ);
            n_read = f_dat.gcount();
            // If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
            // or there was an error. Check if it was the latter
            if((n_read < READ_BUF_SIZ) && !f_dat.eof()){
                    throw std::runtime_error("Error while checking data file " + filename);
            }
            // Process the n_read chars read from the file
            for (read_buf_aux = read_buf_start, read_buf_end = read_buf_start + n_read; read_buf_aux < read_buf_end; read_buf_aux++) {
                col_line++;
                read_buf_uchar = (unsigned char) *read_buf_aux;
                switch (char_to_int[read_buf_uchar]) {
                    case new_line_int:
                        if (N_line != N_expected) {
                            errmsgstream << "In data file, line " << L << ", found " << N_line << " values whereas " << N_expected << " labels were read";
                            throw Exception(errmsgstream.str());
                        }
                        L++;
                        N_line = 0;
                        col_line = 0;
                        break;
                    case val_int:
                        N_line++;
                        break;
                    case sep_int: // skip over any optional separators
                        break;
                    case rest_int: // report error
                        errmsgstream << "In data file, line " << L << ", col " << col_line << ", invalid character '" << read_buf_uchar << "'";
                        throw Exception(errmsgstream.str());
                        break;
                    default: // shouldn't happen
                        errmsgstream << "Error while checking data file " << filename << ", row " << L << ", col " << col_line;
                        throw std::runtime_error(errmsgstream.str());
                        break;
                }
            }
            // Check if the file ended,. If yes, then exit the while loop
            if (f_dat.eof()) {
                L--; // eof == empty line
                break;
            }
        }

        // Close file
        f_dat.close();
    }


    void Genotype::parseEthDataFile(const std::string& filename, unsigned char *data_buf, const std::string& encoding){
        // Just parsing 0 and 1s, ignoring everything else (no format validation)
        ifstream f_dat ;
        int n_read;
        unsigned char char_to_int[256];//Array for converting chars to int fast
        char *read_buf_start, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
        unsigned char read_buf_uchar; unsigned char read_buf_int;
        unsigned char zero = '0'; unsigned char one = '1'; unsigned char two = '2';
        tryOpenFile(filename, f_dat);

        //Try to allocate memory for the buffer, giving an error message if it fails
        std::string read_buf;
        read_buf.resize(READ_BUF_SIZ);
        read_buf_start = &read_buf[0];

        //Initialize the char to int converter
        std::fill_n(char_to_int, 256, 127);
        // We only care about the chars '0', '1' and '2'. Everything else is mapped into the same "bucket" \\LP: dominant/recessive
        if (encoding == "dominant"){
            char_to_int[zero] = 0; char_to_int[one] = 1; char_to_int[two] = 1;
        }
        if (encoding == "recessive"){
            char_to_int[zero] = 0; char_to_int[one] = 0; char_to_int[two] = 1;
        }
        // Read the entire file
        while(1)
        {
            // Try to read READ_BUF_SIZ chars from the file containing the class labels
            f_dat.read(read_buf_start, READ_BUF_SIZ);
            n_read = f_dat.gcount();
            // If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
            // or there was an error. Check if it was the latter
            if((n_read < READ_BUF_SIZ) && !f_dat.eof()){
                    throw std::runtime_error("Error while parsing data file" + filename);
            }
            // Process the n_read chars read from the file
            for(read_buf_aux=read_buf_start,read_buf_end=read_buf_start+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
                read_buf_uchar = (unsigned char) *read_buf_aux;
                read_buf_int = char_to_int[read_buf_uchar];
                //If the character is anything other than '0' or '1' go to process the next char
                if(read_buf_int != 127) {
                    *data_buf++ = read_buf_int;
                }
            }
            // Check if the file ended,. If yes, then exit the while loop
            if(f_dat.eof()) break;
        }

        // Close file
        f_dat.close();
    }

} /* namespace SignificantPattern */
