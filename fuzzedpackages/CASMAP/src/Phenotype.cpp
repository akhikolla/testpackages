/*
 * Phenotype.cpp
 *
 *  Created on: Aug 25, 2016
 *      Author: mabaker
 */

#include "Phenotype.h"
#include "Exception.h"

#include <limits> // numeric_limits
#include <fstream> //ofstream
#include <sstream> //stringstream
#include <string> //stoi



using namespace std;

namespace SignificantPattern
{

/**
 * Deafult constructor
 */
Phenotype::Phenotype ()
{
    setVectorPtr(0);
    setNumObservations(0);
    setNumClasses(0);
}
/// @todo
/// move constructor/operator
/**
 * Copy constructor.
 *
 * Calls copy assignment operator.
 */
Phenotype::Phenotype (const Phenotype& other)
{
    *this = other;
}
/**
 * Copy assignment operator.
 */
Phenotype& Phenotype::operator=(const Phenotype& other)
{
    super::operator=(other);
    if (this != &other) {
        if (!other.isInitialised()) resetNonreusableMemory();
        else copyNonreusableMemory(other);
    }
    return *this;
}
/**
 * Default destructor
 */
Phenotype::~Phenotype ()
{
    cleanupMemory();
    // Note: cleanupMemory() is called in ~ArrayFile()
    resetNonreusableMemory();
}

void Phenotype::setNumClasses(unsigned short d) {
    if (d > MAX_NUMCLASS)
        throw Exception("Unsupported number of labels (too many).");
    nv.resize(d);
    std::fill(nv.begin(),nv.end(),0); // Note: nv.resize(d,0) will set 0 value only for new elements
}
void Phenotype::cleanupMemory() {
    if (isInitialised()) {
        unsigned char *vec = getVectorPtr();
        delete [] vec;
        setVectorPtr(0);
    }
    setNumObservations(0);
}

void Phenotype::resetNonreusableMemory() {
    resetNumClasses();
    resetPlinkFIDAndIIDToLineMap();
}
void Phenotype::copyNonreusableMemory(const Phenotype& other) {
    setNumClasses(other.getNumClasses());
    setNumObservationsInClasses(other.getNumObservationsInClasses());
    setPlinkFIDAndIIDToLineMap(other.getPlinkFIDAndIIDToLineMap());
}

void Phenotype::allocArray(const std::vector<longint>& dimensions) {
    longint N = dimensions[0];
    setVectorPtr(new unsigned char[N]);
    setNumObservations(N);
    super::allocArray(dimensions);
}

void Phenotype::initArray() {
    longint N = getNumObservations();
    #ifdef DEBUG
    fprintf(stderr, "Phenotype::initArray(), N=%ld\n", N);
    #endif
    std::fill_n(getVectorPtr(), N, 0);
    setNumObservationsInClasses({N});
    super::initArray();
}

void Phenotype::initialiseMatrix(longint N)
{
    reallocArray({N});
}

int Phenotype::isInitialised() const {
    return getVectorPtr() ? 1 : 0;
}

void Phenotype::checkNumObservations(const std::string& filename, longint N, longint N_expected) {
    std::stringstream errmsgstream;
    if (N_expected > 0 && N != N_expected) {
        errmsgstream << "Error while checking '"<< filename << "' file: expected " << N_expected << " labels, found " << N;
        throw Exception(errmsgstream.str());
    }
}

void Phenotype::readETHFile(const std::string& filename, longint N_expected)
{
    longint N; unsigned short d;
    checkEthLabelsFile(filename, /*out*/ N, d);
    checkNumObservations(filename, N, N_expected);
    initialiseMatrix(N); setNumClasses(d);
    parseEthLabelsFile(filename, N, /*out*/ getArrayPtr(), getNumObservationsInClassesRef());
}

void Phenotype::readPlinkFamFile(const std::string& filename, longint N_expected) {
    readPlinkLabelsFile(filename, N_expected,
                        NULL, false, PLINK_MINNCOL,
                        1, MAX_NUMCLASS);
}
void Phenotype::readPlinkFamFile(const std::string& filename, const Phenotype& otherOrdering) {
    readPlinkLabelsFile(filename, otherOrdering.getNumObservations(),
                        &otherOrdering, false, PLINK_MINNCOL,
                        1, MAX_NUMCLASS);
}
void Phenotype::readPlinkCovFile(const std::string& filename, longint N_expected) {
    readPlinkLabelsFile(filename, N_expected,
                        NULL, true, PLINK_MINNCOL_SHORT+1, 0, 0);
}
void Phenotype::readPlinkCovFile(const std::string& filename, const Phenotype& otherOrdering) {
    readPlinkLabelsFile(filename, otherOrdering.getNumObservations(),
                        &otherOrdering, true, PLINK_MINNCOL_SHORT+1, 0, 0);
}
void Phenotype::readPlinkLabelsFile(const std::string &filename,
                                 longint N_expected,
                                 const Phenotype *const otherOrderingPtr,
                                 const bool short_format,
                                 const unsigned short value_ncol,
                                 const unsigned short offset_value,
                                 const unsigned short max_value) {
    longint N; unsigned short d;
    std::map<std::string, longint> id2n;

    checkPlinkLabelsFile(filename, otherOrderingPtr, short_format,
                         value_ncol, offset_value, max_value,
                         /*out*/ N, d, id2n);
    checkNumObservations(filename, N, N_expected);

    initialiseMatrix(N); setNumClasses(d); setPlinkFIDAndIIDToLineMap(id2n);
    parsePlinkLabelsFile(filename, N, short_format,
                         value_ncol, offset_value, max_value,
                         /*out*/ getArrayPtr(), getNumObservationsInClassesRef());
}

void Phenotype::checkPlinkLabelsFile(const std::string &filename,
                                     const Phenotype *const otherOrderingPtr,
                                     const bool short_format,
                                     const unsigned short value_ncol,
                                     const unsigned short offset_value,
                                     const unsigned short max_value,
                                     /*out*/ longint &N, unsigned short &d,
                                     std::map<std::string, longint> &id2n) {
    ifstream file;
    tryOpenFile(filename, file);
    std::string line;
    std::string fid, iid, pat, mat, sex, phenotype;
    long j; // target line number (used if otherOrdering is not NULL)
    unsigned short value;
    N = 0; d = 0; id2n.clear();
    longint lineno = 0; // actual line number, incl. empty lines, for error reporting
    std::stringstream errmsgstream;
    do
    {
        std::getline(file, line); lineno++;
        if (line.size() > 0) { // skip empty lines
            // split lines to check for format errors and to read label value
            splitLabelsLine(line, lineno, short_format,
                            value_ncol, offset_value, max_value,
                            fid, iid, pat, mat, sex, phenotype, value);
            if (otherOrderingPtr != NULL) {
                j = otherOrderingPtr->getPlinkLineForFIDAndIID(fid, iid);
                // fprintf(stderr, "[TEMP] Phenotype::checkPlinkFamFile(), key=%s %s => row=%ld\n", fid.c_str(), iid.c_str(), j);
                if (j < 0) {
                    errmsgstream << "In labels file, line " << lineno << ", no matching FID and IID line in other labels file.";
                    throw Exception(errmsgstream.str());
                }
            } else {
                j = N;
            }
            N++;
            // check for duplicate identifiers
            const std::string id_str = getFIDAndIIDKeyStr(fid, iid);
            // fprintf(stderr, "[TEMP] Phenotype::checkPlinkFamFile(), key=%s\n", id_str.c_str());
            auto id2nIter = id2n.find(id_str);
            if (id2nIter != id2n.end()) {
                errmsgstream << "In labels file, line " << lineno << ", duplicate FID and IID: \"" << fid << "\", \"" << iid << "\" .";
                throw Exception(errmsgstream.str());
            }
            // update actual, shifted max value and mapping of id to data row nr
            if (value > d) d = value;
            id2n[id_str] = j;
        }
    }
    while (!file.fail());
    file.close();
    d = d+1; // seen labels in range: 0..d => d+1 labels
}

void Phenotype::parsePlinkLabelsFile(const std::string &filename, longint N,
                                     const bool short_format,
                                     const unsigned short value_ncol,
                                     const unsigned short offset_value,
                                     const unsigned short max_value,
                                     /*out*/ unsigned char *labels_buf,
                                     std::vector<longint> &nv) {
    std::ifstream file;
    tryOpenFile(filename, file);
    std::string line;
    std::string fid, iid, pat, mat, sex, phenotype;
    unsigned short value;
    std::stringstream errmsgstream;
    long lineno=0; long N_read=0;
    do
    {
        std::getline(file, line); lineno++;
        if (line.size() > 0) { // skip empty lines
            // neither splitFamLine nor getPlinkLineForFIDAndIID won't throw
            // errors if checkPlinkFamFile and setPlinkFIDAndIIDToLineMap were
            // run first
            // split lines to check for format errors and to read label value
            splitLabelsLine(line, lineno, short_format,
                            value_ncol, offset_value, max_value,
                            fid, iid, pat, mat, sex, phenotype, value);
            nv[value] = nv[value]+1;
            labels_buf[getPlinkLineForFIDAndIID(fid,iid)] = value;
            N_read++;
        }
    }
    while (!file.fail());
    if (N_read < N) {
        errmsgstream << "Error while parsing labels file: only " << (lineno+1) << " out of " << N << " labels read";
        throw std::runtime_error(errmsgstream.str());
    }
    file.close();
}

void Phenotype::splitLabelsLine(const std::string &line, const long lineno,
                                const bool short_format,
                                const unsigned short value_ncol,
                                const unsigned short offset_value,
                                const unsigned short max_value0,
                                /*out*/ std::string &fid, std::string &iid,
                                std::string &pat, std::string &mat,
                                std::string &sex, std::string &phenotype,
                                unsigned short &value) {
    std::stringstream errmsgstream;
    std::string::size_type pos, lastPos = 0, line_length = line.length();
    bool can_read_more = line_length > 0;
    if (!can_read_more) {
        errmsgstream << "In labels file, empty line " << lineno;
        throw Exception(errmsgstream.str());
    }
    unsigned short min_ncol = short_format ? PLINK_MINNCOL_SHORT : PLINK_MINNCOL;
    unsigned short last_ncol = max(value_ncol, min_ncol);
    unsigned short curr_ncol = 0; // current column number
    unsigned short max_value =
        max_value0
            ? max_value0
            : (std::numeric_limits<unsigned short>::max() + offset_value);
    std::string curr_str;
    // set default values for optional columns
    pat = ""; mat = ""; sex = ""; phenotype = "";
    const std::string delimiters(" \t\r");
    while (can_read_more && curr_ncol < last_ncol) {
        curr_ncol++;
        pos = line.find_first_of(delimiters, lastPos);
        can_read_more = (pos != std::string::npos);

        if (!can_read_more) {
            // error if not enough columns
            if (curr_ncol < last_ncol) {
                errmsgstream << "In labels file, too few columns in line " << lineno;
                throw Exception(errmsgstream.str());
            }
            pos = line_length;
        }
        curr_str = string(line.data() + lastPos, (int) pos - lastPos);

        // set meta values
        switch (curr_ncol) {
            case 1:
                fid = curr_str;
                break;
            case 2:
                iid = curr_str;
                break;
        }
        if (!short_format) {
            switch (curr_ncol) {
                case 3:
                    pat = curr_str;
                    break;
                case 4:
                    mat = curr_str;
                    break;
                case 5:
                    sex = curr_str;
                    break;
                case 6:
                    phenotype = curr_str;
                    break;
            }
        }

        // set target integer value
        if (curr_ncol == value_ncol) {
            value = parseTargetValue(curr_str, lineno, lastPos,
                                     offset_value, max_value);
        }

        lastPos = pos + 1;
    }
}

unsigned short Phenotype::parseTargetValue(const std::string &value_str,
                                           const long lineno, const long colno,
                                           const unsigned short offset_value,
                                           const unsigned short max_value) {
    int value;
    std::stringstream errmsgstream;
    try {
        value = std::stoi(value_str);
    } catch (std::invalid_argument &e) {
        errmsgstream << "In labels file, line " << lineno;
        if (colno > 0) errmsgstream << ", col "  << (colno+1);
        errmsgstream  << ", non-integer target value '" << value_str << "'";
        throw Exception(errmsgstream.str());
    }
    if ((value > max_value || value < offset_value)) {
        errmsgstream << "In labels file, line " << lineno;
        if (colno > 0) errmsgstream << ", col "  << (colno+1);
        errmsgstream << ", target value " << value << " not in {"
                     << offset_value << ", ...," << max_value << "} range.'";
        throw Exception(errmsgstream.str());
    }
    // k,..,d to indexing values 0,..,d-k
    return (unsigned short)(value - offset_value);
}



void Phenotype::writeFileStream(std::ofstream& file) const
{
    unsigned char *ptr = getArrayPtr();
    for (longint i=0; i<getNumObservations(); ++i)
    {
        file << (unsigned short) ptr[i] << endl;
    }
}

void Phenotype::writeETHFile(const std::string& filename) const
{
    writeFile(filename);
}



void Phenotype::checkEthLabelsFile(const std::string &filename,
                                   /*out*/ longint &N, unsigned short &d) {

    ifstream file; //Stream with file containing class labels
    unsigned short curr_value;
    unsigned short max_value = MAX_NUMCLASS;

    tryOpenFile(filename, file);

    // Read the entire file, one int value per line
    std::string line; long lineno = 0;
    N = 0; d = 0; // observations and labels counters
    do
    {
        std::getline(file, line); lineno++;
        if (line.size() > 0) { // skip empty lines
            curr_value = parseTargetValue(line, lineno, 0, 0, max_value);
            N++;
            if (curr_value > d) d = curr_value;
        }
    }
    while (!file.fail());

    //Close the file
    file.close();
    d = d+1; // seen labels in range: 0..d => d+1 labels
}

void Phenotype::parseEthLabelsFile(const std::string &filename, longint N,
                                   /*out*/ unsigned char *labels_buf,
                                   std::vector<longint> &nv) {

    ifstream file; //Stream with file containing class labels
    unsigned short curr_value;
    if (!isSetNumClasses())
        throw Exception("Number of labels (classes) is not set");
    unsigned short max_value = getNumClasses();

    tryOpenFile(filename, file);

    // Read the entire file, one int value per line
    std::string line; long lineno = 0;
    std::stringstream errmsgstream;
    longint N_read = 0; // number of labels actually read
    do
    {
        std::getline(file, line); lineno++;
        if (line.size() > 0) { // skip empty lines
            curr_value = parseTargetValue(line, lineno, 0, 0, max_value);
            N_read++;
            if (N_read > N) { // shouldn't happen
                errmsgstream << "Error while parsing labels file: trying to read more than " << N << " labels";
                throw std::runtime_error(errmsgstream.str());
            }
            nv[curr_value]++;
            *labels_buf++ = curr_value;
        }
    }
    while (!file.fail());

    // Sanity check to see if we successfully read the correct number of labels
    if(N_read < N){ // shouldn't happen
        errmsgstream << "Error while parsing labels file: only " << N_read << " out of " << N << " labels read";
        throw std::runtime_error(errmsgstream.str());
    }

    //Close the file
    file.close();
}

} /* namespace SignificantPattern */
