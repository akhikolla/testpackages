/*
 * conversions.h
 *
 *  Created on: 03.05.2009
 *      Author: daniel
 */

#ifndef CONVERSIONS_H_
#define CONVERSIONS_H_

#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <string>

#include <types.h>


// convert an INTSXP to a std::vector of (unsigned) integers
template <class T>
std::vector<T>
getIntegerVector(SEXP R_input)
{
    // how long is the vector?
    R_len_t input_size = Rf_length(R_input);

    // get the int array data of the R vector
    int* input = INTEGER(R_input);

    // and copy that into the return vector using the
    // iterator constructor of vector
    return std::vector<T>(input, input + input_size);
}

// convert a std::vector of (unsigned) integers to an INTSXP
template <class T>
SEXP
putIntegerVector(std::vector<T> output)
{
    // allocate return INTSXP
    SEXP R_ret;
    Rf_protect(R_ret = Rf_allocVector(INTSXP, output.size()));

    // get the int array data of the R vector
    int* ret = INTEGER(R_ret);

    // copy the std::vector into the int array
    std::copy(output.begin(), output.end(), ret);

    // unprotect and return the R vector
    Rf_unprotect(1);
    return R_ret;
}

// convert a REALSXP to a std::vector of doubles
DoubleVector
getDoubleVector(SEXP R_input);

// convert a std::vector of doubles to a REALSXP
SEXP
putDoubleVector(DoubleVector output);

// convert a std::vector of strings to a STRSXP
SEXP
putStringVector(const std::vector<std::string>& stringVec);

// convert a STRSXP to a std::vector of strings
std::vector<std::string>
getStringVector(SEXP R_input);

// get a list element by name
SEXP
getListElement(SEXP R_list, const std::string &name);

// get real vector element by name
double
getDoubleElement(SEXP R_realVector, const std::string &name);

// R-Matrix to Newmat-Matrix
ReturnMatrix
getMatrix(const SEXP&);

// R-Vector to Newmat-column vector
ColumnVector
vec2col(const SEXP&);

// Newmat-Matrix to R-Matrix
SEXP
putMatrix(const Matrix&);

// convert frequency vector into multiset
Powers
freqvec2multiset(const IntVector&);


#endif /* CONVERSIONS_H_ */
