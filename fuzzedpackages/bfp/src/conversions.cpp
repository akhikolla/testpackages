/*
 * conversions.cpp
 *
 *  Created on: 03.05.2009
 *      Author: daniel
 */

#include <conversions.h>
#include <types.h>

#include <vector>
#include <algorithm>



// convert a REALSXP to a std::vector of doubles
DoubleVector
getDoubleVector(SEXP R_input)
{
    // how long is the vector?
    R_len_t input_size = Rf_length(R_input);

    // get the double array data of the R vector
    double* input = REAL(R_input);

    // and copy that into the return vector using the
    // iterator constructor of vector
    return DoubleVector(input, input + input_size);
}

// convert a std::vector of doubles to a REALSXP
SEXP
putDoubleVector(DoubleVector output)
{
    // count protected SEXP's
    unsigned int nProtect(0);

    // allocate return INTSXP
    SEXP R_ret;
    Rf_protect(R_ret = Rf_allocVector(REALSXP, output.size()));
    nProtect++ ;

    // get the double array data of the R vector
    double* ret = REAL(R_ret);

    // copy the std::vector into the double array
    std::copy(output.begin(), output.end(), ret);

    // unprotect and return the R vector
    Rf_unprotect(nProtect);
    return R_ret;
}

// convert a std::vector of strings to a STRSXP
SEXP
putStringVector(const StringVector& stringVec)
{
        // count protected SEXP's
        unsigned int nProtect(0);

        // allocate return STRSXP
        SEXP R_ret;
        Rf_protect(R_ret = Rf_allocVector(STRSXP, stringVec.size()));
        nProtect++;

        for(std::vector<std::string>::size_type
                        j = 0;
                        j < stringVec.size();
                        j++)
        {
                // insert converted string into return object
                SET_STRING_ELT(R_ret, j, Rf_mkChar(stringVec[j].c_str()));
        }

        Rf_unprotect(nProtect);
        return R_ret;
}

// convert a STRSXP to a std::vector of strings
StringVector
getStringVector(SEXP R_input)
{
        // allocate return vector
        std::vector<std::string> ret;

        // and convert all elements of R_input
        int input_size = Rf_length(R_input);
        for(int j = 0; j < input_size; j++)
        {
                // insert converted string into return object
                ret.push_back(CHAR(STRING_ELT(R_input, j)));
        }

        return ret;
}

// get a list element by name
SEXP
getListElement(SEXP R_list, const std::string &name)
{
        // if it is not found, we return NULL
    SEXP R_elmt = R_NilValue;

    // get the names of the list
    std::vector<std::string> list_names = getStringVector(Rf_getAttrib(R_list, R_NamesSymbol));

    // and search for the given name
    for (std::vector<std::string>::size_type
                i = 0;
                i < list_names.size();
                i++)
    {
        if (list_names[i] == name)
        {
                R_elmt = VECTOR_ELT(R_list, i);
                break;
        }
    }

    return R_elmt;
}

// get real vector element by name
double
getDoubleElement(SEXP R_realVector, const std::string &name)
{
        // if it is not found, we return NA
    double elmt = R_NaReal;

    // extract the names
    std::vector<std::string> vector_names = getStringVector(Rf_getAttrib(R_realVector, R_NamesSymbol));

    // and the array
    double* data = REAL(R_realVector);

    // find the (first) correct name's position
    for (std::vector<std::string>::size_type
                i = 0;
                i < vector_names.size();
                i++)
    {
        if (vector_names[i] == name)
        {
                elmt = data[i];
                break;
        }
    }

    return elmt;
}



ReturnMatrix getMatrix(const SEXP& m) // R-Matrix to Newmat-Matrix
{
        unsigned int mx, my;
        int *mDim;
        double *a;

        a = REAL(m); // read matrix m in array (columnwise)

        mDim = INTEGER(Rf_getAttrib(m, R_DimSymbol)); // get dimensions of m
        mx = mDim[0];
        my = mDim[1];

        Matrix M(my, mx); // copy array into matrix. This is done rowwise, so transpose necessary
        M << a;
        M = M.t();

        M.Release(); return M;
}

ColumnVector vec2col(const SEXP& v) // R-Vector to Newmat-column vector
{
        if (Rf_isMatrix(v))
                Rf_error("Argument of vec2col is a matrix\n");

        ColumnVector V(Rf_length(v));
        double *a = REAL(v);
        V << a;

        return V;
}

SEXP putMatrix(const Matrix& M) // Newmat-Matrix to R-Matrix
{
        unsigned int nProtected = 0, Mx, My, i, j;
        double *a;
        SEXP ret;

        Mx = M.Nrows(); // get dimensions of M
        My = M.Ncols();

        Rf_protect(ret = Rf_allocMatrix(REALSXP, Mx, My)); // allocate return matrix
        ++nProtected;
        a = REAL(ret);

        for (i = 0; i != Mx; i++){      // put values in return matrix
                for (j = 0; j != My; j++){
                        a[i + j * Mx] = M.element(i,j);
                }
        }

        Rf_unprotect(nProtected); // unprotect everything
        return ret;
}


Powers freqvec2multiset(const IntVector& vec) // convert frequency vector into multiset
{
        Powers ret;

        int power = 0;
        for (IntVector::const_iterator
                i = vec.begin();
                i != vec.end();
                ++i, ++power)
        {
            for (int times = 0; times != *i; times++)
            {
                ret.insert(power);
            }
        }

        return ret;
}
