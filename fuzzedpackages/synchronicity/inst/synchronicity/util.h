#ifndef BIGMEMORY_UTIL_HPP
#define BIGMEMORY_UTIL_HPP

#include <vector>
#include <string>

#include <Rcpp.h>

using namespace std;

vector<string> RChar2StringVec( SEXP charVec );

vector<string> RChar2StringVec( SEXP charVec, 
  const vector<unsigned long> &indices );

SEXP String2RChar(const std::string &str);

std::string RChar2String(SEXP str);

SEXP StringVec2RChar( const vector<string> &strVec );

template<typename T>
SEXP StringVec2RChar( const vector<string> &strVec,
  T indices, const unsigned long indicesLength )
{
  if (strVec.empty())
    return R_NilValue;
  SEXP ret = Rf_protect(Rf_allocVector(STRSXP, indicesLength));
  unsigned long i;
  for (i=0; i < indicesLength; ++i)
  {
    SET_STRING_ELT(ret, i, 
      Rf_mkChar(strVec[static_cast<unsigned long>(indices[i])-1].c_str()));
  }
  UNPROTECT(1);
  return ret;
}


template<typename T>
struct NewVec;

template<>
struct NewVec<int>
{SEXP operator()(long n) const {return Rf_allocVector(INTSXP, n);};};

template<>
struct NewVec<double>
{SEXP operator()(long n) const {return Rf_allocVector(REALSXP, n);};};

template<typename T>
struct VecPtr;

template<>
struct VecPtr<int>
{int* operator()(SEXP vec) const {return INTEGER(vec);};};

template<>
struct VecPtr<double>
{double* operator()(SEXP vec) const {return REAL(vec);};};

#endif // BIGMEMORY_UTIL_HPP
