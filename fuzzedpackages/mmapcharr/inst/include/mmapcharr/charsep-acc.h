#ifndef CHAR_SEP_ACC_H
#define CHAR_SEP_ACC_H

#include <mio/mmap.hpp>
#include <system_error> // for std::error_code
#include <Rcpp.h>

using std::size_t;

class charSep {
public:
  charSep(const std::string path, int n, int m, int r): n(n), m(m), r(r) {
    std::error_code error;
    this->ro_mmap.map(path, error);
    if (error) Rcpp::stop("Error when mapping file:\n  %s.\n", error.message());
  }
  
  const unsigned char* matrix() const { return ro_mmap.data(); }
  
  size_t nrow()   const { return n; }
  size_t ncol()   const { return m; }
  size_t nextra() const { return r; }
  
private:
  mio::ummap_source ro_mmap;
  size_t n, m, r;
};


template <typename T, int RTYPE>
class charSepAcc {
public:
  charSepAcc(const charSep * charSepPtr, Rcpp::Vector<RTYPE> code) {
    
    _pMat = charSepPtr->matrix();
    n = charSepPtr->nrow();
    m = charSepPtr->ncol();
    l = (2 * m - 1) + charSepPtr->nextra();
    _code = code;
  };
  
  size_t nrow() const { return n; }
  size_t ncol() const { return m; }
  
  inline T operator() (size_t i, size_t j) {
    return _code[_pMat[i * l + 2 * j]];
  }
  
private:
  const unsigned char* _pMat;
  size_t n;
  size_t m;
  size_t l;
  Rcpp::Vector<RTYPE> _code;
};


template <typename T, int RTYPE>
class charSepAccTranspose : public charSepAcc<T, RTYPE> {
public:
  charSepAccTranspose(const charSep * charSepPtr, Rcpp::Vector<RTYPE> code) 
    : charSepAcc<T, RTYPE>(charSepPtr, code) {}
  
  size_t nrow() const { return charSepAcc<T, RTYPE>::ncol(); }
  size_t ncol() const { return charSepAcc<T, RTYPE>::nrow(); }
  
  inline const unsigned char operator() (size_t i, size_t j) {
    return charSepAcc<T, RTYPE>::operator()(j, i);
  }
};

#endif // CHAR_SEP_ACC_H
