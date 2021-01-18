#ifndef SC_H
#define SC_H

#include <Rcpp.h>

// Called from R
RcppExport SEXP conditionalFNML(SEXP,SEXP);

// Called from R
RcppExport SEXP conditionalQNML(SEXP,SEXP);

// Called from R
RcppExport SEXP indepfNML(SEXP,SEXP,SEXP,SEXP);

// Called from R
RcppExport SEXP indepqNML(SEXP,SEXP,SEXP,SEXP);

// Called from R
RcppExport SEXP indepAsymfNML(SEXP,SEXP,SEXP);

// Called from R
RcppExport SEXP indepAsymqNML(SEXP,SEXP,SEXP);

// Called from R
RcppExport SEXP regret(SEXP,SEXP);

double conditionalNML(SEXP&,SEXP&,bool);

double indepNML(SEXP&,SEXP&,SEXP&,SEXP&,bool);

double indepAsymNML(SEXP&,SEXP&,SEXP&,bool);

std::vector<int> matrixToVector(Rcpp::IntegerMatrix&);

std::vector<int> joinVectors(std::vector<int>,std::vector<int>&);

#endif // SC_H
