#ifndef UTIL_H
#define UTIL_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "SmartSpace.h"
#include "ProductElement.h"
#include "def.h"
#include "ManifoldOptimException.h"

using namespace arma; 
using namespace Rcpp;
using namespace ROPTLIB;

void CopyFrom(SmartSpace* s, const arma::mat& x);
arma::vec ToArmaVec(const SmartSpace* s);
arma::mat ToArmaMat(const SmartSpace* s);
arma::mat ToArmaMat(const ProductElement* s);

void CopyFrom(SmartSpace* s, const NumericMatrix& x);
void CopyFrom(NumericMatrix& x, const SmartSpace* s);

void CopyFrom(SmartSpace* s, const NumericVector& x);
NumericVector ToNumericVector(const SmartSpace* s);

List ExtractElements(const SmartSpace* s);
List ExtractElements(const ProductElement* s);

#endif
