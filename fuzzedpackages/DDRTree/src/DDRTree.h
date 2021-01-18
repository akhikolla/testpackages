#ifndef _DDRTree_DDRTREE_H
#define _DDRTree_DDRTREE_H

#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

void pca_projection_cpp(const MatrixXd& R_C, int dimensions,  MatrixXd& W);
void sq_dist_cpp(const MatrixXd& a, const MatrixXd& b,  MatrixXd& W);

#endif
