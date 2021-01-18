#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

Eigen::SparseMatrix<double> ComputeG(SEXP p, SEXP m, SEXP q_H, Eigen::MatrixXd H_s){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;

  typedef Eigen::Triplet<double> T;

  /*const Map<MatrixXd> X_s(as<Map<MatrixXd> >(X));*/

  int p_s = Rcpp::as<int>(p);
  int m_s = Rcpp::as<int>(m);
  int q_H_s = Rcpp::as<int>(q_H);

  MatrixXd G_jk = MatrixXd(p_s, m_s*p_s).setZero();
  MatrixXd Id = MatrixXd(p_s, p_s).setIdentity();

  std::vector<T> tripletList;
  tripletList.reserve(2*q_H_s*p_s);

  int ind_j = 0;
  int ind_k = 0;

  for(int l = 0; l < q_H_s; l++){

    ind_j = int(H_s(l,0)) - 1;
    ind_k = int(H_s(l,1)) - 1;
    G_jk.block(0, ind_j*p_s, p_s, p_s) = Id;
    G_jk.block(0, ind_k*p_s, p_s, p_s) = -Id;

    for(int l1 = 0; l1 < p_s; l1++){
                   tripletList.push_back(T(l*p_s + l1, ind_j*p_s + l1, 1));
                   tripletList.push_back(T(l*p_s + l1, ind_k*p_s + l1, -1));
                                        }

    G_jk.setZero();
                             }


  Eigen::SparseMatrix<double> G(q_H_s*p_s, m_s*p_s);

  G.setFromTriplets(tripletList.begin(), tripletList.end());

  return G;
                           }
