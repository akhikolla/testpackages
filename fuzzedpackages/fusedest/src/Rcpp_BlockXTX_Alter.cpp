#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

Eigen::SparseMatrix<double> ComputeBlockXTX(Eigen::MatrixXd X_s,
                                               SEXP a,
                                               SEXP ind_strt,
                                               SEXP n_dc){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;
  using Eigen::LLT;

  typedef Eigen::Triplet<double> T;

  const Map<VectorXd> a_s(as<Map<VectorXd> >(a)); /* Vector of group sizes */
  const Map<VectorXd> ind_strt_s(as<Map<VectorXd> >(ind_strt));
  const Map<VectorXd> n_dc_s(as<Map<VectorXd> >(n_dc));

  const int p(X_s.cols());
  const int m(n_dc_s.size());
  int sq_p_m = p*p*m;

  MatrixXd XTX_j = MatrixXd(p, p).setZero();
   MatrixXd Id = MatrixXd(p, p).setIdentity();

  std::vector<T> tripletList;
  tripletList.reserve(sq_p_m);

  int n_dc_j = 0;
  int u = 0;

  for(int j = 0; j < m; j++){

    u = ind_strt_s(j)-1;
    n_dc_j = n_dc_s(j);

    XTX_j = a_s(j)*Id;
    XTX_j += MatrixXd(X_s.block(u, 0, n_dc_j, p).adjoint()*X_s.block(u, 0, n_dc_j, p));

    for(int k = 0; k < p; k++){
      for(int l = 0; l < p; l++){
                   tripletList.push_back(T(p*j + k, p*j + l, XTX_j(k, l)));
                                        }
                                      }
    XTX_j.setZero();

                             }


  Eigen::SparseMatrix<double> XTX_plus_a(m*p, m*p);

  XTX_plus_a.setFromTriplets(tripletList.begin(), tripletList.end());

  return XTX_plus_a;
                           }
