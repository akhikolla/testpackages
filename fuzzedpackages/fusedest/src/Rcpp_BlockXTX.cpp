#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


Eigen::SparseMatrix<double> RcppBlockXTX(Eigen::MatrixXd X_s, SEXP grp_size, SEXP ind_strt){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;

  typedef Eigen::Triplet<double> T;

  /*const Map<MatrixXd> X_s(as<Map<MatrixXd> >(X));*/
  const Map<VectorXd> grp_size_s(as<Map<VectorXd> >(grp_size)); /* Vector of group sizes */
  const Map<VectorXd> ind_strt_s(as<Map<VectorXd> >(ind_strt));

  const int n(X_s.rows());
  const int p(X_s.cols());
  const int m(grp_size_s.size());
  const int max_p_j(grp_size_s.maxCoeff());
  const int sum_sq_grp(grp_size_s.array().square().sum());

  MatrixXd XTX_j = MatrixXd(max_p_j, max_p_j).setZero();

  std::vector<T> tripletList;
  tripletList.reserve(sum_sq_grp);

  int grp_size_j = 0;
  int u = 0;

  for(int j = 0; j < m; j++){

    u = ind_strt_s(j)-1;
    grp_size_j = grp_size_s(j);
    XTX_j.block(0,0, grp_size_j, grp_size_j) = MatrixXd(X_s.block(0, u, n, grp_size_j).adjoint()*X_s.block(0, u, n, grp_size_j));

    for(int k = 0; k < grp_size_j; k++){
      for(int l = 0; l < grp_size_j; l++){
                   tripletList.push_back(T(u + k, u + l, XTX_j(k, l)));
                                        }
                                      }
                             }


  Eigen::SparseMatrix<double> XTX(p, p);

  XTX.setFromTriplets(tripletList.begin(), tripletList.end());

  return XTX;
                           }
