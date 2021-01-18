#ifndef __sqp_misc_eigen_included__        // if include guard for '/misc/sqp_eigen.h' is undefined
#define __sqp_misc_eigen_included__        // define include guard for '/misc/sqp_eigen.h'

#include "sqp.h"
#if SQP_USE_EIGEN == 1
#include<RcppEigen.h>
#include<Eigen/SparseQR>
#include<Eigen/SparseLU> 
#include<Eigen/SparseCholesky>
#endif // SQP_USE_EIGEN == 1

namespace sqp {
namespace misc{




#if SQP_USE_EIGEN == 1

inline Eigen::MatrixXd convert_eigen(arma::mat arma_A) {
  
  Eigen::MatrixXd eigen_B = Eigen::Map<Eigen::MatrixXd>(arma_A.memptr(),
                                                        arma_A.n_rows,
                                                        arma_A.n_cols);
  
  return eigen_B;
}

inline arma::mat convert_eigen(Eigen::MatrixXd eigen_A) {
  
  arma::mat arma_B = arma::mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               false, false);
  
  return arma_B;
}

#endif // SQP_USE_EIGEN == 1 



}
}

#endif                                    // end of include guard for '/misc/sqp_eigen.h'


