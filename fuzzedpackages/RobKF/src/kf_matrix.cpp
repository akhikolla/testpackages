#include <RcppEigen.h>
#include <list>
#include <iostream>
#include "kf.h"

std::list<Eigen::MatrixXd>  kf_matrix(const Eigen::MatrixXd& mu_old,
					const Eigen::MatrixXd& Sigma_old,
					const Eigen::MatrixXd& y,
					const Eigen::MatrixXd& A,
					const Eigen::MatrixXd& b,
					const Eigen::MatrixXd& C,
					const Eigen::MatrixXd& d,
					const Eigen::MatrixXd& R,
					const Eigen::MatrixXd& Q)	
{

  Eigen::MatrixXd m = A*mu_old + b;
  Eigen::MatrixXd P = A*Sigma_old*A.transpose() + Q;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(),A.rows());  
  Eigen::MatrixXd K = P*C.transpose()*(C*P*C.transpose()+R).inverse();
  Eigen::MatrixXd Sigma_new = (I - K*C) * P;
  Eigen::MatrixXd mu_new = m + K*(y - C * m - d);

  std::list<Eigen::MatrixXd> result;
  result.push_back(mu_new);
  result.push_back(Sigma_new);
  return(result);

}
