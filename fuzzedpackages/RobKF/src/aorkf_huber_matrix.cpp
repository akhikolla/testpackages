#include <RcppEigen.h>
#include <list>
#include <iostream>

#include "aorkf_huber.h"

std::list<Eigen::MatrixXd> aorkf_huber_matrix(const Eigen::MatrixXd& mu_old,
					const Eigen::MatrixXd& Sigma_old,
					const Eigen::MatrixXd& y,
					const Eigen::MatrixXd& A,
					const Eigen::MatrixXd& b,
					const Eigen::MatrixXd& C,
					const Eigen::MatrixXd& d,
					const Eigen::MatrixXd& R,
					const Eigen::MatrixXd& Q,
					const double& h)	
{

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(),A.rows());

  Eigen::MatrixXd m = A*mu_old + b;
  Eigen::MatrixXd P = A*Sigma_old*A.transpose() + Q;

  Eigen::MatrixXd Y_GAP  = y - C*m - d;
  Eigen::MatrixXd M      = P*C.transpose()*(C*P*C.transpose()   +   R).inverse();
  Eigen::MatrixXd Update = M*Y_GAP;

  double magnitude = Update.norm();

  if (magnitude > h)
  {
    Update = Update * (h/magnitude);
  }

  Eigen::MatrixXd Sigma_new = (I - M * C ) * P;
  Eigen::MatrixXd mu_new    = m + Update;
  

  std::list<Eigen::MatrixXd> result;
  result.push_back(mu_new);
  result.push_back(Sigma_new);
  return(result);

}
