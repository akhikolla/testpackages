#include <RcppEigen.h>
#include <list>
#include <iostream>
#include "aorkf_t.h"

std::list<Eigen::MatrixXd>  aorkf_t_matrix(const Eigen::MatrixXd& mu_old,
					const Eigen::MatrixXd& Sigma_old,
					const Eigen::MatrixXd& y,
					const Eigen::MatrixXd& A,
					const Eigen::MatrixXd& b,
					const Eigen::MatrixXd& C,
					const Eigen::MatrixXd& d,
					const Eigen::MatrixXd& R,
					const Eigen::MatrixXd& Q,
					const double& s,
					const double& epsilon)	
{

  Eigen::MatrixXd m = A*mu_old + b;
  Eigen::MatrixXd P = A*Sigma_old*A.transpose() + Q;
  Eigen::MatrixXd Gamma = R;
  Eigen::MatrixXd Gamma_new = Gamma;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(A.rows(),A.rows());  
  double squarediff = 0;
  Eigen::MatrixXd Sigma_new;
  Eigen::MatrixXd mu_new;
  do
  {
      Gamma = Gamma_new;
      Eigen::MatrixXd K = (C*P*C.transpose()+Gamma).inverse()*C*P;
      mu_new = m + K.transpose() * (y - C * m - d);
      Sigma_new = K.transpose()*Gamma*K + (I - K.transpose()*C) * P * (I - C.transpose() * K);
      Eigen::MatrixXd delta = y - C*mu_new - d;
      Gamma_new = (s*R + delta*delta.transpose() + C*Sigma_new*C.transpose())/(s+1);
      squarediff = (Gamma_new-Gamma).squaredNorm();
  } while(squarediff > epsilon );

  std::list<Eigen::MatrixXd> result;
  result.push_back(mu_new);
  result.push_back(Sigma_new);
  return(result);
}
