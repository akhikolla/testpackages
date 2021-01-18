#include <RcppEigen.h>
#include <list>
#include <iostream>
#include "kf.h"

#include "exception.h"
#include "user_interupt.h"
#include "check_user_interrupt.h"

// [[Rcpp::export]]
std::list<std::list<Eigen::MatrixXd> > kf_list(const Eigen::MatrixXd& mu_init,
		     const Eigen::MatrixXd& Sigma_init,
		     const std::list<Eigen::MatrixXd>& ys,
		     const Eigen::MatrixXd& A,
		     const Eigen::MatrixXd& b,
		     const Eigen::MatrixXd& C,
		     const Eigen::MatrixXd& d,
		     const Eigen::MatrixXd& R,
		     const Eigen::MatrixXd& Q)	
{

  Eigen::MatrixXd mu = mu_init;
  Eigen::MatrixXd Sigma = Sigma_init;

  std::list<std::list<Eigen::MatrixXd> > result;
  std::list<Eigen::MatrixXd> pair;
  pair.push_back(mu);
  pair.push_back(Sigma);
  result.push_back(pair);
    
  std::list<Eigen::MatrixXd>::const_iterator it = ys.begin();
  while(it != ys.end())
    {

      if(check_user_interrupt())
      {
	  throw_exception("User interrupt");
      }

      pair = kf_matrix(mu,Sigma,*it,A,b,C,d,R,Q);
      std::list<Eigen::MatrixXd>::iterator pair_it = pair.begin();
      mu = *pair_it;
      pair_it++;
      Sigma = *(pair_it);
      result.push_back(pair);
      it++;
    }

  return result;

}
