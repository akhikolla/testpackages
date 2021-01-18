#include <RcppEigen.h>


//equation1 second one
// The model has an analytical solution for the locally D-optimal design. See Dette et al. (2010) for more details.\cr
// The Fisher information matrix does not depend on \eqn{\theta_0}{\theta0}.

// roxygen
//' @title Fisher Information Matrix for the Mixed Inhibition Model
//'
//' @description  It provides the cpp function for the FIM for the model \code{~theta0 + theta1* log(x + theta2)}.
//'
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(theta0, theta1, theta2)}.
//' @return Fisher information matrix.
//' @references Dette, H., Kiss, C., Bevanda, M., & Bretz, F. (2010). Optimal designs for the EMAX, log-linear and exponential models. Biometrika, 97(2), 513-518.
//' @details
//' The FIM of this model does not depend on the parameter \code{theta0}.
//' @export
// [[Rcpp::export]]



Eigen::MatrixXd FIM_loglin(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  theta0, theta1,  theta2, constant;
  theta0 = param[0];
  theta1 = param[1];
  theta2 = param[2];

  theta0 = theta0 + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(3, 3);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < x.size(); i++)
  {

    constant = x[i] + theta2;


    Eigen::MatrixXd f(3, 1);
    f(0, 0) = 1;
    f(1, 0) = log(constant);
    f(2, 0) = theta1/constant;



    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;

  }

  return Fisher_mat;
}
