#include <RcppEigen.h>



// roxygen
//' @title Fisher Information Matrix for the 4-Parameter Logistic Model
//'
//' @description It provides the cpp function for the FIM for the model
//'  \code{~theta1/(1+exp(theta2*x+theta3))+theta4}.
//'  This model is another re-parameterization of the 4-parameter Hill model.
//'   For more details, see Eq. (1) and (2) in Hyun and  Wong (2015).
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(theta1, theta2, theta3, theta4)}.
//' @return Fisher information matrix.
//' @details The fisher information matrix does not depend on \code{theta4}.\cr
//' @references
//' Hyun, S. W., & Wong, W. K. (2015). Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels. The international journal of biostatistics, 11(2), 253-271.
//' @seealso \code{\link{multiple}}
//' @export
//' @examples
//' FIM_logistic_4par(x = c(-6.9, -4.6, -3.9, 6.7 ),
//'                   w = c(0.489, 0.40, 0.061, 0.050),
//'                   param = c(1.563, 1.790, 8.442, 0.137))
// [[Rcpp::export]]



Eigen::MatrixXd FIM_logistic_4par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  theta1, theta2,  theta3, theta4, constant1, constant2;
  theta1 = param[0];
  theta2 = param[1];
  theta3 = param[2];
  theta4 = param[3];

  theta4 = theta4 + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(4, 4);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < x.size(); i++)
  {
    constant1 = exp(x[i]*theta2 + theta3);
    constant2 = 1/(1 + constant1);


    Eigen::MatrixXd f(4, 1);
    f(0, 0) = constant2;
    f(1, 0) = pow(constant2, 2)*-theta1*x[i]*constant1;
    f(2, 0) = pow(constant2, 2)*-theta1*constant1;
    f(3, 0) = 1,


    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;

  }

  return Fisher_mat;
}
