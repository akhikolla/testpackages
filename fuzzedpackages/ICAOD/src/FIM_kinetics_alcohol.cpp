#include <RcppEigen.h>




// roxygen
//' @title Fisher Information Matrix for the Alcohol-Kinetics Model
//' @description It provides the cpp function for FIM for the model \code{~(b3 * x1)/(1 + b1 * x1 + b2 * x2)}
//' @param x1 Vector of design points (first dimension).
//' @param x2 Vector of design points (second dimension).
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(b1, b2, b3)}.
//' @return Fisher information matrix.
//' @export
// [[Rcpp::export]]

Eigen::MatrixXd FIM_kinetics_alcohol(const std::vector<double> x1, const std::vector<double> x2, const std::vector<double> w, const std::vector<double> param)
{
  if(x1.size() != x2.size()){
    Rcpp::stop("'x1' and 'x2' are not of the same length.");
  }
  if(x1.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  beta1, beta2,  beta3, constant;
  beta1 = param[0];
  beta2 = param[1];
  beta3 = param[2];

  Eigen::MatrixXd Fisher_mat(3, 3);
  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < x1.size(); i++)
  {
    constant = (1 + beta1 * x1[i] + beta2 * x2[i]);
    Eigen::MatrixXd f(3, 1);
    f(0, 0) = -beta3 * pow(x1[i], 2)/pow(constant, 2);
    f(1, 0) = -beta3 * (x1[i] *  x2[i])/pow(constant, 2);
    f(2, 0) = x1[i]/constant;
    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;
  }
  return Fisher_mat;
}
