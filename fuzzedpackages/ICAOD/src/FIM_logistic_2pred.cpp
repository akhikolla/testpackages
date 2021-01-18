#include <Rcpp.h>

// roxygen
//' @title Fisher Information Matrix for the Logistic Model with Two Predictors
//' @description It provides the cpp function for FIM for the following model:\cr
//'   \code{~exp(b0+ b1 * x1 + b2 * x2 + b3 * x1 * x2)/(1 + exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))}.
//' @param x1 Vector of design points (for first predictor).
//' @param x2 Vector of design points (for second predictor).
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(b0, b1, b2, b3)}.
//' @return Fisher information matrix.
//' @export
// [[Rcpp::export]]

Rcpp::NumericMatrix FIM_logistic_2pred(const std::vector<double> x1, const std::vector<double> x2, const std::vector<double> w, const std::vector<double> param)
{
  if(x1.size() != x2.size()){
    Rcpp::stop("'x1' and 'x2' are not of the same length.");
  }
  if(x1.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }

  double b0, b1, b2, b3, constant;

  b0 = param[0];
  b1 = param[1];
  b2 = param[2];
  b3 = param[3];
  size_t i;

  Rcpp::NumericMatrix FIM_mat(4, 4);
  for(i=0; i < x1.size(); i++)
  {
    constant = 1/(1+exp(b0 + b1 * x1[i] + b2 * x2[i] + b3 * x1[i] * x2[i]));
    constant = w[i]*constant*(1-constant);
    FIM_mat(0, 0) = constant + FIM_mat(0, 0);
    FIM_mat(1, 0) = x1[i] * constant + FIM_mat(1, 0);
    FIM_mat(2, 0) = x2[i] * constant + FIM_mat(2, 0);
    FIM_mat(3, 0) = x1[i] * x2[i] * constant + FIM_mat(3, 0);

    FIM_mat(1, 1) = pow(x1[i], 2) * constant + FIM_mat(1, 1);
    FIM_mat(2, 1) = x1[i] * x2[i] * constant + FIM_mat(2, 1);
    FIM_mat(3, 1) = pow(x1[i], 2) * x2[i] * constant + FIM_mat(3, 1);

    FIM_mat(2, 2) = pow(x2[i], 2) * constant + FIM_mat(2, 2);
    FIM_mat(3, 2) = pow(x2[i], 2) * x1[i] * constant + FIM_mat(3, 2);

    FIM_mat(3, 3) = pow(x1[i] * x2[i], 2) * constant + FIM_mat(3, 3);
  }

  FIM_mat(0, 1) = FIM_mat(1, 0);
  FIM_mat(0, 2) = FIM_mat(2, 0);
  FIM_mat(0, 3) = FIM_mat(3, 0);

  FIM_mat(1, 2) = FIM_mat(2, 1);
  FIM_mat(1, 3) = FIM_mat(3, 1);

  FIM_mat(2, 3) = FIM_mat(3, 2);

  return FIM_mat;
}


