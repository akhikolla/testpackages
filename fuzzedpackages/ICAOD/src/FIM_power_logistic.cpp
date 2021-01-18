#include <Rcpp.h>


// roxygen
//' @title Fisher Information Matrix for the Power Logistic Model
//' @description It provides the cpp function for FIM for the model  \code{~1/(1 + exp(-b *(x - a)))^s}, but when \code{s} is fixed (a two by two matrix).
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(a, b)}.
//' @param s parameter \code{s}.
//' @return Fisher information matrix.
//' @export
//' @note This matrix is a two by two matrix and not equal to the Fisher information matrix for the power logistic model
//' when the derivative is taken with respect to all the three parameters.
//' This matrix is only given to be used in some illustrative examples.
// [[Rcpp::export]]


Rcpp::NumericMatrix FIM_power_logistic(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param, const double s)
{
    if(x.size() != w.size()){
      Rcpp::stop("'x' and 'w' are not of the same length.");

    }
    double A = 0, B = 0, C = 0, a, b, constant, w_sum;
    a = param[0];
    b = param[1];

    size_t i;

    for(i=0; i < x.size(); i++)
    {
        constant = 1/pow((1+exp(-b*(x[i]-a))), s);
        constant = w[i]*pow(s, 2)*constant*pow(1-pow(constant, 1/s),2)/(1-constant);
        A = pow(b,2)* constant + A;
        B = -b * (x[i] - a) * constant + B;
        C = pow((x[i] - a), 2) * constant + C;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = C;

    return Fisher_mat;
}


