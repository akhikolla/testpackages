#include <Rcpp.h>

// roxygen
//' @title Fisher Information Matrix for the 2-Parameter Logistic (2PL) Model
//' @description It provides the cpp function for FIM for the model  \code{~1/(1 + exp(-b *(x - a)))}.
//' In item response theory (IRT),
//' \eqn{a} is the item difficulty parameter, \eqn{b} is the item discrimination parameter and \eqn{x} is the person ability parameter.
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(a, b)}.
//' @return Fisher information matrix.
//' @export
//' @details
//'  It can be shown that minimax and standardized D-optimal designs for the 2PL model is symmetric around point
//' \eqn{a_M = (a^L + a^U)/2}{aM = (aL + aU)/2} where \eqn{a^L}{aL} and \eqn{a^U}{aU} are the
//' lower bound and upper bound for parameter \eqn{a}, respectively. In \code{\link{ICA.control}},
//'  arguments \code{sym} and \code{sym_point} can be used to specify \eqn{a_M}{aM} and find accurate symmetric optimal designs.
//' @examples
//' FIM_logistic(x = c(1, 2), w = c(.5, .5), param = c(2, 1))
//' @importFrom Rcpp evalCpp
//' @useDynLib ICAOD
// [[Rcpp::export]]

Rcpp::NumericMatrix FIM_logistic(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
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
        constant = 1/(1+exp(-b*(x[i]-a)));
        constant = w[i]*constant*(1-constant);
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


