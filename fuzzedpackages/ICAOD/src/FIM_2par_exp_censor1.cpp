#include <Rcpp.h>


//' @title Fisher Information Matrix for a 2-Parameter Cox Proportional-Hazards Model for Type One Censored Data
//'
//' @description
//' It provides the cpp function for the FIM introduced in  Eq. (3.1) of Schmidt and Schwabe (2015) for type one censored data.
//'
//'
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \eqn{c(\beta_0, \beta_1)}.
//' @param tcensor The experiment is terminated at the fixed time point \code{tcensor}.
//' @return Fisher information matrix.
//' @references Schmidt, D., & Schwabe, R. (2015). On optimal designs for censored data. Metrika, 78(3), 237-257.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix FIM_2par_exp_censor1(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param, const double tcensor)
{


    if(x.size() != w.size()){
      Rcpp::Rcout<<"The length of weights and points is not equal."<<std::endl;

    }
    double A = 0, B = 0, D = 0, alpha, beta, Q, w_sum, theta;
    alpha = param[0];
    beta = param[1];

    size_t i;

    for(i=0; i < x.size(); i++)
    {
      theta = alpha+beta*x[i];
        Q = w[i]*(1-exp(-tcensor*exp(theta)));
        A = 1* Q + A;
        B = x[i] * Q + B;
        D = pow(x[i], 2) * Q + D;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = D;

    return Fisher_mat;
}


