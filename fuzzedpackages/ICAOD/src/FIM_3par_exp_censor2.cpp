#include <Rcpp.h>


//' @title Fisher Information Matrix for a 3-Parameter Cox Proportional-Hazards Model for Random Censored Data
//'
//' @description
//' It provides the cpp function for the FIM introduced in Page 247 of Schmidt and Schwabe (2015) for random censored data (type two censored data).
//'
//'
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \eqn{(\beta_0, \beta_1, \beta_2)}.
//' @param tcensor The experiment is terminated at the fixed time point \code{tcensor}.
//' @return Fisher information matrix.
//' @references Schmidt, D., & Schwabe, R. (2015). On optimal designs for censored data. Metrika, 78(3), 237-257.
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix FIM_3par_exp_censor2(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param, const double tcensor)
{


    if(x.size() != w.size()){
      Rcpp::Rcout<<"The length of weights and points is not equal."<<std::endl;

    }
    double one = 0, x1 = 0, x2 = 0, x3=0, x4=0, beta0, beta1, beta2, theta, Q, w_sum, constant;
    beta0 = param[0];
    beta1 = param[1];
    beta2 = param[2];
    size_t i;

    for(i=0; i < x.size(); i++)
    {
        theta = beta0+beta1*x[i]+beta2*pow(x[i],2);
        constant =  exp(tcensor*exp(theta));
        Q = w[i]*(1+(exp(-constant)-1)/constant);



        one = 1* Q + one;
        x1 = x[i] * Q + x1;
        x2 = pow(x[i], 2) * Q + x2;
        x3 = pow(x[i], 3) * Q + x3;
        x4 = pow(x[i], 4) * Q + x4;

        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(3, 3);
    Fisher_mat(0, 0) = one;
    Fisher_mat(0, 1) = x1;
    Fisher_mat(0, 2) = x2;

    Fisher_mat(1, 0) = x1;
    Fisher_mat(1, 1) = x2;
    Fisher_mat(1, 2) = x3;

    Fisher_mat(2, 0) = x2;
    Fisher_mat(2, 1) = x3;
    Fisher_mat(2, 2) = x4;


    return Fisher_mat;
}


