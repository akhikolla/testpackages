#include <Rcpp.h>

//#include <math.h>
//#include <vector>
//using namespace std;
//using namespace Rcpp;
//ON THE NUMBER OF SUPPORT POINTS OF MAXIMIN AND BAYESIAN OPTIMAL DESIGNS1
// BRAESS AND  DETTE (2007)
// The locally D optimal design is independent of the nominal
// value of \eqn{a} and is equally supported at \eqn{x = 0} and \eqn{x = 1/b}
//  only when \eqn{x \in [0, 1]}{{x belongs to [0, 1]}}. See "Examples".


//' @title Fisher Information Matrix for the 2-Parameter Exponential Model
//'
//' @description
//' It provides the cpp function for FIM for the model  \code{~a + exp(-b*x)}.
//'
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(a, b)}.
//' @return Fisher information matrix.
//' @references Dette, H., & Neugebauer, H. M. (1997). Bayesian D-optimal designs for exponential regression models. Journal of Statistical Planning and Inference, 60(2), 331-349.
//' @details The FIM does not depend on the value of \code{a}.
//' @examples FIM_exp_2par(x = c(1, 2), w = c(.5, .5), param = c(3, 4))
//' @export
// [[Rcpp::export]]




Rcpp::NumericMatrix FIM_exp_2par(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }

  double A = 0, B = 0, C = 0, alpha, beta, constant, w_sum;
  alpha = param[0];
  beta = param[1];


  alpha = alpha + 0; //just to not get a warning


  size_t i;


    for(i=0; i < x.size(); i++)
    {
        constant = exp(-beta*x[i]);

        A = w[i] * 1 + A;
        B = w[i]*(-x[i]) * constant + B;
        C = w[i]*pow(x[i], 2)*exp(-2*beta*x[i]) + C;
        w_sum = w[i] + w_sum;
    }


    Rcpp::NumericMatrix Fisher_mat(2, 2);
    Fisher_mat(0, 0) = A;
    Fisher_mat(0, 1) = B;
    Fisher_mat(1, 0) = B;
    Fisher_mat(1, 1) = C;

    return Fisher_mat;
}
