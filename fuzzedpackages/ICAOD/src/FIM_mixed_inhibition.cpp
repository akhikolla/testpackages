#include <RcppEigen.h>



//Matlab code to produce the derivatives
//syms   I S V K_m K_ic K_iu
//mixed_type = V*S/(K_m*(1+I/K_ic) + S*(1+I/K_iu));
//diff(mixed_type, V)
//diff(mixed_type, K_m)
//diff(mixed_type, K_ic)
//diff(mixed_type, K_iu)
// The model has an analytical solution for the locally D-optimal design. See Bogacka et al. (2011) for details.\cr

// roxygen
//' @title Fisher Information Matrix for the Mixed Inhibition Model.
//'
//' @description
//'  It provides the cpp function for FIM for the model \code{~ V*S/(Km * (1 + I/Kic)+ S * (1 + I/Kiu))}
//'
//' @param S Vector of \code{S} component of design points. \code{S} is the substrate concentration.
//' @param I  Vector of \code{I} component of design points. \code{I} is the inhibitor concentration.
//' @param w Vector of design weight. Its length must be equal to the length of \code{S} and \code{I}, besides  \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(V, Km, Kic, Kiu)}.
//' @return Fisher information matrix of design.
//' @references Bogacka, B., Patan, M., Johnson, P. J., Youdim, K., & Atkinson, A. C. (2011). Optimum design of experiments for enzyme inhibition kinetic models. Journal of biopharmaceutical statistics, 21(3), 555-572.
//' @details
//'  The optimal design does not depend on parameter \eqn{V}.
//' @examples
//' FIM_mixed_inhibition(S = c(30, 3.86, 30, 4.60),
//'                      I = c(0, 0, 5.11, 4.16), w = rep(.25, 4),
//'                      param = c(1.5, 5.2, 3.4, 5.6))
//' @export
// [[Rcpp::export]]


Eigen::MatrixXd FIM_mixed_inhibition(const std::vector<double> S, const std::vector<double> I, const std::vector<double> w, const std::vector<double> param)
{
  if(I.size()/w.size() != S.size()/w.size()){
    Rcpp::stop("The length of 'I' or 'S' is not equal to the length of 'w'.");
  }
  double  V, K_m, K_ic, K_iu, constant;
  V = param[0];
  K_m = param[1];
  K_ic = param[2];
   K_iu = param[3];
  //a = a + 0; //just to not get a warning

  Eigen::MatrixXd Fisher_mat(4, 4);

  size_t i;
  Fisher_mat.setZero();
  for(i=0; i < I.size(); i++)
  {

    Eigen::MatrixXd f(4, 1);
    constant = (K_m*(I[i]/K_ic + 1) + S[i]*(I[i]/K_iu + 1));
    f(0, 0) = S[i]/constant;
    f(1, 0) = -(S[i]*V*(I[i]/K_ic + 1))/pow(constant,2);
    f(2, 0) = (I[i]*K_m*S[i]*V)/(pow(K_ic,2)*pow(constant,2));
    f(3, 0) = (I[i]*pow(S[i],2)*V)/(pow(K_iu,2)*pow(constant,2));
    Fisher_mat = w[i] * f * f.transpose() + Fisher_mat;
  }

  return Fisher_mat;
}
