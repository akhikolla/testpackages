#include <RcppEigen.h>




// roxygen
//' @title Fisher Information Matrix for the Sigmoid Emax Model
//' @description It provides the cpp function for FIM for the model \code{~b1+(b2-b1)*(x^b4)/(x^b4+b3^b4)}
//' @param x Vector of design points.
//' @param w Vector of design weight. Its length must be equal to the length of \code{x} and \code{sum(w) = 1}.
//' @param param Vector of values for the model parameters \code{c(b1, b2, b3, b4)}.
//' The mean of response variable is .
//' @return Fisher information matrix.
//' @export
// [[Rcpp::export]]




Eigen::MatrixXd FIM_sig_emax(const std::vector<double> x, const std::vector<double> w, const std::vector<double> param)//, double sigma)
{
  if(x.size() != w.size()){
    Rcpp::stop("'x' and 'w' are not of the same length.");
  }
  double  theta1, theta2, theta3,  theta4, sigma2, c, sigma;

  theta1 = param[0];
  theta2 = param[1];
  theta3 = param[2];
  theta4 = param[3];
  sigma = 1;
  sigma2 = pow(sigma, 2);


  Eigen::MatrixXd I_mat(4, 4);
  I_mat.setZero();
  Eigen::MatrixXd Fisher_mat(4, 4);
  Fisher_mat.setZero();

  size_t i;

  for(i=0; i < x.size(); i++)
  {
    c = pow((x[i]/theta3), theta4);

    I_mat(0, 0) = pow((1+c), -2);
    I_mat(0, 1) = I_mat(0, 0) * c;
    I_mat(1, 0) = I_mat(0, 1);


    I_mat(0, 2) = I_mat(0, 0) * ((theta1-theta2)*(theta4/theta3) * c)/(1+c);
    I_mat(2, 0) = I_mat(0, 2);



    I_mat(0, 3) = I_mat(0, 0) * ((theta2-theta1) * c * log(x[i]/theta3))/(1+c);
    I_mat(3, 0) = I_mat(0, 3);


    I_mat(1, 2) = I_mat(1, 0) *  I_mat(0, 2)/I_mat(0, 0);
    I_mat(2, 1) = I_mat(1, 2);


    I_mat(1, 3) = I_mat(3, 0) *  I_mat(0, 1)/I_mat(0, 0);
    I_mat(3, 1) = I_mat(1, 3);

    I_mat(2, 3) = I_mat(2, 0) *  I_mat(0, 3)/I_mat(0, 0);
    I_mat(3, 2) = I_mat(2, 3);

    I_mat(1, 1) = pow(I_mat(0, 1), 2)/I_mat(0, 0);
    I_mat(2, 2) = pow(I_mat(0, 2), 2)/I_mat(0, 0);
    I_mat(3, 3) =  pow(I_mat(0, 3), 2)/I_mat(0, 0);


    Fisher_mat = w[i] * I_mat + Fisher_mat;

  }
  Fisher_mat = 1/sigma2 * Fisher_mat;
  return Fisher_mat;
}
