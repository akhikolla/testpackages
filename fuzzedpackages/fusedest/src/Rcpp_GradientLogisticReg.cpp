#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP GradLogisticReg(Eigen::MatrixXd X,
                     Eigen::VectorXd y,
                     SEXP a,
                     Eigen::VectorXd b,
                     Eigen::VectorXd beta_ini,
                     SEXP stepsize,
                     SEXP max_iter,
                     SEXP tol_err){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

            double a_s = Rcpp::as<double>(a);
     double stepsize_s = Rcpp::as<double>(stepsize); /* stepsize */
        int max_iter_s = Rcpp::as<int>(max_iter); /* Maximum number of iterations */
      double tol_err_s = Rcpp::as<double>(tol_err);  /* Tolerance error */

        const int p(X.cols());
      const int n(X.rows());


      VectorXd beta_old = beta_ini;
      VectorXd beta = VectorXd(p).setZero(); /* Initial values for betahat */
      VectorXd exp_eta = VectorXd(n).setZero();
      VectorXd nabla_h_Xbeta = VectorXd(n).setZero();
      VectorXd err_list=VectorXd(max_iter_s).setZero();

        double err = 1e6;

        int l = 0;

        while(err > tol_err_s && l < max_iter_s){

                exp_eta = VectorXd(VectorXd(X*beta_old).array().exp());
          nabla_h_Xbeta = -y + exp_eta.cwiseQuotient(VectorXd(1+ exp_eta.array()));
                   beta = ((1/stepsize_s)*beta_old + b - X.adjoint()*nabla_h_Xbeta)/(1/stepsize_s + a_s);

          err = VectorXd(beta-beta_old).norm()/p;
          err_list(l) = err;
          beta_old = beta;

          l+= 1;
        }


        return  Rcpp::List::create(Rcpp::Named("beta") = Rcpp::wrap(beta),
                                   Rcpp::Named("iter_count") = Rcpp::wrap(l),
                                   Rcpp::Named("err_list") = Rcpp::wrap(err_list.segment(0,l)));
}
