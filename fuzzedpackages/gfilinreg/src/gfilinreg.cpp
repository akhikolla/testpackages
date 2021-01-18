#include <RcppEigen.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
using boost::math::normal;
using boost::math::students_t;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

Eigen::VectorXd qnorm(const Eigen::VectorXd& p) {
  boost::math::normal gaussian;
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++) {
    out(i) = quantile(gaussian, p.coeff(i));
  }
  return out;
}

double logdnorm(const Eigen::VectorXd& x) {
  return -x.cwiseProduct(x).sum() / 2.0;
}

Eigen::VectorXd qcauchy(const Eigen::VectorXd& p) {
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++) {
    out(i) = tan(M_PI * (p.coeff(i) - 0.5));
  }
  return out;
}

double logdcauchy(const Eigen::VectorXd& x) {
  return -x.cwiseProduct(x).array().log1p().sum();
}

Eigen::VectorXd qt(const Eigen::VectorXd& p, const double nu) {
  boost::math::students_t t(nu);
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++) {
    out(i) = quantile(t, p.coeff(i));
  }
  return out;
}

double logdt(const Eigen::VectorXd& x, const double nu) {
  return -(nu + 1.0) / 2.0 * log1p(x.cwiseProduct(x).array() / nu).sum();
}

Eigen::VectorXd qlogistic(const Eigen::VectorXd& p) {
  Eigen::VectorXd out(p.size());
  for(auto i = 0; i < p.size(); i++) {
    out(i) = log(p.coeff(i) / (1.0 - p.coeff(i)));
  }
  return out;
}

double logdlogistic(const Eigen::VectorXd& x) {
  Eigen::VectorXd out(x.size());
  for(auto i = 0; i < x.size(); i++) {
    out(i) = -x.coeff(i) - 2.0 * log1p(exp(-x.coeff(i)));
  }
  return out.sum();
}

// [[Rcpp::export]]
Rcpp::List f_normal(const Eigen::MatrixXd& centers,
                    const Eigen::MatrixXd& XI,
                    const Eigen::MatrixXd& XmI,
                    const Eigen::VectorXd& yI,
                    const Eigen::VectorXd& ymI,
                    const size_t M,
                    const size_t n) {
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  for(size_t m = 0; m < ncenters; m++) {
    Eigen::MatrixXd H(q, q);
    H << XI, qnorm(centers.col(m));
    Eigen::MatrixXd Ht = H.transpose();
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(Ht * H);
    if(lu.isInvertible()) {
      Eigen::VectorXd theta = lu.inverse() * Ht * yI;
      double sigma = theta.coeff(q - 1);
      if(sigma > 0) {
        Eigen::VectorXd v = ymI - XmI * theta.topRows(q - 1);
        J(counter) = logdnorm(v / sigma) - (n - q) * log(sigma);
        Theta.col(counter) = theta;
        counter++;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Theta") = Theta.leftCols(counter).transpose(),
    Rcpp::Named("logWeights") = J.topRows(counter));
}

// [[Rcpp::export]]
Rcpp::List f_cauchy(const Eigen::MatrixXd& centers,
                    const Eigen::MatrixXd& XI,
                    const Eigen::MatrixXd& XmI,
                    const Eigen::VectorXd& yI,
                    const Eigen::VectorXd& ymI,
                    const size_t M,
                    const size_t n) {
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  for(size_t m = 0; m < ncenters; m++) {
    Eigen::MatrixXd H(q, q);
    H << XI, qcauchy(centers.col(m));
    Eigen::MatrixXd Ht = H.transpose();
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(Ht * H);
    if(lu.isInvertible()) {
      Eigen::VectorXd theta = lu.inverse() * Ht * yI;
      double sigma = theta.coeff(q - 1);
      if(sigma > 0) {
        Eigen::VectorXd v = ymI - XmI * theta.topRows(q - 1);
        J(counter) = logdcauchy(v / sigma) - (n - q) * log(sigma);
        Theta.col(counter) = theta;
        counter++;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Theta") = Theta.leftCols(counter).transpose(),
    Rcpp::Named("logWeights") = J.topRows(counter));
}

// [[Rcpp::export]]
Rcpp::List f_student(const Eigen::MatrixXd& centers,
                     const Eigen::MatrixXd& XI,
                     const Eigen::MatrixXd& XmI,
                     const Eigen::VectorXd& yI,
                     const Eigen::VectorXd& ymI,
                     const size_t M,
                     const size_t n,
                     const double nu) {
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  for(size_t m = 0; m < ncenters; m++) {
    Eigen::MatrixXd H(q, q);
    H << XI, qt(centers.col(m), nu);
    Eigen::MatrixXd Ht = H.transpose();
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(Ht * H);
    if(lu.isInvertible()) {
      Eigen::VectorXd theta = lu.inverse() * Ht * yI;
      double sigma = theta.coeff(q - 1);
      if(sigma > 0) {
        Eigen::VectorXd v = ymI - XmI * theta.topRows(q - 1);
        J(counter) = logdt(v / sigma, nu) - (n - q) * log(sigma);
        Theta.col(counter) = theta;
        counter++;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Theta") = Theta.leftCols(counter).transpose(),
    Rcpp::Named("logWeights") = J.topRows(counter));
}

// [[Rcpp::export]]
Rcpp::List f_logistic(const Eigen::MatrixXd& centers,
                      const Eigen::MatrixXd& XI,
                      const Eigen::MatrixXd& XmI,
                      const Eigen::VectorXd& yI,
                      const Eigen::VectorXd& ymI,
                      const size_t M,
                      const size_t n) {
  const size_t ncenters = centers.cols();
  const size_t q = XI.cols() + 1;
  size_t counter = 0;
  Eigen::VectorXd J(M);
  Eigen::MatrixXd Theta(q, M);
  for(size_t m = 0; m < ncenters; m++) {
    Eigen::MatrixXd H(q, q);
    H << XI, qlogistic(centers.col(m));
    Eigen::MatrixXd Ht = H.transpose();
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(Ht * H);
    if(lu.isInvertible()) {
      Eigen::VectorXd theta = lu.inverse() * Ht * yI;
      double sigma = theta.coeff(q - 1);
      if(sigma > 0) {
        Eigen::VectorXd v = ymI - XmI * theta.topRows(q - 1);
        J(counter) = logdlogistic(v / sigma) - (n - q) * log(sigma);
        Theta.col(counter) = theta;
        counter++;
      }
    }
  }
  return Rcpp::List::create(
      Rcpp::Named("Theta") = Theta.leftCols(counter).transpose(),
      Rcpp::Named("logWeights") = J.topRows(counter));
}
