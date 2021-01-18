#include <RcppEigen.h>
#include "lbfgs/lbfgs.c"

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Ref<const VectorXd> cRefVec;
typedef Eigen::Ref<VectorXd> RefVec;
typedef Eigen::Ref<const MatrixXd> cRefMat;
typedef Eigen::Ref<MatrixXd> RefMat;

// seq function
std::vector<double> seq_c(const double start, const double end, const unsigned int length) {
  std::vector<double> v(length);
  if (length <= 1) {
    v[0] = start;
  } else {
    double step = (end - start) / (length-1);
    for (unsigned int i = 0; i < length; i++) {
      v[i] = start + i*step;
    }
  }
  return v;
}

// exp vector function
struct exp_op { double operator() (double d) const { return std::exp(d); } };
std::vector<double> exp_c(const std::vector<double>& x) {
  std::vector<double> y(x.size());
  std::transform(x.begin(), x.end(), y.begin(), exp_op());
  return y;
}

// data struct
struct X_data {
  const cRefMat X;
  const bool scale;
  const unsigned int n;
  const unsigned int p;
  const VectorXd means;
  const VectorXd sds;
  const MatrixXd Xs;

  X_data(const cRefMat& X, const bool scale) :
    X(X), scale(scale), n(X.rows()), p(X.cols()), means(calcMeans()), sds(calcSd()), Xs(scaleX()) {}

  VectorXd calcMeans() const {
    VectorXd means;
    if (scale) {
      means = X.colwise().mean();
    } else {
      means = VectorXd::Zero(p);
    }
    return means;
  }

  VectorXd calcSd() const {
    VectorXd sds;
    if (scale) {
      sds = ((X.rowwise() - means.transpose()).array().square().colwise().sum() / (n-1) ).sqrt();
    } else {
      sds = VectorXd::Ones(p);
    }
    return sds;
  }

  MatrixXd scaleX() const {
    if (scale) {
      return (X.rowwise() - means.transpose()).array().rowwise() / sds.transpose().array();
    } else {
      return X;
    }
  }

  std::vector<double> construct_regpath(const cRefVec& y, const unsigned int nlambda, const double lambda_min_ratio) const {
    const double max_lambda = (( Xs * std::sqrt( (double)n/((double)n-1) ) ).array().colwise() * y.array()).colwise().sum().abs().maxCoeff() / n;
    const double min_lambda = max_lambda * lambda_min_ratio;
    return exp_c(seq_c(std::log(max_lambda), std::log(min_lambda), nlambda));
  }
};

void rescale(RefMat wmat, const cRefVec& means, const cRefVec& sds) {
  RefMat coefs = wmat.block(1,0,wmat.rows()-1,wmat.cols());
  coefs.array().colwise() /= sds.array();
  wmat.row(0).array() -= (coefs.array().colwise() * means.array()).colwise().sum();
}

