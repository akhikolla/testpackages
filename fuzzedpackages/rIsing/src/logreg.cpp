// [[Rcpp::depends(RcppEigen)]]
#include "misc.h"
#include <Rcpp.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::Map<MatrixXd> MapMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Ref<const VectorXd> cRefVec;
typedef Eigen::Ref<VectorXd> RefVec;
typedef Eigen::Ref<const MatrixXd> cRefMat;
typedef Eigen::Ref<MatrixXd> RefMat;

class LogReg
{
private:
  const cRefMat X;
  const cRefVec y;
  const cRefVec b;
  const unsigned int n;
  const unsigned int p1;

protected:
  lbfgsfloatval_t *m_w;

public:
  LogReg(const cRefMat& X, const cRefVec& y, const cRefVec& b) :
  X(X), y(y), b(b), n(b.sum()), p1(X.cols()), m_w(NULL) {}

  virtual ~LogReg()
  {
    if (m_w != NULL) {
      lbfgs_free(m_w);
      m_w = NULL;
    }
  }

  MatrixXd run_regpath(const std::vector<double> lambda) {
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_w = lbfgs_malloc(p1);
    MapVec wvec(m_w,p1);
    wvec.setZero();
    const double y1 = y.sum();
    wvec(0) = std::log(y1/(n-y1));

    const unsigned int nlambda = lambda.size();
    MatrixXd wmat(p1, nlambda);

    // Init parameters
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.orthantwise_start = 1;
    param.orthantwise_c = lambda[0];
    param.linesearch = 2;
    lbfgs(p1, m_w, &fx, _evaluate, NULL, this, &param);
    wmat.col(0) = wvec;

    for (unsigned int i = 1; i < nlambda; i++) {
      param.orthantwise_c = lambda[i];
      lbfgs(p1, m_w, &fx, _evaluate, NULL, this, &param);
      wmat.col(i) = wvec;
    }

    return wmat;
  }

protected:
  static lbfgsfloatval_t _evaluate(
    void *instance,
    const lbfgsfloatval_t *w,
    lbfgsfloatval_t *g,
    const int p1,
    const lbfgsfloatval_t step
  )
  {
    return reinterpret_cast<LogReg*>(instance)->evaluate(w, g, p1, step);
  }

  lbfgsfloatval_t evaluate(
      const lbfgsfloatval_t *w,
      lbfgsfloatval_t *g,
      const int p1,
      const lbfgsfloatval_t step
  )
  {
    Eigen::Map<const VectorXd> wvec(w, p1);
    MapVec gvec(g, p1);

    VectorXd Xw = X*wvec;
    double yXw = Xw.dot(y);
    Xw = Xw.array().exp();
    lbfgsfloatval_t fx = ((b.array() * (Xw.array()+ 1.0).log()).sum() - yXw) / n;

    Xw.array() /= (Xw.array() + 1.0);
    gvec.noalias() = ( X.transpose() * ( ( b.array() * Xw.array() ).matrix() - y ) ) / n;

    return fx;
  }
};

VectorXd loglik(const cRefMat& wmat, const cRefMat& X1, const cRefVec& y, const cRefVec& b) {
  MatrixXd Xwmat = X1 * wmat;
  return ((Xwmat.array().colwise() * y.array()) - ( (1 + Xwmat.array().exp()).log().colwise() * b.array() ) ).colwise().sum().transpose();
}

// Rcpp functions
using namespace Rcpp;

// [[Rcpp::export]]
List logreg_cpp(const NumericMatrix& X_,
                const NumericVector& y_,
                const NumericVector& b_,
                const NumericVector& means,
                const NumericVector& sds,
                const std::vector<double>& lambda) {
  const MapMat X(as<MapMat>(X_));
  const MapVec y(as<MapVec>(y_));
  const MapVec b(as<MapVec>(b_));
  const unsigned int n = b.sum();
  const unsigned int n1 = X.rows();
  const MatrixXd X1 = (MatrixXd(n1, X.cols() + 1) << VectorXd::Ones(n1), X).finished();

  LogReg lr(X1, y, b);
  MatrixXd wmat = lr.run_regpath(lambda);
  VectorXd logliks = loglik(wmat, X1, y, b) / n;
  rescale(wmat, as<MapVec>(means), as<MapVec>(sds));
  return List::create(_["wmat"] = wrap(wmat), _["logliks"] = wrap(logliks));
}

// [[Rcpp::export]]
List logreg_setup(const NumericVector& X,
                  const NumericVector& y,
                  const bool scale,
                  const bool regpath,
                  const unsigned int nlambda,
                  const double lambda_min_ratio) {
  const MapMat X_(as<MapMat>(X));
  const MapVec y_(as<MapVec>(y));
  const X_data X_dat(X_, scale);
  if (regpath) {
    const std::vector<double> lambda = X_dat.construct_regpath(y_, nlambda, lambda_min_ratio);
    return List::create(_["means"] = wrap(X_dat.means),
                        _["sds"] = wrap(X_dat.sds),
                        _["Xs"] = wrap(X_dat.Xs),
                        _["lambda"] = wrap(lambda));
  } else {
    return List::create(_["means"] = wrap(X_dat.means),
                        _["sds"] = wrap(X_dat.sds),
                        _["Xs"] = wrap(X_dat.Xs));
  }
}

// [[Rcpp::export]]
std::vector<double> regpath_ising(const NumericMatrix& Xs_, const NumericVector& y_, const unsigned int nlambda, const double lambda_min_ratio) {
  MapMat Xs(as<MapMat>(Xs_));
  MapVec y(as<MapVec>(y_));
  const unsigned int n = Xs_.rows();
  const double max_lambda = (( Xs * std::sqrt( (double)n/((double)n-1) ) ).array().colwise() * y.array()).colwise().sum().abs().maxCoeff() / n;
  const double min_lambda = max_lambda * lambda_min_ratio;
  return exp_c(seq_c(std::log(max_lambda), std::log(min_lambda), nlambda));
}



// LOGREG TYPE 2
class LogReg2
{
private:
  const cRefMat X;
  const cRefVec y;
  const unsigned int n;
  const unsigned int p1;

protected:
  lbfgsfloatval_t *m_w;

public:
  LogReg2(const cRefMat& X, const cRefVec& y) :
  X(X), y(y), n(X.rows()), p1(X.cols()), m_w(NULL) {}

  virtual ~LogReg2()
  {
    if (m_w != NULL) {
      lbfgs_free(m_w);
      m_w = NULL;
    }
  }

  MatrixXd run_regpath(const std::vector<double> lambda) {
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *m_w = lbfgs_malloc(p1);
    MapVec wvec(m_w,p1);
    wvec.setZero();
    const double y1 = y.sum();
    wvec(0) = std::log(y1/(n-y1));

    const unsigned int nlambda = lambda.size();
    MatrixXd wmat(p1, nlambda);

    // Init parameters
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.orthantwise_start = 1;
    param.orthantwise_c = lambda[0];
    param.linesearch = 2;
    lbfgs(p1, m_w, &fx, _evaluate, NULL, this, &param);
    wmat.col(0) = wvec;

    for (unsigned int i = 1; i < nlambda; i++) {
      param.orthantwise_c = lambda[i];
      lbfgs(p1, m_w, &fx, _evaluate, NULL, this, &param);
      wmat.col(i) = wvec;
    }

    return wmat;
  }

protected:
  static lbfgsfloatval_t _evaluate(
      void *instance,
      const lbfgsfloatval_t *w,
      lbfgsfloatval_t *g,
      const int p1,
      const lbfgsfloatval_t step
  )
  {
    return reinterpret_cast<LogReg2*>(instance)->evaluate(w, g, p1, step);
  }

  lbfgsfloatval_t evaluate(
      const lbfgsfloatval_t *w,
      lbfgsfloatval_t *g,
      const int p1,
      const lbfgsfloatval_t step
  )
  {
    Eigen::Map<const VectorXd> wvec(w, p1);
    MapVec gvec(g, p1);

    VectorXd Xw = X*wvec;
    double yXw = Xw.dot(y);
    Xw = Xw.array().exp();
    lbfgsfloatval_t fx = ((Xw.array()+ 1.0).log().sum() - yXw) / n;

    Xw.array() /= (Xw.array() + 1.0);
    gvec.noalias() = ( X.transpose() * ( Xw - y ) ) / n;

    return fx;
  }
};

VectorXd loglik2(const cRefMat& wmat, const cRefMat& X1, const cRefVec& y) {
  MatrixXd Xwmat = X1 * wmat;
  return ((Xwmat.array().colwise() * y.array()) - (1 + Xwmat.array().exp()).log() ).colwise().sum().transpose();
}

// [[Rcpp::export]]
List logreg_cpp2(const NumericMatrix& X_,
                 const NumericVector& y_,
                 const Nullable<NumericVector>& lambda,
                 const unsigned int nlambda,
                 const double lambda_min_ratio,
                 const bool scale) {
  const MapMat X(as<MapMat>(X_));
  const MapVec y(as<MapVec>(y_));
  const unsigned int n = X.rows();
  const X_data X_dat(X, scale);
  const MatrixXd X1 = (MatrixXd(n, X.cols() + 1) << VectorXd::Ones(n), X_dat.Xs).finished();

  std::vector<double>(lambda_vec);
  if (lambda.isNull()) {
    lambda_vec = X_dat.construct_regpath(y, nlambda, lambda_min_ratio);
  } else {
    lambda_vec = as<std::vector<double> >(lambda);
    std::sort(lambda_vec.rbegin(), lambda_vec.rend());
  }

  LogReg2 lr(X1, y);
  MatrixXd wmat = lr.run_regpath(lambda_vec);
  VectorXd logliks = loglik2(wmat, X1, y) / n;
  rescale(wmat, X_dat.means, X_dat.sds);
  return List::create(_["wmat"] = wrap(wmat), _["logliks"] = wrap(logliks), _["lambda"] = wrap(lambda_vec));
}
