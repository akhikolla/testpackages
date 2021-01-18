#include <RcppArmadilloExtensions/sample.h>
#define NDEBUG 1
#include <RcppNumerical.h>
#include <math.h>
#include <Eigen/LU>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]

typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

//================== Part 1: Class Definitions =========================
//========================================================================

class LogReg_MLE: public Numer::MFuncGrad
{
private:
  const MapMat X;
  const MapVec Y;
public:
  LogReg_MLE(const MapMat x_, const MapVec y_) : X(x_), Y(y_) {}
  
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
  {
    Eigen::VectorXd xbeta;
    xbeta.noalias() = X * beta;
    Eigen::VectorXd eXB = xbeta.array().exp();
    Eigen::VectorXd eXBpOne = eXB.array()+1.0;
    Eigen::VectorXd prob = eXBpOne.array().log();
    const double yxbeta = Y.dot(xbeta);
    const double f = prob.sum()-yxbeta;
    prob = (xbeta - prob).array().exp();
    grad.noalias() = X.transpose() * (prob - Y);
    
    return f;
  }
};
//=======================================================================//

class LogReg_MAP: public Numer::MFuncGrad
{
private:
  const MapVec y;
  const MapMat XX;
  double tau;
  double r;
  const int nlptype;
public:
  LogReg_MAP(const MapMat x_, const MapVec y_, double tau_, double r_, const int nlptype_) :
  y(y_), XX(x_), tau(tau_), r(r_), nlptype(nlptype_) {}
  
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
  {
    Eigen::VectorXd XB;
    XB.noalias() = XX * beta;
    if ((XB.array() > 500).any()){
      for (int i=0; i < XB.size(); i++)
        if (XB(i) >500) XB(i)=500;
    }
    Eigen::VectorXd eXB = XB.array().exp();
    Eigen::VectorXd eXBpOne = eXB.array()+1.0;
    Eigen::VectorXd logstc = eXB.array()/eXBpOne.array();
    Eigen::MatrixXd Xt = XX.transpose();
    const double ypXB = y.dot(XB);
    
    const double negloglik = eXBpOne.array().log().sum()-ypXB;
    Eigen::VectorXd negloglikgrad = Xt*(logstc - y);
    
    double neglogprior;
    Eigen::VectorXd beta2, neglogpriorgrad;
    switch(nlptype){
    case 0:
      beta2 = (r+1)*beta.array().abs().log()+tau/(beta.array().pow(2)); beta2(0) = 0;
      neglogprior = beta2.array().sum();
      neglogpriorgrad = -2*tau/(beta.array().pow(3)) + (r+1)/beta.array();
      neglogpriorgrad(0) = 0;      
      break;
    case 1:
      beta2 = beta.array().pow(2)/(2*tau)-2*r*beta.array().log(); beta2(0) = 0;
      neglogprior = beta2.array().sum();
      neglogpriorgrad = beta.array()/tau - 2*r/beta.array();
      neglogpriorgrad(0) = 0;      
      break;
    default:
      beta2 = (r+1)*beta.array().abs().log()+tau/(beta.array().pow(2)); beta2(0) = 0;
      neglogprior = beta2.array().sum();
      neglogpriorgrad = -2*tau/(beta.array().pow(3)) + (r+1)/beta.array();
      neglogpriorgrad(0) = 0;      
      break;      
    }
    
    const double neglogposterior = negloglik + neglogprior;
    grad.noalias() = negloglikgrad + neglogpriorgrad;
    return neglogposterior;
  }
};
//=======================================================================//

class LogReg_LogMarginal
{
private:
  const MapMat XX;
  Eigen::VectorXd beta;
  const MapVec y;
  double tau;
  double r;
  const int nlptype;
  
public:
  LogReg_LogMarginal(const MapMat x_, Eigen::VectorXd beta_, const MapVec y_, double tau_, double r_, const int nlptype_) :
  XX(x_), beta(beta_), y(y_), tau(tau_), r(r_), nlptype(nlptype_) {}
  
  double marginal_prob() {
    Eigen::VectorXd XB;
    XB.noalias() = XX * beta;
    Eigen::VectorXd eXB = XB.array().exp();
    Eigen::VectorXd eXBpOne = eXB.array()+1.0;
    Eigen::VectorXd B2 = beta.array().square();
    const double ypXB = y.dot(XB);
    int  k = XX.cols();
    double cov;
    double logCond;
    double logPrior;
    double logMarginal;
    Eigen::VectorXd aux;
    
    Eigen::MatrixXd hessian(k,k);
    hessian.setZero();
    Eigen::VectorXd denom = eXBpOne.array().square();
    Eigen::VectorXd SubLogistic = eXB.array()/denom.array();
    
    Eigen::VectorXd diagvec, lp1, lp2, neglogpriorgrad;
    switch(nlptype){
    case 0:
      diagvec = (r+1)/B2.array()-6*tau*beta.array().pow(-4); diagvec(0)=0;
      lp1 = -(r+1)*beta.array().abs().log(); lp1(0) = 0;
      lp2 = tau/B2.array(); lp2(0) = 0;
      logPrior = (k-1)*(r*log(tau)/2-lgamma(r/2))+(lp1-lp2).sum();
      break;
    case 1:
      diagvec = 2*r/B2.array() + 1/tau; diagvec(0)=0;
      lp1 = -B2.array()/(2*tau)+2*r*beta.array().log(); lp1(0) = 0;
      logPrior = -(k-1)*(log(2*M_PI)/2+(r+0.5)*log(tau))+lp1.sum();      
      break;
    default:
      diagvec = (r+1)/B2.array()-6*tau*beta.array().pow(-4); diagvec(0)=0;
      lp1 = -(r+1)*beta.array().abs().log(); lp1(0) = 0;
      lp2 = tau/B2.array(); lp2(0) = 0;
      logPrior = (k-1)*(r*log(tau)/2-lgamma(r/2))+(lp1-lp2).sum();
      break;      
    }
    
    for(int i=0; i < k-1; i++){
      Eigen::VectorXd aux = XX.col(i).array().square();
      hessian(i,i) = diagvec(i)-SubLogistic.dot(aux);
      for(int j=i+1; j < k; j++){
        aux = XX.col(i).array()*XX.col(j).array();
        hessian(i,j) = hessian(j,i) = -SubLogistic.dot(aux);
      }
    }
    aux = XX.col(k-1).array().square(); // set the value for the last element.
    hessian(k-1,k-1) = diagvec(k-1)-SubLogistic.dot(aux);
    cov = 1/(-hessian).determinant();
    logCond = ypXB-eXBpOne.array().log().sum();
    if (cov < 0){
      return (-999999);
    }else{
      logMarginal = 0.5*(k*log(2*M_PI)+log(cov))+logCond+logPrior;
      return(logMarginal);
    }
  }
};
//=======================================================================//

class Cox_MLE: public Numer::MFuncGrad 
{
private:
  const MapMat X;
  const MapVec status;
  const MapVec xbeta_m;
  
public:
  Cox_MLE(const MapMat x_, const MapVec status_, const MapVec xbeta_m_) : X(x_),
  status(status_), xbeta_m(xbeta_m_) {}
  
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
  {
    int n = X.rows();
    int p = X.cols();
    double snewexbval;
    Eigen::MatrixXd XS(p,n);
    XS.setZero();
    Eigen::MatrixXd Xtr = X.transpose();
    Eigen::MatrixXd AX(p,n);
    Eigen::VectorXd sneweXB(n);
    Eigen::VectorXd auxvec(p);
    Eigen::VectorXd xbeta;
    xbeta.noalias() = X * beta;
    if ((xbeta.array() > 500).any()){ 
      for (int i=0; i < xbeta.size(); i++)
        if (xbeta(i) >500) xbeta(i)=500;
    }
    Eigen::VectorXd tmp_xb = xbeta + xbeta_m;
    Eigen::VectorXd enewXB = tmp_xb.array().exp();
    
    snewexbval = 0;
    auxvec.setZero();
    for (int i=n-1; i >= 0; i--){
      snewexbval = enewXB(i) + snewexbval;
      auxvec = enewXB(i)*Xtr.col(i) + auxvec;
      XS.col(i) = auxvec/snewexbval;
      sneweXB(i) = snewexbval;
    }
    
    Eigen::VectorXd lnewseXB = sneweXB.array().log();
    const double f = status.dot(lnewseXB-xbeta);
    grad.noalias() = (XS-Xtr) * status;
    
    return f;
  }
};
//=======================================================================//

class Cox_MAP: public Numer::MFuncGrad 
{
private:
  const MapMat XX;
  const MapVec status;
  double tau;
  double r;
  const int nlptype;
  
public:
  Cox_MAP(const MapMat x_, const MapVec status_, double tau_, double r_, const int nlptype_) :
  XX(x_), status(status_), tau(tau_), r(r_), nlptype(nlptype_) {}
  
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
  {
    int n = XX.rows();
    int p = XX.cols();
    double snewexbval;
    Eigen::MatrixXd XS(p,n);
    XS.setZero();
    Eigen::MatrixXd Xtr = XX.transpose();
    Eigen::MatrixXd AX(p,n);
    Eigen::VectorXd sneweXB(n);
    Eigen::VectorXd auxvec(p);
    Eigen::VectorXd xbeta;
    xbeta.noalias() = XX * beta;
    Eigen::VectorXd enewXB = xbeta.array().exp();
    
    snewexbval = 0;
    auxvec.setZero();
    for (int i=n-1; i >= 0; i--){
      snewexbval = enewXB(i) + snewexbval;
      auxvec = enewXB(i)*Xtr.col(i) + auxvec;
      XS.col(i) = auxvec/snewexbval;
      sneweXB(i) = snewexbval;
    }
    
    Eigen::VectorXd lnewseXB = sneweXB.array().log();
    const double negloglik = status.dot(lnewseXB-xbeta);
    Eigen::VectorXd negloglikgrad  = (XS-Xtr) * status;
    
    double neglogprior;
    Eigen::VectorXd beta2, neglogpriorgrad;
    switch(nlptype){
    case 0:
      beta2 = (r+1)*beta.array().abs().log()+tau/(beta.array().pow(2));
      neglogprior = beta2.array().sum();
      neglogpriorgrad = -2*tau/(beta.array().pow(3)) + (r+1)/beta.array();
      break;
    case 1:
      beta2 = beta.array().pow(2)/(2*tau)-2*r*beta.array().log();
      neglogprior = beta2.array().sum();
      neglogpriorgrad = beta.array()/tau - 2*r/beta.array();
      break;
    default:
      beta2 = (r+1)*beta.array().abs().log()+tau/(beta.array().pow(2));
      neglogprior = beta2.array().sum();
      neglogpriorgrad = -2*tau/(beta.array().pow(3)) + (r+1)/beta.array();
      break;
    }
    
    const double neglogposterior = negloglik + neglogprior;
    grad.noalias() = negloglikgrad + neglogpriorgrad; 
    return neglogposterior;
  }
};
//=======================================================================//

class Cox_LogMarginal
{
private:
  const MapMat XX;
  const MapVec beta;
  const MapVec status;
  double tau;
  double r;
  const int nlptype;
  
public:
  Cox_LogMarginal(const MapMat x_, const MapVec beta_, const MapVec status_, double tau_, double r_, const int nlptype_) :
  XX(x_), beta(beta_), status(status_), tau(tau_), r(r_), nlptype(nlptype_) {}
  
  double marginal_prob(){
    double cov;
    double logCond;
    double logPrior;
    double logMarginal;
    int n = XX.rows();
    int p = XX.cols();
    double snewexbval;
    Eigen::VectorXd xjc, U, W;
    Eigen::MatrixXd hessian(p,p);
    Eigen::MatrixXd temphess(p,p);
    Eigen::MatrixXd XS(p,n);
    XS.setZero();
    Eigen::MatrixXd Xtr = XX.transpose();
    Eigen::MatrixXd AX(p,n);
    Eigen::MatrixXd A(n,p);
    Eigen::MatrixXd tmp_mI, mI, AuxMat, Diagxjc;
    Eigen::VectorXd sneweXB(n);
    Eigen::VectorXd auxvec(p);
    Eigen::VectorXd auxvecH(p);
    Eigen::VectorXd xbeta;
    xbeta.noalias() = XX * beta;
    Eigen::VectorXd enewXB = xbeta.array().exp();
    Eigen::VectorXd B2 = beta.array().square();
    
    snewexbval = 0;
    auxvec.setZero();
    for (int i=n-1; i >= 0; i--){
      snewexbval = enewXB(i) + snewexbval;
      auxvec = enewXB(i)*Xtr.col(i) + auxvec;
      XS.col(i) = auxvec/snewexbval;
      sneweXB(i) = snewexbval;
    }
    
    Eigen::VectorXd lnewseXB = sneweXB.array().log();
    
    Eigen::VectorXd diagvec, lp1, lp2;
    switch(nlptype){
    case 0:
      diagvec = (r+1)/B2.array()-6*tau*beta.array().pow(-4);
      lp1 = -(r+1)*beta.array().abs().log();
      lp2 = tau/B2.array();
      logPrior = p*(r*log(tau)/2-lgamma(r/2))+(lp1-lp2).sum();
      break;
    case 1:
      diagvec = 2*r/B2.array() + 1/tau;
      lp1 = -B2.array()/(2*tau)+2*r*beta.array().log();
      logPrior = -p*(log(2*M_PI)/2+(r+0.5)*log(tau))+lp1.sum();      
      break;
    default:
      diagvec = (r+1)/B2.array()-6*tau*beta.array().pow(-4);
      lp1 = -(r+1)*beta.array().abs().log();
      lp2 = tau/B2.array();
      logPrior = p*(r*log(tau)/2-lgamma(r/2))+(lp1-lp2).sum();
      break;      
    }
    
    for(int s=0; s < p; s++){
      auxvecH.setZero();
      for (int i=n-1; i>=0; i--){
        auxvecH = Xtr.col(i)*XX(i,s)*enewXB(i) + auxvecH;
        U = auxvecH/sneweXB(i);
        W = XS(s,i)*XS.col(i);
        AuxMat = U - W;
        A.row(i) = AuxMat.transpose();
      }
      hessian.row(s) = status.transpose() * A;
    }
    
    Eigen::MatrixXd DiagPriorHess = diagvec.asDiagonal();
    temphess = hessian;
    hessian = temphess - DiagPriorHess;
    
    cov = 1/(hessian).determinant();
    if (cov < 0){
      return(-999999);
    } else {
      logCond = status.dot(xbeta-lnewseXB);
      logMarginal = 0.5*(p*log(2*M_PI)+log(cov))+logCond+logPrior;
      return(logMarginal);
    }
  }
};
//================== Part 2: Function Definitions ========================
//========================================================================

arma::mat cox_order_vecs(const arma::mat exmat) {
  arma::vec time = exmat.col(0);
  arma::vec status = exmat.col(1);
  int n = exmat.n_rows;
  int p = exmat.n_cols;
  arma::mat outmat(n,p);
  arma::uvec idx;
  idx.zeros(n); // = no_init(n);
  auto begin = std::begin(idx), end = std::end(idx);
  std::iota(begin, end, static_cast<size_t>(0));
  auto comparator = [&](const size_t & a, const size_t & b){ return time[a] < time[b]; };
  std::sort(begin, end, comparator);
  for (int i=0; i < n; i++){
    outmat.row(i) = exmat.row(idx[i]);
  }
  return(outmat);
}
//========================================================================

arma::uvec order_c(arma::vec in_vec, int d) {
  int n = in_vec.size();
  arma::uvec idx;
  arma::uvec out_vec(d);
  idx.zeros(n);
  auto begin = std::begin(idx), end = std::end(idx);
  std::iota(begin, end, static_cast<size_t>(0));
  auto comparator = [&](const size_t & a, const size_t & b){ return in_vec[a] < in_vec[b]; };
  std::sort(begin, end, comparator);
  std::reverse(begin, end);
  for (int i=0; i < d; i++){
    out_vec(i) = idx[i];
  }
  return(out_vec);
}
//=======================================================================//

arma::uvec c_setdiff(arma::uvec a, arma::uvec b){
  arma::uvec bb = arma::sort(b);
  std::vector<int> outvec;
  std::set_difference(a.begin(),a.end(),bb.begin(),bb.end(),
                      std::inserter(outvec,outvec.begin()));
  arma::uvec out = arma::conv_to<arma::uvec>::from(outvec);
  return(out);
}
//==================================================//

arma::uvec c_union(arma::uvec a, arma::uvec b){
  int s1 = a.size();
  int s2 = b.size();
  int k = s1+s2;
  arma::uvec out(k);
  for (int i=0; i< s1; i++) out(i) = a(i);
  for (int i=s1; i < k; i++) out(i) = b(i-s1);
  arma::uvec out_ind = find_unique(out);
  arma::uvec fout = out(out_ind);
  return(fout);
}
//==================================================//

double log_msize_prob(int p,int k,int a, int b){
  double z;
  z = lgamma(a+k)+lgamma(b+p-k)+lgamma(a+b)-lgamma(a)-lgamma(b)-lgamma(a+b+p);
  return(z);
}
//==================================================//

arma::uvec push_begin(arma::uvec a, int val){
  arma::uvec av(1); av(0) = val;
  a.insert_rows(0,av.row(0));
  return(a);
}
//==================================================//

arma::uvec rm_begin(arma::uvec a){
  arma::uvec b; b = a;
  b.shed_row(0);
  return(b);
}
//==================================================//
arma::uvec cox_c_sample(arma::uvec vector, int n, arma::vec probs, bool replace){
  arma::uvec yy = Rcpp::RcppArmadillo::sample(vector,n,replace,probs);
  return(yy);
}
//==================================================//

int my_find(arma::vec a, double num){
  std::vector<double> x = arma::conv_to< std::vector<double> >::from(a);
  std::vector<double>::iterator it;
  it = std::find(x.begin(), x.end(), num);
  if (it != x.end())
    return(it-x.begin());
  else
    return(-1);
}
//==================================================//

double calc_key(arma::uvec b){
  double out = 0.0;
  arma::uvec a = arma::sort(b);
  int s = a.size();
  double val;
  for(int i=0; i < s; i++){
    val = a(i) + 1.0;
    out = out + pow(2,log(val)) + M_PI*log(val);
  }
  return(out);
}
//==================================================//

arma::vec repzero(arma::vec a, double val){
  arma::vec b = a;
  arma::uvec c = find(b==0);
  b(c) = val*arma::ones<arma::vec>(c.size());
  return(b);
}
//==================================================//

bool isequal(arma::uvec a, arma::uvec b){
  bool y;
  (std::equal(a.begin(), a.end(), b.begin())) ? y=true : y=false;
  return(y);
}
//==================================================//

arma::uvec seq_gen(int k){
  arma::uvec idx(k);
  std::iota(std::begin(idx), std::end(idx), static_cast<size_t>(0));
  return(idx);
}
//==================================================//

double LogReg_Model_Prob(const MapMat XX, const MapVec y, double tau, double r, int a, int b, int n, int p, const int nlptype)
{
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  int k = XX.cols();
  Eigen::VectorXd glm_coef(k);
  glm_coef.setZero();
  
  LogReg_MLE nll(XX,y);
  double fopt_glm;
  int glm_status = Numer::optim_lbfgs(nll, glm_coef, fopt_glm, maxit, eps_f, eps_g);
  if (glm_status < 0) return(-999999);
  Eigen::VectorXd beta = glm_coef;
  
  LogReg_MAP CM(XX,y,tau,r,nlptype);
  double fopt;
  int status = Numer::optim_lbfgs(CM, beta, fopt, maxit, eps_f, eps_g);
  if (status < 0) return(-999999);
  Eigen::VectorXd betavec = beta;
  
  LogReg_LogMarginal calclm(XX,betavec,y,tau,r,nlptype);
  double logmarg = calclm.marginal_prob();
  double model_prob = logmarg + log_msize_prob(p,k-1,a,b);
  if (std::isnan(model_prob) || std::isinf(model_prob)){
    return(-999999);
  }else{
    return(model_prob);
  }
}
//==================================================//

double Cox_Model_Prob(const MapMat XX, const MapVec status, arma::vec beta, double tau, double r, int a, int b, int p, const int nlptype)
{
  int k = XX.cols();
  const MapVec betavec(beta.memptr(),k);
  Cox_LogMarginal calclm(XX,betavec,status,tau,r,nlptype);
  
  double logmarg = calclm.marginal_prob();
  double model_prob = logmarg + log_msize_prob(p,k,a,b);
  if (std::isnan(model_prob) || std::isinf(model_prob)){
    return(-999999);
  }else{
    return(model_prob);
  }
}
//==================================================//

Rcpp::NumericVector cox_beta_est(arma::mat cur_model, const MapVec Status, double tau, double r, const int nlptype){
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  int k = cur_model.n_cols;
  int n = cur_model.n_rows;
  arma::vec cur_xbeta = arma::zeros<arma::vec>(n);
  const MapVec Cur_XBeta(cur_xbeta.memptr(), n);
  MapMat fxx(cur_model.memptr(), n, k);
  Eigen::VectorXd coxph_coef(k);
  coxph_coef.setZero();
  
  Cox_MLE nll(fxx,Status,Cur_XBeta);
  double fopt_glm;
  int mle_status = Numer::optim_lbfgs(nll, coxph_coef, fopt_glm, maxit, eps_f, eps_g);
  if (mle_status < 0){
    return(Rcpp::wrap(-999999));
  } else {
    Eigen::VectorXd betahat = coxph_coef;
    Cox_MAP FM(fxx,Status,tau,r,nlptype);
    double fopt;
    int map_status = Numer::optim_lbfgs(FM, betahat, fopt, maxit, eps_f, eps_g);
    if (map_status < 0){
      return(Rcpp::wrap(-999999));
    } else {
      return(Rcpp::wrap(betahat));
    }
  }
}

//==================================================//
arma::uvec find_max_utils(const arma::mat ordexmat,arma::uvec kc_cols, const MapVec Status, const MapVec Cur_XBeta, int d)
{
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  arma::mat x_kc = ordexmat.cols(kc_cols);
  int k = kc_cols.size();
  arma::vec utl_vec(k);
  for (int i=0; i<k; i++){
    MapMat fdes(x_kc.colptr(i), x_kc.n_rows, 1);
    Eigen::VectorXd mle_coef(1);
    mle_coef.setZero();
    Cox_MLE cph(fdes,Status,Cur_XBeta);
    double fopt_coxph;
    int coxph_status = optim_lbfgs(cph, mle_coef, fopt_coxph, maxit, eps_f, eps_g);
    if (coxph_status < 0) utl_vec(i) = -999999; else utl_vec(i) = -fopt_coxph;
  }
  arma::uvec sel_inds = order_c(utl_vec,d);
  return(sel_inds);
}

//================== Part 3: The Main Functions =========================
//=======================================================================

// [[Rcpp::export]]
Rcpp::NumericVector null_mle_lreg(arma::mat& XX, int n, int p, int cons, double a, double b,
                                  double sprob, int niters)
{
  int ncol,k;
  double probs = 1/p;
  arma::uvec indices = seq_gen(p)+1;
  arma::uvec sel_cols, cur_cols;
  Eigen::VectorXd out_null(2000000), betahat;
  int counts = 0;
  double eps_f = 1.0e-08; 
  double eps_g = 1.0e-05;
  int maxit = 300;
  Rcpp::NumericVector yy,rbinout,prob_msize;
  arma::mat final_x;
  yy = Rcpp::rbinom(n,1,sprob);
  MapVec y = Rcpp::as<MapVec>(yy);
  
  for (int i=0; i < niters; i++){
    ncol = 0;
    while (ncol== 0 || ncol > cons){
      prob_msize = Rcpp::rbeta(1,a,b);
      rbinout = Rcpp::rbinom(1,p,prob_msize[0]);
      ncol = rbinout[0];
    }
    
    sel_cols = Rcpp::RcppArmadillo::sample(indices,ncol,false,probs);
    
    k = sel_cols.n_elem+1;
    arma::uvec cur_cols = push_begin(sel_cols,0);
    final_x = XX.cols(cur_cols);
    MapMat fxx(final_x.memptr(), n, k);
    Eigen::VectorXd glm_coef(k);
    glm_coef.setZero();
    
    LogReg_MLE nll(fxx,y);
    double fopt_glm;
    int glm_status = Numer::optim_lbfgs(nll, glm_coef, fopt_glm, maxit, eps_f, eps_g);
    if (glm_status >= 0){
      int t = (glm_coef.array().abs() > 7).sum();
      if (t==0){
        for (int j=0; j < k-1; j++){
          out_null(counts+j) = glm_coef(j+1);
        }
        counts = counts+(k-1);
      }
    }
  }
  betahat=out_null.head(counts-1);
  return(Rcpp::wrap(betahat));
}
//==================================================//

//' Non-parallel version of Bayesian variable selector for logistic regression
//' data using nonlocal priors
//' @description This function performs Bayesian variable selection for
//' logistic regression data in a non-parallel fashion. It does not contain any
//' pre-processing step or variable initialization. Moreover it does not have
//' the feature to be run in parallel for performing the coupling algorithm.
//' Therefore in general, it is NOT recommended to be used unless the user
//' knows how to initialize all the parameters. However, this function is
//' called by \code{\link{bvs}} function, the recommended way to run Bayesian
//' variable selection for such datasets.
//' @param exmat An extended matrix where the first column is binary resonse
//' vector and the rest is the design matrix which has its first column all 1
//' to account for the intercept in the model and is the output of
//' \code{PreProcess} code where the fixed columns are moved to the beginning.
//' @param chain1 The first chain or initial model where the MCMC algorithm
//' starts from. Note that the first \code{nf+1} elements are \code{1} where
//' \code{nf} is the number of fixed covariates that do not enter the selection
//' procedure and \code{1} is for the intercept.
//' @param nf The number of fixed covariates that do not enter the selection
//' procedure.
//' @param tau The paramter \code{tau} of the iMOM prior.
//' @param r The paramter \code{r} of the iMOM prior.
//' @param nlptype Determines the type of nonlocal prior that is used in the
//' analyses. \code{0} is for piMOM and \code{1} is for pMOM.
//' @param a The first parameter in beta distribution used as prior on model
//' size. This parameter is equal to 1 when uinform-binomial prior is used.
//' @param b The second paramter in beta distribution used as prior on model
//' size. This parameter is equal to 1 when uinform-binomial prior is used.
//' @param in_cons The average model size. This value under certain conditions
//' and when \code{p} is large, is equal to parameter \code{a} of the
//' beta-binomial prior.
//' @param loopcnt Number of iterations for MCMC procedure.
//' @param cplng A boolean variable indicating the coupling algorithm to be
//' performed or not.
//' @param chain2 Second chain or model for starting the MCMC procedure. This
//' parameter is only used when \code{cplng=TRUE}. Thus, it could be simply
//' set to \code{chain1} when it is not used.
//'
//' @return It returns a list containing following objects:
//' \item{max_chain}{A \code{1} by \code{p+1} binary vector showing the selected
//' model with maximum probability. \code{1} means a specific variable is
//' selected. The first variable is always the intercept.}
//' \item{beta_hat}{The coefficient vector for the selected model. The first one
//'  is always for the intercept.}
//' \item{max_prop}{The unnormalized probability of the model with highest
//' posterior probability.}
//'  \item{num_iterations}{The number of MCMC iterations that are executed.
//'  This is used when \code{cplng=TRUE} to check whether the total designated
//'  MCMC iterations were used or two chains are coupled sooner than that.}
//' \item{cplng_flag}{This is used when \code{cplng=TRUE} and indicates whether
//' two chains are coupled or not.}
//' \item{num_vis_models}{Number of visited models in search for the highest
//' probability model. This contains redundant models too and is not the number
//' of unique models.}
//' \item{hash_key}{This is only used when \code{cplng = FALSE}. This is a
//' vector containing real numbers uniquely assigned to each model for
//' distinguishing them.}
//' \item{hash_prob}{This is only used when \code{cplng = FALSE}. This is a
//' vector of probabilities for each visited model.}
//' \item{vis_covs}{This is only used when \code{cplng = FALSE}. This is a
//' list where each element contains indices of covariates for each visited
//' model.}
//' @author Amir Nikooienejad
//' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2016). Bayesian
//' variable selection for binary outcomes in high dimensional genomic studies
//' using nonlocal priors. Bioinformatics, 32(9), 1338-1345.\cr\cr
//' Nikooienejad, A., Wang, W., and Johnson, V. E. (2017). Bayesian Variable
//' Selection in High Dimensional Survival Time Cancer Genomic Datasets using
//' Nonlocal Priors. arXiv preprint, arXiv:1712.02964.\cr\cr
//' Johnson, V. E., and Rossell, D. (2010). On the use of non-local prior
//' densities in Bayesian hypothesis tests. Journal of the Royal Statistical
//' Society: Series B (Statistical Methodology), 72(2), 143-170.\cr\cr
//' Johnson, V. E. (1998). A coupling-regeneration scheme for
//' diagnosing convergence in Markov chain Monte Carlo algorithms. Journal of
//' the American Statistical Association, 93(441), 238-248.
//' @seealso \code{\link{bvs}}
//' @examples
//' ### Initializing parameters
//' n <- 200
//' p <- 40
//' set.seed(123)
//' Sigma <- diag(p)
//' full <- matrix(c(rep(0.5, p*p)), ncol=p)
//' Sigma <- full + 0.5*Sigma
//' cholS <- chol(Sigma)
//' Beta <- c(-1.7,1.8,2.5)
//' X <- matrix(rnorm(n*p), ncol=p)
//' X <- X%*%cholS
//' colnames(X) <- paste("gene_",c(1:p),sep="")
//' beta <- numeric(p)
//' beta[c(1:length(Beta))] <- Beta
//' XB <- X%*%beta
//' probs <- as.vector(exp(XB)/(1+exp(XB)))
//' y <- rbinom(n,1,probs)
//' exmat <- cbind(y,X)
//' tau <- 0.5; r <- 1; a <- 3; b <- p-a; in_cons <- a;
//' loopcnt <- 100; cplng <- FALSE;
//' initProb <- in_cons/p
//'
//' ### Initializing Chains
//' schain <- p
//' while (schain > in_cons || schain == 0) {
//' chain1 <- rbinom(p, 1, initProb)
//'  schain <- sum(chain1)
//' }
//' chain1 <- as.numeric(c(1, chain1))
//' chain2 <- chain1
//' nlptype <- 0 ## PiMOM nonlocal prior
//' nf <- 0 ### No fixed columns
//'
//' ### Running the function
//' bvsout <- logreg_bvs(exmat,chain1,nf,tau,r,nlptype,a,b,in_cons,loopcnt,cplng,chain2)
//'
//' ### Number of visited models for this specific run:
//' bvsout$num_vis_models
//'
//' ### The selected model:
//' which(bvsout$max_chain > 0)
//'
//' ### Estimated coefficients:
//' bvsout$beta_hat
//'
//' ### The unnormalized probability of the selected model:
//' bvsout$max_prob
// [[Rcpp::export]]
Rcpp::List logreg_bvs(const arma::mat& exmat, arma::uvec chain1, int nf, double tau, double r, const int nlptype, int a, int b,
                      int in_cons, int loopcnt, bool cplng, arma::uvec chain2)
{
  arma::mat XX = exmat;
  arma::vec y = XX.col(0);
  XX.shed_col(0);
  int n = XX.n_rows;
  int p = XX.n_cols - 1;
  arma::uvec indices = seq_gen(p)+1;
  indices = indices.tail(p-nf);
  const MapVec Y(y.memptr(), n);
  arma::mat cur_x1, cur_x2, final_x;
  arma::uvec cur_ind, cand1, cand2, maxchain, cur_cols1, cur_cols2;
  int cons = in_cons+1;
  bool getout = false;
  int i, elsum1, elsum2, lcnt = 0, cpt, k1, k2;
  bool cflag = false;
  double logprob1, logprob2, logcand1, logcand2, maxprob=-1000, u1, u2;
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  
  Rcpp::List hash_vis_covs(10000000);
  arma::vec hash_key = arma::zeros<arma::vec>(10000000);
  arma::vec hash_prob = arma::zeros<arma::vec>(10000000);
  arma::vec hash_prob_sub, hash_key_sub;
  Rcpp::List vis_covs;
  
  int count = 0, k;
  double thiskey;
  
  cur_cols1 = arma::find(chain1);
  thiskey = calc_key(cur_cols1);
  cur_x1 = XX.cols(cur_cols1);
  k1 = cur_x1.n_cols;
  const MapMat X1(cur_x1.memptr(), n, k1);
  logprob1 = LogReg_Model_Prob(X1,Y,tau,r,a,b,n,p,nlptype);
  
  hash_key(count) = thiskey;
  hash_prob(count) = logprob1;
  hash_vis_covs[count] = cur_cols1;
  count++;
  
  if (cplng){
    hash_key(0) = 0;
    hash_prob(0) = 0;
    hash_vis_covs[0] = 0;
    cur_cols2 = arma::find(chain2);
    cur_x2 = XX.cols(cur_cols2);
    k2 = cur_x2.n_cols;
    const MapMat X2(cur_x2.memptr(), n, k2);
    logprob2 = LogReg_Model_Prob(X2,Y,tau,r,a,b,n,p,nlptype);
    count++;
    
    maxprob = std::max(logprob1,logprob2);
    
    while(!getout && lcnt < loopcnt){
      cur_ind = indices;
      cur_ind = arma::shuffle(cur_ind);
      for (int j=0; j < p-nf; j++){
        cand1 = chain1;
        cand2 = chain2;
        i = cur_ind(j);
        (cand1(i)==0)?cand1(i)=1:cand1(i)=0;
        (cand2(i)==0)?cand2(i)=1:cand2(i)=0;
        elsum1 = arma::sum(cand1);
        elsum2 = arma::sum(cand2);
        if (elsum1 > cons) logcand1 = -9999;
        else{
          cur_cols1 = arma::find(cand1);
          cur_x1 = XX.cols(cur_cols1);
          const MapMat X1(cur_x1.memptr(), n, elsum1);
          logcand1 = LogReg_Model_Prob(X1,Y,tau,r,a,b,n,p,nlptype);
          count++;
        }
        if (elsum2 > cons) logcand2 = -9999;
        else{
          cur_cols2 = arma::find(cand2);
          cur_x2 = XX.cols(cur_cols2);
          const MapMat X2(cur_x2.memptr(), n, elsum2);
          logcand2 = LogReg_Model_Prob(X2,Y,tau,r,a,b,n,p,nlptype);
          count++;
        }
        
        u1 = unif_rand();
        (chain1(i)==chain2(i)) ? u2=u1: u2=1-u1;
        
        if(logcand1-logprob1 > log(u1)){
          chain1 = cand1;
          logprob1 = logcand1;
        }
        if(logcand2-logprob2 > log(u2)){
          chain2 = cand2;
          logprob2 = logcand2;
        }
        if( logprob1 > maxprob ){
          maxprob = logprob1;
          maxchain = chain1;
        }
        if( logprob2 > maxprob ){
          maxprob = logprob2;
          maxchain = chain2;
        }
      }
      if( isequal(chain1, chain2) ){
        cflag = true;
        getout = true;
      }
      lcnt++ ;
    }
    (cflag)?cpt=lcnt:cpt=loopcnt;
    
  } else {
    maxprob = logprob1;
    while(lcnt < loopcnt){
      cur_ind = indices;
      cur_ind = arma::shuffle(cur_ind);
      for (int j=0; j < p-nf; j++){
        cand1 = chain1;
        i = cur_ind(j);
        (cand1(i)==0)?cand1(i)=1:cand1(i)=0;
        elsum1 = arma::sum(cand1);
        if (elsum1 > cons) logcand1 = -999999;
        else{
          cur_cols1 = arma::find(cand1);
          thiskey = calc_key(cur_cols1);
          cur_x1 = XX.cols(cur_cols1);
          const MapMat X1(cur_x1.memptr(), n, elsum1);
          logcand1 = LogReg_Model_Prob(X1,Y,tau,r,a,b,n,p,nlptype);
          
          hash_key(count) = thiskey;
          hash_prob(count) = logcand1;
          hash_vis_covs[count] = cur_cols1;
          count++;
        }
        u1 = unif_rand();
        if(logcand1-logprob1 > log(u1)){
          chain1 = cand1;
          logprob1 = logcand1;
        }
        if( logprob1 > maxprob ){
          maxprob = logprob1;
          maxchain = chain1;
        }
      }
      lcnt++ ;
    }
    cpt=lcnt;
    
    hash_prob_sub = hash_prob.subvec(0,count-1);
    hash_key_sub = hash_key.subvec(0,count-1);
    Rcpp::IntegerVector Lcov_idx(count);
    std::iota(Lcov_idx.begin(), Lcov_idx.end(), 0); vis_covs = hash_vis_covs[Lcov_idx];
  }
  
  cur_cols1 = arma::find(maxchain);
  final_x = XX.cols(cur_cols1);
  k = final_x.n_cols;
  const MapMat Fxx(final_x.memptr(), n, k);
  Eigen::VectorXd glm_coef(k);
  glm_coef.setZero();
  
  LogReg_MLE nll(Fxx,Y);
  double fopt_glm;
  int glm_status = Numer::optim_lbfgs(nll, glm_coef, fopt_glm, maxit, eps_f, eps_g);
  if (glm_status < 0) Rcpp::stop("surprisingly the optimization for finding coefficients of the selected model did not converge!!");
  Eigen::VectorXd betahat = glm_coef;
  LogReg_MAP FM(Fxx,Y,tau,r,nlptype);
  double fopt;
  int status = Numer::optim_lbfgs(FM, betahat, fopt, maxit, eps_f, eps_g);
  if (status < 0) Rcpp::stop("surprisingly the optimization for finding coefficients of the selected model did not converge!!");
  return(Rcpp::List::create(Rcpp::Named("max_chain")=maxchain, Rcpp::Named("beta_hat")=betahat,
                            Rcpp::Named("num_iterations")=cpt, Rcpp::Named("max_prob")=maxprob, Rcpp::Named("cplng_flag")=cflag,
                            Rcpp::Named("num_vis_models")=count, Rcpp::Named("hash_key")=hash_key_sub, Rcpp::Named("hash_prob")=hash_prob_sub,
                            Rcpp::Named("vis_covs")=vis_covs));
}
//==================================================//

// [[Rcpp::export]]
Rcpp::NumericVector lreg_coef_est(const arma::mat& exmat, arma::uvec mod_cols, double tau, double r, const int nlptype)
{
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  arma::mat XX = exmat;
  arma::vec y = XX.col(0);
  XX.shed_col(0);
  int n = XX.n_rows;
  arma::uvec cur_cols = push_begin(mod_cols,0);
  arma::mat fx = XX.cols(cur_cols);
  int k = fx.n_cols;
  MapMat Fxx(fx.memptr(), n, k);
  const MapVec Y(y.memptr(), n);
  
  Eigen::VectorXd glm_coef(k);
  glm_coef.setZero();
  LogReg_MLE nll(Fxx,Y);
  double fopt_glm;
  int glm_status = Numer::optim_lbfgs(nll, glm_coef, fopt_glm, maxit, eps_f, eps_g);
  if (glm_status < 0) {Rcpp::stop("The optimization function to estimate coefficients did not converge!");}
  Eigen::VectorXd beta = glm_coef;
  
  LogReg_MAP CM(Fxx,Y,tau,r,nlptype);
  double fopt;
  int status = Numer::optim_lbfgs(CM, beta, fopt, maxit, eps_f, eps_g);
  if(status < 0) {Rcpp::stop("The optimization function to estimate coefficients did not converge!");}
  Eigen::VectorXd betavec = beta;
  return(Rcpp::wrap(betavec));
}
//==================================================//

// [[Rcpp::export]]
double lreg_mod_prob(const arma::mat& exmat, arma::uvec mod_cols, double tau, double r, int a, int b, const int nlptype)
{
  arma::mat XX = exmat;
  arma::vec y = XX.col(0);
  XX.shed_col(0);
  int n = XX.n_rows;
  int p = XX.n_cols;
  arma::uvec cur_cols = mod_cols - 1;
  arma::mat fx = XX.cols(cur_cols);
  int k = fx.n_cols;
  MapMat Fxx(fx.memptr(), n, k);
  const MapVec Y(y.memptr(), n);
  double modprob = LogReg_Model_Prob(Fxx,Y,tau,r,a,b,n,p,nlptype);
  return(modprob);
}
//==================================================//

// [[Rcpp::export]]
Rcpp::NumericVector null_mle_cox(arma::mat& XX, int n, int p, int cons,
                                 double a, double b, double csr, int niters)
{
  int ncol,k;
  double probs = 1/p;
  double lambda = 1;
  double lambda2 = csr*lambda/(1-csr);
  arma::uvec sel_cols;
  arma::uvec indices = seq_gen(p);
  Eigen::VectorXd out_null(2000000), betahat;
  int counts = 0;
  double eps_f = 1.0e-08;
  double eps_g = 1.0e-05;
  int maxit = 300;
  Rcpp::NumericVector rbinout, prob_msize;
  
  arma::vec exb; exb.ones(n);
  Rcpp::NumericVector rvec = Rcpp::runif(n,0,1);
  arma::vec ruvec = Rcpp::as<arma::vec>(rvec);
  arma::vec times3 = -log(ruvec)/(lambda*exb);
  Rcpp::NumericVector times_nv = Rcpp::NumericVector(times3.begin(),times3.end());
  Rcpp::NumericVector cst_nv = Rcpp::rexp(n,lambda2);
  arma::vec cstvec = Rcpp::as<arma::vec>(cst_nv);
  Rcpp::NumericVector times2 = Rcpp::pmin(times_nv,cst_nv);
  arma::vec times = Rcpp::as<arma::vec>(times2);
  arma::uvec status2 = (times3 <= cstvec);
  arma::vec status_org = arma::conv_to<arma::vec>::from(status2);
  
  arma::vec cur_xbeta = arma::zeros<arma::vec>(n);
  const MapVec Cur_XBeta(cur_xbeta.memptr(), n);
  arma::mat TS = arma::join_rows(times,status_org);
  arma::mat exmat = arma::join_rows(TS,XX);
  arma::mat ord_temp = cox_order_vecs(exmat);
  arma::vec status = ord_temp.col(1);
  const MapVec Status(status.memptr(), n);
  ord_temp.shed_cols(0,1);
  const arma::mat ordexmat = ord_temp;
  
  for (int i=0; i < niters; i++){
    ncol = 0;
    while (ncol== 0 || ncol > cons){
      prob_msize = Rcpp::rbeta(1,a,b);
      rbinout = Rcpp::rbinom(1,p,prob_msize[0]);
      ncol = rbinout[0];
    }
    k = ncol;
    sel_cols = Rcpp::RcppArmadillo::sample(indices,ncol,false,probs);
    arma::mat cur_model = ordexmat.cols(sel_cols);
    MapMat fxx(cur_model.memptr(), n, k);
    Eigen::VectorXd coxph_coef(k);
    coxph_coef.setZero();
    Cox_MLE nll(fxx,Status,Cur_XBeta);
    double fopt_glm;
    int mle_status = Numer::optim_lbfgs(nll, coxph_coef, fopt_glm, maxit, eps_f, eps_g);
    if (mle_status >= 0){
      int t = (coxph_coef.array().abs() > 5).sum();
      if (t==0){
        for (int j=0; j<k; j++) out_null(counts+j)=coxph_coef(j);
        counts = counts+k;
      }
    }
  }
  betahat=out_null.head(counts-1);
  return(Rcpp::wrap(betahat));
}
//==================================================//

//' Non-parallel version of Bayesian variable selector for survival data using
//' nonlocal priors
//' @description This function performs Bayesian variable selection for
//' survival data in a non-parallel fashion. It runs modified S5 algorithm to
//' search the model space but since this is only on one CPU, the number of
//' visited models will not be large and therefore is NOT recommended for
//' high dimensional datasets. This function is called by \code{\link{bvs}}
//' function in a parllel fashion and therefore that function is recommended
//' to be used.
//'
//' @param exmat An extended matrix where the first two columns are survival
//' times and status, respectively and the rest is the design matrix which is
//' produced by \code{PreProcess} function.
//' @param cur_cols A vector containing indices of the initial model for
//' variable selector to start the S5 algorithm from. Note that the first
//' \code{nf} indices are \code{1} to \code{nf} where \code{nf} is the number
//' of fixed covariates that do not enter the selection procedure.
//' @param nf The number of fixed covariates that do not enter the selection
//' procedure.
//' @param tau The paramter \code{tau} of the iMOM prior.
//' @param r The paramter \code{r} of the iMOM prior.
//' @param nlptype Determines the type of nonlocal prior that is used in the
//' analyses. \code{0} is for piMOM and \code{1} is for pMOM. 
//' @param a The first parameter in beta distribution used as prior on model
//' size. This parameter is equal to 1 when uinform-binomial prior is used.
//' @param b The second paramter in beta distribution used as prior on model
//' size. This parameter is equal to 1 when uinform-binomial prior is used.
//' @param d This is the number of candidate covariates picked from top
//' variables with highest utility function value and used in S5 algorithm.
//' @param L Number of temperatures in S5 algorithm.
//' @param J Number of iterations at each temperature in S5 algorithm.
//' @param temps Vector of temperatuers used in S5 algorithm.
//'
//' @return It returns a list containing following objects:
//' \item{max_model}{A \code{1} by \code{p} binary vector showing the selected
//' model with maximum probability. \code{1} means a specific variable is
//' selected.}
//' \item{hash_key}{A column vector indicating the generated key for each model
//' that is used to track visited models and growing dictionary.}
//' \item{max_prob}{The unnormalized probability of the model with highest
//' posterior probability.}
//' \item{all_probs}{A vector containing unnormalized probabilities of all
//' visited models.}
//' \item{vis_covs_list}{A list containing the covariates in each visited model
//' in the stochastic search process.}
//' @author Amir Nikooienejad
//' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2017). Bayesian
//' Variable Selection in High Dimensional Survival Time Cancer Genomic
//' Datasets using Nonlocal Priors. arXiv preprint, arXiv:1712.02964.\cr\cr
//' Shin, M., Bhattacharya, A., and Johnson, V. E. (2017). Scalable
//' Bayesian variable selection using nonlocal prior densities in ultrahigh
//' dimensional settings. Statistica Sinica.\cr\cr
//' Johnson, V. E., and Rossell, D. (2010). On the use of non-local prior
//' densities in Bayesian hypothesis tests. Journal of the Royal Statistical
//' Society: Series B (Statistical Methodology), 72(2), 143-170.
//' @seealso \code{\link{bvs}}
//' @examples
//' ### Initializing the parameters
//' n <- 100
//' p <- 40
//' set.seed(123)
//' Sigma <- diag(p)
//' full <- matrix(c(rep(0.5, p*p)), ncol=p)
//' Sigma <- full + 0.5*Sigma
//' cholS <- chol(Sigma)
//' Beta <- c(-1.8, 1.2, -1.7, 1.4, -1.4, 1.3)
//' X = matrix(rnorm(n*p), ncol=p)
//' X = X%*%cholS
//' X <- scale(X)
//' beta <- numeric(p)
//' beta[c(1:length(Beta))] <- Beta
//' XB <- X%*%beta
//' sur_times <- rexp(n,exp(XB))
//' cens_times <- rexp(n,0.2)
//' times <- pmin(sur_times,cens_times)
//' status <- as.numeric(sur_times <= cens_times)
//' exmat <- cbind(times,status,X)
//' L <- 10; J <- 10
//' d <- 2 * ceiling(log(p))
//' temps <- seq(3, 1, length.out = L)
//' tau <- 0.5; r <- 1; a <- 6; b <- p-a
//' nlptype <- 0 ### PiMOM nonlocal prior
//' cur_cols <- c(1,2,3) ### Starting model for the search algorithm
//' nf <- 0 ### No fixed columns
//' 
//'### Running the Function
//' coxout <- cox_bvs(exmat,cur_cols,nf,tau,r,nlptype,a,b,d,L,J,temps)
//' 
//' ### The number of visited model for this specific run:
//' length(coxout$hash_key)
//' 
//'
//' ### The selected model:
//' which(coxout$max_model>0)
//'
//' ### The unnormalized probability of the selected model:
//' coxout$max_prob
//' 
// [[Rcpp::export]]
Rcpp::List cox_bvs(const arma::mat& exmat, arma::uvec cur_cols, int nf, double tau, double r, const int nlptype, int a,
                   int b, int d, int L, int J, arma::vec temps)
{
  int n = exmat.n_rows;
  arma::mat ord_temp = cox_order_vecs(exmat);
  arma::mat fdes, cur_model, kc_model, tmp_des, mcout;
  arma::uvec idx, kc_cols, max_ind, cand_add_cols, new_model_cols, idx_add, idx_neg, idx_fin = seq_gen(2), add_cand, neg_cand, final_model;
  arma::uvec tsel_add, tsel_neg, fsel, tmp_cur_cols;
  arma::vec add_probs, neg_probs, add_s_probs, neg_s_probs, f_prob(2);
  arma::vec status = ord_temp.col(1), tmp_add_s_probs, tmp_neg_s_probs;
  arma::vec betahat;
  arma::urowvec tmp_chain, max_chain;
  arma::vec cxbeta; Rcpp::NumericVector bhat;
  Rcpp::List vis_covs;
  ord_temp.shed_cols(0,1);
  
  const MapVec Status(status.memptr(), n);
  double thiskey, thiskey_cand, tmp_prob, max_prob;
  const arma::mat ordexmat = ord_temp;
  int p = ordexmat.n_cols;
  int cur_cols_size = cur_cols.size();
  cur_cols = cur_cols-1;
  idx = seq_gen(p); idx = idx.tail(p-nf);
  
  arma::mat xx;
  arma::vec hash_key = arma::zeros<arma::vec>(100000);
  arma::vec hash_key_cand = arma::zeros<arma::vec>(100000);
  Rcpp::List hash_vis_covs(100000);
  Rcpp::List hash_cand_covs(100000);
  arma::vec hash_prob = arma::zeros<arma::vec>(100000);
  arma::vec hash_prob_sub, hash_key_sub;
  int count = 0, k, add_size, neg_size, k1, d1;
  int count_cand = 0;
  int  shifter, madd, mneg, foundkey, foundkey_cand;
  bool no_add_flag = false;
  
  thiskey = calc_key(cur_cols);
  xx = ordexmat.cols(cur_cols);
  tmp_chain = arma::zeros<arma::urowvec>(p);
  tmp_chain(cur_cols) = arma::ones<arma::urowvec>(cur_cols_size);
  max_chain = tmp_chain;
  bhat = cox_beta_est(xx, Status,tau,r,nlptype);
  if (bhat[0]==-999999){
    count++;
    hash_key(0) = thiskey; hash_prob(0) = -999999999; hash_vis_covs[0] = arma::find(tmp_chain);
    hash_prob_sub = hash_prob.subvec(0,count-1); hash_key_sub = hash_key.subvec(0,count-1);
    Rcpp::IntegerVector Lcov_idx(count);
    std::iota(Lcov_idx.begin(), Lcov_idx.end(), 0); vis_covs = hash_vis_covs[Lcov_idx];
    
    return(Rcpp::List::create(Rcpp::Named("max_model")=max_chain, Rcpp::Named("hash_key")=hash_key_sub,
                              Rcpp::Named("max_prob")= -999999999, Rcpp::Named("all_probs")=hash_prob_sub,
                              Rcpp::Named("vis_covs_list")=vis_covs));
  }
  betahat = Rcpp::as<arma::vec>(bhat);
  MapMat XX(xx.memptr(), n, cur_cols_size);
  tmp_prob = Cox_Model_Prob(XX,Status,betahat,tau,r,a,b,p,nlptype);
  max_prob = tmp_prob;
  hash_key(count) = thiskey;
  hash_prob(count) = tmp_prob;
  hash_vis_covs[count] = arma::find(tmp_chain);
  count++;
  
  for (int l=0; l < L; l++){
    for (int j=0; j < J; j++){
      k = cur_cols.size();
      thiskey_cand = calc_key(cur_cols);
      foundkey_cand = my_find(hash_key_cand,thiskey_cand);
      if(foundkey_cand == -1){
        hash_key_cand(count_cand) = thiskey_cand;
        cur_model = ordexmat.cols(cur_cols);
        bhat = cox_beta_est(cur_model, Status,tau,r,nlptype);
        betahat = Rcpp::as<arma::vec>(bhat);
        cxbeta = cur_model*betahat;
        const MapVec Cur_XBeta(cxbeta.memptr(),n);
        tmp_cur_cols = cur_cols.tail(k-nf);
        kc_cols = c_setdiff(idx,tmp_cur_cols);
        k1 = kc_cols.size();
        d1 = std::min(k1,d);
        max_ind = find_max_utils(ordexmat,kc_cols,Status,Cur_XBeta,d1);
        cand_add_cols = kc_cols(max_ind);
        hash_cand_covs[count_cand] = cand_add_cols;
        count_cand++;
      } else {
        cand_add_cols = Rcpp::as<arma::uvec>(hash_cand_covs[foundkey_cand]);
      }
      
      // //=========== From Here.....
      // cur_model = ordexmat.cols(cur_cols);
      // bhat = cox_beta_est(cur_model, Status,tau,r,nlptype);
      // betahat = Rcpp::as<arma::vec>(bhat);
      // cxbeta = cur_model*betahat;
      // const MapVec Cur_XBeta(cxbeta.memptr(),n);
      // tmp_cur_cols = cur_cols.tail(k-nf);
      // kc_cols = c_setdiff(idx,tmp_cur_cols);
      // max_ind = find_max_utils(ordexmat,kc_cols,Status,Cur_XBeta,d);
      // cand_add_cols = kc_cols(max_ind);
      // //============ Up to here....
      
      add_size = cand_add_cols.size();
      neg_size = k-nf;
      idx_add = seq_gen(add_size);
      idx_neg = seq_gen(neg_size);
      add_cand = cur_cols; add_cand.resize(k+1);
      add_probs.set_size(add_size), neg_probs.set_size(neg_size);
      
      if (add_size == 0){
        add_probs.set_size(1);
        add_probs(0) = -999999999;
        no_add_flag=true;
      } else {
        for (int i=0; i < add_size; i++){
          add_cand(k) = cand_add_cols(i);
          thiskey = calc_key(add_cand);
          foundkey = my_find(hash_key,thiskey);
          if(foundkey == -1){
            hash_key(count) = thiskey;
            xx = ordexmat.cols(add_cand);
            bhat = cox_beta_est(xx, Status,tau,r,nlptype);
            if (bhat[0]==-999999){
              tmp_prob = -999999999;
            } else {
              betahat = Rcpp::as<arma::vec>(bhat);
              MapMat XX(xx.memptr(), n, k+1);
              tmp_prob = Cox_Model_Prob(XX,Status,betahat,tau,r,a,b,p,nlptype);
            }
            hash_prob(count) = tmp_prob;
            tmp_chain = arma::zeros<arma::urowvec>(p);
            tmp_chain(add_cand) = arma::ones<arma::urowvec>(k+1);
            hash_vis_covs[count] = arma::find(tmp_chain);
            
            if (tmp_prob > max_prob ){
              max_prob = tmp_prob;
              max_chain = tmp_chain;
            }
            count++;
          }else{
            tmp_prob = hash_prob(foundkey);
          }
          add_probs(i) = tmp_prob;
        }
      }
      
      if (neg_size > 1) {
        for (int i=0; i < neg_size; i++){
          neg_cand = cur_cols;
          neg_cand.shed_row(nf+i);
          thiskey = calc_key(neg_cand);
          foundkey = my_find(hash_key,thiskey);
          if(foundkey == -1){
            hash_key(count) = thiskey;
            xx = ordexmat.cols(neg_cand);
            bhat = cox_beta_est(xx, Status,tau,r,nlptype);
            
            if (bhat[0]==-999999){
              tmp_prob = -999999999;
            } else {
              betahat = Rcpp::as<arma::vec>(bhat);
              MapMat XX(xx.memptr(), n, k-1);
              tmp_prob =  Cox_Model_Prob(XX,Status,betahat,tau,r,a,b,p,nlptype);
            }
            hash_prob(count) = tmp_prob;
            tmp_chain = arma::zeros<arma::urowvec>(p);
            tmp_chain(neg_cand) = arma::ones<arma::urowvec>(k-1);
            hash_vis_covs[count] = arma::find(tmp_chain);
            if (tmp_prob > max_prob ){
              max_prob = tmp_prob;
              max_chain = tmp_chain;
            }
            count++;
          }else{
            tmp_prob = hash_prob(foundkey);
          }
          neg_probs(i) = tmp_prob;
        }
        madd = ceil(max(add_probs)); mneg = ceil(max(neg_probs));
        (madd >= mneg) ? shifter = madd : shifter = mneg;
        tmp_neg_s_probs = pow(exp(neg_probs - shifter),1/temps(l));
        if (sum(tmp_neg_s_probs)==0){
          f_prob(1) = 0;
        }else{
          neg_s_probs = tmp_neg_s_probs/sum(tmp_neg_s_probs);
          tsel_neg = cox_c_sample(idx_neg,1,neg_s_probs,false);
          f_prob(1) = tmp_neg_s_probs(tsel_neg(0));
        }
      }else{
        shifter = ceil(max(add_probs));
        f_prob(1) = 0;
      }
      
      tmp_add_s_probs = pow(exp(add_probs - shifter),1/temps(l));
      if(sum(tmp_add_s_probs) == 0 || no_add_flag){
        f_prob(0)=0;
        no_add_flag = false;
      }else{
        add_s_probs = tmp_add_s_probs/sum(tmp_add_s_probs);
        tsel_add = cox_c_sample(idx_add,1,add_s_probs,false);
        f_prob(0) = tmp_add_s_probs(tsel_add(0));
      }
      
      if (sum(f_prob)==0){
        f_prob(0)=0.5; f_prob(1)=0.5;
      }else{
        f_prob = f_prob/sum(f_prob);
      }
      
      fsel = cox_c_sample(idx_fin,1,f_prob,false);
      switch (fsel(0)){
      case 0:
        new_model_cols = add_cand; new_model_cols(k) = cand_add_cols(tsel_add(0));
        break;
      case 1:
        new_model_cols = cur_cols;
        if (k>1) new_model_cols.shed_row(tsel_neg(0)+nf);
        break;
      default:
        Rcpp::stop("Non-acceptable value for the next step in stochastic search algorithm");
      }
      cur_cols = new_model_cols;
    }
  }
  hash_prob_sub = hash_prob.subvec(0,count-1);
  hash_key_sub = hash_key.subvec(0,count-1);
  Rcpp::IntegerVector Lcov_idx(count);
  std::iota(Lcov_idx.begin(), Lcov_idx.end(), 0); vis_covs = hash_vis_covs[Lcov_idx];
  
  return(Rcpp::List::create(Rcpp::Named("max_model")=max_chain, Rcpp::Named("hash_key")=hash_key_sub,
                            Rcpp::Named("max_prob")=max_prob, Rcpp::Named("all_probs")=hash_prob_sub,
                            Rcpp::Named("vis_covs_list")=vis_covs));
}
//==================================================//

// [[Rcpp::export]]
arma::vec inc_prob_calc(arma::vec all_probs, Rcpp::List vis_covs, int p)
{
  int n = vis_covs.size();
  arma::vec inc_probs = arma::zeros<arma::vec>(p);
  double shifter = ceil(max(all_probs));
  arma::vec all_new_probs = exp(all_probs - shifter);
  arma::vec norm_probs = all_new_probs/sum(all_new_probs);
  arma::uvec mod_inds;
  arma::uvec out_inds;
  
  for (int i=0; i < n; i++){
    Rcpp::IntegerVector tmp_ind(vis_covs[i]);
    mod_inds = Rcpp::as<arma::uvec>(tmp_ind);
    inc_probs.elem(mod_inds) += norm_probs(i);
  }
  return(inc_probs);
}
//==================================================//

// [[Rcpp::export]]
Rcpp::NumericVector cox_coef_est(const arma::mat& exmat, arma::uvec mod_cols, double tau, double r, const int nlptype)
{
  int n = exmat.n_rows;
  arma::mat ord_temp = cox_order_vecs(exmat);
  arma::vec status = ord_temp.col(1);
  const MapVec Status(status.memptr(), n);
  ord_temp.shed_cols(0,1);
  const arma::mat ordexmat = ord_temp;
  mod_cols = mod_cols - 1;
  arma::mat cur_model = ordexmat.cols(mod_cols);
  
  Rcpp::NumericVector bhat = cox_beta_est(cur_model, Status,tau,r,nlptype);
  if (bhat[0]==-999999) {Rcpp::stop("The optimization function to estimate coefficients did not converge!");}
  return(bhat);
}
//==================================================//

// [[Rcpp::export]]
double cox_mod_prob(const arma::mat& exmat, arma::uvec mod_cols, double tau, double r, int a, int b, const int nlptype)
{
  double prob;
  mod_cols = mod_cols-1;
  int n = exmat.n_rows;
  int k = mod_cols.size();
  arma::mat ord_temp = cox_order_vecs(exmat);
  arma::vec status = ord_temp.col(1);
  ord_temp.shed_cols(0,1);
  const arma::mat ordexmat = ord_temp;
  int p = ordexmat.n_cols;
  const MapVec Status(status.memptr(), n);
  arma::mat cur_model = ordexmat.cols(mod_cols);
  MapMat XX(cur_model.memptr(), n, k);
  Rcpp::NumericVector bhat = cox_beta_est(cur_model, Status,tau,r,nlptype);
  if (bhat[0]==-999999) {Rcpp::stop("The optimization function to estimate coefficients did not converge!");}
  arma::vec betahat = Rcpp::as<arma::vec>(bhat);
  prob = Cox_Model_Prob(XX,Status,betahat,tau,r,a,b,p,nlptype);
  return(prob);
}

//================== Part 4: BMA Predictive Accuracy Measurements =========================
//=========================================================================================

arma::mat sort_TS(const arma::mat& TS) {
  arma::vec times = TS.col(0);
  arma::vec status = TS.col(1);
  int n = TS.n_rows;
  arma::mat outTS(n,2);
  arma::uvec idx;
  idx.zeros(n);
  auto begin = std::begin(idx), end = std::end(idx);
  std::iota(begin, end, static_cast<size_t>(0));
  auto comparator = [&](const size_t & a, const size_t & b){ return times[a] < times[b]; };
  std::sort(begin, end, comparator);
  for (int i=0; i < n; i++){
    outTS.row(i) = TS.row(idx[i]);
  }
  return(outTS);
}
//==================================================//

arma::mat KMestimate(const arma::mat& TS){
  arma::vec times = TS.col(0);
  arma::vec status = 1-TS.col(1);
  int n = TS.n_rows;
  arma::mat kmes(n,2);
  int d, ar;
  
  double current=1.0;
  for(int i=0; i < n; i++){
    d=0; ar=0;
    for(int j=0; j < n; j++){
      ar += (times(i) <= times(j));
      d += (times(i) == times(j)) && (status(i));
    }
    current = current * (1.0 - (double) d / (double) ar);
    kmes(i,0) = current;
  }
  kmes.col(1)=times;
  return(kmes);
}

//==================================================//

arma::vec g_est(const arma::mat& TS_te, const arma::mat& kmes_tr)
{
  arma::vec ntime = TS_te.col(0);
  arma::vec pretime = kmes_tr.col(1);
  arma::uvec idx;
  int maxind;
  int n_new = TS_te.n_rows;
  arma::vec gout(n_new);
  
  for (int i=0; i < n_new; i++){
    idx = arma::find(pretime <= ntime(i));
    if (idx.n_elem){
      maxind = max(idx);
      gout(i) = kmes_tr(maxind,0);
    } else {
      gout(i) = 1;
    }
  }
  return(gout);
}
//==================================================//

arma::vec calc_marker(arma::mat xcols_tr, arma::mat xcols_te, arma::vec coefs)
{
  arma::vec marker;
  arma::mat mtr = mean(xcols_tr);
  double m2 = arma::dot(mtr,coefs);
  arma::mat m1 = xcols_te * coefs;
  marker = m1-m2;
  return(marker);
}
//==================================================//

// [[Rcpp::export]]
Rcpp::List aucBMA_logistic(const arma::mat& X_tr, const arma::vec& y_tr, const arma::mat& X_te, const arma::vec& y_te,
                           double tau, double r, const int nlptype, arma::vec probs, Rcpp::ListOf<Rcpp::IntegerVector> models, int k)
{
  arma::uvec mod_cols, mod1;
  arma::vec cfs, mark, mark1, mark2, mark3;
  int n_new = y_te.n_elem;
  Rcpp::NumericVector tmp_cfs;
  arma::mat tmp_Markers(n_new,k), xcols_te, m1;
  int n_tr = y_tr.n_elem; int p = X_tr.n_cols;
  arma::mat exmat_tr(n_tr,p+1); exmat_tr.col(0) = y_tr; exmat_tr.cols(1,p) = X_tr;
  
  for (int i=0; i < k; i++){
    mod1 = Rcpp::as<arma::uvec>(models[i]); mod1 = mod1 - 1;
    mod_cols = rm_begin(mod1);
    tmp_cfs = lreg_coef_est(exmat_tr, mod_cols, tau, r,nlptype);
    cfs = Rcpp::as<arma::vec>(tmp_cfs);
    xcols_te = X_te.cols(mod1);
    m1 = xcols_te*cfs; mark = m1.col(0);
    mark1 = arma::exp(mark); mark2 = mark1+1; mark3 = mark1/mark2;
    tmp_Markers.col(i) = mark3;
  }
  arma::mat Markers = tmp_Markers;
  
  arma::vec thresh = arma::unique(Markers.col(0));
  int n_th = thresh.n_elem;
  
  int P = sum(y_te); int N = n_new-P;
  arma::vec tmp_tp, tmp_fp;
  arma::vec tpr(n_th), fpr(n_th);
  
  for(int j=0; j < n_th; j++){
    tmp_tp.zeros(k); tmp_fp.zeros(k);
    for(int i=0; i< n_new; i++){
      for (int r=0; r < k; r++){
        if (Markers(i,r) > thresh(j)){
          tmp_tp(r) += y_te(i)*probs(r);
          tmp_fp(r) += (1-y_te(i))*probs(r);
        }
      }
    }
    tpr(j) = sum(tmp_tp)/P; fpr(j) = sum(tmp_fp)/N;
  }
  
  arma::mat roc_mat(n_th+2,2);
  roc_mat(0,0) = roc_mat(0,1) = 1; 
  roc_mat(n_th+1,0) = roc_mat(n_th+1,1) = 0; 
  roc_mat(arma::span(1,n_th),0) = fpr;  roc_mat(arma::span(1,n_th),1) = tpr;
  
  double base, height, auc = 0;
  for (int i=0; i < n_th+1; i++){
    base = roc_mat(i,0) - roc_mat(i+1,0);
    height = 0.5*(roc_mat(i,1) + roc_mat(i+1,1));
    auc += base * height;
  }
  return(Rcpp::List::create(Rcpp::Named("auc")=auc, Rcpp::Named("roc")=roc_mat));
}
//==================================================//

// [[Rcpp::export]]
arma::vec aucBMA_survival(const arma::mat& X_tr, const arma::mat& TS_tr, const arma::mat& X_te, const arma::mat& TS_te,
                          double tau, double r, const int nlptype, arma::vec times, arma::vec probs,
                          Rcpp::ListOf<Rcpp::IntegerVector> models, int k)
{
  arma::uvec mod_cols, mod1;
  int n_tr = TS_tr.n_rows; int p = X_tr.n_cols;
  arma::mat exmat_tr(n_tr,p+2); exmat_tr.cols(0,1) = TS_tr; exmat_tr.cols(2,p+1) = X_tr;
  arma::mat sTS_tr = sort_TS(TS_tr);
  arma::vec cfs, mark, mark1, mark2, mark3;
  
  Rcpp::NumericVector tmp_cfs;
  int n_t = times.n_elem;
  int n_new = TS_te.n_rows;
  arma::vec auc; auc.zeros(n_t);
  
  arma::mat Markers(n_new,k), xcols_tr, xcols_te;
  for (int i=0; i < k; i++){
    mod1 = Rcpp::as<arma::uvec>(models[i]); mod_cols = mod1; mod1 = mod1 - 1;
    tmp_cfs = cox_coef_est(exmat_tr, mod_cols, tau, r, nlptype);
    cfs = Rcpp::as<arma::vec>(tmp_cfs);
    xcols_tr = X_tr.cols(mod1); xcols_te = X_te.cols(mod1);
    mark = calc_marker(xcols_tr,xcols_te,cfs);
    mark1 = arma::exp(mark); mark2 = mark1+1; mark3 = mark1/mark2;
    Markers.col(i) = mark3; 
  }
  
  arma::vec thresh = arma::sort(arma::unique(Markers.col(0)));
  int n_th = thresh.n_elem;  
  
  arma::mat kmes = KMestimate(sTS_tr);
  arma::vec G = g_est(TS_te, kmes);
  double senumer, sedenom, spnumer, spdenom, se_bma_cur, sp_bma_cur, tmp_var=0.0;
  arma::vec tmp_spec(k); arma::vec tmp_sens(k);
  arma::vec surv_new = TS_te.col(0);
  arma::vec st_new = TS_te.col(1);
  
  arma::vec sens; sens.ones(n_t*(n_th+1));
  arma::vec spec; spec.zeros(n_t*(n_th+1));
  for (int r = 1; r < n_th + 1; r++){
    for (int j = 0; j < n_t; j++){
      senumer= 0.0, sedenom=0.0, spnumer=0.0, spdenom=0.0;
      for (int i = 0; i < n_new; i++){
        
        se_bma_cur = 0.0;
        if(times(j) >= surv_new(i)){
          for (int z=0; z < k; z++){
            se_bma_cur += (Markers(i,z) > thresh(r-1)) * probs(z);
          }
          senumer += se_bma_cur*st_new(i) / G(i);
          sedenom += st_new(i) / G(i);
        }
        
        sp_bma_cur = 0.0;
        tmp_var = times(j) < surv_new(i);
        for (int z=0; z < k; z++){
          sp_bma_cur += (Markers(i,z) <= thresh(r-1)) * tmp_var * probs(z);
        }
        spnumer += sp_bma_cur;
        spdenom += tmp_var;
      }
      
      if(sedenom > FLT_EPSILON){
        sens(r*n_t+j) = senumer/sedenom;
      }else{
        sens(r*n_t+j) = 0.0;
      }
      if(spdenom > FLT_EPSILON){
        spec(r*n_t+j) = spnumer/spdenom;
      }else{
        spec(r*n_t+j) = 0.0;
      }
    }
  }
  
  for (int i = 0; i < n_t; i++){
    for (int j = 0; j < n_th; j++){
      auc(i) += ((sens(i+n_t*j) + sens(i+n_t*(1+j)))/2.0) * fabs((1.0-spec(i+n_t*j)) - (1.0-spec(i+n_t*(1+j))));
    }
  }
  return(auc);
}

