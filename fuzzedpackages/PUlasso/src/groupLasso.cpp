#include "groupLasso.h"
using namespace Eigen;

//Constructor
template <class TX>
groupLassoFit<TX>::groupLassoFit(TX & X_,
                                 VectorXd & y_,
                                 VectorXd & icoef_,
                                 ArrayXd & gsize_,
                                 ArrayXd & pen_,
                                 ArrayXd & lambdaseq_,
                                 bool isUserLambdaseq_,
                                 int pathLength_,
                                 double lambdaMinRatio_,
                                 int maxit_,
                                 VectorXd & wei_,
                                 bool weiOption_,
                                 double tol_,
                                 bool verbose_,
                                 int trace_)
  :X(X_),y(y_), gsize(gsize_), pen(pen_),lambdaseq(lambdaseq_), isUserLambdaseq(isUserLambdaseq_),
   pathLength(pathLength_),lambdaMinRatio(lambdaMinRatio_),maxit(maxit_),wei(wei_),weiOption(weiOption_),
   tol(tol_),verbose(verbose_),trace(trace_), iter(0),resid(y_),converged_CD(false),converged_KKT(false)
{
  
  checkDesignMatrix(X);
  
  N = static_cast<int>(X.rows());
  p = static_cast<int>(X.cols())+1;
  J = static_cast<int>(gsize.size());
  K = isUserLambdaseq?(static_cast<int>(lambdaseq.size())):(pathLength);
  centerFlag = false;
  grpSIdx=ArrayXi::Zero(J);
  
  for(int ii=2;ii<J;++ii){
    grpSIdx(ii)=grpSIdx(ii-1)+gsize(ii-1);
  }
  
  iters = ArrayXi::Zero(K);
  coefficients = MatrixXd::Zero(p, K);
  std_coefficients = MatrixXd::Zero(p, K);
  
  //Normalize weights
  if(weiOption){wei=(N*wei)/wei.sum();}
  
  //Calculate Xcenter, Rinvs
  //For a dense class X, X = P0X, Sparse or Map class, no change in X
  
  Xcenter = VectorXd::Ones(p-1);
  Rinvs.resize(J);
  //if(verbose){Rcpp::Rcout<<"QR decompositions\n";}
  Rinvs_X();
  
  //Initialize beta
  beta =org_to_std(icoef_);
  
  //Initialize gradient
  g.resize(J);
  
  //Initialize active/inactive set
  //inactiveSet={1,..,J-1}
  activeSet.clear();
  for(int j=1;j<J;j++)
  {
    if(beta[j]==0) {
      inactiveSet.insert(j);
    } else {
      activeSet.insert(j);
    }
  }
  inactiveSet1.clear();
  inactiveSet2.clear();
  convFlag.resize(K);
  convFlag.setZero();
  
  //if(verbose){Rcpp::Rcout<<"end of construction\n";}
}

//Getters
template <typename TX>
MatrixXd groupLassoFit<TX>::getCoefficients(){return coefficients;}

template <typename TX>
MatrixXd groupLassoFit<TX>::getStdCoefficients(){return std_coefficients;}

template <typename TX>
ArrayXi groupLassoFit<TX>::getIters(){return iters;}

template <typename TX>
ArrayXd groupLassoFit<TX>::getLambdaSequence(){return lambdaseq;}

template <typename TX>
ArrayXi groupLassoFit<TX>::getconvFlag(){return convFlag;}

//Misc functions
template <typename TX>
VectorXd groupLassoFit<TX>::back_to_org(const VectorXd & beta)
{
  VectorXd gamma(beta);
  
  for(int j=1;j<J;++j)
  {
    gamma.segment(grpSIdx(j)+1,gsize(j))= Rinvs[j]*beta.segment(grpSIdx(j)+1,gsize(j));
  }
  gamma(0) = beta(0)-gamma.segment(1,p-1).adjoint()*Xcenter;
  
  return gamma;
}

template <typename TX>
VectorXd groupLassoFit<TX>::org_to_std(const VectorXd & gamma)
{
  VectorXd beta(gamma);
  
  for(int j=1;j<J;++j)
  {
    beta.segment(grpSIdx(j)+1,gsize(j)) =
      Rinvs[j].inverse()*gamma.segment(grpSIdx(j)+1,gsize(j));
  }
  
  beta(0)= gamma(0)+ gamma.segment(1,p-1).adjoint()*Xcenter ;
  
  return beta;
}

//linpred = beta0+Q1beta1+...+Qpbetap
template <typename TX>
VectorXd groupLassoFit<TX>::linpred(const VectorXd & beta)
{
  VectorXd lpred(N);
  lpred.setZero();
  
  for (int j=1; j<J; ++j)
  {
    int sind = grpSIdx(j);
    lpred+=X.block(0,sind,N,gsize(j))*(Rinvs[j]*beta.segment(sind+1,gsize(j)));
  }
  // beta0 + P0 (sum Qj_beta_j)
  // if X already centered (with weight), wmean = 0
  double wmean(0);
  if(weiOption){wmean = (lpred.array()*wei.array()).mean();}else{wmean = lpred.mean();}
  //    std::cout<<"wmean"<<wmean<<std::endl;
  lpred = lpred.array()-(wmean-beta(0));
  return lpred;
}

template <typename TX>
VectorXd groupLassoFit<TX>::linpred_update(const VectorXd & new_resid, const VectorXd & old_resid, const VectorXd & old_lpred)
{
  return old_resid-new_resid+old_lpred;
}

//(R^{-1})'(P0Xj)'Wresid/N
// P0Xj = Xj when dense matrix
template <class TX>
VectorXd groupLassoFit<TX>::gr(const TX & X, int j, const VectorXd & resid){
  
  int sind = grpSIdx(j);
  VectorXd wresid(resid),gr;
  
  if(weiOption){wresid = wei.array()*resid.array();}
  gr = Rinvs[j].adjoint()*((X.block(0,sind,N,gsize(j)).adjoint()*wresid))/N;
  
  return gr;
}

template<>
VectorXd groupLassoFit<SparseMatrix<double> >::gr(const SparseMatrix<double> & X, int j, const VectorXd & resid)
{
  int sind = grpSIdx(j);
  VectorXd wresid(resid),gr;
  
  if(weiOption){wresid = wei.array()*resid.array();}
  gr = Rinvs[j].adjoint()*((X.block(0,sind,N,gsize(j)).adjoint()*wresid))/N
    - Rinvs[j].adjoint()*Xcenter.segment(sind,gsize(j))*wresid.mean();
    
    return gr;
}


template <class TX>
ArrayXd groupLassoFit<TX>::computeLambdaSequence(const VectorXd & y)
{
  ArrayXd lambda_path(pathLength);
  VectorXd gradnorm(J);
  double lammax(1);
  double TOLERANCE(1e-08);
  
  gradnorm.setZero();
  
  double w_ymean(0);
  //w_ymean = 1'Wy/n
  if(weiOption){w_ymean = (wei.array()*y.array()).mean();}else{w_ymean = y.mean();}
  
  // resid0 = y-beta0;
  VectorXd resid0 = y.array()-w_ymean;
  for (int j=1; j<J;++j)
  {
    g[j] = gr(X, j, resid0);
    gradnorm(j) = g[j].norm()/pen(j);
    
  }
  
  lammax = gradnorm.maxCoeff()+TOLERANCE;
  
  double logDiff=std::log(lammax)-std::log(lambdaMinRatio*lammax);
  double ratio=std::exp(-logDiff/(pathLength-1));
  
  lambda_path(0)=lammax;
  for (int i=1; i<pathLength;++i)
  {
    lambda_path(i)=lambda_path(i-1)*ratio;
  }
  
  return lambda_path;
}


template <typename TX>
bool groupLassoFit<TX>::checkKKT_j(int j, const VectorXd & resid, const ArrayXd & lambda_k)
{
  bool kkt_at_j(true);
  int sind = grpSIdx(j);
  Map<VectorXd> bj(&beta.coeffRef(sind+1),gsize(j));
  VectorXd bj_old = bj;
  VectorXd zj;
  if(j>0){
    g[j] = gr(X, j, resid);
    zj = g[j]+bj_old;
  }
  
  kkt_at_j = zj.norm()<lambda_k(j);
  return kkt_at_j;
}

template <typename TX>
bool groupLassoFit<TX>::KKT(const VectorXd & resid, const ArrayXd & lambda_k, int setidx)
{
  bool violation(false);
  std::set<int> set;
  if(setidx==1){set = inactiveSet1;}
  else if (setidx==2){set = inactiveSet2;}
  
  if(set.size()>0)
  {
    std::set<int>::const_iterator it;
    it = set.begin();
    std::set<int> newActiveSet;
    
    while(it!=set.end())
    {
      if(!checkKKT_j(*it, resid, lambda_k))
        newActiveSet.insert(*it);
      it++;
    }
    
    if(newActiveSet.size()>0)
    {
      violation = true;
      if(setidx==1)
      {
        
        for(std::set<int>::const_iterator it_nAS = newActiveSet.begin();
            it_nAS!= newActiveSet.end();it_nAS++)
        {
          inactiveSet.erase(*it_nAS);
          inactiveSet1.erase(*it_nAS);
          activeSet.insert(*it_nAS);
        }
      }
      else if(setidx==2)
      {
        for(std::set<int>::const_iterator it_nAS = newActiveSet.begin();
            it_nAS!= newActiveSet.end();it_nAS++)
        {
          inactiveSet.erase(*it_nAS);
          inactiveSet2.erase(*it_nAS);
          activeSet.insert(*it_nAS);
        }
      }
    }//if newActiveSet.size > 0
    iter++;
  }//if set.size > 0
  
  return violation;
}
template <class TX>
void groupLassoFit<TX>::coordinateDescent_0(VectorXd & resid)
{
  double w_resid_mean(0);
  Map<VectorXd> bj(&beta.coeffRef(0),gsize(0));
  VectorXd bj_old = bj;
  
  if(weiOption){w_resid_mean = (wei.array()*resid.array()).mean();}else{w_resid_mean = resid.mean();}
  bj = w_resid_mean+bj_old.array();//bj = mean(y)
  resid = resid.array()-(bj-bj_old).coeff(0,0);
  
  //iter++;
}
//Do BCD in active set until active set coefficients converge or number of cycle reaches iter_max
//Do coordinate descent in activeset
template <typename TX>
void groupLassoFit<TX>::blockCoordinateDescent(VectorXd & resid, const ArrayXd & lambda_k, const double tol)
{
  std::set<int>::const_iterator it;
  VectorXd beta_old(p);
  VectorXd diff(p);
  double error(1);
  converged_CD = false;
  
  while(!converged_CD&&iter<maxit)
  {
    beta_old = beta;
    
    if(activeSet.size()>0){
      it = activeSet.begin();
      while(it!=activeSet.end())
      {
        coordinateDescent_j(*it, resid, lambda_k);
        ++it;
      }//end of one cycle
      
      iter++;
    }
    diff = beta-beta_old;
    error = diff.cwiseAbs().maxCoeff();
    
    if(error<tol)
    {
      converged_CD=true;
    }
    
  }//end of while loop
}

template <typename TX>
bool groupLassoFit<TX>::quadraticBCD(VectorXd & resid, const ArrayXd & lambda_k, const double tol)
{
  converged_CD = false;
  converged_KKT = false;
  bool violation(false);
  bool violation1(false);
  
  while(iter<maxit)
  {
    while(iter<maxit)
    {
      //do BCD until alg converged or iters reach maxit in active set
      blockCoordinateDescent(resid, lambda_k, tol);
      //KKT in strong set
      violation1 = KKT(resid,lambda_k,1);
      //converged and no violation in inactiveSet1
      if(converged_CD&&!violation1)
        break;
    }//do BCD on A, check KKT on inactiveSet1
    
    violation = KKT(resid,lambda_k,2);
    
    if(converged_CD&&!violation){break;}
  }//inner loop + KKT for inactiveSet2
  
  converged_KKT = !violation1 && !violation;
  return converged_CD && converged_KKT;
}

//Rinvs
template <typename TX>
void groupLassoFit<TX>::Rinvs_X()
{
  MatrixXd w_Xj;
  VectorXd weisqrt(wei.array().sqrt());
  VectorXd w_Xl;
  
  // weighted mean
  for (int l=0;l<(p-1);++l)
  {
    w_Xl = X.col(l);
    if(weiOption){w_Xl = wei.array()*w_Xl.array();}
    Xcenter(l) = w_Xl.mean();
    X.col(l) = X.col(l).array()-Xcenter(l);
  }
  centerFlag = true;
  
  for(int j=1;j<J;++j){
    int sind = grpSIdx(j);
    if(gsize(j)>1)
    {
      //Do a QR decomposition
      // AP = QR <-> A = QRP^{-1}
      
      w_Xj = X.block(0,sind,N,gsize(j));
      if(weiOption){w_Xj = weisqrt.asDiagonal()*w_Xj;}
      
      ColPivHouseholderQR<MatrixXd> qr(w_Xj);
      if(qr.rank() < gsize(j)){throw std::invalid_argument("X(j) does not have full column rank");}
      
      MatrixXd R = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).triangularView<Upper>();
      MatrixXd P = qr.colsPermutation();
      R=R*P.inverse()/std::sqrt(N);
      
      Rinvs.at(j)= R.inverse();// QtQ = NIn. R' = R/sqrt(N)
    }
    else
    {
      w_Xj = X.block(0,sind,N,gsize(j));
      if(weiOption){w_Xj = weisqrt.asDiagonal()*w_Xj;}
      Rinvs.at(j) = w_Xj.adjoint()*w_Xj/N;
      Rinvs.at(j) = Rinvs.at(j).array().sqrt().inverse();
    }
  }
}

template <class TX>
void groupLassoFit<TX>::decenterX()
{
  for (int l=0;l<(p-1);++l)
  {
    X.col(l) = X.col(l).array()+Xcenter(l);
  }
  centerFlag = false;
}

template <class TX>
void groupLassoFit<TX>::centerX()
{
  for (int l=0;l<(p-1);++l)
  {
    X.col(l) = X.col(l).array()-Xcenter(l);
  }
}

template <class TX>
void groupLassoFit<TX>::coordinateDescent_j(int j, VectorXd & resid, const ArrayXd & lambda_k)
{
  int sind = grpSIdx(j);
  Map<VectorXd> bj(&beta.coeffRef(sind+1),gsize(j));
  VectorXd bj_old = bj;
  VectorXd zj;
  VectorXd update;
  double zjnorm(0);
  
  g[j] = gr(X,j,resid);
  zj = g[j]+bj_old;
  zjnorm = zj.norm();
  bj = ((zjnorm>lambda_k(j))?(1-(lambda_k(j)/zjnorm)):0)*zj;
  update =X.block(0,sind,N,gsize(j))*(Rinvs[j]*(bj-bj_old));
  resid -= update;
  
}
template <>
void groupLassoFit<SparseMatrix<double> >::coordinateDescent_j(int j, VectorXd & resid, const ArrayXd & lambda_k)
{
  int sind = grpSIdx(j);
  Map<VectorXd> bj(&beta.coeffRef(sind+1),gsize(j));
  VectorXd bj_old = bj;
  VectorXd zj;
  VectorXd update;
  double zjnorm(0);
  
  g[j] = gr(X,j,resid);
  zj = g[j] + bj_old;
  zjnorm = zj.norm();
  bj = ((zjnorm>lambda_k(j))?(1-(lambda_k(j)/zjnorm)):0)*zj;
  
  update = X.block(0,sind,N,gsize(j))*(Rinvs[j]*(bj-bj_old));
  double wupdate_mean(0);
  if(weiOption){wupdate_mean = (update.array()*wei.array()).mean(); }else{wupdate_mean = update.mean();}
  update = update.array()-wupdate_mean;
  resid -= update;
  
}
///////////////////////////////////////////////////////
//Specialization Sparse
///////////////////////////////////////////////////////
template <>
void groupLassoFit<SparseMatrix<double> >::Rinvs_X()
{
  
  MatrixXd w_Xj;
  VectorXd weisqrt(wei.array().sqrt());
  VectorXd w_Xl;
  
  for (int l=0;l<(p-1);++l)
  {
    w_Xl = X.col(l);
    if(weiOption){w_Xl = wei.array()*w_Xl.array();}
    Xcenter(l) = w_Xl.mean();
  }
  
  for(int j=1;j<J;++j)
  {
    int sind = grpSIdx(j);
    w_Xj = X.block(0,sind,N,gsize(j));
    
    int k(0);
    for(int l=sind; l<(sind+gsize(j)); ++l)
    {
      w_Xl = X.col(l);
      k = l-sind;
      w_Xj.col(k) = w_Xl.array()-Xcenter(l);
    }
    if(gsize(j)>1)
    {
      //Do a QR decomposition
      if(weiOption){w_Xj = weisqrt.asDiagonal()*w_Xj;}
      
      ColPivHouseholderQR<MatrixXd> qr(w_Xj);
      
      if(qr.rank() < gsize(j)){throw std::invalid_argument("X(j) does not have full column rank");}
      
      MatrixXd R = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).triangularView<Upper>();
      MatrixXd P =qr.colsPermutation();
      R=R*P.inverse()/std::sqrt(N);
      
      Rinvs.at(j)= R.inverse();// QtQ = NIn. R' = R/sqrt(N)
    }
    else
    {
      if(weiOption){w_Xj = weisqrt.asDiagonal()*w_Xj;}
      Rinvs.at(j) = w_Xj.adjoint()*w_Xj/N;
      Rinvs.at(j) = Rinvs.at(j).array().sqrt().inverse();
    }
    
  }
  
}
template <>
void groupLassoFit<SparseMatrix<double> >::decenterX()
{
  //do nothing
}

template <>
void groupLassoFit<SparseMatrix<double> >::centerX()
{
  //do nothing
}


template<class TX>
void groupLassoFit<TX>::checkDesignMatrix(const TX & X)
{
  for(int j=0;j<X.cols();j++)
  {
    if((X.col(j).array()==0).all()){throw std::invalid_argument("each column should have at least one non-zero element");}
  }
}


template<>
void groupLassoFit<SparseMatrix<double> >::checkDesignMatrix(const SparseMatrix<double> & X)
{
  for(int j=0;j<X.cols();j++)
  {
    if(X.col(j).nonZeros()==0){throw std::invalid_argument("each column should have at least one non-zero element");}
  }
}

//Explicit Instantiation
template class groupLassoFit<MatrixXd>;
template class groupLassoFit<SparseMatrix<double> >;
template class groupLassoFit<Map<MatrixXd> >;
