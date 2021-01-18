#include <algorithm>
#include "sampleint.h"
#include "AdMatAlgFoo.h"
#include "tle.h"

void trialstep(const NumericMatrix& X, const unsigned n,const unsigned p,
  const int maxrefstps, const bool Poolm,const int m, const int Cnf, const double limlnk2,
  double c0, const unsigned kdblstar, vector<int>& Set);

double refinementstep(const NumericMatrix& X, const unsigned n, const unsigned p,
  const int Cnf, const int maxnsteps, const double limlnk2, double c0, const unsigned k,
  const vector<int>& iSet, vector<int>& Set1, Estimate& tmpsol,
	bool ClctSt, int rep, int* nstepspt, NumericMatrix* StpLogLikpt);

void highstobsllik(const NumericMatrix& X,const int n,const int p,const int Cnf,
  double c0,const int k,const vector<int>& Set0,vector<int>& Set1,Estimate& tmpsol,
  const double limlnk2, vector<double>& obsllik);

RcppExport
SEXP Cfasttle(SEXP X_s, SEXP n_s, SEXP p_s, SEXP Poolm_s, SEXP m_s, SEXP kdblstar_s, SEXP k_s, SEXP nrep_s,
	SEXP Cnf_s, SEXP c0_s, SEXP maxrefstps_s, SEXP limlnk2_s, SEXP ClctSt_s)
{
  RNGScope scope;
  #define SCOPE

  int n(as<int>(n_s)), p(as<int>(p_s)), kdblstar(as<int>(kdblstar_s)), k(as<int>(k_s)), m(as<int>(m_s));
  bool Poolm(as<int>(Poolm_s)!=0), ClctSt(as<int>(ClctSt_s)!=0); 
  int nrep(as<int>(nrep_s)), Cnf(as<int>(Cnf_s)), maxrefstps(as<int>(maxrefstps_s));
  double c0(as<double>(c0_s)), limlnk2(as<double>(limlnk2_s));
  double tmploglik,bestloglik(-Inf);
  NumericMatrix X(X_s);
  Estimate tmpsol(p);
  vector<int> trlsolSet, rsolSet(k), bestSet(k);
  if (!Poolm) trlsolSet.resize(kdblstar);
  else {
    int trsolSetsize(0);
    for (int i=0; i<m; ++i) {
      int iind(i*n/m),nxtiind((i+1)*n/m);
      trsolSetsize += (nxtiind-iind+p+1)/2;
    }
    trlsolSet.resize(trsolSetsize);
  }
 
  NumericVector RepLogLik(nrep);
  IntegerVector RepSteps(nrep);
  NumericMatrix StpLogLik(nrep,maxrefstps), *StpLogLikpt;
  if (ClctSt)  {
    std::fill(StpLogLik.begin(),StpLogLik.end(),NumericVector::get_na());
    StpLogLikpt = &StpLogLik;
  }
  else StpLogLikpt = NULL;
  int nsteps;

  for (int rep=0;rep<nrep;rep++) {
    trialstep(X,n,p,maxrefstps,Poolm,m,Cnf,limlnk2,c0,kdblstar,trlsolSet);
    tmploglik = refinementstep(X,n,p,Cnf,maxrefstps,limlnk2,c0,k,trlsolSet,rsolSet,tmpsol,ClctSt,rep,&nsteps,StpLogLikpt);
    if (tmploglik > bestloglik) {
      bestSet = rsolSet;
      bestloglik = tmploglik;
    }
    if (ClctSt) {
      RepSteps[rep] = nsteps-1;
      RepLogLik[rep] = tmploglik;
    }
  }

  if (!isfinite(bestloglik) || bestloglik==-Inf) {
    if (!ClctSt)  return List::create(
      Named("LogLik") = -Inf,
      Named("Set") = R_NilValue,
      Named("raw.cov") = R_NilValue
    );
    else return List::create(
      Named("LogLik") = -Inf,
      Named("Set") = R_NilValue,
      Named("raw.cov") = R_NilValue,
      Named("RepSteps") = RepSteps,
      Named("RepLogLik") = RepLogLik,
      Named("StpLogLik") = StpLogLik
    );
  } else {
    if (!ClctSt)  return List::create(
      Named("LogLik") = bestloglik,
      Named("Set") = bestSet,
      Named("raw.cov") = tmpsol.SigmaE()
    );
    else  return List::create(
      Named("LogLik") = bestloglik,
      Named("Set") = bestSet,
      Named("raw.cov") = tmpsol.SigmaE(),
      Named("RepSteps") = RepSteps,
      Named("RepLogLik") = RepLogLik,
      Named("StpLogLik") = StpLogLik
    );
  }
}

void trialstep(const NumericMatrix& X, const unsigned n, const unsigned p,
  const int maxrefstps, const bool Poolm, const int m, const int Cnf, const double limlnk2,
  double c0, const unsigned kdblstar, vector<int>& Set)
{
  if (!Poolm) sampleint(n,kdblstar,Set,true);	
  else {
    NumericMatrix Xsbst(n/m+1,p);
    static vector<int> vx,trlsolSet,rsolSet;  
    if (vx.size()!=n) vx.resize(n);           
    if (trlsolSet.size()!=kdblstar) trlsolSet.resize(kdblstar);
    static vector<double> rndnmbs;
    if (rndnmbs.size()!=n) rndnmbs.resize(n);
    Estimate tmpsol(p);

    rndnmbs = as<vector<double> >(runif(n,0.,1.));
    for(unsigned i=0; i<n; ++i) vx[i] = i;
    sort(vx.begin(), vx.end(), Comp<double>(rndnmbs));  
    for (int i=0,setind=0; i<m; ++i) {
      unsigned iind(i*n/m),lind((i+1)*n/m-1),nind(lind-iind+1),kstar((nind+p+1)/2);
      if (rsolSet.size()!=kstar) rsolSet.resize(kstar);
      for (unsigned r=0;r<nind;++r) for (unsigned c=0;c<p;++c) Xsbst(r,c) = X(vx[iind+r],c);
      sampleint(nind,kdblstar,trlsolSet,true);
      refinementstep(Xsbst,nind,p,Cnf,maxrefstps,limlnk2,c0,kstar,trlsolSet,rsolSet,tmpsol,false,0,NULL,NULL);
      copy(rsolSet.begin(),rsolSet.end(),Set.begin()+setind);
      setind += kstar; 
    }
  }
  return;
}

double refinementstep(const NumericMatrix& X, const unsigned n, const unsigned p, const int Cnf, const int maxnsteps, const double limlnk2,
  double c0, const unsigned k, const vector<int>& iSet, vector<int>& Set1, Estimate& tmpsol,
  bool ClctSt, int rep, int* nstepspt, NumericMatrix* StpLogLikpt)
{
    static vector<int> Set0;
    if (Set0.size()!=k) Set0.resize(k);
    static vector<double> obsllik;
    if (obsllik.size()!=n) obsllik.resize(n);
    double loglik;
         
    highstobsllik(X,n,p,Cnf,c0,k,iSet,Set0,tmpsol,limlnk2,obsllik);
    int nrfinsteps(0);
    bool RSet01(true);
    while  ( nrfinsteps==0 || (!equal(Set0.begin(),Set0.end(),Set1.begin()) && nrfinsteps<maxnsteps) )
    {
      if (RSet01) highstobsllik(X,n,p,Cnf,c0,k,Set0,Set1,tmpsol,limlnk2,obsllik);
      else highstobsllik(X,n,p,Cnf,c0,k,Set1,Set0,tmpsol,limlnk2,obsllik);
      if (ClctSt) {
        loglik = 0.;
        if (RSet01) for (unsigned i=0;i<k;i++) loglik += obsllik[Set1[i]];
        else for (unsigned i=0;i<k;i++) loglik += obsllik[Set0[i]];
        (*StpLogLikpt)(rep,nrfinsteps) = loglik;
      }
      RSet01 = !RSet01;
      ++nrfinsteps;
    }
    if (nrfinsteps==maxnsteps && RSet01) Set1 = Set0;
    loglik = 0.;
    for (unsigned i=0;i<k;i++) loglik += obsllik[Set1[i]];
    if (ClctSt) *nstepspt = nrfinsteps;
    return loglik;
}

void parcolmeans(const NumericMatrix& X,const vector<int>& Set,arma::vec& res);

void parcovll(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
  const double c0, const double limlnk2, vector<double>& obsllik);

void parcovll3(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
  const double c0, vector<double>& obsllik);

void parcovll4(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
  const double c0, const double limlnk2, vector<double>& obsllik);

void parcovll5(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
  const double c0, vector<double>& obsllik);

void highstobsllik(const NumericMatrix& X, const int n, const int p, const int Cnf,
  double c0, const int k, const vector<int>& Set0, vector<int>& Set1, Estimate& tmpsol,
  const double limlnk2, vector<double>& obsllik)
{
  static arma::mat Xdev;
  if (Xdev.n_rows!=n || Xdev.n_cols!=p) Xdev.set_size(n,p);
  for (int i=0;i<n;i++) obsllik[i] = 0.;

  parcolmeans(X,Set0,tmpsol.muE());
  double muEj;
  for (int j=0;j<p;j++) {
    muEj = tmpsol.muE()(j);
    for (int i=0;i<n;i++) Xdev(i,j) = X(i,j) - muEj;
  }
  switch (Cnf)  {
    case 1:
      parcovll(Xdev,Set0,tmpsol.SigmaE(),c0,limlnk2,obsllik);
      break;    
    case 3:
      parcovll3(Xdev,Set0,tmpsol.SigmaE(),c0,obsllik);
      break;    
    case 4:
      parcovll4(Xdev,Set0,tmpsol.SigmaE(),c0,limlnk2,obsllik);
      break;    
    case 5:
      parcovll5(Xdev,Set0,tmpsol.SigmaE(),c0,obsllik);
      break;    
    } 
    highestval<double>(n,k,obsllik,Set1,true);
    return;
}

void parcolmeans(const NumericMatrix& X,const vector<int>& Set,arma::vec& res)
{
  int n(Set.size());
  double psum;
  for (int j=0;j<X.ncol();j++) {
    psum = 0.;
    for (int i=0;i<n;i++) psum += X(Set[i],j);
    res(j) =psum/n;
  }
  return;
}

void parcovll(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
	const double c0, const double limlnk2, vector<double>& obsllik)
{
  int n(Xdev.n_rows), p(Xdev.n_cols), k(Set.size());
  double psum;

  for (int j1=0;j1<p;j1++) for (int j2=j1;j2<p;j2++) {
    psum = 0.;
    for (int i=0;i<k;i++) psum += Xdev(Set[i],j1)*Xdev(Set[i],j2);
    Sigma(j1,j2) = psum/k;
    if (j2>j1) Sigma(j2,j1) = Sigma(j1,j2);
  }

  static arma::mat SigmaI;
  if (SigmaI.n_rows!=p || SigmaI.n_cols!=p) SigmaI.set_size(p,p);
  double lnSigdet,viol,minegv,maxegv;
  if ( !safepdsolve(Sigma,SigmaI,lnSigdet,viol,minegv,maxegv,MINLNDET,limlnk2,true) ) {
    for (int i=0;i<n;i++) obsllik[i] = -Inf;
    return;
  }
  double lnSigSrIdet(-lnSigdet/2);
  double c1(c0+lnSigSrIdet);
  for (int i=0;i<n;i++) {
    obsllik[i] += c1;
    for (int j1=0;j1<p;j1++)  {
      double Xdevij1 = Xdev(i,j1);
      obsllik[i] -= SigmaI(j1,j1)*Xdevij1*Xdevij1 / 2;
      for (int j2=0;j2<j1;j2++) obsllik[i] -= SigmaI(j1,j2)*Xdevij1*Xdev(i,j2);
    }
  }
  return;
}

void parcovll3(const arma::mat& Xdev,const vector<int>& Set,arma::mat& Sigma,
	const double c0,vector<double>& obsllik)
{
  int n(Xdev.n_rows), p(Xdev.n_cols), q(p/2), k(Set.size());

  Sigma.zeros(p,p);
  double psumMP,psumLR,psumMPLG;
  for (int j=0;j<q;j++) {
    psumMP = psumLR = psumMPLG = 0.;
    for (int i=0;i<k;i++) {
      int i1(Set[i]);
      psumMP += Xdev(i1,j)*Xdev(i1,j);
      psumLR += Xdev(i1,q+j)*Xdev(i1,q+j);
      psumMPLG += Xdev(i1,j)*Xdev(i1,q+j);
    }
    Sigma(j,j) = psumMP/k;
    Sigma(q+j,q+j) = psumLR/k;
    Sigma(j,q+j) = Sigma(q+j,j) = psumMPLG/k;
  }

  double a,b,c,d,hlflnd,c1;
  for (int j=0;j<q;j++) {
    a = Sigma(j,j);
    b = Sigma(j,q+j);
    c = Sigma(q+j,q+j);
    hlflnd = log(d=a*c-b*b)/2; 
    if (j==0) c1 = c0 - hlflnd; 
    else c1 = -hlflnd; 
    for (int i=0;i<n;i++) obsllik[i] += 
      c1 - ( (c*Xdev(i,j)*Xdev(i,j) + a*Xdev(i,q+j)*Xdev(i,q+j)) / 2 - b*Xdev(i,j)*Xdev(i,q+j) ) / d;
  }

  return;
}

void parcovll4(const arma::mat& Xdev, const vector<int>& Set, arma::mat& Sigma,
	const double c0, const double limlnk2, vector<double>& obsllik)
{
  int n(Xdev.n_rows), p(Xdev.n_cols), q(p/2);
  static arma::mat pSigma;
  if (pSigma.n_rows!=q || pSigma.n_cols!=q) pSigma.set_size(q,q);
  Sigma.zeros(p,p);
  parcovll(Xdev.submat(0,0,n-1,q-1),Set,pSigma,c0/2,limlnk2,obsllik);
  Sigma.submat(0,0,q-1,q-1) = pSigma;
  parcovll(Xdev.submat(0,q,n-1,p-1),Set,pSigma,c0/2,limlnk2,obsllik);
  Sigma.submat(q,q,p-1,p-1) = pSigma;

  return;
}

void parcovll5(const arma::mat& Xdev,const vector<int>& Set,arma::mat& Sigma,
	const double c0,vector<double>& obsllik)
{
  int n(Xdev.n_rows), p(Xdev.n_cols), k(Set.size());

  Sigma.zeros(p,p);
  double psum,Xdevij;
  for (int j=0;j<p;j++) {
    psum = 0.;
    for (int i=0;i<k;i++) {
      Xdevij = Xdev(Set[i],j);
      psum += Xdevij*Xdevij;
    }
    Sigma(j,j) = psum/k;
  }

  double SigIjj,lnSigIjj;
  for (int j=0;j<p;j++) {
    lnSigIjj = log(SigIjj=1./Sigma(j,j));
    for (int i=0;i<n;i++) 
      if (j==0) obsllik[i] = c0 + ( lnSigIjj - SigIjj*Xdev(i,j)*Xdev(i,j) )/2;
      else obsllik[i] += (lnSigIjj - SigIjj*Xdev(i,j)*Xdev(i,j) )/2;
  }

  return;
}

