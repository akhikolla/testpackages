#include "msnCP_dev.h"
#include "msnCP_Aux.h"
#include <limits>

#include "AdMatAlgFoo.h"
#include "AuxFoo.h"

const double ln2pi = log(2.*PI);

void cov2cor(const int p,const mat& S,mat& R)
{
  for (int r=0;r<p;r++)  {
    R(r,r) = 1.;
    for (int c=0;c<r;c++) {
      R(r,c) = S(r,c)/sqrt(S(r,r)*S(c,c));
      R(c,r) = R(r,c);
    }
  }
  return;
} 

void cnvCPtoDP(const int p,const NumericVector mu,const mat& Sigma,const NumericVector gamma1,
		vec& ksi, mat& Omega, vec& alpha, mat& Omegabar, vec& delta,
		double* c2, bool* admissible, double* viol, const double c2tol, const double ldRtol, const double limlnk2, const bool FixedArrays)
{
  static vec c,muz,omega,mu0,sigmaz,tmp;
  static mat mu0OtP;
  double logDet,singviol,minegv,maxegv;

  if (!FixedArrays)  {
    if (tmp.size()!=p) tmp.set_size(p);
    if (c.size()!=p) c.set_size(p);
    if (muz.size()!=p) muz.set_size(p);
    if (omega.size()!=p) omega.set_size(p);
    if (mu0.size()!=p) mu0.set_size(p);
    if (sigmaz.size()!=p) sigmaz.set_size(p);
    if (mu0OtP.n_rows!=p || mu0OtP.n_cols!=p)  mu0OtP.set_size(p,p);
  }
  for (int i=0;i<p;i++) { 
    c(i) = pow(2*fabs(gamma1(i))/(4.-PI),1./3);
    if (gamma1(i)<0) c(i) *= -1;
    muz(i) = c(i)/sqrt(1.+c(i)*c(i));
    sigmaz(i) = sqrt(1.-muz(i)*muz(i));
    delta(i) = muz(i)/b;
    omega(i) = sqrt(Sigma(i,i))/sigmaz(i);
    mu0(i) = omega(i)*muz(i);
    ksi(i) = mu(i)-mu0(i);
  }
  outerprod(p,mu0,mu0OtP);
  Omega = Sigma+mu0OtP;
  cov2cor(p,Omega,Omegabar);
  if ( !safepdsolve(Omegabar,delta,tmp,logDet,singviol,minegv,maxegv,ldRtol,limlnk2,false) )
  {
    for (int i=0;i<p;i++) alpha(i) = (double)NAN;
    *c2 = (double)NAN;
    *admissible = false;
    *viol = -det(Omegabar);
    return;
  }
  *c2 = 1.- dot(delta,tmp);
  if (*c2 < c2tol)  {
    alpha = tmp/c2tol;
    *admissible = false;
    *viol = -*c2;
  }  else  {
    alpha = tmp/sqrt(*c2);
    *admissible = true;
    *viol = (double)NAN;
  }
  return;
} 

double msnCP_dev1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind, 
		const int Config, const int n, const int p, const int k, const double limlnk2, 
		const bool trace, const double c2tol, const double ldRtol, 
		const double PenF, const double PenC, const bool nopenalty,
//		const double MachineEPS, const bool FixedArrays)
		const double MachineEPS, const bool FixedArrays, const bool Srpar)
{
  double dbltmp,penalty,DPc2,DPviol;
  bool DPadmissible;
  int q(p/2), nSigmapar(ncovp(Config,q,p));
  NumericVector::iterator mu1ptr(param.begin());
  NumericVector::iterator beta2kptr(mu1ptr+p);
  NumericVector::iterator Sigmaptr(mu1ptr+k*p);
  NumericVector::iterator gamma1ptr(Sigmaptr+nSigmapar);

  static mat Sigma,OmegaInv,DPOmega, DPOmegabar;
  static mat y0;
  static NumericMatrix beta2k;
  static vec omega,alphoveromg,DPksi1,DPalpha,DPdelta;     
  static rowvec y0i;

  if (y0.n_rows!=n || y0.n_cols!=p) y0.set_size(n,p);
  if (!FixedArrays)  {
    if (Sigma.n_rows!=p || Sigma.n_cols!=p)  Sigma.set_size(p,p);
    if (OmegaInv.n_rows!=p || OmegaInv.n_cols!=p)  OmegaInv.set_size(p,p);
    if (omega.size()!=p) omega.set_size(p);
    if (alphoveromg.size()!=p) alphoveromg.set_size(p);
    if (y0i.size()!=p) y0i.set_size(p);
    if (DPksi1.size()!=p) DPksi1.set_size(p);
    if (DPOmega.n_rows!=p || DPOmega.n_cols!=p)  DPOmega.set_size(p,p);
    if (DPalpha.size()!=p) DPalpha.set_size(p);
    if (DPOmegabar.n_rows!=p || DPOmegabar.n_cols!=p)  DPOmegabar.set_size(p,p);
    if (DPdelta.size()!=p) DPdelta.set_size(p);
  }

//  Sigma = RestCov(q,Sigmaptr,Config,FixedArrays);
  Sigma = RestCov(q,Sigmaptr,Config,FixedArrays,Srpar);

// Rprintf("Srpar = %u Sigma =\n",Srpar);
// RprintM<mat>(p,p,Sigma);

  double ldR(log_det(Sigma).real());
  if ( !(nopenalty) && ldR < ldRtol ) {
    dbltmp = ldRtol-ldR;
    penalty = PenF * (PenC+dbltmp*dbltmp);
  }
  else penalty = 0.;
  NumericVector mu1(mu1ptr,mu1ptr+p);
  if (k>1)  beta2k = NumericMatrix(k-1,p,beta2kptr);
  NumericVector gamma1(gamma1ptr,gamma1ptr+p);

  cnvCPtoDP(p,mu1,Sigma,gamma1,DPksi1,DPOmega,DPalpha,DPOmegabar,DPdelta,&DPc2,&DPadmissible,&DPviol,
    c2tol,ldRtol,limlnk2,FixedArrays);

  if (!nopenalty && DPc2 < c2tol) {
    dbltmp = c2tol-DPc2; 
    penalty += PenF * (PenC + dbltmp*dbltmp); 
  }
  if (!DPadmissible)  {
    if (nopenalty) return INFINITY;  
    else return penalty;
  }       

  double logDet,viol,minegv,maxegv;
  if (!safepdsolve(DPOmega,OmegaInv,logDet,viol,minegv,maxegv,ldRtol,limlnk2,true)) return INFINITY; 
  for (int i=0;i<p;i++) {
    omega(i) = sqrt(DPOmega(i,i));
    alphoveromg(i) = DPalpha(i)/omega(i);
  }
  for (int r=0;r<n;r++) for(int c=0;c<p;c++) 
    if (grpind(r)<0) y0(r,c) = y(r,c) - DPksi1(c); 
    else y0(r,c) = y(r,c) - DPksi1(c) - beta2k(grpind(r),c); 

  double dev(n*(p*ln2pi+logDet));
  for (int obs=0;obs<n;obs++) {
    y0i = y0.row(obs);
    dev += dot(y0i,OmegaInv*y0i.t()) - 2*zeta(0,dot(y0i,alphoveromg));
  }

  if(trace) { 
    Rprintf("msnCP.dev %f\n",dev);
    Rprintf("Centred parameters:\n");
    Rprintf("mu1 = ") ; 
    if (k>1) { Rprintf("beta2k =") ; RprintM<NumericMatrix>(k-1,p,beta2k); }
    Rprintf("gamma1 = ") ; Rprintv<NumericVector>(p,gamma1);
    Rprintf("Sigma =\n") ; RprintM<mat>(p,p,Sigma);
    Rprintf("Direct parameters:\n");
    Rprintf("ksi1 = ") ; Rprintv<vec>(p,DPksi1);
    Rprintf("alpha = ") ; Rprintv<vec>(p,DPalpha);
    Rprintf("Omega =\n") ; RprintM<mat>(p,p,DPOmega);
  }

  return dev+penalty;
}

