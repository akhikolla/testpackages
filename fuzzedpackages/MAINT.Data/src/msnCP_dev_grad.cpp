#include <limits>
#include "AuxFoo.h"
#include "AdMatAlgFoo.h"
#include "RestCovGrad.cpp"
#include "msnCP_dev.h"
#include "msnCP_dev_grad.h"

const double b0 = 2./(4.-PI);
const double cubrootb0 = pow(b0,(1./3));

void msnCP_ll_grad(const NumericVector& mu1, const NumericMatrix& beta2k, const mat& Sigma, 
  const NumericVector& gamma1, const NumericMatrix& y, const IntegerVector& grpind, 
  const int n, const int p, const int k, const int nvcovpar, const double limlnk2, const bool trace,
  const double PenF, const double ldRtol, const double c2tol, const double beta0tol,
  vec& CPgrad, const double Machinetol, const bool FixedArrays)
{
  int p2(p*p);

  static vec sigma,mu0,SigmaImu0,mu0bar,ksi1;
  static rowvec y0sum,sumz1eta,y0i;
  static vec z1eta;
  static vec penaltygrad;
  static mat SigmaI,Ip,D33,Dtld33,tmpM,tmpM2; // See if there is a diagonal matrix -- for Dtld33!!!
  static mat y0,yg0sum,D23,Dtld32,D32;     	// Try sparse matrix representations for Dtld32  !!!
  static int Ipdim(0);
  static ivec dvecind;
  static int dvecinddim(0);
  static mat Omega,OmegaInv,mu0OtP,sclMatI,OmgbI; 
  static vec omega,bomega,delta,OmgbIdelta,tmpv,eta;
  static vec ksi1grad;
  static rowvec etagrad,gamma1grad;
  static mat beta2kgrad;
  static rowvec Omegagrad,OmegaInvgrad1,Sigmagrad;
  static mat sumgz1eta,OmegaInvgrad2;
  static mat nOmgminusS0;

  if (dvecinddim!=p2)  {
    dvecind.set_size(p2);
    for (int r=0,i1=0,i2=0;r<p;r++) for (int c=0;c<p;c++,i1++)
      if (c<r) dvecind(i1) = (int)NAN;
      else dvecind(i1) = i2++;
      dvecinddim = p2;
  }

  if (y0.n_rows!=n || y0.n_cols!=p) y0.set_size(n,p);
  SetZero(D23,nvcovpar,p,true);
  SetZero(D32,p,nvcovpar,true);
  SetZero(Dtld32,p,nvcovpar,true);
  SetZero(Dtld33,p,p,!FixedArrays);
  if (Omegagrad.size()!=nvcovpar)  Omegagrad.set_size(nvcovpar);
  if (OmegaInvgrad1.size()!=nvcovpar)  OmegaInvgrad1.set_size(nvcovpar);
  if (OmegaInvgrad2.n_rows!=nvcovpar || OmegaInvgrad2.n_cols!=nvcovpar)  OmegaInvgrad2.set_size(nvcovpar,nvcovpar);
  if (Sigmagrad.size()!=nvcovpar)  Sigmagrad.set_size(nvcovpar);
  if (y0sum.size()!=p) y0sum.set_size(p);
  SetZero(yg0sum,k,p,true);
  SetZero(sumgz1eta,k,p,true);

  if (!FixedArrays)  {
    if (sigma.size()!=p) sigma.set_size(p);
    if (mu0.size()!=p) mu0.set_size(p);
    if (SigmaI.n_rows!=p || SigmaI.n_cols!=p)  SigmaI.set_size(p,p);
    if (SigmaImu0.size()!=p)  SigmaImu0.set_size(p);
    if (mu0bar.size()!=p)  mu0bar.set_size(p);
    if (Dtld33.n_rows!=p || Dtld33.n_cols!=p)  Dtld33.set_size(p,p);
    if (ksi1.size()!=p)  ksi1.set_size(p);
    if (Omega.n_rows!=p || Omega.n_cols!=p)  Omega.set_size(p,p);
    if (OmegaInv.n_rows!=p || OmegaInv.n_cols!=p)  OmegaInv.set_size(p,p);
    if (mu0OtP.n_rows!=p || mu0OtP.n_cols!=p)  mu0OtP.set_size(p,p);
    if (omega.size()!=p)  omega.set_size(p);
    if (bomega.size()!=p)  bomega.set_size(p);
    if (delta.size()!=p)  delta.set_size(p);
    if (OmgbIdelta.size()!=p)  OmgbIdelta.set_size(p);
    if (sclMatI.n_rows!=p || sclMatI.n_cols!=p)  sclMatI.set_size(p,p);
    if (OmgbI.n_rows!=p || OmgbI.n_cols!=p)  OmgbI.set_size(p,p);
    if (tmpM.n_rows!=p || tmpM.n_cols!=p)  tmpM.set_size(p,p);
    if (tmpM2.n_rows!=p || tmpM2.n_cols!=p)  tmpM2.set_size(p,p);
    if (tmpv.size()!=p)  tmpv.set_size(p);
    if (eta.size()!=p)  eta.set_size(p);
    if (D33.n_rows!=p || D33.n_cols!=p)  D33.set_size(p,p);
    if (ksi1grad.size()!=p)  ksi1grad.set_size(p);
    if (etagrad.size()!=p)  etagrad.set_size(p);
    if (beta2kgrad.n_rows!=k-1 || beta2kgrad.n_cols!=p)  beta2kgrad.set_size(k-1,p);
    if (nOmgminusS0.n_rows!=p || nOmgminusS0.n_cols!=p)  nOmgminusS0.set_size(p,p);
    if (y0i.size()!=p) y0i.set_size(p);
    if (sumz1eta.size()!=p) sumz1eta.set_size(p);
    if (gamma1grad.size()!=p) gamma1grad.set_size(p);
    if (z1eta.size()!=p) z1eta.set_size(p);
  }

  int ngradpar((k+1)*p+nvcovpar);
  SetZero(penaltygrad,ngradpar,!FixedArrays);
  if (Ipdim!=p) Ip = eye(p,p);

  for (int i=0;i<p;i++) { 
    sigma(i) = sqrt(Sigma(i,i));
    mu0(i) = sigma(i)*pow(b0*fabs(gamma1(i)),1./3);
    if (gamma1(i)<0) mu0(i) *= -1;
  }

  double lRdet,viol,minegv,maxegv;
  if ( !safepdsolve(Sigma,SigmaI,lRdet,viol,minegv,maxegv,ldRtol,limlnk2,true) || lRdet < -std::numeric_limits<double>::max() ) 
  {
    for (int i=0;i<ngradpar;i++) CPgrad(i) = 0.;
    return;
  }

  if (lRdet < ldRtol) {
    double cnst(PenF*(ldRtol-lRdet)); 
    for (int c=0,i=k*p;c<p;c++) {
      penaltygrad(i) = SigmaI(c,c) - 1./Sigma(c,c); 
      for (int r=c;r<p;r++,i++) {
        if (r!=c) penaltygrad(i) = 2*SigmaI(r,c); 
        penaltygrad(i) *= cnst;
      }
    }
  }

  SigmaImu0 = SigmaI * mu0;
  double beta02(dot(mu0,SigmaImu0));
  double beta0(sqrt(beta02));
  if (beta0>beta0tol) { 
    for (int i=0;i<p;i++) mu0bar(i) = mu0(i)/(sigma(i)*beta0);
    for (int c=0;c<p;c++) for(int j=0;j<p;j++)
    if (j==c) D23(dvecind(c*p+j),c) = 2*mu0(j);
    else if (j<c) {
      D23(dvecind(j*p+c),c) = mu0(j);
    } else {
      D23(dvecind(c*p+j),c) = mu0(j);
    }
    for (int i=0;i<p;i++)  {
      double mu0bi(mu0bar(i));
      Dtld32(i,dvecind(i*(p+1))) = beta0*mu0bi/(2*sigma(i));
      Dtld33(i,i) = (b0/(3*beta02))*sigma(i)/(mu0bi*mu0bi);
    }
  }
  else {
    SetZero(Dtld32,p,nvcovpar,true);
    SetZero(Dtld33,p,p,true);
  }

  for (int i=0;i<p;i++) ksi1(i) = mu1(i)-mu0(i);
  outerprod(p,mu0,mu0OtP);
  Omega = Sigma+mu0OtP;
  double logDet;
  if ( !safepdsolve(Omega,OmegaInv,logDet,viol,minegv,maxegv,ldRtol,limlnk2,true) ) {
    for (int i=0;i<p;i++) CPgrad(i) = 0.;
    return;
  }
  for (int i=0;i<p;i++)  {
    omega(i) = sqrt(Omega(i,i));
    bomega(i) = b*omega(i);
    delta(i) = mu0(i)/(bomega(i));
  }
  outerprod(p,omega,sclMatI);
  OmgbI = sclMatI % OmegaInv; 
  OmgbIdelta = OmgbI * delta;
  double c2(1.-dot(delta,OmgbIdelta));

  if (c2 < c2tol) {    

    static mat dOmgbdOmg,tmpextv,dOmgdgam1;     // Try sparse matrix representations for dOmgbdOmg !!!
    static rowvec vtmpM,vdOmgb;
    static mat Omgb; 
    static arma::ivec diagind;
    static int diaginddim(0);
    static vec dAomgm1,Adomgm1,gamma1grp1;
    static rowvec gamma1grp2;
    double tmpc,penc(PenF*(c2-c2tol));

    if (diaginddim!=p)  {
      diagind.set_size(p);
      diagind(0) = 0;
      int cnst(p*(p+1));
      for (int i0=1,i1=p-1;i0<p;i0++,i1--) diagind(i0) = (cnst-i1*(i1+1))/2;
      diaginddim = p;
    }

    if (vtmpM.size()!=nvcovpar) vtmpM.set_size(nvcovpar);  
    if (vdOmgb.size()!=nvcovpar) vdOmgb.set_size(nvcovpar);  
    if (tmpextv.n_rows!=1 || tmpextv.n_cols!=nvcovpar) tmpextv.set_size(1,nvcovpar);  
    SetZero(dOmgbdOmg,nvcovpar,nvcovpar,true); 
    if (dOmgdgam1.n_rows!=nvcovpar || dOmgdgam1.n_cols!=p) dOmgdgam1.set_size(nvcovpar,p);  
    if (!FixedArrays)  {
      if (Omgb.n_rows!=p || Omgb.n_cols!=p)  Omgb.set_size(p,p);
      if (dAomgm1.size()!=p)  dAomgm1.set_size(p);
      if (Adomgm1.size()!=p)  Adomgm1.set_size(p);
      if (gamma1grp1.size()!=p)  gamma1grp1.set_size(p);
      if (gamma1grp2.size()!=p)  gamma1grp2.set_size(p);
    }
    Omgb = Omega / sclMatI; 
    for (int c=0,i=1;c<p-1;c++,i++) for (int r=c+1;r<p;r++,i++) {
      dOmgbdOmg(i,i) = 1./sclMatI(r,c);
      dOmgbdOmg(i,diagind(r)) = -Omgb(r,c)/(2*Omega(r,r));
      dOmgbdOmg(i,diagind(c)) = -Omgb(r,c)/(2*Omega(c,c));
    }
    outerprod(p,OmgbIdelta,tmpM);
    for (int i=0,c=0;c<p;c++) for (int r=c;r<p;r++,i++)
      if (r==c) vtmpM(i) = tmpM(r,c);
      else vtmpM(i) = 2*tmpM(r,c);
      vdOmgb = vtmpM * dOmgbdOmg;
      tmpextv = vdOmgb * D23 * Dtld32;
      for (int i=0;i<nvcovpar;i++) penaltygrad(k*p+i) -= penc*(vdOmgb(i)+tmpextv(i));
 
      dOmgdgam1 = D23 * Dtld33;		
      for (int i=0;i<p;i++) {
        tmpc = cubrootb0*sigma(i);
        dAomgm1(i) = tmpc / (3*bomega(i)*pow(gamma1(i)*gamma1(i),1./3));
        Adomgm1(i) = -( tmpc*pow(fabs(gamma1(i)),1./3) / (bomega(i)*omega(i)) ) * dOmgdgam1(diagind(i),i)/(2*omega(i)) ;
        if (gamma1(i) < 0.) Adomgm1(i) *= -1; 
      }
      gamma1grp1 = -2*OmgbIdelta % (dAomgm1+Adomgm1);    
      gamma1grp2 = vdOmgb*dOmgdgam1; 
      for (int i=0;i<p;i++) penaltygrad(k*p+nvcovpar+i) -= penc*(gamma1grp1(i)+gamma1grp2(i)); 
      if (c2 < Machinetol) {
        CPgrad = penaltygrad;
        return;
      }
    }

    double b2(b*b);
    double c1(sqrt((b2-(1.-b2)*beta02)/(1.+beta02)));
    double q1(1./(c1*(1.+beta02)));    
    double q2(q1*(2*c1-q1)/2);
    double q1q2(q1*q2);
    tmpv = sqrt(fabs(q1q2)) * SigmaImu0;
    outerprod(p,tmpv,tmpM);
    if (q1q2 < 0.) tmpM *= -1;
    D33 = q1*SigmaI-tmpM;    
    for (int i=0;i<p;i++) for (int c=0,j=0;c<p;c++) for (int r=c;r<p;r++,j++) 
    if (r==c) D32(i,j) = -SigmaImu0(c) * D33(i,c) ;
    else D32(i,j) = -SigmaImu0(c)*D33(i,r) - SigmaImu0(r)*D33(i,c) ;
    D33 -= tmpM;    
    eta = q1*SigmaImu0;

    for (int r=0;r<n;r++) for(int c=0;c<p;c++) { 
      if (grpind(r)<0) y0(r,c) = y(r,c) - ksi1(c); 
      else y0(r,c) = y(r,c) - ksi1(c) - beta2k(grpind(r),c);
    } 

    nOmgminusS0 = n*Omega - y0.t()*y0;
    for (int c=0,j=0;c<p;c++) for (int r=c;r<p;r++,j++)  { 
    if (r==c) {
      OmegaInvgrad1(j) = nOmgminusS0(c,c)/2 ;
      outerprod(p,OmegaInv.col(c),tmpM);
      for (int c1=0,j1=0;c1<p;c1++) for (int r1=c1;r1<p;r1++,j1++)  OmegaInvgrad2(j1,j) = -tmpM(r1,c1);
    } else {
      OmegaInvgrad1(j) = nOmgminusS0(r,c) ;
      outerprod(p,OmegaInv.col(c),OmegaInv.col(r),tmpM);
      outerprod(p,OmegaInv.col(r),OmegaInv.col(c),tmpM2);
      for (int c1=0,j1=0;c1<p;c1++) for (int r1=c1;r1<p;r1++,j1++)  OmegaInvgrad2(j1,j) = -tmpM(r1,c1)-tmpM2(r1,c1) ;
    }
  }

  Omegagrad = OmegaInvgrad1 * OmegaInvgrad2;
  yg0sum.row(grpind(0)+1) = y0i = y0.row(0);   
  double z1(zeta(1,dot(y0i,eta)));
  z1eta = z1*eta;
  sumgz1eta.row(grpind(0)+1) = z1eta.t();
  etagrad = z1*y0i;
  for (int obs=1;obs<n;obs++) {
    y0i = y0.row(obs);
    yg0sum.row(grpind(obs)+1) += y0i;
    z1 = zeta(1,dot(y0i,eta));
    z1eta = z1*eta;
    sumgz1eta.row(grpind(obs)+1) += z1eta.t();
    etagrad += z1*y0i;
  }
  y0sum = yg0sum.row(0);
  sumz1eta = sumgz1eta.row(0);
  for (int g=1;g<k;g++) {
    y0sum += yg0sum.row(g);
    sumz1eta += sumgz1eta.row(g);
  }

  ksi1grad = OmegaInv*y0sum.t() - sumz1eta.t();
  for (int g=1;g<k;g++) beta2kgrad.row(g-1) = yg0sum.row(g)*OmegaInv - sumgz1eta.row(g);
  for (int i=0;i<p;i++) CPgrad(i) = penaltygrad(i) + ksi1grad(i);
  for (int i=p,j=0;j<p;j++) for (int g=1;g<k;g++,i++) CPgrad(i) = penaltygrad(i) + beta2kgrad(g-1,j);
  Sigmagrad = -ksi1grad.t()*Dtld32 + Omegagrad + Omegagrad*D23*Dtld32 + etagrad*(D32+D33*Dtld32) ;
  for (int i0=0,i=k*p;i0<nvcovpar;i0++,i++) CPgrad(i) = penaltygrad(i) + Sigmagrad(i0);
  gamma1grad = -ksi1grad.t()*Dtld33 + Omegagrad*D23*Dtld33 + etagrad*D33*Dtld33 ;
  for (int i0=0,i=k*p+nvcovpar;i0<p;i0++,i++) CPgrad(i) = penaltygrad(i)+gamma1grad(i0);

  if(trace) {         
    Rprintf("msnCP.ll.grad -- penalty ="); Rprintv<vec>(2*p+nvcovpar,penaltygrad);
    Rprintf("\nCentred parameters:");
    Rprintf("mu1 = ") ; Rprintv<NumericVector>(p,mu1);
    if (k>1) { Rprintf("beta2k = ") ; RprintM<NumericMatrix>(k-1,p,beta2k); }
    Rprintf("gamma1 = ") ; Rprintv<NumericVector>(p,gamma1);
    Rprintf("Sigma =\n") ; RprintM<mat>(p,p,Sigma);
    Rprintf("Direct parameters:\n");
    Rprintf("ksi1 = ") ; Rprintv<vec>(p,ksi1);
    Rprintf("alpha = ") ; Rprintv<vec>(p,omega % eta);
    Rprintf("Omega =\n") ; RprintM<mat>(p,p,Omega);
  }

  return;
}

NumericVector msnCP_dev_grad1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind,
  const int Config, const int n, const int p, const int k, const double limlnk2,
  const bool trace, const double c2tol, const double ldRtol, const double beta0tol, 
//  const double PenF, const double MachineEPS, const bool FixedArrays)
  const double PenF, const double MachineEPS, const bool FixedArrays, const bool Srpar)
{
  int q(p/2), nvcovpar(p*(p+1)/2), nvcovsrpar(ncovp(Config,q,p));
  int CPgradl((k+1)*p+nvcovpar), gradl((k+1)*p+nvcovsrpar);
  NumericVector::iterator mu1ptr(param.begin());
  NumericVector::iterator beta2kptr(mu1ptr+p);
  NumericVector::iterator Sigmaptr(mu1ptr+k*p);
  NumericVector::iterator gamma1ptr(Sigmaptr+nvcovsrpar);
  double Machinetol(sqrt(MachineEPS));

  static mat Sigma;
  static NumericMatrix beta2k;
  static vec CPgrad;
  static rowvec tmpv;
  static rowvec SigmaSrgrad;
  static mat DSigSigSr,tmpM;   // Try sparse matrix representation for DSigSigSr !!! 

  std::vector<int> SigtoSrindmap(nvcovsrpar);
  for (int i=0;i<p;i++) SigtoSrindmap[i] = i*(p+1)-i*(i+1)/2;
  for (int i0=0,i=p,c=0;c<p;c++) for (int r=c;r<p;r++,i0++) if (c!=r) {
    if ( (c<q && r<q) || (c>=q && r>=q) ) {
      if (Config==1 || Config==4) SigtoSrindmap[i++] = i0;
    } else { 
      if ( Config==1 || (Config==3 && r==c+q) ) SigtoSrindmap[i++] = i0;
    }
  }

  static NumericVector paramgrad;
  if (paramgrad.size()!=gradl) paramgrad = clone(param); 
/*
for (int i=0;i<gradl;i++) paramgrad(i) = 0.;
return paramgrad;
*/

  SetZero(DSigSigSr,nvcovpar,nvcovsrpar,true); 
  if (SigmaSrgrad.size()!=nvcovsrpar)  SigmaSrgrad.set_size(nvcovsrpar);
  if (tmpM.n_rows!=nvcovpar || tmpM.n_cols!=nvcovsrpar) tmpM.set_size(nvcovpar,nvcovsrpar);
  if ( !FixedArrays)  {
    if (CPgrad.size()!=CPgradl)  CPgrad.set_size(CPgradl);
    if (Sigma.n_rows!=p || Sigma.n_cols!=p)  Sigma.set_size(p,p);
    if (tmpv.size()!=nvcovpar)  tmpv.set_size(nvcovpar);
  }
  
  NumericVector mu1(mu1ptr,mu1ptr+p);
  if (k>1)  beta2k = NumericMatrix(k-1,p,beta2kptr);
  Sigma = RestCov(q,Sigmaptr,Config,FixedArrays,Srpar);
  NumericVector gamma1(gamma1ptr,gamma1ptr+p);

  msnCP_ll_grad( mu1,beta2k,Sigma,gamma1,y,grpind,n,p,k,nvcovpar,limlnk2,trace,
    PenF,ldRtol,c2tol,beta0tol,CPgrad,MachineEPS,FixedArrays );

  bool zerograd(true);
  for (int i=0;zerograd&&i<CPgradl;i++)  if (fabs(CPgrad(i)) > Machinetol) zerograd = false;
  if (zerograd) {
    for (int i=0;i<gradl;i++) paramgrad(i) = 0.;
    return paramgrad;
  }
  for (int i=0;i<k*p;i++) paramgrad(i) = -2*CPgrad(i);
  for (int i0=k*p+nvcovpar,i1=k*p+nvcovsrpar;i0<CPgradl;i0++,i1++) paramgrad(i1) = -2*CPgrad(i0);

  if (Srpar) {
    RestCov_grad<mat,mat,vec>(p,q,nvcovpar,Config,param.begin()+k*p,FixedArrays,DSigSigSr);
    for (int c=0,ind=0;c<p;c++) for (int r=c;r<p;r++,ind++) {
      tmpv(ind) = CPgrad(p*k+ind);
      tmpM.row(ind) = DSigSigSr.row(utind1(c,r));
    }
    SigmaSrgrad = tmpv * tmpM; 
    for (int i0=0,i1=k*p;i0<nvcovsrpar;i0++,i1++) paramgrad(i1) = -2*SigmaSrgrad(i0);
  } else { 
      for (int i0=0,i1=k*p;i0<nvcovsrpar;i0++,i1++) paramgrad(i1) = -2*CPgrad(k*p+SigtoSrindmap[i0]);
  }

  return paramgrad;

}

