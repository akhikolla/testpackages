#include<vector>
#include <numeric>
#include <limits>
#include "AdMatAlgFoo.h"
#include "CEMGauss.h"
#include "MDataGaussLogLik.h"
#include "MinMaxEgval.h"
#include "RcppArmadillo.h"

#include "msnCP_dev.h"

using namespace Rcpp;

const double MINLIKSUM = 1e-300;
const double MINLIKEXP = -200.;
const double INF = 1e99;

//bool FillParm(const mat& X, const NumericMatrix& z, const int Cf, const bool Homoc, const double tautol, 
bool FillParm(const mat& X, const NumericMatrix& z, const int Cf, const bool Homoc, 
              const double tautol, const double MaxVarGRt, 
              const int n, const int p, const int k, std::vector<double>& nk, mat& Wk, mat& wdev,
              NumericVector& tau, NumericMatrix& muk, mat& Sigma, cube& Sigmak);

void SetRestCov(mat& Sigma,const int Cf);

RcppExport
SEXP CEMGauss(SEXP X_s, SEXP k_s, SEXP Cf_s, SEXP Homoc_s, SEXP maxiter_s, SEXP tautol_s, SEXP convtol_s,
    SEXP k2max_s, SEXP MaxSctEgvlRt_s,
    SEXP InitSolz_s, SEXP InitSoltau_s, SEXP InitSolmuk_s, SEXP InitSolSigma_s, SEXP InitSolSigmak_s, 
    SEXP InitSolLnLik_s, SEXP startwithM_s, SEXP SctEgvCnstr_s, SEXP MaxVarGRt_s)        
{

   static const double PenF = 1e6;	           // penalty factor for violation of maximal ln covariance eigenvalue ratio 	
   static const double EPSILON = std::numeric_limits<double>::epsilon();
   const double MINLNDET = log(std::numeric_limits<double>::min());
   bool validsol(true);

   NumericMatrix Xr(X_s);
   int n(Xr.nrow()),p(Xr.ncol());
   mat X(Xr.begin(),n,p,false);
   int k(as<int>(k_s)), Cf(as<int>(Cf_s)), maxiter(as<int>(maxiter_s));
   bool Homoc(static_cast<bool>(as<int>(Homoc_s))), startwithM(static_cast<bool>(as<int>(startwithM_s))); 
   bool SctEgvCnstr(static_cast<bool>(as<int>(SctEgvCnstr_s)));
   double tautol(as<double>(tautol_s)), convtol(as<double>(convtol_s));
   double maxlnk2(log(as<double>(k2max_s))), MaxSctlnEgvlRt(log(as<double>(MaxSctEgvlRt_s))), MaxVarGRt(as<double>(MaxVarGRt_s));
   NumericMatrix muk0,muk1;
   vec nrmfct(n);
   mat Likk(n,k),LikExp(n,k),Sigma0,Sigma1(p,p),SigmaInv;
   cube Sigmak0(p,p,k),Sigmak1(p,p,k);
   NumericMatrix z0,z1;
   NumericVector tau0,tau1;
   std::vector<double> nk(k); 
   double SiglnDet,viol,Likall,LnLik1;
   double mincorregv,maxcorregv;

   if (Homoc) SigmaInv.set_size(p,p);
   
   NumericMatrix InitSolz,InitSolmuk,InitSolSigma;
   NumericVector InitSoltau,InitSolSigmak;

   if (!Rf_isNull(InitSolz_s)) InitSolz = InitSolz_s;
   if (!Rf_isNull(InitSolmuk_s)) InitSolmuk = InitSolmuk_s;
   if (!Rf_isNull(InitSolSigma_s))  InitSolSigma = InitSolSigma_s;
   if (!Rf_isNull(InitSoltau_s))  InitSoltau = InitSoltau_s;
   if (!Rf_isNull(InitSolSigmak_s))  {
     InitSolSigmak = InitSolSigmak_s;
     InitSolSigmak.attr("dim") = IntegerVector::create(p,p,k);
   }  

   double LnLik(as<double>(InitSolLnLik_s));

   mat Wk(p,p),wdev(n,p); 
   static std::vector<double> tmpvct;
   if (tmpvct.size()!=n) tmpvct.resize(n);

   z0 = clone(InitSolz);
   z1 = clone(InitSolz);
   tau0 = clone(InitSoltau);
   tau1 = clone(InitSoltau);
   muk0 = clone(InitSolmuk);
   muk1 = clone(InitSolmuk);   
    
   if (Homoc) Sigma0 = as<mat>(clone(InitSolSigma));
   else {
     int p2 = p*p;
     for (int g=0,ind0=0;g<k;++g,ind0+=p2) Sigmak0.slice(g) = mat(&(InitSolSigmak(ind0)),p,p);
   }   
   if (Cf!=1) {
     if (Homoc) SetRestCov(Sigma0,Cf);
     else for (int g=0;g<k;++g) SetRestCov(Sigmak0.slice(g),Cf); 
   }

   NumericMatrix *zpin=&z0,*zpout=&z1;
   NumericVector *taupin=&tau0,*taupout=&tau1;
   NumericMatrix *mukpin=&muk0,*mukpout=&muk1;
   mat *Sigmapin=&Sigma0,*Sigmapout=&Sigma1;
   cube *Sigmakpin=&Sigmak0,*Sigmakpout=&Sigmak1;
   bool in0out1 = true;   

  bool converg(false);
  if (startwithM) {
     validsol = FillParm(X,z0,Cf,Homoc,tautol,MaxVarGRt,n,p,k,nk,Wk,wdev,tau0,muk0,Sigma0,Sigmak0);
     if (!validsol) return R_NilValue;
  }


  int iter=0;   
  while (!converg && iter<maxiter)  {
      
    iter++;
       
  //  E-step

    if (Homoc) validsol = safepdsolve(*Sigmapin,SigmaInv,SiglnDet,viol,mincorregv,maxcorregv,MINLNDET,maxlnk2,true); 
    for (int g=0; validsol && g<k; ++g)  {
        vec mukg(mukpin->begin()+g*p,p); 
      if (Homoc) {
        MDataGaussLogLik(n,p,Cf,X,mukg,0,&SigmaInv,&SiglnDet,tmpvct,validsol,maxlnk2,false);
      } else {
        mat Sigmakg = Sigmakpin->slice(g);           
        MDataGaussLogLik(n,p,Cf,X,mukg,&Sigmakg,0,&SiglnDet,tmpvct,validsol,maxlnk2,true);
      }          
      for (int obs=0;obs<n;++obs) LikExp(obs,g) = tmpvct[obs];
    }
    if (validsol) { 
      for (int obs=0;obs<n;++obs) {
        double maxLikExp = LikExp(obs,0);
        for (int g=1;g<k;++g) maxLikExp = fmax(LikExp(obs,g),maxLikExp);
        if (maxLikExp < MINLIKEXP) { 
          for (int g=0;g<k;++g) LikExp(obs,g) -= maxLikExp;
        }
        for (int g=0;g<k;++g) Likk(obs,g) = (*taupin)[g] * exp(LikExp(obs,g));
      }     
  
      for (int obs=0;obs<n;++obs) {
        Likall = 0.;
        for (int g=0;g<k;++g) Likall += Likk(obs,g);
          if (Likall>MINLIKSUM) for (int g=0;g<k;++g) (*zpout)(obs,g) = Likk(obs,g)/Likall;
          else for (int g=0;g<k;++g) (*zpout)(obs,g) = 1./k;
      }
    }

  //  M-step
  
    if (validsol) 
      validsol = FillParm(X,*zpout,Cf,Homoc,tautol,MaxVarGRt,n,p,k,nk,Wk,wdev,*taupout,*mukpout,*Sigmapout,*Sigmakpout);
    if (!validsol) LnLik1=-INF;
    else {
      LnLik1 = 0.;
      double maxSctlnMidPEgvl(0.),minSctlnMidPEgvl(INF),maxSctlnLogREgvl(0.),minSctlnLogREgvl(INF);
      double minmpegval,maxmpegval,minlregval,maxlregval,singviol;
      
      if (Homoc) {
        validsol = safepdsolve(*Sigmapout,SigmaInv,SiglnDet,viol,mincorregv,maxcorregv,MINLNDET,maxlnk2,true); 
        if (SctEgvCnstr) {
          bool PDCovar = MinMaxEgval(*Sigmapout,Cf,singviol,minmpegval,maxmpegval,minmpegval,maxmpegval,EPSILON,MINLNDET,maxlnk2);
          if (!PDCovar) validsol = false;
          minSctlnMidPEgvl = fmin(minSctlnMidPEgvl,minmpegval);
          maxSctlnMidPEgvl = fmax(maxSctlnMidPEgvl,maxmpegval);
          minSctlnLogREgvl = fmin(minSctlnLogREgvl,minlregval);
          maxSctlnLogREgvl = fmax(maxSctlnLogREgvl,maxlregval);
        }
      }
      for (int g=0; validsol && g<k; ++g)  {
        vec mukg(mukpout->begin()+g*p,p); 
        if (Homoc) {
          MDataGaussLogLik(n,p,Cf,X,mukg,0,&SigmaInv,&SiglnDet,tmpvct,validsol,maxlnk2,false);
        } else {
          mat Sigmakg = Sigmakpout->slice(g);           
          if (SctEgvCnstr) {
            bool PDCovar = MinMaxEgval(Sigmakg,Cf,singviol,minmpegval,maxmpegval,minmpegval,maxmpegval,EPSILON,MINLNDET,maxlnk2);
            if (!PDCovar) {
              validsol = false;
              break;
            }
            minSctlnMidPEgvl = fmin(minSctlnMidPEgvl,minmpegval);
            maxSctlnMidPEgvl = fmax(maxSctlnMidPEgvl,maxmpegval);
            minSctlnLogREgvl = fmin(minSctlnLogREgvl,minlregval);
            maxSctlnLogREgvl = fmax(maxSctlnLogREgvl,maxlregval);
          }
          MDataGaussLogLik(n,p,Cf,X,mukg,&Sigmakg,0,&SiglnDet,tmpvct,validsol,maxlnk2,true);
        }

        for (int obs=0;obs<n;++obs) LikExp(obs,g) = tmpvct[obs];

      }
      if (SctEgvCnstr) {
        double lnSctlnEgvlRt =fmax(log(maxSctlnMidPEgvl/minSctlnMidPEgvl),log(maxSctlnLogREgvl/minSctlnLogREgvl));
        if (lnSctlnEgvlRt > MaxSctlnEgvlRt) {
          LnLik1 -= PenF*(lnSctlnEgvlRt-MaxSctlnEgvlRt);
          validsol = false;
        }
      }

      if (validsol) {
        for (int obs=0;obs<n;++obs) {
          nrmfct(obs) = 0.;
          double maxLikExp = LikExp(obs,0);
          for (int g=1;g<k;++g) maxLikExp = fmax(LikExp(obs,g),maxLikExp);
          if (maxLikExp < MINLIKEXP) { 
            for (int g=0;g<k;++g) LikExp(obs,g) -= maxLikExp;
            nrmfct(obs) = maxLikExp;
          }
          for (int g=0;g<k;++g) Likk(obs,g) =  (*taupout)(g) * exp(LikExp(obs,g));
        }     
        for (int obs=0;obs<n;++obs) {
          Likall = 0.;
          for (int g=0;g<k;++g) Likall += Likk(obs,g);
          LnLik1 += nrmfct(obs) + log(Likall);
        }
      }

    }   

    if (!validsol) LnLik1 = -INF;
    if (!std::isfinite(LnLik1) || LnLik1==-INF) {
      converg = true;
    } else {
      if (LnLik1 < LnLik+convtol) converg = true; 
      if (LnLik1>LnLik) LnLik = LnLik1;
    }    
    
    if (!converg || !std::isfinite(LnLik1) || LnLik1 < LnLik)
    {
      if (in0out1) { 
        zpin = &z1;  zpout = &z0;
        taupin = &tau1; taupout = &tau0;
        mukpin = &muk1; mukpout = &muk0;
        Sigmapin = &Sigma1; Sigmapout = &Sigma0;
        Sigmakpin = &Sigmak1; Sigmakpout = &Sigmak0;                                   
      } else { 
        zpin = &z0; zpout = &z1;
        taupin = &tau0; taupout = &tau1;
        mukpin = &muk0; mukpout = &muk1;
        Sigmapin = &Sigma0; Sigmapout = &Sigma1;
        Sigmakpin = &Sigmak0; Sigmakpout = &Sigmak1;         
      } 
      in0out1 = !in0out1;
    } 
  }   

  if (LnLik==-INF) return R_NilValue;      
    
  IntegerVector grp(n);
  for (int obs=0;obs<n;++obs) {
    double maxz = 0.;
    for (int g=0;g<k;++g) if ((*zpout)(obs,g) > maxz) {
      grp[obs] = g+1;
      maxz = (*zpout)(obs,g);
    }
  }  
  
  int npar,Sigmapar;
  if (Cf==1) Sigmapar = p*(p+1)/2;	
  else if (Cf==2) Sigmapar = 3*p/2;	
  else if (Cf==3) Sigmapar = p*(p/2+1)/2;	
  else if (Cf==4) Sigmapar = p;
  if (Homoc) npar = p*k + Sigmapar + k-1;		
  else npar = (p+Sigmapar)*k + k-1;
//  double BIC = -2*LnLik + log(n)*npar;	 		
  double BIC = -2*LnLik + log(static_cast<double>(n))*npar;	 		
  double AIC = -2*LnLik + 2*npar;	 		    
      
  if (Homoc) return List::create(
    Named("tau")=*taupout,
    Named("muk")=*mukpout,
    Named("Sigma")=*Sigmapout,
    Named("z")=*zpout,
    Named("clusters")=grp,
    Named("LnLik")=LnLik,
    Named("npar")=npar,
    Named("BIC")=BIC,
    Named("AIC")=AIC
  );
  else return List::create(
    Named("tau")=*taupout,
    Named("muk")=*mukpout,
    Named("Sigmak")=*Sigmakpout,
    Named("z")=*zpout,
    Named("clusters")=grp,
    Named("LnLik")=LnLik,
    Named("npar")=npar,
    Named("BIC")=BIC,
    Named("AIC")=AIC
  );
 
}

bool FillParm(const mat& X, const NumericMatrix& z, const int Cf, const bool Homoc, 
              const double tautol, const double MaxVarGRt, 
              const int n, const int p, const int k, std::vector<double>& nk, mat& Wk, mat& wdev,
              NumericVector& tau, NumericMatrix& muk, mat& Sigma, cube& Sigmak)
{
  bool stopcycle(false); 
  double nktol(tautol/n);  
  static std::vector<double> minVar,maxVar;

  for (int g=0; !stopcycle && g<k; ++g) {
    nk[g] = 0;
    for (int obs=0;obs<n;++obs) nk[g] += z(obs,g);
    if (nk[g] < nktol) stopcycle = true;
    tau[g] = fmax(nk[g]/n,0.);
  }
  if (stopcycle) return false;
  if (Homoc) Sigma.zeros();
  else {
   if (minVar.size()!=p) minVar.resize(p);
   if (maxVar.size()!=p) maxVar.resize(p);
  }

  for (int g=0; g<k; ++g)  {
    for (int j=0; j<p; ++j) {
      double tmpsum = 0;
      for (int obs=0; obs<n; ++obs) tmpsum += z(obs,g) * X(obs,j);
      muk(j,g) = tmpsum/nk[g];
    }    
    for (int obs=0;obs<n;++obs) {
      double zweight = sqrt(z(obs,g)); 
      for (int j=0; j<p; ++j) wdev(obs,j) = zweight * (X(obs,j)-muk(j,g));
    }   
    
    Wk.zeros();
    for (int j1=0;j1<p;++j1) {
      double psum = 0.;
      for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j1);
      Wk(j1,j1) = psum;
      if (Cf==1) {
        for (int j2=0;j2<j1;++j2) {
          double psum = 0.;
          for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
          Wk(j1,j2) = Wk(j2,j1) = psum;
        }
      } 
    }
    if (Cf!=1 && Cf!=4) {
      int q = p/2;
      if (Cf==2) {
        for (int j1=0;j1<q;++j1) {
          double psum = 0.;
          for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,q+j1);
          Wk(j1,q+j1) = Wk(q+j1,j1) = psum;
        }
      } else if (Cf==3) {
        for (int j1=0;j1<q;++j1) for (int j2=0;j2<j1;++j2) {
          double psum = 0.;
          for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
          Wk(j1,j2) = Wk(j2,j1) = psum;
        }
        for (int j1=q;j1<p;++j1) for (int j2=q;j2<j1;++j2) {
          double psum = 0.;
          for (int obs=0;obs<n;obs++) psum += wdev(obs,j1)*wdev(obs,j2);
          Wk(j1,j2) = Wk(j2,j1) = psum;
        }
      } 
    }
    if (!Homoc) {
      for (int j1=0;j1<p;++j1) {
        double varj1g;
        varj1g = Sigmak(j1,j1,g) = Wk(j1,j1)/nk[g];
        if (g==0) minVar[j1] = maxVar[j1] = varj1g;
        else {
          minVar[j1] = fmin(minVar[j1],varj1g);
          maxVar[j1] = fmax(maxVar[j1],varj1g);
        }
        for (int j2=0;j2<j1;++j2)
          Sigmak(j1,j2,g) = Sigmak(j2,j1,g) = Wk(j1,j2)/nk[g];
      } 
    }          
    if (Homoc) {
      for (int j1=0;j1<p;++j1) {
        Sigma(j1,j1) += Wk(j1,j1)/n;
        for (int j2=0;j2<j1;++j2) {
          Sigma(j1,j2) += Wk(j1,j2)/n;
          Sigma(j2,j1) = Sigma(j1,j2);
        }
      }
    }
  }

  if (!Homoc) {
    double VarGRt = maxVar[0]/minVar[0];
    for (int j=1;j<p;++j) VarGRt = fmax(VarGRt,maxVar[j]/minVar[j]); 
    if (VarGRt>MaxVarGRt) return false;
  }
 
  return true;
}

void SetRestCov(mat& Sigma,const int Cf)
{
  int p = Sigma.n_rows,q=p/2;

  if (Cf==2) {
    for (int v1=0;v1<q;++v1) for (int v2=0;v2<q;++v2)
      if (v1!=v2)  
        Sigma(v1,v2) = Sigma(v2,v1) = Sigma(q+v1,q+v2) = Sigma(q+v2,q+v1) =
        Sigma(v1,q+v2) = Sigma(q+v2,v1) = Sigma(q+v1,v2) = Sigma(v2,q+v1) = 0.;
  } else if (Cf==3) { 
    for (int v1=0;v1<q;++v1) for (int v2=0;v2<q;++v2)
      Sigma(v1,q+v2) = Sigma(q+v2,v1) = Sigma(q+v1,v2) = Sigma(v2,q+v1) = 0.;
  } else if (Cf==4) { 
    for (int j1=0;j1<p;++j1) for (int j2=0;j2<p;++j2)
      if (j1!=j2)  
        Sigma(j1,j2) = Sigma(j2,j1)= 0.;
  }
}

