#ifndef _Minmaxegval_cpp
#define _Minmaxegval_cpp

#include <limits>
#include "AdMatAlgFoo.h"
#include "MinMaxEgval.h"

const double INF = std::numeric_limits<double>::max();

bool MinMaxEgval(const arma::mat& Sigma, const int Cf, double& lnviol,
  double& minmppegv, double& maxmpegv, double& minlregv, double& maxlregv,  
  const double eps, const double minlndet, const double maxlnk2)
{
  double minegv,maxegv,mincorregv,maxcorregv;
  int p(Sigma.n_rows);
  int q = p/2;
  lnviol = 0.;

  if (Cf==1) {
    arma::vec egval(p);
    double logDet; 
    bool pdmat = chcksing(Sigma,logDet,lnviol,mincorregv,maxcorregv,minlndet,maxlnk2,true);
    if (!pdmat) {
      return false;
    }
    pdmat = eig_sym(egval,Sigma); 
    minegv = egval(0);
    maxegv = egval(p-1);
    if (!pdmat || minegv <= 0. ) {
      lnviol = INF;
      return false;
    }
    if (minegv < eps) {
      lnviol = log(eps/minegv); 
      return false;
    } 
    pdmat = eig_sym(egval,Sigma.submat(0,0,q-1,q-1));  // MidPoints 
    if (!pdmat) return false;
    minmppegv = egval(0);
    maxmpegv = egval(q-1);
    pdmat = eig_sym(egval,Sigma.submat(q,q,p-1,p-1));  // LogRanges 
    if (!pdmat) return false;
    minlregv = egval(0);
    maxlregv = egval(q-1);
    return true;

  }  else if (Cf==2) {
    double a,b,c,det,tr,d,srtqd;
    double corr,corr2,Rdet,Rtr,bind;
    minegv=minmppegv=minlregv=mincorregv=INF;
    maxegv=maxmpegv=maxlregv=maxcorregv=0.;
    for (int j=0;j<q;++j) {
      a = Sigma(j,j);  
      b = Sigma(q+j,q+j);  
      minmppegv = fmin(minmppegv,a);
      maxmpegv = fmax(maxmpegv,a);
      minlregv = fmin(minlregv,b);
      maxlregv = fmax(maxlregv,b);
      c = Sigma(j,q+j);
      det = a*b - c*c;
      tr = a+b;
      d = tr*tr-4*det;
      if (d < 0.) {
        lnviol = -d;
      } else { 
        srtqd = sqrt(d);
        minegv = fmin(minegv,(tr-srtqd)/2);
        maxegv = fmax(maxegv,(tr+srtqd)/2);
        if (minegv < eps) lnviol = fmax(lnviol,log(eps/minegv));
      }
      corr2 = c*c/(a*b); 
      corr = sqrt(corr2); 
      Rdet = 1.-corr2; 
      Rtr = 2.;       
      bind = sqrt(4.-4*Rdet);
      maxcorregv = fmax(maxcorregv,Rtr+bind); 
      mincorregv = fmin(mincorregv,Rtr-bind);   
    }
    if (log(maxcorregv/mincorregv) > maxlnk2) lnviol = INF;        
    if (lnviol > 0.) return false;
    else return true;

  }  else if (Cf==3) {
    arma::vec egval(q);
    double minmpcorregv,maxmpcorregv,minlrcorregv,maxlrcorregv,logDet; 
    arma::mat MidPCov = Sigma.submat(0,0,q-1,q-1); 
    bool pdmat = chcksing(MidPCov,logDet,minmpcorregv,maxmpcorregv,lnviol,minlndet,maxlnk2,true);
    if (!pdmat) return false;
    pdmat = eig_sym(egval,MidPCov);   //  MidPoints
    minmppegv = egval(0);
    maxmpegv = egval(q-1);
    if (!pdmat || minmppegv <= 0.) {
      lnviol = INF;
      return false;
    }
    arma::mat LogRCov = Sigma.submat(q,q,p-1,p-1); 
    pdmat = chcksing(LogRCov,logDet,minlrcorregv,maxlrcorregv,lnviol,minlndet,maxlnk2,true);
    mincorregv = fmin(minmpcorregv,minlrcorregv);
    maxcorregv = fmax(maxmpcorregv,maxlrcorregv);
    if (!pdmat || log(maxcorregv/mincorregv) > maxlnk2) {
      lnviol = INF;
      return false;
    }
    pdmat = eig_sym(egval,LogRCov);   //  LogRanges 
    minlregv = egval(0);
    maxlregv = egval(q-1);
    minegv = fmin(minmppegv,minlregv);
    maxegv = fmax(maxmpegv,maxlregv);
    if (!pdmat || minegv <= 0. ) {
      lnviol = INF;
      return false;
    }
    if (minegv < eps) {
      lnviol = log(eps/minegv); 
      return false;
    } 
    else return true;

  } else if (Cf==4) {
    minmppegv = maxmpegv = Sigma(0,0);
    minlregv = maxlregv = Sigma(q,q);
    for (int j=1;j<q;++j) {
      double midpel = Sigma(j,j);
      minmppegv = fmin(minmppegv,midpel);
      maxmpegv = fmax(maxmpegv,midpel);
      double logrel = Sigma(q+j,q+j);
      minlregv = fmin(minlregv,logrel);
      maxlregv = fmax(maxlregv,logrel);
      if (midpel < eps) lnviol = fmax(lnviol,log(eps/midpel));
      if (logrel < eps) lnviol = fmax(lnviol,log(eps/logrel));
    }
    if (lnviol > 0.) return false;
    else return true;
  }

  return false;  // Cf is not any of the valid values  
}


#endif

