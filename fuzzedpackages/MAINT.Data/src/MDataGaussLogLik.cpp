#include<cmath>
#include "AdMatAlgFoo.h"
#include "MDataGaussLogLik.h"

#include "msnCP_dev.h"

const double INF = std::numeric_limits<double>::infinity();

void MDataGaussLogLik(const int n, const int p, const int Cf, const arma::mat& X, const arma::vec& u, 
                      arma::mat* Sigmap, arma::mat* SigmaInvp, double* lndetSigp, std::vector<double>& res, bool& validsol, 
//                     const double maxlnk2, const double MaxSctlnEgvlRt, const bool chksing)
                      const double maxlnk2, const bool chksing)
{
  const double ZERO = std::numeric_limits<double>::min();  //  small value that is considered to be numerically identical to zero
  static const double PenF = 1e6;	  // penalty factor for numerically singular covariance matrices	
  double mincorregv(INF),maxcorregv(0.),Singviol,penalty;
  validsol = true;

  if (Cf==1) {

    arma::mat dev,tmp;
    if (dev.n_rows!=p || dev.n_cols!=n)  dev.set_size(p,n);
    if (tmp.n_rows!=p || tmp.n_cols!=n)  tmp.set_size(p,n);
    for (int obs=0;obs<n;++obs) for (int j=0;j<p;++j) dev(j,obs) = X(obs,j) - u(j);

    if (Sigmap) {
      if (chksing) validsol = safepdsolve(*Sigmap, dev, tmp, *lndetSigp, Singviol, mincorregv, maxcorregv, MINLNDET, maxlnk2, true);
      else validsol = pdsolve(*Sigmap, dev, tmp, lndetSigp);  
    } else if (SigmaInvp) {
      tmp = *SigmaInvp * dev;
    } else {
      Rprintf("Error on MDataGaussLogLik: Sigmap and SigmaInvp cannot both be NULL pointers\n");
      return;
    }
    if (!validsol) {
      penalty = -PenF*Singviol;
      for (int obs=0;obs<n;++obs) res[obs] = penalty;
      return;
    }  

    double c0 = -(p*LN2PI + *lndetSigp)/2;
    for (int obs=0;obs<n;++obs) res[obs] = c0 - dot(dev.col(obs),tmp.col(obs))/2;
    return;
  }
 
  else if (Cf==2) {

    int q = p/2;
    double a,b,c,det,lndet,abscorr,lnk2,dev1,dev2;

    double c0 = -(p*LN2PI)/2;
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    for (int j=0;j<q;++j) {

      if (Sigmap) {

        a = (*Sigmap)(j,j);  
        b = (*Sigmap)(q+j,q+j);  
        c = (*Sigmap)(j,q+j);
        if (a <= ZERO || b <= ZERO) {  
          validsol = false;
          penalty = -INF;
          for (int obs=0;obs<n;++obs) res[obs] = penalty;
          return;
        }
        det = a*b - c*c;
        if (det <= ZERO) {  
          validsol = false;
          penalty = -INF;
          for (int obs=0;obs<n;++obs) res[obs] = penalty;
          return;
        }
        lndet = log(det); 
        if (chksing) {
          if (lndet < MINLNDET) {
            validsol = false;
            penalty = -PenF*(MINLNDET-lndet);
            for (int obs=0;obs<n;++obs) res[obs] = penalty;
            return;
          }
          abscorr = fabs(c)/sqrt(a*b); 
          maxcorregv = fmax(maxcorregv,1+abscorr);
          mincorregv = fmin(mincorregv,1-abscorr);
        }

      } else {
        a = (*SigmaInvp)(j,j);  
        b = (*SigmaInvp)(q+j,q+j);  
        c = (*SigmaInvp)(j,q+j);
      }  

      for (int obs=0;obs<n;++obs) {
        dev1 = X(obs,j) - u(j);  
        dev2 = X(obs,q+j) - u(q+j);  
        if (Sigmap) res[obs] -= ( lndet + (a*dev2*dev2 + b*dev1*dev1 -2*c*dev1*dev2)/det ) / 2;
        else res[obs] -= (a*dev1*dev1 + b*dev2*dev2 +2*c*dev1*dev2) / 2;
      }
    }

    if (chksing) {
      if (mincorregv <= 0.) { 
        validsol = false;
        penalty = -PenF*(-mincorregv);
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }
      lnk2 = log(maxcorregv/mincorregv);
      if (lnk2 > maxlnk2) { 
        validsol = false;
        penalty = -PenF*(lnk2-maxlnk2);
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }
    }

    if (!Sigmap) {
      double hlflndetSigp = *lndetSigp/2;
      for (int obs=0;obs<n;++obs) res[obs] -= hlflndetSigp;
    }
    return;    
  }

  else if (Cf==3) {

    arma::mat dev,tmp;
    int q = p/2;
    if (dev.n_rows!=q || dev.n_cols!=n)  dev.set_size(q,n);
    if (tmp.n_rows!=q || tmp.n_cols!=n)  tmp.set_size(q,n);
    bool singularMat;
    double minmpcorreg,maxmpcorreg,minlrcorreg,maxlrcorreg,lnk2,c0(-p*LN2PI/2);
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    for (int block=0;block<2;++block) {   // block==0: MidPoints ; block==1: LogRanges  
      if (block==0) {
        for (int obs=0;obs<n;++obs) for (int j=0;j<q;++j) dev(j,obs) = X(obs,j) - u(j);
        if (Sigmap) {
          if (chksing) {
            singularMat = !safepdsolve(Sigmap->submat(0,0,q-1,q-1), dev, tmp, *lndetSigp, Singviol, minmpcorreg, maxmpcorreg, MINLNDET, maxlnk2, true);
          }  else singularMat = !pdsolve(Sigmap->submat(0,0,q-1,q-1), dev, tmp, lndetSigp);
        } 
        else tmp = SigmaInvp->submat(0,0,q-1,q-1) * dev;
      } else {
        for (int obs=0;obs<n;++obs) for (int j=0;j<q;++j) dev(j,obs) = X(obs,q+j) - u(q+j);
        if (Sigmap) {
          if (chksing) {
            singularMat = !safepdsolve(Sigmap->submat(q,q,p-1,p-1), dev, tmp, *lndetSigp, Singviol, minlrcorreg, maxlrcorreg, MINLNDET, maxlnk2, true);
          }  else singularMat = !pdsolve(Sigmap->submat(q,q,p-1,p-1), dev, tmp, lndetSigp);
        } 
        else tmp = SigmaInvp->submat(q,q,p-1,p-1) * dev;
      } 

      if (singularMat) {
        validsol = false;
        penalty = -PenF*Singviol;
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return;
      }
  
      if (Sigmap) for (int obs=0;obs<n;++obs) res[obs] -= ( *lndetSigp + dot(dev.col(obs),tmp.col(obs)) ) / 2;
      else for (int obs=0;obs<n;++obs) res[obs] -= dot(dev.col(obs),tmp.col(obs)) / 2;
    }

    if (chksing) {
      maxcorregv = fmax(maxmpcorreg,maxlrcorreg);
      mincorregv = fmin(minmpcorreg,minlrcorreg);
      lnk2 = log(maxcorregv/mincorregv);
      if (lnk2 > maxlnk2) { 
        validsol = false;
        penalty = -PenF*(lnk2-maxlnk2);
        for (int obs=0;obs<n;++obs) res[obs] = penalty;
        return; 
      }
    }    
    if (!Sigmap) {
      double hlflndetSigp = *lndetSigp/2;
      for (int obs=0;obs<n;++obs) res[obs] -= hlflndetSigp;
    }

    return;
  }

  else if (Cf==4) {

    double detj,lndetj,dev,lnk2;
    double c0 = -p*LN2PI/2;
    for (int obs=0;obs<n;++obs) res[obs] = c0;

    if (Sigmap) { 

      for (int j=0;j<p;++j) {

        detj = (*Sigmap)(j,j);

        if (detj <= ZERO) {  
          validsol = false;
          penalty = -INF;
          for (int obs=0;obs<n;++obs) res[obs] = penalty;
          return;
        }
        lndetj = log(detj); 

        if (chksing) {
          if (lndetj < MINLNDET) {
            validsol = false;
            penalty = -PenF*(MINLNDET-lndetj);
            for (int obs=0;obs<n;++obs) res[obs] = penalty;
            return;
          }
          maxcorregv = fmax(maxcorregv,detj); 
          mincorregv = fmin(mincorregv,detj);   
        }

        for (int obs=0;obs<n;++obs) {
          dev = X(obs,j) - u(j);  
          res[obs] -= (lndetj + dev*dev/detj) / 2;
        }
      }
    } else { 

      double hlflndetSigp = *lndetSigp/2;
      for (int obs=0;obs<n;++obs) {
        for (int j=0;j<p;++j) {
          dev = X(obs,j) - u(j);  
          res[obs] -= dev*dev*(*SigmaInvp)(j,j) / 2;
        }
        res[obs] -= hlflndetSigp;
      }
    }

    return;    
  }

}  


