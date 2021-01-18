#include <Rmath.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma; 


RcppExport SEXP FHUBNB2 (SEXP m1_, SEXP m2_, SEXP od1_, SEXP od2_){
    
  double m1 = as<double>(m1_);
  double m2 = as<double>(m2_);
  double od1 = as<double>(od1_);
  double od2 = as<double>(od2_);

  int N1 = 0; int N2 = 0; int i; int j; double corr;
  double SZ1 = 1/od1; double PB1 = SZ1/(SZ1+m1);
  double SZ2 = 1/od2; double PB2 = SZ2/(SZ2+m2);
  while(R::pnbinom(N1, SZ1, PB1, 1, 0) < 1)   N1 += 1;
  while(R::pnbinom(N2, SZ2, PB2, 1, 0) < 1)   N2 += 1;
  if(N1 > 9000) return(wrap(100));
  if(N2 > 9000) return(wrap(100));
  arma::mat Exytmp(N1, N2);
  for (i=0; i<N1; i++){
    for (j=0; j<N2; j++){
    Exytmp(i,j) = 1-R::fmax2(R::pnbinom(i, SZ1, PB1, 1, 0), 
                             R::pnbinom(j, SZ2, PB2, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1+od1*m1)*(1+od2*m2)));
  return wrap(corr);
}
  
  
RcppExport SEXP FHUBZIP (SEXP m1_, SEXP m2_, SEXP od1_, SEXP od2_){
    
    double m1 = as<double>(m1_);
    double m2 = as<double>(m2_);
    double od1 = as<double>(od1_);
    double od2 = as<double>(od2_);
    
    int N1 = 0; int N2 = 0; int i; int j; double corr;
    double pi1 = od1/(1+od1); double lambda1 = (1+od1)*m1;
    double pi2 = od2/(1+od2); double lambda2 = (1+od2)*m2;
    while(pi1+(1-pi1)*R::ppois(N1, lambda1, 1, 0) < 1)  N1 += 1;
    while(pi2+(1-pi2)*R::ppois(N2, lambda2, 1, 0) < 1)  N2 += 1;
    if(N1 > 9000) return(wrap(100));
    if(N2 > 9000) return(wrap(100));
    arma::mat Exytmp(N1, N2);
    for (i=0; i<N1; i++){
        for (j=0; j<N2; j++){
            Exytmp(i,j) = 1-R::fmax2(pi1+(1-pi1)*R::ppois(i, lambda1, 1, 0), 
                                     pi2+(1-pi2)*R::ppois(j, lambda2, 1, 0));
        }
    }
    corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1+od1*m1)*(1+od2*m2)));
    return wrap(corr);
}



RcppExport SEXP FHUBbinom (SEXP m1_, SEXP m2_, SEXP n1_, SEXP n2_){
    
    double m1 = as<double>(m1_);
    double m2 = as<double>(m2_);
    double n1 = as<double>(n1_);
    double n2 = as<double>(n2_);
    
   int i; int j; double corr;
   double p1 = m1/n1; double p2 = m2/n2;
   if(n1 > 9000) return(wrap(100));
   if(n2 > 9000) return(wrap(100));
   arma::mat Exytmp(n1, n2);
   for (i=0; i<n1; i++){
     for (j=0; j<n2; j++){
       Exytmp(i,j) = 1-R::fmax2(R::pbinom(i, n1, p1, 1, 0), 
                                R::pbinom(j, n2, p2, 1, 0));
    }
  }
   corr = (arma::accu(Exytmp)-m1*m2)/(sqrt(m1*m2*(1-p1)*(1-p2)));
   return wrap(corr);
}


RcppExport SEXP FHUBZIPNB2 (SEXP zipmu_, SEXP nbmu_, SEXP zipod_, SEXP nbod_){
    
    double zipmu = as<double>(zipmu_);
    double nbmu = as<double>(nbmu_);
    double zipod = as<double>(zipod_);
    double nbod = as<double>(nbod_);
    
  int N1 = 0; int N2 = 0; int i; int j; double corr;
  double pi = zipod/(1+zipod); double lambda = (1+zipod)*zipmu;
  double SZ = 1/nbod; double PB = SZ/(SZ+nbmu);
  while(pi+(1-pi)*R::ppois(N1, lambda, 1, 0) < 1)  N1 += 1;
  while(R::pnbinom(N2, SZ, PB, 1, 0) < 1)   N2 += 1;
  if(N1 > 9000) return(wrap(100));
  if(N2 > 9000) return(wrap(100));
  arma::mat Exytmp(N1, N2);
  
  for (i=0; i<N1; i++){
    for (j=0; j<N2; j++){
      Exytmp(i,j) = 1-R::fmax2( pi+(1-pi)*R::ppois(i, lambda, 1, 0), 
                                R::pnbinom(j, SZ, PB, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-zipmu*nbmu)/(sqrt(zipmu*nbmu*(1+zipod*zipmu)*(1+nbod*nbmu)));
  return wrap(corr);
}


RcppExport SEXP FHUBZIPbinomial (SEXP zipmu_, SEXP bmu_, SEXP zipod_, SEXP bn_){
    
    double zipmu = as<double>(zipmu_);
    double bmu = as<double>(bmu_);
    double zipod = as<double>(zipod_);
    double bn = as<double>(bn_);
    
  int N1 = 0; int i; int j; double corr;
  double pi = zipod/(1+zipod); double lambda = (1+zipod)*zipmu;
  double p = bmu/bn;
  while(pi+(1-pi)*R::ppois(N1, lambda, 1, 0) < 1)  N1 += 1;
  if(N1 > 9000) return(wrap(100));
  if(bn > 9000) return(wrap(100));
  arma::mat Exytmp(N1, bn);
  
  for (i=0; i<N1; i++){
    for (j=0; j<bn; j++){
      Exytmp(i,j) = 1-R::fmax2( pi+(1-pi)*R::ppois(i, lambda, 1, 0), 
                                R::pbinom(j, bn, p, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-zipmu*bmu)/(sqrt(zipmu*bmu*(1+zipod*zipmu)*(1-p)));
  return wrap(corr);
}



RcppExport SEXP FHUBNB2binomial (SEXP nbmu_, SEXP bmu_, SEXP nbod_, SEXP bn_){
    
    double nbmu = as<double>(nbmu_);
    double bmu = as<double>(bmu_);
    double nbod = as<double>(nbod_);
    double bn = as<double>(bn_);
    
  int N1 = 0; int i; int j; double corr;
  double SZ = 1/nbod; double PB = SZ/(SZ+nbmu);
  double p = bmu/bn;
  while(R::pnbinom(N1, SZ, PB, 1, 0) < 1)   N1 += 1;
  if(N1 > 9000) return(wrap(100));
  if(bn > 9000) return(wrap(100));
  arma::mat Exytmp(N1, bn);
  for (i=0; i<N1; i++){
    for (j=0; j<bn; j++){
      Exytmp(i,j) = 1-R::fmax2( R::pnbinom(i, SZ, PB, 1, 0), 
                                R::pbinom(j, bn, p, 1, 0));
    }
  }
  corr = (arma::accu(Exytmp)-nbmu*bmu)/(sqrt(nbmu*bmu*(1+nbod*nbmu)*(1-p)));
  return wrap(corr);
}

