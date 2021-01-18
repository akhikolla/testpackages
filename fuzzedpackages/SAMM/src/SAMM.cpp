
#include <RcppArmadillo.h>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace arma;
using namespace Rcpp;
const double log2pi = std::log(2.0 * M_PI);
#define FALSE  0
#define TRUE   1
const double ALPHA =  1.0;       /* reflection coefficient */
const double BETA   =     .2;       /* contraction coefficient */
const double GAMMA    =   2;       /* expansion coefficient */



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat nearPD(const arma::mat & x, const double & epsilon) {
  if( (rcond(x)<epsilon )||(rcond(x)>1/epsilon )) {
    arma::mat Eigenvecs;
    arma::vec Eigenval;
    arma::eig_sym(Eigenval , Eigenvecs ,x+std::log10(static_cast<double>(x.n_cols))*eye(x.n_cols,x.n_cols));
    Eigenval=Eigenval-std::log10(static_cast<double>(x.n_cols));
    arma::mat xnew=x;
    for (unsigned int i = 0; i < Eigenval.n_elem; i++) {
      if (Eigenval(i)<=0) {
        Eigenval(i)=epsilon;
      }
    }
    xnew=Eigenvecs*diagmat(Eigenval)*Eigenvecs.t();
    return xnew;
  }
  else {
    return x;
  }
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat inv_sympdsamm(const arma::mat & x, const double & epsilon) {
  double epsilon0=epsilon;
  bool invsuccess=false;
  arma::mat xinv;
  int counteri=0;
  do {
    try {
      if (counteri==0) {
        counteri=counteri+1;
        xinv=inv_sympd(x+epsilon0*eye(x.n_cols,x.n_cols));
        invsuccess=true;
      }
      else {
        xinv=pinv(x+epsilon0*eye(x.n_cols,x.n_cols));
        invsuccess=true;
      }
    }
    catch(...) {
      epsilon0=1.05*epsilon0;
      //Rcpp::Rcout  <<  epsilon0 << std::endl;
      invsuccess=false;
      counteri=counteri+1;
    }
  } while (!invsuccess &&  counteri<100);
  
  if (!invsuccess) {
    xinv=diagmat(1/diagvec(x+epsilon));
  }
  return xinv;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat cholsammlower(const arma::mat & x, const double & epsilon) {
  arma::mat xchol=diagmat(pow(diagvec(x+epsilon),.5));
  try {
    chol(xchol,x,"lower");
  }
  catch(...) {
    try {
      chol(xchol,x+epsilon*eye(x.n_cols,x.n_cols),"lower");
    }
    catch(...) {
      xchol=diagmat(pow(diagvec(x+epsilon),.5));
    }
  }
  return xchol;
  
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat cholsammupper(const arma::mat & x, const double & epsilon) {
  arma::mat xchol=diagmat(pow(diagvec(x+epsilon),.5));
  try {
    chol(xchol,x,"upper");
    
  }
  catch(...) {
    try {
      chol(xchol,x+epsilon*eye(x.n_cols,x.n_cols),"upper");
    }
    catch(...) {
      xchol=diagmat(pow(diagvec(x+epsilon),.5));
    }
  }
  return xchol;
  
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::vec Mahalanobis(const arma::mat & x, const arma::vec & center, const arma::mat & cov, bool covinv=false) {
  
  int n = x.n_cols;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (unsigned int i=0; i < n; i++) {
    x_cen.col(i) = x.col(i) - center;
  }
  if (covinv) {
    return sum((x_cen.t() * cov) % x_cen.t(), 1);
  }
  else {
    return sum((x_cen.t() * inv_sympdsamm(cov,1e-10)) % x_cen.t(), 1);
  }
}





// [[Rcpp::depends("RcppArmadillo")]]

double dmvnorm_arma2(const arma::mat & x, const arma::vec & mean, const arma::mat & sigma, bool logd = false, bool covinv=false) {
  
  double offset = std::log10(static_cast<double>(sigma.n_cols));
  double logdet;
  arma::mat sigmaandoffset = sigma + offset * eye(sigma.n_cols,sigma.n_cols);
  arma::vec distval = Mahalanobis(x,  mean, sigma, covinv);
  arma::vec eigvals=arma::eig_sym(sigmaandoffset)-offset;
  if (covinv) {
    logdet = sum(arma::log(1/eigvals));
  }
  else {
    logdet = sum(arma::log(eigvals));
  }
  arma::vec logretval = -( (x.n_rows * log2pi + logdet + distval)/2  ) ;
  
  if (logd) {
    return(as_scalar(logretval));
  } else {
    return(exp(as_scalar(logretval)));
  }
}



double dmvnorm_arma(const arma::mat & x, const arma::vec & mean, const arma::mat & sigma, bool logd = false, bool covinv=false) {
  double out;
  
  if (!covinv) {
    int xdim = x.n_rows;
    arma::mat rooti(xdim,xdim);
    if (covinv) {
      rooti = trimatu(cholsammlower(sigma, 1e-10));
    }
    else {
      rooti = arma::inv(trimatu(cholsammlower(sigma, 1e-10)));
    }
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    
    arma::vec z = rooti * (x - mean) ;
    out = constants - 0.5 * arma::sum(z%z) + rootisum;
    
    
    if (logd == false) {
      out = exp(out);
    }
    return out;
  }
  else {
    out=dmvnorm_arma2(x,mean,sigma,logd,covinv);
    return out;
  }
}





////////////////////////////////////////COVFUNCTIONS
// [[Rcpp::depends("RcppArmadillo")]]

arma::mat disteucarma(arma::mat Ar, arma::mat Br) {
  int m = Ar.n_rows,
    n = Br.n_rows,
    k = Ar.n_cols;
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false);
  
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  
  return sqrt(C);
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat ar1cov_cpp(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho = (2/M_PI)*atan(params(0));
  arma::mat H = disteucarma(times, times);
  arma::mat V(data(0,0),data(0,0));
  for (unsigned int i = 0; i < data(0,0); i++) {
    for (unsigned int j = 0; j < data(0,0); j++) {
      V(i,j)=powf(rho,H(i,j));
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
arma::mat ar1hetcov_cpp(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho = (2/M_PI)*atan(params(0));
  arma::mat H = disteucarma(times, times);
  arma::mat V(data(0,0),data(0,0));
  for (unsigned int i = 0; i < data(0,0); i++) {
    for (unsigned int j = 0; j < data(0,0); j++) {
      V(i,j)=powf(rho,H(i,j));
    }
  }
  arma::vec vvec(data(0,0));
  vvec(0)=1;
  for (unsigned int i = 1; i < data(0,0); i++) {
    vvec(i)=exp(params(i));
  }
  return diagmat(vvec)*V*diagmat(vvec);
}




// [[Rcpp::depends("RcppArmadillo")]]
arma::mat arma11cov_cpp(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho =(2/M_PI)*atan(params(0));
  double lambda =(2/M_PI)*atan(params(1));
  
  arma::mat H = disteucarma(times, times);
  
  arma::mat V(data(0,0),data(0,0));
  for (int i = 0; i < data(0,0); i++) {
    for (int j = 0; j < data(0,0); j++) {
      if (abs(i-j)>1) {
        V(i,j)=lambda*powf(rho,H(i,j));
      } else {
        if (i==j) {
          V(i,j)=1;
        } else {
          V(i,j)=lambda;
        }
      }
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmcov_cpp(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =data(0,0);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < data(0,0); i++) {
    V(i,i)=1;
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmhetcov_cpp(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =data(0,0);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < data(0,0); i++) {
    V(i,i)=1;
  }
  arma::vec vvec(data(0,0));
  vvec(0)=1;
  for (unsigned int i = 1; i < data(0,0); i++) {
    vvec(i)=exp(params(i));
  }
  return diagmat(vvec)*V*diagmat(vvec);
}




// [[Rcpp::depends("RcppArmadillo")]]

arma::mat lincombcov_cpp(const arma::vec & params, const arma::mat  & data) {
  int k=params.n_elem;
  arma::vec weights=params-min(params);
  weights=weights/sum(weights);
  arma::mat V(data.n_rows,data.n_rows);
  for (unsigned int i = 0; i < k; i++) {
    V=V+weights(i)*data.submat(0,i*data.n_rows,data.n_rows-1,(i+1)*data.n_rows-1);
  }
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat unstrcov_cpp(const arma::vec & params, const arma::mat  & data) {
  int d1=data(0,0);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat V=diagmat(d2vec)*D1*diagmat(d2vec);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat diagcov_cpp(const arma::vec & params, const arma::mat & data) {
  int k=params.n_elem;
  arma::vec v=ones(k+1);
  for (unsigned int i = 1; i < (k+1); ++i) {
    v(i)=exp(params(i-1));
  }
  arma::mat V=diagmat(v);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat unstrKronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(diagmat(d2vec)*D1*diagmat(d2vec),D2);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat diagKronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1=zeros(d1,d1);
  D1.diag()=join_cols(ones(1,1),exp(params));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat ar1KronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat ar1hetKronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1hetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat arma11KronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=arma11cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmKronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmhetKronKcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmhetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}






// [[Rcpp::depends("RcppArmadillo")]]

arma::mat UnstrKronUnstrcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d1vec(d1);
  d1vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d1vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2(d2,d2);
  
  for (unsigned int i = 0; i < d2; i++) {
    for (unsigned int j = i; j < d2; j++) {
      if (i==j) {
        D2(i,j)=1;
      }
      else {
        D2(i,j)=(2/M_PI)*atan(params(ii));
        D2(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d2);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d2; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  
  arma::mat V=kron(diagmat(d1vec)*D1*diagmat(d1vec),diagmat(d2vec)*D2*diagmat(d2vec));
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat rbfcov_cpp(const arma::vec & params, const arma::mat & data) {
  arma::mat C=disteucarma(data,data);
  return exp(-exp(params(0))*pow(C, 2));
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat expcov_cpp(const arma::vec & params, const arma::mat & data) {
  arma::mat C=disteucarma(data,data);
  return exp(-(exp(params(0)))*C);
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat rbfdistcov_cpp(const arma::vec & params, const arma::mat & data) {
  return exp(-exp(params(0))*pow(data, 2));
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat expdistcov_cpp(const arma::vec & params, const arma::mat & data) {
  return exp(-(exp(params(0)))*data);
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat relmatcov_cpp(const arma::vec & params, const arma::mat & data) {
  arma::vec pks(data.n_cols);
  double c=0;
  for (unsigned int iter = 0; iter < data.n_cols; ++iter) {
    pks(iter)=sum(data.col(iter)+ones(data.n_rows))/(2*data.n_rows);
    c=c+2*pks(iter)*(1-pks(iter));
  }
  arma::mat W=data+1-2*ones(data.n_rows,1)*pks.t();
  arma::mat Amat=(1/c)*W*W.t();
  return Amat+params(0)*eye(Amat.n_cols,Amat.n_cols);
}
// [[Rcpp::depends("RcppArmadillo")]]

arma::mat ConstMatcov_cpp(const arma::vec & params, const arma::mat & data) {
  return data;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKronunstrcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,diagmat(d2vec)*D1*diagmat(d2vec));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKrondiagcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1=zeros(d1,d1);
  D1.diag()=join_cols(ones(1,1),exp(params));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKronar1cov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKronar1hetcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1hetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKronarma11cov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=arma11cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKroncompsymmcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat KKroncompsymmhetcov_cpp(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmhetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat sppowcov_cpp(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      V(i,j)=pow(rho, data(i,j));
    }
  }
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat splincov_cpp(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      if (rho*data(i,j)<=1) {
        V(i,j)=(1-rho*data(i,j));
      } else {
        V(i,j)=0;
      }
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat splinlogcov_cpp(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      if (log(rho*data(i,j))<=1) {
        V(i,j)=(1-rho*log(data(i,j)));
      } else {
        V(i,j)=0;
      }
    }
  }
  return V;
}




typedef arma::mat (*funcPtr)(const vec & params,const arma::mat & data);

XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "ar1")
    return(XPtr<funcPtr>(new funcPtr(&ar1cov_cpp)));
  else if (fstr == "ar1het")
    return(XPtr<funcPtr>(new funcPtr(&ar1cov_cpp)));
  else if (fstr == "rbf")
    return(XPtr<funcPtr>(new funcPtr(&rbfcov_cpp)));
  else if (fstr == "exp")
    return(XPtr<funcPtr>(new funcPtr(&expcov_cpp)));
  else if (fstr == "rbfdist")
    return(XPtr<funcPtr>(new funcPtr(&rbfdistcov_cpp)));
  else if (fstr == "expdist")
    return(XPtr<funcPtr>(new funcPtr(&expdistcov_cpp)));
  else if (fstr == "arma11")
    return(XPtr<funcPtr>(new funcPtr(&arma11cov_cpp)));
  else if (fstr == "compsymm")
    return(XPtr<funcPtr>(new funcPtr(&compsymmcov_cpp)));
  else if (fstr == "compsymmhet")
    return(XPtr<funcPtr>(new funcPtr(&compsymmhetcov_cpp)));
  else if (fstr == "lincomb")
    return(XPtr<funcPtr>(new funcPtr(&lincombcov_cpp)));
  else if (fstr == "unstr")
    return(XPtr<funcPtr>(new funcPtr(&unstrcov_cpp)));
  else if (fstr == "RelMat")
    return(XPtr<funcPtr>(new funcPtr(&relmatcov_cpp)));
  else if (fstr == "Diag")
    return(XPtr<funcPtr>(new funcPtr(&diagcov_cpp)));
  else if (fstr == "UnstrKronK")
    return(XPtr<funcPtr>(new funcPtr(&unstrKronKcov_cpp)));
  else if (fstr == "ar1KronK")
    return(XPtr<funcPtr>(new funcPtr(&ar1KronKcov_cpp)));
  else if (fstr == "ar1hetKronK")
    return(XPtr<funcPtr>(new funcPtr(&ar1hetKronKcov_cpp)));
  else if (fstr == "arma11KronK")
    return(XPtr<funcPtr>(new funcPtr(&arma11KronKcov_cpp)));
  else if (fstr == "compsymmKronK")
    return(XPtr<funcPtr>(new funcPtr(&compsymmKronKcov_cpp)));
  else if (fstr == "compsymmhetKronK")
    return(XPtr<funcPtr>(new funcPtr(&compsymmhetKronKcov_cpp)));
  else if (fstr == "DiagKronK")
    return(XPtr<funcPtr>(new funcPtr(&diagKronKcov_cpp)));
  else if (fstr == "KKronUnstr")
    return(XPtr<funcPtr>(new funcPtr(&KKronunstrcov_cpp)));
  else if (fstr == "KKronar1")
    return(XPtr<funcPtr>(new funcPtr(&KKronar1cov_cpp)));
  else if (fstr == "KKronar1het")
    return(XPtr<funcPtr>(new funcPtr(&KKronar1hetcov_cpp)));
  else if (fstr == "KKronarma11")
    return(XPtr<funcPtr>(new funcPtr(&KKronarma11cov_cpp)));
  else if (fstr == "KKroncompsymm")
    return(XPtr<funcPtr>(new funcPtr(&KKroncompsymmcov_cpp)));
  else if (fstr == "KKroncompsymmhet")
    return(XPtr<funcPtr>(new funcPtr(&KKroncompsymmhetcov_cpp)));
  else if (fstr == "KKronDiag")
    return(XPtr<funcPtr>(new funcPtr(&KKrondiagcov_cpp)));
  else if (fstr == "UnstrKronUnstr")
    return(XPtr<funcPtr>(new funcPtr(&UnstrKronUnstrcov_cpp)));
  else if (fstr == "Const")
    return(XPtr<funcPtr>(new funcPtr(&ConstMatcov_cpp)));
  else if (fstr == "sppow")
    return(XPtr<funcPtr>(new funcPtr(&sppowcov_cpp)));
  else if (fstr == "splin")
    return(XPtr<funcPtr>(new funcPtr(&splincov_cpp)));
  else if (fstr == "splinlog")
    return(XPtr<funcPtr>(new funcPtr(&splinlogcov_cpp)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}




arma::mat callViaString(const arma::vec params, const arma::mat data, std::string funname) {
  XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
  funcPtr fun = *xpfun;
  arma::mat y = fun(params, data);
  return (y);
}





// [[Rcpp::depends("RcppArmadillo")]]

arma::mat symMroot(const arma::mat & M) {
  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, M+eye(M.n_cols,M.n_cols));
  eigval=eigval-1;
  
  for (unsigned int i=0; i<eigval.size(); i++) {
    if(eigval(i)<1e-8) {
      eigval(i)=1e-8;
    }
  }
  return eigvec*diagmat(sqrt(eigval))*eigvec.t();
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec loglikfuncmmmkmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const Rcpp::List & sigmahatlist,const arma::mat & B,const arma::mat & W,const arma::mat & R ) {
  int k=Zlist.size();
  arma::mat V=kron(W*R*W.t(),Rcpp::as<arma::mat>(sigmahatlist(k)));
  arma::mat Vgt=Rcpp::as<arma::mat>(sigmahatlist(k));
  arma::mat Z;
  arma::mat K;
  arma::mat ZKZt;
  arma::vec loglik;
  for (unsigned int i = 0; i < k; i++) {
    K=Rcpp::as<arma::mat>(Klist(i));
    Z=Rcpp::as<arma::mat>(Zlist(i));
    Vgt=Rcpp::as<arma::mat>(sigmahatlist(i));
    ZKZt = Z*K*Z.t();
    V=V+kron(ZKZt,Vgt);
  }
  loglik=dmvnorm_arma(arma::vectorise(Y),arma::vectorise(X*B),V,true, false);
  return loglik;
}





//////////////////////////////SIGMA Functions

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat diagSig_cpp(const arma::vec & params, const arma::mat & data) {
  arma::mat V=diagmat(exp(params));
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]

arma::mat IdentSig_cpp(const arma::vec & params, const arma::mat & data) {
  arma::mat V=exp(params(0))*eye(as_scalar(data),as_scalar(data));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat UnstrKronIdentSig_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d1vec(d1);
  for (unsigned int i = 0; i < d1; i++) {
    d1vec(i)=exp(params(ii));
    ii=ii+1;
  }
  arma::mat V=kron(diagmat(d1vec)*D1*diagmat(d1vec),eye(d2,d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat IdentKronUnstrSig_cpp(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D2(d2,d2);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d2; i++) {
    for (unsigned int j = i; j < d2; j++) {
      if (i==j) {
        D2(i,j)=1;
      }
      else {
        D2(i,j)=(2/M_PI)*atan(params(ii));
        D2(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d2);
  for (unsigned int i = 0; i < d2; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat V=kron(eye(d1,d1),diagmat(d2vec)*D2*diagmat(d2vec));
  return V;
}




// [[Rcpp::depends("RcppArmadillo")]]

arma::mat FA1hetSig_cpp(const arma::vec & params, const arma::mat & data) {
  int dim1=as_scalar(data);
  arma::vec d1=params.subvec(0,(dim1-1));
  arma::vec d2=params.subvec((dim1),(2*dim1-1));
  arma::mat V=d1*d1.t()+diagmat(exp(d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat FA1homSig_cpp(const arma::vec & params, const arma::mat & data) {
  int dim1=as_scalar(data);
  arma::vec d1=params.subvec(0,(dim1-1));
  double d2=params(dim1);
  arma::mat V=d1*d1.t()+exp(d2)*eye(dim1,dim1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmhomSig_cpp(const arma::vec & params, const arma::mat  & data) {
  double sigma = exp(params(0));
  double rho = (2/M_PI)*atan(params(1));
  int dim =as_scalar(data);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < as_scalar(data); i++) {
    V(i,i)=1;
  }
  return sigma*V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat compsymmhetSig_cpp(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =as_scalar(data);
  arma::mat V=rho*ones(dim,dim);
  arma::vec d1=params.subvec(1,(dim));
  return diagmat(d1)*V*diagmat(d1);
}





// [[Rcpp::depends("RcppArmadillo")]]

arma::mat FAhetSig_cpp(const arma::vec & params, const arma::mat  & data) {
  arma::mat Lambda=zeros(data(0,0),data(0,1));
  unsigned int ii=0;
  for (unsigned int i = 0; i < data(0,1); i++) {
    for (unsigned int j = i; j < data(0,0); j++) {
      Lambda(j,i)=params(ii);
      ii=ii+1;
    }
  }
  arma::vec d2(data(0,0));
  for (unsigned int j = 0; j < data(0,0); j++) {
    d2(j)=params(ii);
    ii=ii+1;
  }
  arma::mat V=Lambda*Lambda.t()+diagmat(exp(d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]

arma::mat FAhomSig_cpp(const arma::vec & params, const arma::mat  & data) {
  arma::mat Lambda=zeros(data(0,0),data(0,1));
  unsigned int ii=0;
  for (unsigned int i = 0; i < data(0,1); i++) {
    for (unsigned int j = i; j < data(0,0); j++) {
      Lambda(j,i)=params(ii);
      ii=ii+1;
    }
  }
  double d2=params(ii);
  arma::mat V=Lambda*Lambda.t()+exp(d2)*eye(data(0,0),data(0,0));
  return V;
}



XPtr<funcPtr> putFunPtrInXPtrSigma(std::string fstr) {
  if (fstr == "compsymmhom")
    return(XPtr<funcPtr>(new funcPtr(&compsymmhomSig_cpp)));
  else if (fstr == "compsymmhet")
    return(XPtr<funcPtr>(new funcPtr(&compsymmhetSig_cpp)));
  else if (fstr == "diag")
    return(XPtr<funcPtr>(new funcPtr(&diagSig_cpp)));
  else if (fstr == "Ident")
    return(XPtr<funcPtr>(new funcPtr(&IdentSig_cpp)));
  else if (fstr == "UnstrKronIdent")
    return(XPtr<funcPtr>(new funcPtr(&UnstrKronIdentSig_cpp)));
  else if (fstr == "IdentKronUnstr")
    return(XPtr<funcPtr>(new funcPtr(&IdentKronUnstrSig_cpp)));
  else if (fstr == "FA1het")
    return(XPtr<funcPtr>(new funcPtr(&FA1hetSig_cpp)));
  else if (fstr == "FA1hom")
    return(XPtr<funcPtr>(new funcPtr(&FA1homSig_cpp)));
  else if (fstr == "FAhet")
    return(XPtr<funcPtr>(new funcPtr(&FAhetSig_cpp)));
  else if (fstr == "FAhom")
    return(XPtr<funcPtr>(new funcPtr(&FAhomSig_cpp)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}



arma::mat callViaStringSigma(const arma::vec params, const arma::mat data, std::string funname) {
  XPtr<funcPtr> xpfun = putFunPtrInXPtrSigma(funname);
  funcPtr fun = *xpfun;
  arma::mat y = fun(params, data);
  return (y);
}





/////////////////////////////////////////////////////////////


// [[Rcpp::depends("RcppArmadillo")]]


double minimfuncreml(double delta, const arma::vec & eta, const arma::vec & lambda,const int & n,const int & q) {
  int df = n - q;
  double reml=df* log(sum(pow(eta,2)/(lambda + exp(delta)))) + sum(log(lambda + exp(delta)));
  //Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  // Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  
  return -reml;
}



// [[Rcpp::depends("RcppArmadillo")]]


double minimfuncremlfn1(double x, void * params)
{
  Rcpp::List * paramsk = (Rcpp::List*)(params);
  arma::vec  eta =as<arma::vec>(paramsk[0]);
  arma::vec  lambda=as<arma::vec>(paramsk[1]);
  int  n=as<int>(paramsk[2]);
  int  q=as<int>(paramsk[3]);
  int df = n - q;
  double reml=df*log(sum(pow(eta,2)/(lambda + exp(x)))) + sum(log(lambda + exp(x)));
  //Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  // Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  
  return -reml;
 
}






double zerofuncreml(double delta, arma::vec eta, arma::vec lambda, int n, int q) {
  int df = n - q;
  double remlzero=.5*df * (sum(pow(eta,2)/pow(lambda + exp(delta),2)))/(sum(pow(eta,2)/(lambda + exp(delta)))) - .5*sum(1/(lambda + exp(delta)));
  
  return remlzero;
}

double minimfuncml(double delta, arma::vec eta, arma::vec lambda, arma::vec epsilon, int n, int q) {
  double loglik=n * log(sum(pow(eta,2)/(lambda + exp(delta)))) + sum(log(epsilon + exp(delta)));
  loglik = -0.5 * (loglik + n + n * log(2 * M_PI/n));
  // Rcpp::Rcout  << "loglik:" << std::endl <<loglik << std::endl;
  
  return -loglik;
}


double zerofuncml(double delta, arma::vec eta, arma::vec lambda, arma::vec epsilon, int n, int q) {
  double loglikzero=.5*(n )* (sum(pow(eta,2)/pow(lambda + exp(delta),2)))/(sum(pow(eta,2)/(lambda + exp(delta)))) - .5*sum(1/(epsilon + exp(delta)));
  // Rcpp::Rcout  << "loglikzero:" << std::endl <<loglikzero << std::endl;
  
  return loglikzero;
}





// [[Rcpp::depends("RcppArmadillo")]]

arma::vec glominremlbrent(double lower, double upper, double tolparconv, unsigned int maxiter, arma::vec eta, arma::vec lambda, int n, int q)
{
  double a = lower;
  double b = upper;
  double fa =zerofuncreml(a, eta,  lambda, n, q);	// calculated now to save function calls
  double fb = zerofuncreml(b, eta,  lambda, n, q);	// calculated now to save function calls
  double fs = 0;		// initialize
  
  if (!(fa * fb < 0))
  {
    arma::vec yandx(2);
    yandx(0)=-11;
    yandx(1)=-11;
    
    return yandx;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < maxiter; ++iter)
  {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tolparconv)
    {
      arma::vec yandx(2);
      yandx(0)=fs;
      yandx(1)=s;
      
      return yandx;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tolparconv) ) ||
         ( !mflag && (std::abs(c-d) < tolparconv))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = zerofuncreml(s, eta,  lambda, n, q);	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s)
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for
  
  arma::vec yandx(2);
  yandx(0)=fs;
  yandx(1)=s;
  
  return yandx;
} // end brents_fun

// [[Rcpp::depends("RcppArmadillo")]]

List emmreml_arma_reml(const arma::colvec & y,const arma::mat & X,arma::mat  Z, arma::mat & K,  const double & tolparconv=1e-7, const int maxiter=1000, bool geterrors=false, bool Hinv=false) {
  
  int q = X.n_cols;
  int n = X.n_rows;
  arma::mat spI = eye(n,n);
  
  arma::mat S = spI - X* (X.t()* X).i()*X.t();
  double offset = log(static_cast<double>(n));
  arma::mat ZK = Z * K;
  arma::mat ZKZt = ZK*Z.t();
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  arma::vec yandx(2);
  double m;
  
  yandx=glominremlbrent(-(1e+2)/3,(1e+2)/3,tolparconv, maxiter,eta, lambda,  n, q);
  m=yandx(1);
  //  Rcpp::Rcout  << "m:" << std::endl << m << std::endl;
  
  arma::mat Hinvhat =(ZKZt + exp(m) * spI).i();
  arma::mat XtHinvhat = X.t()*Hinvhat;
  arma::mat betahat =solve(XtHinvhat*X, XtHinvhat*y);
  arma::vec ehat =y -X*betahat;
  arma::mat Hinvhatehat = Hinvhat * ehat;
  double sigmausqhat = sum(pow(eta,2)/(lambda + exp(m)))/(n - q);
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  double sigmaesqhat = exp(m) * sigmausqhat;
  arma::vec uhat = ZK.t()*Hinvhatehat;
  int df = n - q;
  double reml=(n - q) * log(sum(pow(eta,2)/(lambda + exp(m)))) + sum(log(lambda + exp(m)));
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Vinv*X);
    arma::vec Varuhat=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    arma::vec PEV=diagvec(sigmausqhat*K)-Varuhat;
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
  } else {
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("Hinv") = Hinvhat);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml);
      return outaa;
    }
  }
}




// [[Rcpp::depends("RcppArmadillo")]]

arma::vec glominmlbrent(double lower, double upper, double tolparconv, unsigned int maxiter, arma::vec eta, arma::vec lambda, arma::vec epsilon,int n, int q)
{
  double a = lower;
  double b = upper;
  double fa =zerofuncml(a, eta,  lambda, epsilon,n, q);	// calculated now to save function calls
  double fb = zerofuncml(b, eta,  lambda, epsilon,n, q);// calculated now to save function calls
  double fs = 0;		// initialize
  
  if (!(fa * fb < 0))
  {
    arma::vec yandx(2);
    yandx(0)=-11;
    yandx(1)=-11;
    
    return yandx;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < maxiter; ++iter)
  {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tolparconv)
    {
      arma::vec yandx(2);
      yandx(0)=fs;
      yandx(1)=s;
      
      return yandx;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tolparconv) ) ||
         ( !mflag && (std::abs(c-d) < tolparconv))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = zerofuncml(s, eta,  lambda, epsilon,n, q);	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s)
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for
  
  arma::vec yandx(2);
  yandx(0)=fs;
  yandx(1)=s;
  
  return yandx;
} // end brents_fun

// [[Rcpp::depends("RcppArmadillo")]]


List emmreml_arma_ml(const arma::colvec & y,const arma::mat & X, arma::mat  Z,  arma::mat  & K,  const double & tolparconv=1e-7, const int maxiter=1000, bool geterrors=false,  bool Hinv=false) {
  
  int q = X.n_cols;
  int n = X.n_rows;
  arma::mat spI = eye(n,n);
  
  arma::mat S = spI - X* (X.t()* X).i()*X.t();
  double offset = log(static_cast<double>(n));
  arma::mat ZK = Z * K;
  arma::mat ZKZt = ZK*Z.t();
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::vec epsilon;
  arma::eig_sym(epsilon,ZKZtandoffset);
  epsilon=epsilon-offset;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  arma::vec yandx(2);
  double m;
  yandx=glominmlbrent(-(1e+2)/3,(1e+2)/3,tolparconv, maxiter,eta, lambda, epsilon, n, q);
  m=yandx(1);
  arma::mat Hinvhat=(ZKZt + exp(m) * spI).i();
  arma::mat XtHinvhat = X.t()*Hinvhat;
  arma::mat betahat =solve(XtHinvhat*X, XtHinvhat*y);
  arma::vec ehat =y -X*betahat;
  arma::mat Hinvhatehat = Hinvhat * ehat;
  double sigmausqhat = sum(pow(eta,2)/(lambda + exp(m)))/(n - q);
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  double sigmaesqhat = exp(m) * sigmausqhat;
  arma::vec uhat = ZK.t()*Hinvhatehat;
  double loglik=n * log(sum(pow(eta,2)/(lambda + exp(m)))) + sum(log(epsilon + exp(m)));
  
  loglik = -0.5 * (loglik + n + n * log(2 * M_PI/n));
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Vinv*X);
    arma::vec Varuhat=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    arma::vec PEV=diagvec(sigmausqhat*K)-Varuhat;
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV);
      return outaa;
    }
  } else {
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik);
      return outaa;
    }
  }
}



// [[Rcpp::depends("RcppArmadillo")]]


List emm(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const bool & REML=true, const double & tolparconv=1e-7, const int maxiter=1000, bool geterrors=false, bool Hinv=false) {
  Rcpp::List out;
  arma::mat Ztemp;
  Ztemp=as<arma::mat>(Zlist(0));
  arma::mat Ktemp=Rcpp::as<arma::mat>(Klist(0));
  
  if (REML) {
    out=emmreml_arma_reml(y,X,Ztemp ,Ktemp,tolparconv, maxiter, geterrors, Hinv);
  }
  else {
    out=emmreml_arma_ml(y,X,Ztemp,Ktemp, tolparconv, maxiter,geterrors, Hinv);
  }
  
  return out;
}
//////////////////////////////////////////////

// [[Rcpp::depends("RcppArmadillo")]]

double dFunc(const arma::vec & params,const arma::mat & data, const std::string & funname,const  arma::vec & y, const arma::mat & X, const arma::mat & Zt, const arma::mat & Whalf,  int n,int q, int p,const arma::mat & ZtW , const arma::mat & XtWX,const arma::vec & XtWy, const arma::mat & ZtWX, const arma::vec & ZtWy, bool REML, double tolparinv=1e-10) {
  arma::mat K=exp(params(0))*callViaString(params(span(1, params.n_elem-1)),data,funname);
  arma::mat Lambdat=diagmat(pow(diagvec(K)+tolparinv,.5));
  Lambdat=cholsammupper(K,tolparinv);
  
  arma::vec b(q);                 // conditional mode of random effects
  arma::vec beta(p);              // conditional estimate of fixed-effects
  arma::vec cu(q);                // intermediate solution
  arma::mat DD=XtWX;                      // down-dated XtWX
  // stored here b/c x slot will be updated
  arma::mat L=eye(Lambdat.n_cols,Lambdat.n_cols);
  L=cholsammlower(Lambdat*ZtW*(Lambdat*ZtW).t()+eye(Lambdat.n_cols,Lambdat.n_cols),tolparinv);
  arma::vec mu(n);                // conditional mean of response
  arma::vec u(q);                 // conditional mode of spherical random effects
  cu=vectorise(solve(trimatl(L), Lambdat * ZtWy));
  
  arma::mat RZX=(solve(trimatl(L), Lambdat * ZtWX));
  DD = (XtWX - RZX.t()*RZX);
  
  // conditional estimate of fixed-effects coefficients (solve eqn. 33)
  beta=vectorise(solve(DD, XtWy - RZX.t()* cu));
  
  // conditional mode of the spherical random-effects coefficients (eqn. 34)
  u= vectorise(solve(trimatu(L.t()), cu - RZX* beta));
  
  b= vectorise(Lambdat.t()*u);
  
  
  // conditional mean of the response
  mu=vectorise(Zt.t()*b + X * beta);
  arma::vec  wtres= Whalf*(y-mu);      // weighted residuals
  double pwrss = as_scalar(sum(pow(wtres,2)) + sum(pow(u,2))); // penalized, weighted residual sum-of-squares
  int fn = mu.n_elem;
  double ld;
  double sign;
  
  log_det(ld, sign, L);
  ld=2*ld; // log determinant
  
  if (REML) {
    double ld2;
    log_det(ld2, sign, DD);
    ld=ld + ld2;
    fn =fn - beta.n_elem;
  }
  // profiled deviance or REML criterion
  ld=ld + fn*(1 + log(2*M_PI*pwrss) - log(static_cast<double>(fn)));
  return ld;
}

// [[Rcpp::depends("RcppArmadillo")]]

List outFunc(const arma::vec & params,const arma::mat & data, const std::string & funname, const arma::vec & y, const arma::mat & X, const arma::mat & Zt, const arma::mat & Whalf,  int n,int q, int p,const arma::mat & ZtW , const arma::mat & XtWX,const arma::vec & XtWy, const arma::mat & ZtWX,const arma::vec & ZtWy, bool REML, double tolparinv=1e-10, bool geterrors=false, bool Hinv=false ) {
//params(0) is lambda=sigmau/sigmae
  
  arma::mat K=exp(params(0))*callViaString(params(span(1, params.n_elem-1)),data,funname);
  arma::mat Lambdat=diagmat(pow(diagvec(K)+tolparinv,.5));
  Lambdat=cholsammupper(K,tolparinv);
  
  arma::vec b(q);                 // conditional mode of random effects
  arma::vec beta(p);              // conditional estimate of fixed-effects
  arma::vec cu(q);                // intermediate solution
  arma::mat DD=XtWX;                      // down-dated XtWX
  // stored here b/c x slot will be updated
  arma::mat L=eye(Lambdat.n_cols,Lambdat.n_cols);
  L=cholsammlower(Lambdat*ZtW*(Lambdat*ZtW).t()+eye(Lambdat.n_cols,Lambdat.n_cols),tolparinv);
  arma::vec mu(n);                // conditional mean of response
  arma::vec u(q);                 // conditional mode of spherical random effects
  cu=vectorise(solve(trimatl(L), Lambdat * ZtWy));
  
  arma::mat RZX=(solve(trimatl(L), Lambdat * ZtWX));
  DD = (XtWX - RZX.t()*RZX);
  
  // conditional estimate of fixed-effects coefficients (solve eqn. 33)
  beta=vectorise(solve(DD, XtWy - RZX.t()* cu));
  
  // conditional mode of the spherical random-effects coefficients (eqn. 34)
  u= vectorise(solve(trimatu(L.t()), cu - RZX* beta));
  
  b= vectorise(Lambdat.t()*u);
  
  
  // conditional mean of the response
  mu=vectorise(Zt.t()*b + X * beta);
  arma::vec  wtres= Whalf*(y-mu);      // weighted residuals
  double pwrss = as_scalar(sum(pow(wtres,2)) + sum(pow(u,2))); // penalized, weighted residual sum-of-squares
  int fn = mu.n_elem;
  double ld;
  double sign;
  
  log_det(ld, sign, L);
  ld=2*ld; // log determinant
  
  if (REML) {
    double ld2;
    log_det(ld2, sign, DD);
    ld=ld + ld2;
    fn =fn - beta.n_elem;
  }
  // profiled deviance or REML criterion
  ld=ld + fn*(1 + log(2*M_PI*pwrss) - log(static_cast<double>(fn)));
  
  ////////////////////////////////
  
  double sigmaesqhat =  pwrss/n;
  double sigmausqhat =exp(params(0))*sigmaesqhat;
  arma::mat ZK=Zt.t()*K/sigmausqhat;
  arma::mat Hhat=(ZK*Zt/sigmausqhat+(sigmaesqhat/sigmausqhat)*Whalf*Whalf.t());
  arma::mat Hinvhat=Hhat.i();
  arma::vec ehat =y-X*beta;
  arma::mat Vhat = (sigmausqhat) * Hhat;
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  arma::vec paramsK=params.subvec(span(1, params.n_elem-1));
  Rcpp::List corfuncparamslist=List::create(paramsK,NULL);
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Vinv*X);
    arma::vec Varuhat=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    arma::vec PEV=diagvec(K)-Varuhat;
    double loglik=dmvnorm_arma(y, X*beta,Vinv, true, true);
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = beta,Rcpp::Named("uhat") = b, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV, Rcpp::Named("corfuncparamslist") =corfuncparamslist);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = beta,Rcpp::Named("uhat") = b, Rcpp::Named("loglik") = loglik, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV, Rcpp::Named("corfuncparamslist") =corfuncparamslist);
      return outaa;
    }
  } else {
    double loglik=dmvnorm_arma(y, X*beta,Vinv,true, true);
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = beta,Rcpp::Named("uhat") = b, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("corfuncparamslist") =corfuncparamslist);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = beta,Rcpp::Named("uhat") = b, Rcpp::Named("loglik") = loglik, Rcpp::Named("corfuncparamslist") =corfuncparamslist);
      return outaa;
    }
  }
}



// [[Rcpp::depends("RcppArmadillo")]]

// /List simplex_minfuncforcormv(arma::vec start,const double tolparconv,const  int  maxiter, const arma::mat data, std::string funname, arma::mat Omega, arma::mat Z, arma::mat ZKZtlisti,arma::vec e, double sigmasq)
List simplex_dfunc(arma::vec start,const double  & tolparconv,const  int & maxiter, const arma::mat & data, const std::string & funname, const arma::vec & y, const arma::mat & X, const arma::mat & Zt, const arma::mat & Whalf,  int n,int q, int p,const arma::mat & ZtW ,const  arma::mat & XtWX,const arma::vec & XtWy, const arma::mat & ZtWX,const arma::vec & ZtWy, bool REML, double tolparinv=1e-10)
{
  
  const  double scale=1;
  int npar=start.n_elem;
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    
    // minfuncforcor(const arma::vec  params,const arma::mat data, std::string funname, arma::mat Omegainv, arma::mat Z,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambda)
    //minfuncforcor(as_scalar(vvec),,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda);
    if (funname=="Const") {
      f[j] =dFunc(join_cols(vvec,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
    }
    else {
      f[j] =dFunc(vvec,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
    }
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    
    if (funname=="Const") {
      fr =dFunc(join_cols(vr,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
    }
    else {
      fr =dFunc(vr,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
    }
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      
      
      if (funname=="Const") {
        fe =dFunc(join_cols(ve,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
      }
      else {
        fe =dFunc(ve,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
      }
      
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        
        
        if (funname=="Const") {
          fc =dFunc(join_cols(vc,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        else {
          fc =dFunc(vc,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        
        if (funname=="Const") {
          fc =dFunc(join_cols(vc,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        else {
          fc =dFunc(vc,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        
        if (funname=="Const") {
          f[vg] =dFunc(join_cols(vvec,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        else {
          f[vg] =dFunc(vvec,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        ksimplex++;
        
        
        vvec=(v.row(vh)).t();
        if (funname=="Const") {
          f[vh] =dFunc(join_cols(vvec,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        else {
          f[vh] =dFunc(vvec,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
        }
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  if (funname=="Const") {
    minval =dFunc(join_cols(vvec,zeros(1,1)),data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
  }
  else {
    minval =dFunc(vvec,data,funname,y, X, Zt, Whalf,n,q,  p,ZtW , XtWX, XtWy,ZtWX,ZtWy,REML,tolparinv);
  }
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}



// [[Rcpp::depends("RcppArmadillo")]]

List dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false, bool Hinv=false) {
  arma::mat Zt=Z.t();
  int n= y.n_elem;
  int p =X.n_cols;
  int q =Zt.n_rows;
  arma::mat Whalf=symMroot(W*R*W.t());
  arma::mat  WX = Whalf*X;
  arma::vec Wy = Whalf *y;
  arma::mat  ZtW = Zt * Whalf;
  arma::mat  XtWX = WX.t()*WX;
  arma::vec XtWy = WX.t()*Wy;
  arma::mat ZtWX= ZtW * WX;
  arma::vec ZtWy = ZtW * Wy;
  std::string Kfunctiontemp= Rcpp::as<std::string>(Kfunc(0));
  arma::vec Kparamstemp= Rcpp::as<arma::vec>(Kfunc(1));
  arma::vec corfuncparamslist=join_cols(zeros(1,1),Kparamstemp);
  arma::mat Kdatatemp= Rcpp::as<arma::mat>(Kfunc(2));
  List out;
  if (Kfunctiontemp=="Const") {
    out=simplex_dfunc(corfuncparamslist(span(0,corfuncparamslist.n_elem-2)),tolparconv,maxiter,Kdatatemp, Kfunctiontemp, y, X, Zt, Whalf, n,q,p, ZtW , XtWX,XtWy,ZtWX,ZtWy, REML, tolparinv);
  } else {
    out=simplex_dfunc(corfuncparamslist,tolparconv,maxiter,Kdatatemp, Kfunctiontemp, y, X, Zt, Whalf, n,q,p, ZtW , XtWX,XtWy,ZtWX,ZtWy, REML, tolparinv);
  }
  arma::vec out1 =Rcpp::as<arma::vec>(out(1));
  
  if (Kfunctiontemp=="Const") {
    out1=join_cols(out1,zeros(1,1));
  }
  
  List outaaa=outFunc(out1,Kdatatemp, Kfunctiontemp, y, X,Zt, Whalf,n,q,p,ZtW ,XtWX,XtWy,ZtWX,ZtWy,REML,tolparinv,geterrors, Hinv);
  return outaaa;
}



////////////////////////////


// [[Rcpp::depends("RcppArmadillo")]]


List mm(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool Hinv=false) {
  /////////////
  int k=2;
  arma::vec dimvec(2);
  int dimen=0;
  int q = X.n_cols;
  arma::vec betahat(q);
  
  int n = X.n_rows;
  arma::vec ehat(n);
  arma::vec ytrans(n);
  arma::mat Xtrans(n,q);
  arma::mat  Ktemp;
  arma::mat multiplier;
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n,n);
  arma::mat Omegainv(n,n);
  arma::vec sigmasqhatvec=ones(k);
  Rcpp::List ZKZtlist(k);
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
  }
  
  arma::mat A = as<arma::mat>(ZKZtlist(0));
  arma::mat B= as<arma::mat>(ZKZtlist(1));
  arma::mat C= cholsammlower(B, tolparinv);
  arma::mat G=C.i()*A*(C.i()).t();
  
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, G+std::log10(static_cast<double>(G.n_cols))*eye(G.n_cols,G.n_cols));
  eigval=eigval-std::log10(static_cast<double>(G.n_cols));
  arma::mat S=diagmat(sqrt(eigval))*eigvec.t()*G.i();
  arma::mat U=S.t();
  vec UBUdiag=1/sqrt(diagvec(S*B*U));
  U=(U*diagmat(UBUdiag));
  ytrans=U.t()*y;
  Xtrans=U.t()*X;
  
  
  arma::vec wvec(n);
  unsigned int itercountouter=0;
  do {
    itercountouter=itercountouter+1;
    for(unsigned int i=0; i<n; i++) {
      wvec(i)=1/(sigmasqhatvec(0)*as_scalar(eigval(i))+sigmasqhatvec(1));
    }
    Omegainv=diagmat(wvec);
    
    betahat=solve(Xtrans.t()*Omegainv*Xtrans, Xtrans.t()*(wvec%ytrans));
    ehat=ytrans-Xtrans*betahat;
    multiplier=(ehat.t()*diagmat(wvec%eigval%wvec)*ehat)/sum(wvec%eigval);
    sigmasqhatvec(0)=sigmasqhatvec(0)*powf(multiplier(0),.5);
    multiplier=(ehat.t()*diagmat(wvec%wvec)*ehat)/sum(wvec);
    sigmasqhatvec(1)=sigmasqhatvec(1)*powf(multiplier(0),.5);
    
    
  }  while (itercountouter<maxiter && powf(powf(1-sqrt(multiplier(0)), 2),.5)>tolparconv);
  
  
  ////////////////////////////////
  double sigmausqhat =sigmasqhatvec(0);
  double sigmaesqhat = sigmasqhatvec(1);
  arma::mat Hhat=(as<arma::mat>(ZKZtlist(0))+(sigmaesqhat/sigmausqhat)*as<arma::mat>(ZKZtlist(1)));
  arma::mat Hinvhat=Hhat.i();
  ehat =y-X*betahat;
  arma::mat Hinvhatehat = Hinvhat * ehat;
  arma::mat Vhat = (sigmausqhat) * Hhat;
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  arma::mat ZK=as<arma::mat>(Zlist(0))*as<arma::mat>(Klist(0));
  arma::vec uhat = ZK.t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv_sympdsamm(X.t()*Vinv*X, tolparinv);
    arma::vec Varuhat=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    arma::vec PEV=diagvec(sigmausqhat*as<arma::mat>(Klist(0)))-Varuhat;
    double loglik=dmvnorm_arma(y, X*betahat,Vinv, true, true);
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Covbetahat") = Covbetahat,Rcpp::Named("PEV") = PEV);
      return outaa;
    }
  } else {
    double loglik=dmvnorm_arma(y, X*betahat,Vinv,true, true);
    if (Hinv) {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat);
      return outaa;
    }
    else {
      List outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik);
      return outaa;
    }
  }
  
  ///////////
}

/////////////////////////////////////////////////////////////


double minfuncforcor(const arma::vec & params, const arma::mat & data, const std::string & funname, const arma::mat & Omega, const arma::mat & Z, const arma::mat & ZKZtlisti,const arma::vec & e, const double & sigmasq) {
  arma::mat Omegatemp=Omega-sigmasq*ZKZtlisti;
  arma::mat K=callViaString(params,data, funname);
  //arma::mat Kinv=callViaStringinv(paramstemp,data, funname);
  //arma::mat Kdiff=callViaStringdiff(paramstemp,data, funname);
  Omegatemp=Omegatemp+sigmasq*Z*K*Z.t();
  //Rcpp::Rcout  << det(Omegatemp) << std::endl;
  double funcmin;
  
  
  funcmin=-dmvnorm_arma(e, zeros(e.n_elem,1), Omegatemp, true, false);
  
  //Rcpp::Rcout  << funcmin << std::endl;
  return funcmin;
}







// [[Rcpp::depends("RcppArmadillo")]]


List simplex_minfuncforcormv(arma::vec  start,const double & tolparconv,const  int  & maxiter, const arma::mat & data, const std::string & funname, const arma::mat & Omega, const arma::mat & Z, const arma::mat & ZKZtlisti,const arma::vec & e, const double & sigmasq)
{
  double scale=1;
  int npar=start.n_elem;
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    
    
    f[j] =  minfuncforcor(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    
    fr = minfuncforcor(vr,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      
      
      fe = minfuncforcor(ve,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        
        
        fc= minfuncforcor(vc,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        
        fc= minfuncforcor(vc,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        
        f[vg] =minfuncforcor(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
        
        
        vvec=(v.row(vh)).t();
        f[vh] =minfuncforcor(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  
  minval =minfuncforcor(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}




// [[Rcpp::depends("RcppArmadillo")]]


List dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0, bool Hinv=false) {
  /////////////
  
  int k=Zlist.size();
  List corfuncparamslist(k+1);
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  unsigned int ii=0;
  unsigned int iirow=0;
  unsigned int iicol=0;
  arma::colvec betahat(q);
  int n = X.n_rows;
  arma::vec ehat(n);
  arma::mat  Ktempinv;
  arma::mat  Ktempdiff;
  arma::vec Kparamstemp;
  
  std::string Kfunctiontemp;
  
  arma::mat Kdatatemp;
  arma::mat  Ktemp;
  arma::mat  Ztemp;
  arma::mat Ztemp0;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n,n);
  arma::mat Omegainv(n,n);
  arma::mat P(n,n);
  arma::vec logsigmahatvec(k+1);
  double logsigma=log(var(y)/(k+1));
  Rcpp::List ZKZtlist(k+1);
  //Rcpp::Rcout  << 9999992 << std::endl;
  for(unsigned int i=0; i<k; i++) {
    if (corfunc(i)) {
      // Rcpp::Rcout  << 1999999 << std::endl;
      Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
      Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(Klist(i))(1));
      corfuncparamslist(i)=Kparamstemp;
      Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
      // Rcpp::Rcout  << 9999993 << std::endl;
      
      ////Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
    } else {
      Ktemp=as<arma::mat>(Klist(i));
    }
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp0=as<arma::mat>(Zlist(i));
    Ztemp=Ztemp0;
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    //  Rcpp::Rcout  << 1999999 << std::endl;
    logsigmahatvec(i)=logsigma;
    //  Rcpp::Rcout  << 9999994 << std::endl;
  }
  
  if (corfunc(k)) {
    //  Rcpp::Rcout  << 9999995 << std::endl;
    Kfunctiontemp= Rcpp::as<std::string>(R(0));
    Kparamstemp= Rcpp::as<arma::vec>(R(1));
    corfuncparamslist(k)=Kparamstemp;
    Kdatatemp= Rcpp::as<arma::mat>(R(2));
    Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
    //  Rcpp::Rcout  << 99996 << std::endl;
    
  } else {
    Ktemp=as<arma::mat>(R(0));
  }
  Ztemp=W;
  ZKZtlist(k)=Ztemp*Ktemp*Ztemp.t();
  
  logsigmahatvec(k)=logsigma;
  arma::vec logsigmahatvec0(k);
  // Rcpp::Rcout  << 997 << std::endl;
  
  arma::vec dlogLvec(k+1);
  arma::mat AImat(k+1,k+1);
  
  Rcpp::List ZKZtdifflist(k+1);
  arma::mat ZKZtdiff(0,n);
  List brentout(2);
  unsigned int itercountouter=0;
  do {
    logsigmahatvec0=logsigmahatvec;
    itercountouter=itercountouter+1;
    for(unsigned int i=0; i<k; i++) {
      
      ZKZtdiff=arma::mat(0,n);
      if (corfunc(i)) {
        
        Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
        Kparamstemp =  Rcpp::as<arma::vec>(corfuncparamslist(i));
        Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
        Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
        Ztemp0=as<arma::mat>(Zlist(i));
        Ztemp=Ztemp0;
        ZKZtdiff=join_cols(ZKZtdiff,exp(logsigmahatvec(i))*Ztemp*Ktemp*Ztemp.t());
      } else {
        Ztemp0=as<arma::mat>(Zlist(i));
        Ztemp=Ztemp0;
        Ktemp=as<arma::mat>(Klist(i));
        ZKZtdiff=join_cols(ZKZtdiff,exp(logsigmahatvec(i))*Ztemp*Ktemp*Ztemp.t());
      }
      ZKZtdifflist(i)= ZKZtdiff;
    }
    ZKZtdiff=arma::mat(0,n);
    
    if (corfunc(k)) {
      
      Kfunctiontemp= Rcpp::as<std::string>(R(0));
      Kparamstemp =  Rcpp::as<arma::vec>(corfuncparamslist(k));
      Kdatatemp= Rcpp::as<arma::mat>(R(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
      
      Ztemp=W;
      ZKZtdiff=join_cols(ZKZtdiff,exp(logsigmahatvec(k))*Ztemp*Ktemp*Ztemp.t());
    } else {
      Ktemp=as<arma::mat>(R(0));
      Ztemp=W;
      ZKZtdiff=join_cols(ZKZtdiff,exp(logsigmahatvec(k))*Ztemp*Ktemp*Ztemp.t());
    }
    
    
    ZKZtdifflist(k)= ZKZtdiff;
    
    
    
    Omega=zeros(n,n);
    for(unsigned int i=0; i<(k+1); i++) {
      Omega=Omega+exp(logsigmahatvec(i))*as<arma::mat>(ZKZtlist(i));
    }
    Omegainv=inv_sympdsamm(Omega,1e-10);
    P=Omegainv-Omegainv*X*inv_sympdsamm(X.t()*Omegainv*X, 1e-10)*X.t()*Omegainv;
    
    //dlogLvec
    //
    ii=0;
    for(unsigned int i=0; i<(k+1); i++) {
      dlogLvec(ii)=-.5*trace(P*Rcpp::as<arma::mat>(ZKZtdifflist(i)).rows(0,n-1));
      dlogLvec(ii)=dlogLvec(ii)+.5*as_scalar(y.t()*P*Rcpp::as<arma::mat>(ZKZtdifflist(i)).rows(0,n-1)*P*y);
      ii=ii+1;
    }
    if (methodopt!=0) {
      //AImat
      iirow=0;
      for(unsigned int irow=0; irow<(k+1); irow++) {
        for(unsigned int jrow=0; jrow<(Rcpp::as<arma::mat>(ZKZtdifflist(irow)).n_rows/n); jrow++) {
          iicol=0;
          for(unsigned int icol=0; icol<(k+1); icol++) {
            for(unsigned int jcol=0; jcol<(Rcpp::as<arma::mat>(ZKZtdifflist(icol)).n_rows/n); jcol++) {
              if (methodopt==1) {
                AImat(iirow, iicol)=.5*trace(P*Rcpp::as<arma::mat>(ZKZtdifflist(irow)).rows(jrow*n,(jrow+1)*n-1)*P*Rcpp::as<arma::mat>(ZKZtdifflist(icol)).rows(jcol*n,(jcol+1)*n-1));
              }
              if (methodopt==2) {
                AImat(iirow, iicol)=.5*as_scalar(y.t()*P*Rcpp::as<arma::mat>(ZKZtdifflist(irow)).rows(jrow*n,(jrow+1)*n-1)*P*Rcpp::as<arma::mat>(ZKZtdifflist(icol)).rows(jcol*n,(jcol+1)*n-1)*P*y);
              }
              iicol=iicol+1;
            }
            
          }
          iirow=iirow+1;
        }
        
      }
      
    }
    
    ///update current parameter vector
    if (methodopt==0) {
      logsigmahatvec=logsigmahatvec+dlogLvec;
    }
    else {
      logsigmahatvec=logsigmahatvec+solve(AImat+eye(AImat.n_cols,AImat.n_cols),dlogLvec);
    }
    
    
    
    if (itercountouter % 10==0) {
      for(unsigned int i=0; i<(k+1); i++) {
        
        if (exp(logsigmahatvec(i))>1e-10) {
          ZKZt=as<arma::mat>(ZKZtlist(i));
          
          //////////////////////Here estimate the covfunc parameter, recalculate things that depend on K
          //Rcpp::Rcout  << sigmasqhatvec(i) << std::endl;
          ////   // Rcpp::Rcout  << "check:" << std::endl << "4.1" << std::endl;
          
          Omegainv=inv_sympdsamm(Omega,tolparinv);
          betahat=solve(X.t()*Omegainv*X, X.t()*Omegainv*y);
          ehat=y-X*betahat;
          if (corfunc(i)) {
            if (!corfuncfixed(i)) {
              if (exp(logsigmahatvec(i))>1e-10) {
                if (i<k) {
                  Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
                  Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
                  Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
                } else {
                  Kfunctiontemp= Rcpp::as<std::string>(R(0));
                  Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
                  Kdatatemp= Rcpp::as<arma::mat>(R(2));
                }
                ////   // Rcpp::Rcout  << "check:" << std::endl << "5" << std::endl;
                
                //Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
                
                if (i<k) {
                  Ztemp0=as<arma::mat>(Zlist(i));
                  Ztemp=Ztemp0;
                  ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
                } else {
                  Ztemp=W;
                  ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
                }
                
                ////   // Rcpp::Rcout  << "check:" << std::endl << "5" << std::endl;
                
                brentout=simplex_minfuncforcormv(Kparamstemp,tolparconv,maxiter,Kdatatemp, Kfunctiontemp, Omega,Ztemp,ZKZt,ehat, exp(logsigmahatvec(i)));
                Kparamstemp=Rcpp::as<arma::vec>(brentout(1));
                corfuncparamslist(i)=Kparamstemp;
                
                //////////////recalculate parts etc,,,,
                ZKZtlist(i)= Ztemp*callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp)*Ztemp.t();
                
              }
            }
          }
        }
      }
    }
    
    //Rcpp::Rcout  << dlogLvec<< std::endl;
    //Rcpp::Rcout  << paramsvec<< std::endl;
  }  while (itercountouter<maxiter && sum(abs(logsigmahatvec0-logsigmahatvec))>tolparconv);
  //Rcpp::Rcout  << AImat<< std::endl;
  
  
  /////////
  arma::vec weights=exp(logsigmahatvec.subvec(0,k-1))/sum(exp(logsigmahatvec.subvec(0,k-1)));
  arma::mat Hinvhat =Omegainv*sum(exp(logsigmahatvec.subvec(0,k-1)));
  arma::mat XtHinvhat = X.t()*Hinvhat;
  
  ehat=y-X*betahat;
  
  arma::mat Hinvhatehat = Hinvhat * ehat;
  
  
  
  List ZKList(k);
  
  for(unsigned int i=0; i<k; i++) {
    if (corfunc(i)) {
      if (!corfuncfixed(i)) {
        Ztemp0=as<arma::mat>(Zlist(i));
        Ztemp=Ztemp0;
        Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
        Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
        Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
        Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
        ZKList(i)=Ztemp*Ktemp;
      } else {
        Ztemp0=as<arma::mat>(Zlist(i));
        Ztemp=Ztemp0;
        Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
        Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(Klist(i))(1));
        Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
        Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
        ZKList(i)=Ztemp*Ktemp;
      }
    }
    else {
      Ktemp=as<arma::mat>(Klist(i));
      Ztemp0=as<arma::mat>(Zlist(i));
      Ztemp=Ztemp0;
      ZKList(i)=Ztemp*Ktemp;
    }
  }
  
  arma::mat ZK(n,dimen, fill::zeros);
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=as<arma::mat>(ZKList(i));
    idx = idx + dimvec(i) ;
  }
  
  
  idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=weights(i)*ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 );
    idx = idx + dimvec(i) ;
  }
  
  
  arma::vec uhat = (ZK).t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Omegainv*X);
    arma::vec PEV=diagvec(pow(sum(exp(logsigmahatvec.subvec(0,k-1))),2)*ZK.t()*Omegainv*(ZK-X*Covbetahat*X.t()*Omegainv*ZK));
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("Ve") =exp(logsigmahatvec(k)),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = exp(logsigmahatvec.subvec(0,k-1))/sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("corfuncparamslist")=corfuncparamslist,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("Ve") =exp(logsigmahatvec(k)),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = exp(logsigmahatvec.subvec(0,k-1))/sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("corfuncparamslist")=corfuncparamslist,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      
    }
    double sigmasqu=outaa["Vu"];
    arma::vec weights=as<arma::vec>(outaa["weights"]);
    
    arma::vec PEVinsert;
    for (unsigned int i=0; i<k; i++) {
      PEVinsert=join_cols(PEVinsert,sigmasqu*weights(i)*diagvec(as<arma::mat>(Klist(i))));
      
    }
    
    outaa["PEV"] = PEVinsert-PEV;
    
    
    return outaa;
  } else {
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu")=sum(exp(logsigmahatvec.subvec(0,k-1))) ,Rcpp::Named("Ve") =exp(logsigmahatvec(k)),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = exp(logsigmahatvec.subvec(0,k-1))/sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("corfuncparamslist")=corfuncparamslist);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu")=sum(exp(logsigmahatvec.subvec(0,k-1))) ,Rcpp::Named("Ve") =exp(logsigmahatvec(k)),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = exp(logsigmahatvec.subvec(0,k-1))/sum(exp(logsigmahatvec.subvec(0,k-1))),Rcpp::Named("corfuncparamslist")=corfuncparamslist);
      return outaa;
    }
  }////////////
}

///////////


/////////////////////////////////////////////////////////////


// [[Rcpp::depends("RcppArmadillo")]]


List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist,const  arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10 ,const double & tolparinv=1e-10,const int & maxiter=200, bool geterrors=false, const double & lambda=0, bool Hinv=false) {
  /////////////
  
  int k=Zlist.size();
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  arma::vec betahat(q);
  
  int n = X.n_rows;
  arma::vec ehat(n);
  arma::mat  Ktemp;
  arma::mat multiplier;
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n,n);
  double Omega0tr;
  arma::mat Omegainv(n,n);
  arma::vec Omegainvehat(n);
  arma::vec sigmasqhatvec=ones(k+1);
  Rcpp::List ZKZtlist(k+1);
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    // Rcpp::Rcout  << "ZKZtlist:" << std::endl << Ztemp*Ktemp*Ztemp.t() << std::endl;
    
  }
  
  ZKZtlist(k)=W*R*W.t();
  
  //Rcpp::Rcout  << "ZKZtlist:" << std::endl << as<arma::mat>(ZKZtlist(0)) << std::endl;
  //Rcpp::Rcout  << "ZKZtlist:" << std::endl << as<arma::mat>(ZKZtlist(1)) << std::endl;
  unsigned int itercountouter=0;
  do {
    itercountouter=itercountouter+1;
    Omega=zeros(n,n);
    for(unsigned int i=0; i<(k+1); i++) {
      Omega=Omega+sigmasqhatvec(i)*as<arma::mat>(ZKZtlist(i));
    }
    Omegainv=inv_sympdsamm(Omega,1e-10);
    betahat=solve(X.t()*Omegainv*X, X.t()*Omegainv*y);
    ehat=y-X*betahat;
    Omegainvehat=Omegainv*ehat;
    for(unsigned int i=0; i<(k+1); i++) {
      if (sigmasqhatvec(i)>1e-10) {
        ZKZt=as<arma::mat>(ZKZtlist(i));
        Omega0tr=sum(vectorise(Omegainv%ZKZt));
        multiplier=Omegainvehat.t()*ZKZt*Omegainvehat;
        
        if (i<k) {
          //Rcpp::Rcout  << lambda << std::endl;
          
          multiplier=multiplier/(Omega0tr+(-std::log10(1-lambda)));
          //  Rcpp::Rcout  <<-std::log10(1-lambda) << std::endl;
          
        } else {
          multiplier=multiplier/(Omega0tr+1e-10);
        }
        
        
        if (multiplier(0)-1<0) {
          multiplier(0)=1-1.5*(1-multiplier(0));
        }
        if (multiplier(0)-1>0) {
          multiplier(0)=1+1.5*(multiplier(0)-1);
        }
        
        sigmasqhatvec(i)=sigmasqhatvec(i)*powf(multiplier(0),.5);
      }
    }
    //Rcpp::Rcout  << "sigmasqhatvec:" << std::endl << sigmasqhatvec << std::endl;
    
  }      while (itercountouter<maxiter && powf(powf(1-sqrt(multiplier(0)), 2),.5)>tolparconv);
  

  
  /////////
  arma::vec weights=sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1));
  arma::mat Hinvhat =Omegainv*sum(sigmasqhatvec.subvec(0,k-1));
  arma::mat XtHinvhat = X.t()*Hinvhat;
  
  ehat=y-X*betahat;
  
  arma::mat Hinvhatehat = Hinvhat * ehat;
  
  
  
  List ZKList(k);
  
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    Ztemp=as<arma::mat>(Zlist(i));
    ZKList(i)=Ztemp*Ktemp;
  }
  
  arma::mat ZK(n,dimen, fill::zeros);
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=as<arma::mat>(ZKList(i));
    idx = idx + dimvec(i) ;
  }
  
  
  idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=weights(i)*ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 );
    idx = idx + dimvec(i) ;
  }
  
  
  arma::vec uhat = (ZK).t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Omegainv*X);
    arma::vec PEV=diagvec(pow(sum(sigmasqhatvec.subvec(0,k-1)),2)*ZK.t()*Omegainv*(ZK-X*Covbetahat*X.t()*Omegainv*ZK));
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Ve")=sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Ve")=sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      
    }
    double sigmasqu=outaa["Vu"];
    arma::vec weights=as<arma::vec>(outaa["weights"]);
    
    arma::vec PEVinsert;
    for (unsigned int i=0; i<k; i++) {
      PEVinsert=join_cols(PEVinsert,sigmasqu*weights(i)*diagvec(as<arma::mat>(Klist(i))));
      
    }
    
    outaa["PEV"] = PEVinsert-PEV;
    
    
    return outaa;
  } else {
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu")=sum(sigmasqhatvec.subvec(0,k-1)) ,Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)));
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu")=sum(sigmasqhatvec.subvec(0,k-1)) ,Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)));
      return outaa;
    }
  }////////////
  
}

/////////////////////////////////////////////////////////////


double minimfuncremlemmmk(double delta, arma::vec eta, arma::vec lambda, int n, int q) {
  int df = n - q;
  double reml=df*log(sum(pow(eta,2)/(lambda + exp(delta)))) + sum(log(lambda + exp(delta)));
  //Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  // Rcpp::Rcout  << "reml:" << std::endl <<reml << std::endl;
  
  return -reml;
}

double zerofuncremlemmmk(double delta, arma::vec eta, arma::vec lambda, int n, int q) {
  int df = n - q;
  double remlzero=.5*df * (sum(pow(eta,2)/pow(lambda + exp(delta),2)))/(sum(pow(eta,2)/(lambda + exp(delta)))) - .5*sum(1/(lambda + exp(delta)));
  
  return remlzero;
}

double minimfuncmlemmmk(double delta, arma::vec eta, arma::vec lambda, arma::vec epsilon, int n, int q) {
  double loglik=n * log(sum(pow(eta,2)/(lambda + exp(delta)))) + sum(log(epsilon + exp(delta)));
  loglik = -0.5 * (loglik + n + n * log(2 * M_PI/n));
  // Rcpp::Rcout  << "loglik:" << std::endl <<loglik << std::endl;
  
  return -loglik;
}


double zerofuncmlemmmk(double delta, arma::vec eta, arma::vec lambda, arma::vec epsilon, int n, int q) {
  double loglikzero=.5*(n ) * (sum(pow(eta,2)/pow(lambda + exp(delta),2)))/(sum(pow(eta,2)/(lambda + exp(delta)))) - .5*sum(1/(epsilon + exp(delta)));
  // Rcpp::Rcout  << "loglikzero:" << std::endl <<loglikzero << std::endl;
  
  return loglikzero;
}



// [[Rcpp::depends("RcppArmadillo")]]

arma::vec glominremlbrentemmmk(double lower, double upper, double tolparconv, unsigned int maxiter, arma::vec eta, arma::vec lambda, int n, int q)
{
  double a = lower;
  double b = upper;
  double fa =zerofuncremlemmmk(a, eta,  lambda, n, q);	// calculated now to save function calls
  double fb = zerofuncremlemmmk(b, eta,  lambda, n, q);	// calculated now to save function calls
  double fs = 0;		// initialize
  
  if (!(fa * fb < 0))
  {
    arma::vec yandx(2);
    yandx(0)=-11;
    yandx(1)=-11;
    
    return yandx;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < maxiter; ++iter)
  {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tolparconv)
    {
      arma::vec yandx(2);
      yandx(0)=fs;
      yandx(1)=s;
      
      return yandx;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tolparconv) ) ||
         ( !mflag && (std::abs(c-d) < tolparconv))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = zerofuncremlemmmk(s, eta,  lambda, n, q);	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s)
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for
  
  arma::vec yandx(2);
  yandx(0)=fs;
  yandx(1)=s;
  
  return yandx;
} // end brents_fun
// [[Rcpp::depends("RcppArmadillo")]]

double minimfunctionouter_reml(arma::vec  weights,const arma::vec  &y,
                               int  k,arma::vec  dimvec,int  dimen,int  q,
                               int  n,arma::mat  & spI,arma::mat  &S,double  offset,
                               List  &ZKZtList,List  &ZKList,
                               const double  tolparconv=1e-9, const int  maxiter=1000) {
  
  
  arma::mat ZKZt(n,n, fill::zeros);
  
  
  for(unsigned int i=0; i<k; i++) {
    ZKZt=ZKZt+weights(i)*as<arma::mat>(ZKZtList(i));
  }
  
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  
  arma::vec yandx(2);
  double m;
  yandx=glominremlbrentemmmk(-(1e+2)/3,(1e+2)/3, tolparconv, maxiter,eta, lambda,  n, q);
  m=yandx(1);
  int df = n - q;
  double reml=(n - q) * log(sum(pow(eta,2)/(lambda + exp(m)))) + sum(log(lambda + exp(m)));
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  
  return -reml;
}


// [[Rcpp::depends("RcppArmadillo")]]

List minimfunctionouter_reml2(arma::vec  weights,const arma::vec  &y,const arma::mat  &X,
                              int  k,arma::vec  dimvec,int  dimen,int  q,
                              int  n,arma::mat  & spI,arma::mat  &S,double  offset,
                              List  &ZKZtList,arma::mat  ZK,
                              const double  tolparconv=1e-9, const int  maxiter=1000, bool geterrors=false, bool Hinv=false) {
  //Rcpp::Rcout  << "locy" << std::endl <<&y << std::endl;
  
  arma::mat ZKZt(n,n, fill::zeros);
  
  
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZKZt=ZKZt+weights(i)*as<arma::mat>(ZKZtList(i));
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=weights(i)*ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 );
    idx = idx + dimvec(i) ;
  }
  
  
  
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  
  arma::vec yandx(2);
  double m;
  yandx=glominremlbrentemmmk(-(1e+2)/3,(1e+2)/3, tolparconv, maxiter,eta, lambda,  n, q);
  m=yandx(1);
  int df = n - q;
  double reml=(n - q) * log(sum(pow(eta,2)/(lambda + exp(m)))) + sum(log(lambda + exp(m)));
  reml = -0.5 * (reml + df + df * log(2 * M_PI/df));
  
  arma::mat Hinvhat =(ZKZt + exp(m) * spI).i();
  arma::mat XtHinvhat = X.t()*Hinvhat;
  arma::mat betahat =solve(XtHinvhat*X, XtHinvhat*y);
  
  arma::vec ehat=y-X*betahat;
  
  arma::mat Hinvhatehat = Hinvhat * ehat;
  double sigmausqhat = sum(pow(eta,2)/(lambda + exp(m)))/(n - q);
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  double sigmaesqhat = exp(m) * sigmausqhat;
  arma::vec uhat = (ZK).t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Vinv*X);
    arma::vec PEV=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = weights,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("weights") = weights,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
  } else {
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = weights);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("reml") = reml, Rcpp::Named("weights") = weights);
      return outaa;
    }
  }
}


// [[Rcpp::depends("RcppArmadillo")]]


List simplex_reml(arma::vec start,int npar, double scale,const arma::vec  &y,int  k,arma::vec  dimvec,int  dimen,int  q,
                  int  n,arma::mat  &spI,arma::mat  &S,double  offset,
                  List  &ZKZtList,List  &ZKList,
                  const double  tolparconv=1e-9, const int  maxiter=1000, const bool constrain=true )
{
  
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  if (constrain== true) {
    if (min(v.row(j))<0) {
      v.row(j)=(v.row(j)-min(v.row(j)))/sum((v.row(j)-min(v.row(j))));
    } else {
      v.row(j)=(v.row(j))/sum((v.row(j)));
    }
  }
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    f[j] = minimfunctionouter_reml(vvec,y,
                                   k,dimvec,dimen,q,
                                   n,spI,S,offset,
                                   ZKZtList, ZKList,
                                   tolparconv, maxiter);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    if (constrain== true) {
      if (min(vr)<0) {
        vr=(vr-min(vr))/sum((vr-min(vr)));
      } else {
        vr=(vr)/sum((vr));
      }
    }
    
    fr = minimfunctionouter_reml(vr,y,
                                 k,dimvec,dimen,q,
                                 n,spI,S,offset,
                                 ZKZtList, ZKList,
                                 tolparconv, maxiter);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      if (constrain== true) {
        if (min(ve)<0) {
          ve=(ve-min(ve))/sum((ve-min(ve)));
        } else {
          ve=(ve)/sum((ve));
        }
      }
      
      fe = minimfunctionouter_reml(ve,y,
                                   k,dimvec,dimen,q,
                                   n,spI,S,offset,
                                   ZKZtList, ZKList,
                                   tolparconv, maxiter);
      
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        if (constrain== true) {
          if (min(vc)<0) {
            vc=(vc-min(vc))/sum((vc-min(vc)));
          } else {
            vc=(vc)/sum((vc));
          }
        }
        
        fc= minimfunctionouter_reml(vc,y,
                                    k,dimvec,dimen,q,
                                    n,spI,S,offset,
                                    ZKZtList, ZKList,
                                    tolparconv, maxiter);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        if (constrain== true) {
          if (min(vc)<0) {
            vc=(vc-min(vc))/sum((vc-min(vc)));
          } else {
            vc=(vc)/sum((vc));
          }
        }
        fc= minimfunctionouter_reml(vc,y,
                                    k,dimvec,dimen,q,
                                    n,spI,S,offset,
                                    ZKZtList, ZKList,
                                    tolparconv, maxiter);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        if (constrain== true) {
          if (min(v.row(vg))<0) {
            v.row(vg)=(v.row(vg)-min(v.row(vg)))/sum((v.row(vg)-min(v.row(vg))));
          } else {
            v.row(vg)=(v.row(vg))/sum((v.row(vg)));
          }
        }
        f[vg] =minimfunctionouter_reml(vvec,y,
                                       k,dimvec,dimen,q,
                                       n,spI,S,offset,
                                       ZKZtList, ZKList
                                         ,tolparconv, maxiter);
        ksimplex++;
        
        if (constrain== true) {
          if (min(v.row(vh))<0) {
            v.row(vh)=(v.row(vh)-min(v.row(vh)))/sum((v.row(vh)-min(v.row(vh))));
          } else {
            v.row(vh)=(v.row(vh))/sum((v.row(vh)));
          }
        }
        vvec=(v.row(vh)).t();
        f[vh] =minimfunctionouter_reml(vvec,y,
                                       k,dimvec,dimen,q,
                                       n,spI,S,offset,
                                       ZKZtList, ZKList
                                         ,tolparconv, maxiter);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  if (sum(vvec)<=0 || std::isnan(vvec(0))) {
    vvec=ones(k)/k;
  }
  minval =minimfunctionouter_reml(vvec,y,
                                  k,dimvec,dimen,q,
                                  n,spI,S,offset,
                                  ZKZtList, ZKList
                                    ,tolparconv, maxiter);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}




// [[Rcpp::depends("RcppArmadillo")]]

List emmremlmk_arma(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool geterrors=false, bool Hinv=false) {
  //Rcpp::Rcout  << "locy" << std::endl <<&y << std::endl;
  int k=Zlist.length();
  arma::vec dimvec(k);
  int dimen=0;
  arma::mat  Ktemp;
  arma::mat  Ztemp;
  int q = X.n_cols;
  int n = X.n_rows;
  arma::mat spI = eye(n,n);
  
  arma::mat S = spI - X* (X.t()* X).i()*X.t();
  double offset = std::log(static_cast<double>(n));
  List ZKZtList(k);
  List ZKList(k);
  ///////////
  
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtList(i)=Ztemp*Ktemp*Ztemp.t();
    ZKList(i)=Ztemp*Ktemp;
  }
  
  arma::mat ZK(n,dimen, fill::zeros);
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=as<arma::mat>(ZKList(i));
    idx = idx + dimvec(i) ;
  }
  
  arma::vec startvec=ones(k)/k;
  List minimout=simplex_reml(startvec,k,.1,y,
                             k,dimvec,dimen,q,
                             n,spI,S,offset,
                             ZKZtList, ZKList,
                             tolparconv, maxiter,true);
  arma::vec vvec=minimout(1);
  List outaa=minimfunctionouter_reml2(vvec,y,X,
                                      k,dimvec,dimen,q,
                                      n,spI,S,offset,
                                      ZKZtList,  ZK,
                                      tolparconv,maxiter, geterrors, Hinv);
  
  if (geterrors==true) {
    double sigmasqu=outaa["Vu"];
    arma::vec weights=as<arma::vec>(outaa["weights"]);
    arma::vec PEV=as<arma::vec>(outaa["PEV"]);
    
    arma::vec PEVinsert;
    for (unsigned int i=0; i<k; i++) {
      PEVinsert=join_cols(PEVinsert,sigmasqu*weights(i)*diagvec(as<arma::mat>(Klist(i))));
      
    }
    
    
    outaa["PEV"] = PEVinsert-PEV;
    return outaa;
  } else {
    return outaa;
  }
}



//////////////////////////////////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]

arma::vec glominmlbrentemmmk(double lower, double upper, double tolparconv, unsigned int maxiter, arma::vec eta, arma::vec lambda, arma::vec epsilon, int n, int q)
{
  double a = lower;
  double b = upper;
  double fa =zerofuncmlemmmk(a, eta,  lambda, epsilon, n, q);	// calculated now to save function calls
  double fb = zerofuncmlemmmk(b, eta,  lambda, epsilon,  n, q);	// calculated now to save function calls
  double fs = 0;		// initialize
  
  if (!(fa * fb < 0))
  {
    arma::vec yandx(2);
    yandx(0)=-11;
    yandx(1)=-11;
    
    return yandx;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < maxiter; ++iter)
  {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tolparconv)
    {
      arma::vec yandx(2);
      yandx(0)=fs;
      yandx(1)=s;
      
      return yandx;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tolparconv) ) ||
         ( !mflag && (std::abs(c-d) < tolparconv))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = zerofuncmlemmmk(s, eta,  lambda, epsilon, n, q);	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s)
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for
  
  arma::vec yandx(2);
  yandx(0)=fs;
  yandx(1)=s;
  
  return yandx;
} // end brents_fun
// [[Rcpp::depends("RcppArmadillo")]]

double minimfunctionouter_ml(arma::vec  weights,const arma::vec  &y,
                             int  k,arma::vec  dimvec,int  dimen,int  q,
                             int  n,arma::mat  &spI,arma::mat  &S,double  offset,
                             List  &ZKZtList,List  &ZKList,
                             const double  tolparconv=1e-9, const int  maxiter=1000) {
  
  
  arma::mat ZKZt(n,n, fill::zeros);
  
  
  for(unsigned int i=0; i<k; i++) {
    ZKZt=ZKZt+weights(i)*as<arma::mat>(ZKZtList(i));
  }
  
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::vec epsilon;
  arma::eig_sym(epsilon,ZKZtandoffset);
  epsilon=epsilon-offset;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  
  arma::vec yandx(2);
  double m;
  yandx=glominmlbrentemmmk(-(1e+2)/3,(1e+2)/3, tolparconv, maxiter,eta, lambda,epsilon,  n, q);
  m=yandx(1);
  double ml=minimfuncmlemmmk(m,eta, lambda,epsilon, n,q);
  return ml;
}


// [[Rcpp::depends("RcppArmadillo")]]

List minimfunctionouter_ml2(arma::vec  weights,const arma::vec  &y,const arma::mat  &X,
                            int  k,arma::vec  dimvec,int  dimen,int  q,
                            int  n,arma::mat  &spI,arma::mat  &S,double  offset,
                            List  &ZKZtList,arma::mat  ZK,
                            const double  tolparconv=1e-9, const int  maxiter=1000, bool geterrors=false, bool Hinv=false) {
  
  arma::mat ZKZt(n,n, fill::zeros);
  
  
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZKZt=ZKZt+weights(i)*as<arma::mat>(ZKZtList(i));
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=weights(i)*ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 );
    idx = idx + dimvec(i) ;
  }
  
  
  
  arma::mat ZKZtandoffset = ZKZt + offset * spI;
  arma::vec epsilon;
  arma::eig_sym(epsilon,ZKZtandoffset);
  epsilon=epsilon-offset;
  arma::mat SZKZtSandoffset = (S*ZKZtandoffset)*S;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,SZKZtSandoffset);
  arma::mat Ur = Eigenvecs.cols(q, (n -1));
  arma::vec lambda = Eigenval.rows(q,(n -1)) - offset;
  arma::vec eta = Ur.t()*y;
  
  arma::vec yandx(2);
  double m;
  yandx=glominmlbrentemmmk(-(1e+2)/3,(1e+2)/3, tolparconv, maxiter,eta, lambda,  epsilon, n, q);
  m=yandx(1);
  
  
  double loglik=n * log(sum(pow(eta,2)/(lambda + exp(m)))) + sum(log(epsilon + exp(m)));
  
  loglik = -0.5 * (loglik + n + n * log(2 * M_PI/n));
  arma::mat Hinvhat =(ZKZt + exp(m) * spI).i();
  arma::mat XtHinvhat = X.t()*Hinvhat;
  arma::mat betahat =solve(XtHinvhat*X, XtHinvhat*y);
  
  arma::vec ehat=y-X*betahat;
  
  arma::mat Hinvhatehat = Hinvhat * ehat;
  double sigmausqhat = sum(pow(eta,2)/(lambda + exp(m)))/(n - q);
  arma::mat Vinv = (1/sigmausqhat) * Hinvhat;
  double sigmaesqhat = exp(m) * sigmausqhat;
  arma::vec uhat = (ZK).t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Vinv*X);
    arma::vec PEV=diagvec(pow(sigmausqhat,2)*ZK.t()*Vinv*(ZK-X*Covbetahat*X.t()*Vinv*ZK));
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = weights,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = weights,Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      return outaa;
    }
  } else {
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = weights);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sigmausqhat ,Rcpp::Named("Ve") =sigmaesqhat,Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = weights);
      return outaa;
    }
    
  }
}



// [[Rcpp::depends("RcppArmadillo")]]


List simplex_ml(arma::vec start,int npar, double scale,const arma::vec  &y,int  k,arma::vec  dimvec,int  dimen,int  q,
                int  n,arma::mat  &spI,arma::mat  &S,double  offset,
                List  &ZKZtList,List  &ZKList,
                const double  tolparconv=1e-9, const int  maxiter=1000, const bool constrain=true )
{
  
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  if (constrain== true) {
    if (min(v.row(j))<0) {
      v.row(j)=(v.row(j)-min(v.row(j)))/sum((v.row(j)-min(v.row(j))));
    } else {
      v.row(j)=(v.row(j))/sum((v.row(j)));
    }
  }
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    f[j] = minimfunctionouter_ml(vvec,y,
                                 k,dimvec,dimen,q,
                                 n,spI,S,offset,
                                 ZKZtList, ZKList,
                                 tolparconv, maxiter);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    if (constrain== true) {
      if (min(vr)<0) {
        vr=(vr-min(vr))/sum((vr-min(vr)));
      } else {
        vr=(vr)/sum((vr));
      }
    }
    
    fr = minimfunctionouter_ml(vr,y,
                               k,dimvec,dimen,q,
                               n,spI,S,offset,
                               ZKZtList, ZKList,
                               tolparconv, maxiter);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      if (constrain== true) {
        if (min(ve)<0) {
          ve=(ve-min(ve))/sum((ve-min(ve)));
        } else {
          ve=(ve)/sum((ve));
        }
      }
      
      fe = minimfunctionouter_ml(ve,y,
                                 k,dimvec,dimen,q,
                                 n,spI,S,offset,
                                 ZKZtList, ZKList,
                                 tolparconv, maxiter);
      
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        if (constrain== true) {
          if (min(vc)<0) {
            vc=(vc-min(vc))/sum((vc-min(vc)));
          } else {
            vc=(vc)/sum((vc));
          }
        }
        
        fc= minimfunctionouter_ml(vc,y,
                                  k,dimvec,dimen,q,
                                  n,spI,S,offset,
                                  ZKZtList, ZKList,
                                  tolparconv, maxiter);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        if (constrain== true) {
          if (min(vc)<0) {
            vc=(vc-min(vc))/sum((vc-min(vc)));
          } else {
            vc=(vc)/sum((vc));
          }
        }
        fc= minimfunctionouter_ml(vc,y,
                                  k,dimvec,dimen,q,
                                  n,spI,S,offset,
                                  ZKZtList, ZKList,
                                  tolparconv, maxiter);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        if (constrain== true) {
          if (min(v.row(vg))<0) {
            v.row(vg)=(v.row(vg)-min(v.row(vg)))/sum((v.row(vg)-min(v.row(vg))));
          } else {
            v.row(vg)=(v.row(vg))/sum((v.row(vg)));
          }
        }
        f[vg] =minimfunctionouter_ml(vvec,y,
                                     k,dimvec,dimen,q,
                                     n,spI,S,offset,
                                     ZKZtList, ZKList
                                       ,tolparconv, maxiter);
        ksimplex++;
        
        if (constrain== true) {
          if (min(v.row(vh))<0) {
            v.row(vh)=(v.row(vh)-min(v.row(vh)))/sum((v.row(vh)-min(v.row(vh))));
          } else {
            v.row(vh)=(v.row(vh))/sum((v.row(vh)));
          }
        }
        vvec=(v.row(vh)).t();
        f[vh] =minimfunctionouter_ml(vvec,y,
                                     k,dimvec,dimen,q,
                                     n,spI,S,offset,
                                     ZKZtList, ZKList
                                       ,tolparconv, maxiter);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  if (sum(vvec)<=0 || std::isnan(vvec(0))) {
    vvec=ones(k)/k;
  }
  minval =minimfunctionouter_ml(vvec,y,
                                k,dimvec,dimen,q,
                                n,spI,S,offset,
                                ZKZtList, ZKList
                                  ,tolparconv, maxiter);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}





// [[Rcpp::depends("RcppArmadillo")]]

List emmmlmk_arma(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool geterrors=false, bool Hinv=false) {
  int k=Zlist.length();
  arma::vec dimvec(k);
  int dimen=0;
  arma::mat  Ktemp;
  arma::mat  Ztemp;
  int q = X.n_cols;
  int n = X.n_rows;
  arma::mat spI = eye(n,n);
  
  arma::mat S = spI - X* (X.t()* X).i()*X.t();
  double offset = log(static_cast<double>(n));
  List ZKZtList(k);
  List ZKList(k);
  ///////////
  
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtList(i)=Ztemp*Ktemp*Ztemp.t();
    ZKList(i)=Ztemp*Ktemp;
  }
  
  arma::mat ZK(n,dimen, fill::zeros);
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=as<arma::mat>(ZKList(i));
    idx = idx + dimvec(i) ;
  }
  
  arma::vec startvec=ones(k)/k;
  List minimout=simplex_ml(startvec,k,.1,y,
                           k,dimvec,dimen,q,
                           n,spI,S,offset,
                           ZKZtList, ZKList,
                           tolparconv, maxiter,true);
  arma::vec vvec=minimout(1);
  List outaa=minimfunctionouter_ml2(vvec,y,X,
                                    k,dimvec,dimen,q,
                                    n,spI,S,offset,
                                    ZKZtList,  ZK,
                                    tolparconv,maxiter, geterrors, Hinv);
  if (geterrors==true) {
    double sigmasqu=outaa["Vu"];
    arma::vec weights=as<arma::vec>(outaa["weights"]);
    arma::vec PEV=as<arma::vec>(outaa["PEV"]);
    
    arma::vec PEVinsert;
    for (unsigned int i=0; i<k; i++) {
      PEVinsert=join_cols(PEVinsert,sigmasqu*weights(i)*diagvec(as<arma::mat>(Klist(i))));
      
    }
    
    outaa["PEV"] = PEVinsert-PEV;
    
    return outaa;
  } else {
    return outaa;
  }
  
}





// [[Rcpp::depends("RcppArmadillo")]]


List emmmk(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool REML=true, bool geterrors=false, bool Hinv=false) {
  Rcpp::List out;
  
  if (REML==true) {
    out=emmremlmk_arma(y,X, Zlist,Klist, tolparconv, maxiter, geterrors, Hinv);
    
  }
  else {
    out=emmmlmk_arma(y,X, Zlist,Klist, tolparconv, maxiter, geterrors, Hinv);
  }
  
  return out;
}


/////////////////////////////
// [[Rcpp::depends("RcppArmadillo")]]

arma::vec loglikfuncemmremlmv_arma(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const Rcpp::List & Vgtlist,const arma::mat & Vet,const arma::mat & B) {
  int k=Zlist.size();
  int n = Y.n_rows;
  arma::mat V=kron(eye(n,n),Vet);
  arma::mat Vgt=Vet;
  arma::mat Z;
  arma::mat K;
  arma::mat ZKZt;
  arma::vec loglik;
  for (unsigned int i = 0; i < k; i++) {
    K=Rcpp::as<arma::mat>(Klist(i));
    Z=Rcpp::as<arma::mat>(Zlist(i));
    Vgt=Rcpp::as<arma::mat>(Vgtlist(i));
    ZKZt = Z*K*Z.t();
    V=V+kron(ZKZt,Vgt);
  }
  loglik=dmvnorm_arma(arma::vectorise(Y),arma::vectorise(X*B),V,true, false);
  return loglik;
}


// [[Rcpp::depends("RcppArmadillo")]]


List emmmv(const arma::mat & Yt,const arma::mat & Xt,const arma::mat & Z, const arma::mat & K,const double & tolparconv=1e-8,const double tolparinv=1e-8,const int & maxiter=500,const bool & geterrors=false, bool Hinv=false) {
  arma::mat Y=Yt.t();
  arma::mat X=Xt.t();
  arma::mat  KZt = K*Z.t();
  arma::mat  ZKZt = Z*KZt;
  int N = K.n_rows;
  arma::mat Eigenvecs;
  arma::vec Eigenval;
  arma::eig_sym(Eigenval , Eigenvecs ,ZKZt+std::log10(static_cast<double>(ZKZt.n_cols))*eye(ZKZt.n_cols,ZKZt.n_cols));
  Eigenval=Eigenval-std::log10(static_cast<double>(ZKZt.n_cols));
  
  int n = ZKZt.n_rows;
  
  arma::mat Eigenvecs2(n,n);
  arma::vec Eigenval2(n);
  for (unsigned int i = 0; i < n; i++) {
    Eigenvecs2.col(i)=Eigenvecs.col(n-i-1);
    Eigenval2(i)=Eigenval(n-i-1);
  }
  int d = Y.n_rows;
  int q = X.n_rows;
  //   Rcpp::Rcout  << "convnum:" << std::endl << Eigenvecs2 << std::endl;
  
  arma::mat Ytrans = Y * Eigenvecs2;
  arma::mat Xtrans = X * Eigenvecs2;
  arma::mat Vgt = cov(Y.t())/2;
  arma::mat Vet = cov(Y.t())/2;
  arma::mat XttinvXtXtt = Xtrans.t() * inv_sympdsamm(Xtrans*Xtrans.t(), tolparinv);
  arma::mat Bt = Ytrans * XttinvXtXtt;
  arma::mat Vetm1(d,d);
  arma::mat Gt(d,n);
  arma::mat ytl(d,1);
  arma::mat xtl(q,1);
  arma::mat Vlt(d,d);
  arma::mat Sigmalt(d,d);
  arma::mat etl(d,1);
  double convnum=1e+20;
  arma::mat Vetsum(d,d);
  arma::mat Vetsumhalf(d,d);
  arma::mat Vgtsum(d,d);
  Rcpp::List LL(n);
  Rcpp::List LLi;
  arma::mat Sigmalti;
  int counter=0;
  
  do {
    counter=counter+1;
    Vetm1=Vet;
    Vetsum.zeros();
    Gt.zeros();
    for (unsigned int i = 0; i < n; i++) {
      Vlt=Eigenval2(i) * Vgt + Vet;
      Gt.col(i)=Eigenval2(i) * Vgt * arma::solve(Vlt+tolparinv*eye(Vlt.n_cols,Vlt.n_cols) , (Ytrans.col(i) - Bt * Xtrans.col(i)));
      Sigmalt=Eigenval2(i) * Vgt;
      Sigmalt=Sigmalt- (Sigmalt) * arma::solve(Vlt+tolparinv*eye(Vlt.n_cols,Vlt.n_cols) ,Sigmalt);
      LL(i)=Rcpp::List::create(Sigmalt);
    }
    
    Bt = (Ytrans - Gt) * XttinvXtXtt;
    Vgtsum.zeros();
    Vgt.zeros();
    for (unsigned int i = 0; i < n; i++) {
      LLi=Rcpp::as<Rcpp::List>(LL(i));
      Sigmalti=Rcpp::as<arma::mat>(LLi(0));
      Vgtsum = Vgtsum+(1/(Eigenval2(i)))*(Gt.col(i)*(Gt.col(i)).t()+Sigmalti);
    }
    Vgt=Vgtsum/n;
    
    Vetsum.zeros();
    Vet.zeros();
    //  Rcpp::Rcout  << "Vetsum:" << std::endl << Vetsum << std::endl;
    
    for (unsigned int i = 0; i < n; i++) {
      LLi=Rcpp::as<Rcpp::List>(LL(i));
      Sigmalti=Rcpp::as<arma::mat>(LLi(0));
      etl=Ytrans.col(i)- Bt * Xtrans.col(i)- Gt.col(i);
      Vetsumhalf=etl*etl.t();
      Vetsum=Vetsum+Vetsumhalf + Sigmalti;
    }
    Vet=Vetsum/n;
    //  Rcpp::Rcout  << "Vet:" << std::endl << Vet << std::endl;
    convnum=fabs(sum(diagvec(Vet - Vetm1)))/fabs(1+sum(diagvec(Vetm1)));
    //Rcpp::Rcout  << "convnum:" << std::endl << convnum << std::endl;
  } while ((convnum >tolparconv) && (counter<maxiter));
  Bt = (Ytrans - Gt) * XttinvXtXtt;
  unsigned int iterinv=0;
  do {
    iterinv=iterinv+1;
    Vet=Vet+tolparinv*eye(Vet.n_cols,Vet.n_cols);
  } while ((det(Vet)==0) && (iterinv<100)) ;
  
  iterinv=0;
  do {
    iterinv=iterinv+1;
    Vgt=Vgt+tolparinv*eye(Vgt.n_cols,Vgt.n_cols);
  } while ((det(Vgt)==0) && (iterinv<100)) ;
  
  arma::mat Hobs=kron(ZKZt, Vgt) + kron(eye(n,n),Vet);
  iterinv=0;
  do {
    iterinv=iterinv+1;
    Hobs=Hobs+tolparinv*eye(n*d,n*d);
  } while ((det(Hobs)==0)&& (iterinv<10));
  arma::mat HobsInv = inv_sympdsamm(Hobs,tolparinv);
  vec ehat = arma::vectorise(Y - Bt * X);
  
  arma::mat HobsInve = HobsInv * ehat;
  arma::mat varvecG = kron(K, Vgt);
  arma::mat gpred(d*n,1);
  arma::mat Zforvec=kron(Z, eye(d,d));
  
  gpred = varvecG * (Zforvec).t() * HobsInve;
  arma::mat ZKforvec=Zforvec * varvecG;
  arma::mat Gpred(d,N);
  Gpred= reshape(gpred, d, N);
  if (geterrors==true) {
    arma::mat Xforvec=kron(X.t(),eye(d,d));
    arma::mat varBhat = Xforvec.t()*HobsInv*Xforvec;
    iterinv=0;
    do {
      iterinv=iterinv+1;
      varBhat+tolparinv*eye(varBhat.n_cols,varBhat.n_cols);
    } while ((det(varBhat)==0)&& (iterinv<10));
    varBhat =inv_sympdsamm(varBhat, tolparinv);
    
    arma::mat P= HobsInv - HobsInv *Xforvec* varBhat*Xforvec.t()*HobsInv;
    
    arma::vec varGhat = diagvec(ZKforvec.t()* P* ZKforvec);
    
    arma::vec PEV= diagvec(varvecG) - varGhat;
    arma::vec loglik=loglikfuncemmremlmv_arma(Y.t(),X.t(),List::create(Rcpp::Named("Z") =Z), List::create(Rcpp::Named("K") =K),List::create(Rcpp::Named("Vgt") =Vgt),Vet,Bt.t());
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =Vgt,Rcpp::Named("Ve") =Vet,Rcpp::Named("betahat") =Bt, Rcpp::Named("Gpred")=Gpred,Rcpp::Named("loglik")=loglik,  Rcpp::Named("varBhat")=varBhat, Rcpp::Named("PEV")=PEV, Rcpp::Named("Hinv")=HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =Vgt,Rcpp::Named("Ve") =Vet,Rcpp::Named("betahat") =Bt, Rcpp::Named("Gpred")=Gpred,Rcpp::Named("loglik")=loglik,  Rcpp::Named("varBhat")=varBhat, Rcpp::Named("PEV")=PEV);
      return outaa;
    }
  } else {
    arma::vec loglik=loglikfuncemmremlmv_arma(Y.t(),X.t(),List::create(Rcpp::Named("Z") =Z), List::create(Rcpp::Named("K") =K),List::create(Rcpp::Named("Vgt") =Vgt),Vet,Bt.t());
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =Vgt,Rcpp::Named("Ve") =Vet,Rcpp::Named("betahat") =Bt, Rcpp::Named("Gpred")=Gpred,Rcpp::Named("loglik")=loglik, Rcpp::Named("Hinv")=HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =Vgt,Rcpp::Named("Ve") =Vet,Rcpp::Named("betahat") =Bt, Rcpp::Named("Gpred")=Gpred,Rcpp::Named("loglik")=loglik);
      return outaa;
    }
  }
}




///////////////////////////////////




// [[Rcpp::depends("RcppArmadillo")]]


List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200, bool geterrors=false, bool Hinv=false) {
  /////////////
  int k=Zlist.size();
  
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  int d = Y.n_cols;
  arma::colvec betahat(q*d);
  arma::mat Betahat(q, d);
  int n = X.n_rows;
  arma::vec ehat(n*d);
  arma::mat Ehat(n,d);
  arma::mat Rmat(n,d);
  arma::mat bigX=kron(eye(d,d), X);
  arma::mat  Ktemp;
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n*d,n*d);
  arma::mat Omegainv(n*d,n*d);
  arma::vec L1diag(d);
  arma::mat LL1t(d,d);
  arma::mat L1t(d,d);
  arma::mat L1tinv(d,d);
  arma::mat N1t;
  arma::vec L2diag(d);
  arma::mat LL2t(d,d);
  arma::mat L2t(d,d);
  arma::mat L2tinv(d,d);
  arma::mat N2t;
  List sigmahatlist(k+1);
  arma::mat sigma(d,d);
  sigma=cov(Y)/(k+1);
  arma::mat Onesn=ones(n,1);
  arma::mat IdOnesn=kron(eye(d,d),Onesn);
  
  Rcpp::List ZKZtlist(k+1);
  
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    
    sigmahatlist(i)=sigma;
  }
  
  ZKZtlist(k)=W*R*W.t();
  sigmahatlist(k)=sigma;
  
  //////////////////////
  ////insert1
  arma::mat A = as<arma::mat>(ZKZtlist(0));
  arma::mat B= as<arma::mat>(ZKZtlist(1));
  arma::mat C= cholsammlower(B, tolparinv);
  arma::mat G=C.i()*A*(C.i()).t();
  
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, G);
  arma::mat S=diagmat(sqrt(eigval))*eigvec.t()*G.i();
  arma::mat U=S.t();
  vec UBUdiag=1/sqrt(diagvec(S*B*U));
  U=(U*diagmat(UBUdiag));
  arma::mat Ytrans=U.t()*Y;
  arma::mat Xtrans=U.t()*X;
  arma::vec Dvec=eigval;
  
  
  ///insert1 till here
  arma::mat Ytrans2;
  arma::mat BigXtrans;
  arma::mat Gamma1(d,d);
  arma::mat Gamma2(d,d);
  arma::mat Psi;
  arma::mat Psit;
  arma::vec Lambdavec;
  
  vec PsitGamma2Psidiag;
  unsigned int itercountouter=0;
  
  do {
    sigma=as<arma::mat>(sigmahatlist(1));
    itercountouter=itercountouter+1;
    ///insert2
    Gamma1 =as<arma::mat>(sigmahatlist(0));
    Gamma2= as<arma::mat>(sigmahatlist(1));
    C=cholsammlower(Gamma2+tolparinv*eye(Gamma2.n_cols,Gamma2.n_cols), tolparinv);
    G=C.i()*Gamma1*(C.i()).t();
    eig_sym(eigval, eigvec, G);
    Psit=diagmat(sqrt(eigval))*eigvec.t()*G.i();
    Psi=Psit.t();
    PsitGamma2Psidiag=1/sqrt(diagvec(Psit*Gamma2*Psi));
    
    Psi=(Psi*diagmat(PsitGamma2Psidiag));
    
    Lambdavec=eigval;
    
    Ytrans2=Ytrans*Psi;
    BigXtrans=kron(Psi.t(),Xtrans);
    
    
    ///////////This part can be done much faster
    Omega=diagmat(kron(Lambdavec,Dvec)+1);
    
    Omegainv=diagmat(1/(kron(Lambdavec,Dvec)+1));
    
    betahat=solve(BigXtrans.t()*Omegainv*BigXtrans, BigXtrans.t()*Omegainv*reshape(Ytrans2,n*d,1));
    Betahat=reshape(betahat,q, d);
    
    
    
    
    for(unsigned int i=0; i<(d); i++) {
      L1diag(i)=sum(Dvec%(1/(eigval(i)*Dvec+1)));
      L2diag(i)=sum((1/(eigval(i)*Dvec+1)));
    }
    
    
    LL1t=Psi*diagmat(L1diag)*Psi.t();
    L1t=diagmat(pow(diagvec(LL1t)+tolparinv,.5));
    L1t=cholsammupper(LL1t,tolparinv);
    L1tinv=inv(trimatu(L1t+tolparinv*eye(LL1t.n_rows,LL1t.n_rows)));
    LL2t=Psi*diagmat(L2diag)*Psi.t();
    L2t=diagmat(pow(diagvec(LL2t)+tolparinv,.5));
    L2t=cholsammupper(LL2t,tolparinv);
    L2tinv=inv(trimatu(L2t+tolparinv*eye(LL2t.n_rows,LL2t.n_rows)));
    
    N1t=diagmat(sqrt(Dvec))*(((Ytrans-Xtrans*Betahat)*Psi)/(Dvec*Lambdavec.t()+1))*diagmat(Lambdavec)*Psi.i();
    N2t=(((Ytrans-Xtrans*Betahat)*Psi)/(Dvec*Lambdavec.t()+1))*Psi.i();
    sigmahatlist(0)=L1tinv*symMroot(L1t.t()*N1t.t()*N1t*L1t)*L1tinv.t();
    sigmahatlist(1)=L2tinv*symMroot(L2t.t()*N2t.t()*N2t*L2t)*L2tinv.t();
    // Rcpp::Rcout  << "check" << std::endl << 10<< std::endl;
    
    /////inset2 till here
  } while (itercountouter<maxiter && fabs(sum(diagvec(as<arma::mat>(sigmahatlist(k))-sigma)))/(fabs(1+sum(diagvec(sigma))))>tolparconv/n);
  
  
  arma::mat HobsInv;
  arma::vec eehat;
  arma::mat Omegainve;
  arma::mat varvecG;
  arma::mat gpred;
  List Gpredlist(k);
  List PEVlist(k);
  arma::mat Vet=as<arma::mat>(sigmahatlist(k));
  Omega=zeros(n*d,n*d);
  
  for(unsigned int i=0; i<(k+1); i++) {
    Omega=Omega+kron(as<arma::mat>(ZKZtlist(i)),as<arma::mat>(sigmahatlist(i)));
  }
  Omegainv=inv_sympdsamm(Omega,tolparinv);
  
  arma::mat Xforvec=kron(X,eye(d,d));
  arma::mat varBhat =inv(Xforvec.t()*Omegainv*Xforvec);
  
  int nk;
  for(unsigned int i=0; i<(k); i++) {
    
    nk=  as<arma::mat>(Klist(i)).n_cols;
    arma::vec eehat =arma::vectorise((Y-X*Betahat).t());
    Omegainve = Omegainv * eehat;
    varvecG = kron(as<arma::mat>(Klist(i)),as<arma::mat>(sigmahatlist(i)));
    gpred =  varvecG*(kron(as<arma::mat>(Zlist(i)).t(),eye(d,d)))*Omegainve;
    Gpredlist(i)= (reshape(gpred, d, nk)).t();
    
  }
  if (geterrors==true) {
    for(unsigned int i=0; i<(k); i++) {
      arma::mat P= Omegainv - Omegainv *Xforvec* varBhat*Xforvec.t()*Omegainv;
      arma::mat Zforvec=kron(as<arma::mat>(Zlist(i)), eye(d,d));
      arma::mat ZKforvec=Zforvec * varvecG;
      arma::vec varGhat = diagvec(ZKforvec.t()* P* ZKforvec);
      PEVlist(i)= diagvec(varvecG) - varGhat;
    }
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klist,sigmahatlist,Betahat, W, R);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("Betahat") =Betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("Hinv")=HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("Betahat") =Betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik);
      return outaa;
    }
  } else {
    
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klist,sigmahatlist,Betahat, W, R);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("Betahat") =Betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("Hinv")=HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("Betahat") =Betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik);
      return outaa;
    }
    
  }
  
}

/////////////////////////////////////////////////////////////




// [[Rcpp::depends("RcppArmadillo")]]


List mmmkmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200, bool geterrors=false, bool Hinv=false) {
  /////////////
  int k=Zlist.size();
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  int d = Y.n_cols;
  arma::colvec betahat(q*d);
  arma::mat Betahat(q, d);
  int n = X.n_rows;
  arma::vec ehat(n*d);
  arma::mat Ehat(n,d);
  arma::mat Rmat(n,d);
  arma::mat bigX=kron(eye(d,d), X);
  arma::mat  Ktemp;
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n*d,n*d);
  arma::mat Omegainv(n*d,n*d);
  arma::mat LLt(d,d);
  arma::mat Lt(d,d);
  arma::mat Ltinv(d,d);
  List sigmahatlist(k+1);
  arma::mat sigma(d,d);
  sigma=cov(Y)/(k+1);
  arma::mat Onesn=ones(n,1);
  arma::mat IdOnesn=kron(eye(d,d),Onesn);
  
  Rcpp::List ZKZtlist(k+1);
  
  ///////////
  for(unsigned int i=0; i<k; i++) {
    Ktemp=as<arma::mat>(Klist(i));
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    sigmahatlist(i)=sigma;
  }
  
  ZKZtlist(k)=W*R*W.t();
  sigmahatlist(k)=sigma;
  
  unsigned int itercountouter=0;
  do {
    sigma=as<arma::mat>(sigmahatlist(k));
    itercountouter=itercountouter+1;
    Omega=zeros(n*d,n*d);
    for(unsigned int i=0; i<(k+1); i++) {
      Omega=Omega+kron(as<arma::mat>(sigmahatlist(i)),as<arma::mat>(ZKZtlist(i)));
    }
    
    Omegainv=inv_sympdsamm(Omega,tolparinv);
    betahat=solve(bigX.t()*Omegainv*bigX, bigX.t()*Omegainv*reshape(Y,n*d,1));
    
    
    ehat=reshape(Y, n*d,1)-bigX*betahat;
    
    Ehat=reshape(ehat,n,d);
    
    Rmat=reshape(Omegainv*ehat,n,d);
    
    for(unsigned int i=0; i<(k+1); i++) {
      ZKZt=as<arma::mat>(ZKZtlist(i));
      LLt=IdOnesn.t()*(kron(ones(d,d),ZKZt)%Omegainv)*IdOnesn;
      Lt=diagmat(pow(diagvec(LLt)+tolparinv,.5));
      Lt=cholsammupper(LLt,tolparinv);
      Ltinv=inv(trimatu(Lt+tolparinv*eye(LLt.n_rows,LLt.n_rows)));
      // sigmahatlist(i)=(Ltinv*Lt*(as<arma::mat>(sigmahatlist(i))*Rmat.t()*ZKZt*Rmat*as<arma::mat>(sigmahatlist(i)))*Lt.t()*Ltinv.t());
      sigmahatlist(i)=Ltinv*symMroot(Lt*(as<arma::mat>(sigmahatlist(i))*Rmat.t()*ZKZt*Rmat*as<arma::mat>(sigmahatlist(i)))*Lt.t())*Ltinv.t();
    }
    
  }  while (itercountouter<maxiter && fabs(sum(diagvec(as<arma::mat>(sigmahatlist(k))-sigma)))/(fabs(1+sum(diagvec(sigma))))>tolparconv/n);
  
  Betahat=reshape(betahat,q, d);
  arma::mat HobsInv;
  arma::vec eehat;
  arma::mat Omegainve;
  arma::mat varvecG;
  arma::mat gpred;
  List Gpredlist(k);
  List PEVlist(k);
  arma::mat Vet=as<arma::mat>(sigmahatlist(k));
  Omega=zeros(n*d,n*d);
  
  for(unsigned int i=0; i<(k+1); i++) {
    Omega=Omega+kron(as<arma::mat>(ZKZtlist(i)),as<arma::mat>(sigmahatlist(i)));
  }
  Omegainv=inv_sympdsamm(Omega,tolparinv);
  
  arma::mat Xforvec=kron(X,eye(d,d));
  arma::mat varBhat =inv(Xforvec.t()*Omegainv*Xforvec);
  
  int nk;
  for(unsigned int i=0; i<(k); i++) {
    
    nk=  as<arma::mat>(Klist(i)).n_cols;
    arma::vec eehat =arma::vectorise((Y-X*Betahat).t());
    Omegainve = Omegainv * eehat;
    varvecG = kron(as<arma::mat>(Klist(i)),as<arma::mat>(sigmahatlist(i)));
    gpred =  varvecG*(kron(as<arma::mat>(Zlist(i)).t(),eye(d,d)))*Omegainve;
    Gpredlist(i)= (reshape(gpred, d, nk)).t();
    
  }
  
  if (geterrors==true) {
    for(unsigned int i=0; i<(k); i++) {
      arma::mat P= Omegainv - Omegainv *Xforvec* varBhat*Xforvec.t()*Omegainv;
      arma::mat Zforvec=kron(as<arma::mat>(Zlist(i)), eye(d,d));
      arma::mat ZKforvec=Zforvec * varvecG;
      arma::vec varGhat = diagvec(ZKforvec.t()* P* ZKforvec);
      PEVlist(i)= diagvec(varvecG) - varGhat;
    }
    
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klist,sigmahatlist,betahat.t(), W, R);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("Hinv") =HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik);
      return outaa;
    }
  } else {
    
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klist,sigmahatlist,betahat.t(), W, R);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("Hinv") =HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik);
      return outaa;
    }
    
  }
  
}

/////////////////////////////////////////////////////////////


double minfuncforcormmmk(const arma::vec  params, const arma::mat data, std::string funname, arma::mat Omega, arma::mat Z, arma::mat ZKZtlisti,arma::vec e, double sigmasq) {
  arma::mat Omegatemp=Omega-sigmasq*ZKZtlisti;
  arma::mat K=callViaString(params,data, funname);
  //arma::mat Kinv=callViaStringinv(paramstemp,data, funname);
  //arma::mat Kdiff=callViaStringdiff(paramstemp,data, funname);
  Omegatemp=Omegatemp+sigmasq*Z*K*Z.t();
  //Rcpp::Rcout  <<"check1"<< std::endl;
  double funcmin;
  
  //double dmvnorm_arma(const arma::mat & x, const arma::vec & mean, const arma::mat & sigma, bool log = false, bool covinv=false) {
  
  funcmin=-dmvnorm_arma(e, zeros(e.n_elem,1), Omegatemp, true, false);
  
  //Rcpp::Rcout  << funcmin << std::endl;
  return funcmin;
}







// [[Rcpp::depends("RcppArmadillo")]]


List simplex_minfuncforcormmmkmv(arma::vec start,const double tolparconv,const  int  maxiter, const arma::mat data, std::string funname, arma::mat Omega, arma::mat Z, arma::mat ZKZtlisti,arma::vec e, double sigmasq)
{
  double scale=1;
  int npar=start.n_elem;
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    
    
    f[j] =  minfuncforcormmmk(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    
    fr = minfuncforcormmmk(vr,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      
      
      fe = minfuncforcormmmk(ve,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        
        
        fc= minfuncforcormmmk(vc,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        
        fc= minfuncforcormmmk(vc,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        
        f[vg] =minfuncforcormmmk(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
        
        
        vvec=(v.row(vh)).t();
        f[vh] =minfuncforcormmmk(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  
  minval =minfuncforcormmmk(vvec,data,funname,Omega,Z,ZKZtlisti,e,sigmasq);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}


/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



// [[Rcpp::depends("RcppArmadillo")]]


List mmmkcorfuncmvopt(const arma::colvec & y,const arma::mat & X, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R,const double & tolparconv=1e-10,const double & tolparinv=1e-10,const int & maxiter=200, bool geterrors=false, double lambda=0, bool Hinv=false) {
  /////////////
  
  int k=Zlist.size();
  //Rcpp::Rcout  << 1 << std::endl;
  List corfuncparamslist(k+1);
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  arma::vec betahat(q);
  
  int n = X.n_rows;
  arma::vec ehat(n);
  arma::mat  Ktemp;
  arma::mat  Ktempinv;
  arma::mat  Ktempdiff;
  arma::vec Kparamstemp;
  std::string Kfunctiontemp;
  
  arma::mat Kdatatemp;
  arma::mat multiplier;
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat ZtOmegainv;
  arma::mat Omega(n,n);
  
  arma::mat Omegainv(n,n);
  arma::vec sigmasqhatvec=ones(k+1);
  Rcpp::List ZKZtlist(k+1);
  arma::vec r;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    if (corfunc(i)) {
      Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
      Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(Klist(i))(1));
      corfuncparamslist(i)=Kparamstemp;
      Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
      //Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
    } else {
      Ktemp=as<arma::mat>(Klist(i));
    }
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    //   //Rcpp::Rcout  << "check:" << std::endl << "1" << std::endl;
    
  }
  //Rcpp::Rcout  <<2 << std::endl;
  
  if (corfunc(k)) {
    Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(R(0))(0));
    Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(R(0))(1));
    corfuncparamslist(k)=Kparamstemp;
    Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(R(0))(2));
    Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
    //Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
  } else {
    Ktemp=as<arma::mat>(R(0));
  }
  
  ZKZtlist(k)=W*Ktemp*W.t();
  
  //Rcpp::Rcout  << "check:" << std::endl << "2" << std::endl;
  
  //Rcpp::Rcout  << "ZKZtlist:" << std::endl << as<arma::mat>(ZKZtlist(0)) << std::endl;
  //Rcpp::Rcout  << "ZKZtlist:" << std::endl << as<arma::mat>(ZKZtlist(1)) << std::endl;
  unsigned int itercountouter=0;
  List brentout(2);
  //Rcpp::Rcout  << 3 << std::endl;
  
  do {
    itercountouter=itercountouter+1;
    Omega=zeros(n,n);
    for(unsigned int i=0; i<(k+1); i++) {
      Omega=Omega+sigmasqhatvec(i)*as<arma::mat>(ZKZtlist(i));
    }
    //   // Rcpp::Rcout  << "itercountouter:" << std::endl << itercountouter << std::endl;
    
    try {
      Omegainv=inv_sympdsamm(Omega,tolparinv);
    }
    catch(...) {
      Omegainv=inv_sympdsamm(Omega,1e-7);
      
    }
    
    
    betahat=solve(X.t()*Omegainv*X, X.t()*Omegainv*y);
    ehat=y-X*betahat;
    for(unsigned int i=0; i<(k+1); i++) {
      if (sigmasqhatvec(i)>1e-10) {
        ZKZt=as<arma::mat>(ZKZtlist(i));
        multiplier=ehat.t()*Omegainv*ZKZt*Omegainv*ehat;
        
        if (i<k) {
          if (lambda>=1) {
            lambda=1-1e-10;
          }
          if (lambda<=0) {
            lambda=0+1e-10;
          }
          multiplier=multiplier/(trace(Omegainv*ZKZt)+(-std::log10(1-lambda))+1e-10);
        } else {
          multiplier=multiplier/(trace(Omegainv*ZKZt)+1e-10);
        }
        
        //   // Rcpp::Rcout  << "check:" << std::endl << "4" << std::endl;
        
     //   if (multiplier(0)-1<0) {
      //    multiplier(0)=1-1.5*(1-multiplier(0));
      //  }
      //  if (multiplier(0)-1>0) {
      //    multiplier(0)=1+1.5*(multiplier(0)-1);
      //  }
        
        sigmasqhatvec(i)=sigmasqhatvec(i)*powf(multiplier(0),.5);
        if (sigmasqhatvec(i)<=0) {
          sigmasqhatvec(i)=0;
        }
        //////////////////////Here estimate the covfunc parameter, recalculate things that depend on K
        //Rcpp::Rcout  << sigmasqhatvec(i) << std::endl;
        ////   // Rcpp::Rcout  << "check:" << std::endl << "4.1" << std::endl;
        
        if (itercountouter % 10==0) {
          if (corfunc(i)) {
            if (!corfuncfixed(i)) {
              if (sigmasqhatvec(i)>1e-10) {
                if (i<k) {
                  Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
                  Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
                  Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
                } else {
                  Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(R(0))(0));
                  Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
                  Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(R(0))(2));
                }
                ////   // Rcpp::Rcout  << "check:" << std::endl << "5" << std::endl;
                
                //Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
                
                if (i<k) {
                  Ztemp=Rcpp::as<arma::mat>(Zlist(i));
                  ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
                } else {
                  Ztemp=W;
                  ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
                }
                
                ////   // Rcpp::Rcout  << "check:" << std::endl << "5" << std::endl;
                
                brentout=simplex_minfuncforcormmmkmv(Kparamstemp,tolparconv,1000,Kdatatemp, Kfunctiontemp, Omega,Ztemp,ZKZt,ehat, sigmasqhatvec(i));
                Kparamstemp=Rcpp::as<arma::vec>(brentout(1));
                corfuncparamslist(i)=Kparamstemp;
                
                //////////////recalculate parts etc,,,,
                ZKZtlist(i)= Ztemp*callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp)*Ztemp.t();
                
              }
            }
          }
        }
      }
    }
    // Rcpp::Rcout  << "sigmasqhatvec:" << std::endl << multiplier(0) << std::endl;
    
  }  while (itercountouter<maxiter && powf(powf(1-sqrt(multiplier(0)), 2),.5)>tolparconv);
  
  
  
  /////////
  arma::vec weights=sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1));
  arma::mat Hinvhat =Omegainv*sum(sigmasqhatvec.subvec(0,k-1));
  arma::mat XtHinvhat = X.t()*Hinvhat;
  
  ehat=y-X*betahat;
  
  arma::mat Hinvhatehat = Hinvhat * ehat;
  
  
  List ZKList(k);
  
  for(unsigned int i=0; i<k; i++) {
    if (corfunc(i)) {
      Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
      Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
      Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
    } else {
      Ktemp=as<arma::mat>(Klist(i));
    }
    Ztemp=as<arma::mat>(Zlist(i));
    ZKList(i)=Ztemp*Ktemp;
  }
  
  arma::mat ZK(n,dimen, fill::zeros);
  unsigned int idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=as<arma::mat>(ZKList(i));
    idx = idx + dimvec(i) ;
  }
  idx=0;
  ///////////
  for(unsigned int i=0; i<k; i++) {
    ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 )=weights(i)*ZK.submat(0, idx, n-1,idx+ dimvec(i) - 1 );
    idx = idx + dimvec(i) ;
  }
  
  
  arma::vec uhat = (ZK).t()*Hinvhatehat;
  if (geterrors==true) {
    arma::vec Covbetahat=inv(X.t()*Omegainv*X);
    arma::vec PEV=diagvec(pow(sum(sigmasqhatvec.subvec(0,k-1)),2)*ZK.t()*Omegainv*(ZK-X*Covbetahat*X.t()*Omegainv*ZK));
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu") =sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
    }
    else {
      outaa=List::create(Rcpp::Named("Vu") =sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("Covbetahat") = Covbetahat, Rcpp::Named("PEV") = PEV);
      
    }
    double sigmasqu=outaa["Vu"];
    arma::vec weights=as<arma::vec>(outaa["weights"]);
    
    arma::vec PEVinsert;
    for (unsigned int i=0; i<k; i++) {
      PEVinsert=join_cols(PEVinsert,sigmasqu*weights(i)*diagvec(as<arma::mat>(Klist(i))));
      
    }
    
    outaa["PEV"] = PEVinsert-PEV;
    
    
    return outaa;
  } else {
    double loglik=dmvnorm_arma(y, X*betahat,Omegainv, true, true);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("Vu")=sum(sigmasqhatvec.subvec(0,k-1)) ,Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("Hinv") = Hinvhat, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("corfuncparamslist")=corfuncparamslist);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("Vu")=sum(sigmasqhatvec.subvec(0,k-1)) ,Rcpp::Named("Ve") =sigmasqhatvec(k),Rcpp::Named("betahat") = betahat,Rcpp::Named("uhat") = uhat, Rcpp::Named("loglik") = loglik, Rcpp::Named("weights") = sigmasqhatvec.subvec(0,k-1)/sum(sigmasqhatvec.subvec(0,k-1)),Rcpp::Named("corfuncparamslist")=corfuncparamslist);
      return outaa;
    }
    
  }////////////
  
}

/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//Omegainv,Ztemp,ZKZt,Rmat, as<arma::mat>(sigmahatlist(i))
double minfuncforcor(const arma::vec  params,const arma::mat data, std::string funname, arma::mat Omegainv, arma::mat Z,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambda, const double tolparinv) {
  
  arma::vec paramstemp=params;
  arma::mat K=callViaString(paramstemp,data, funname);
  
  ////Rcpp::Rcout  << det(Omegatemp) << std::endl;
  double funcmin;
  arma::mat ZKZt=Z*K*Z.t();
  funcmin=trace(Omegainv*kron(Lambda, ZKZt));
  funcmin=funcmin+as_scalar(vectorise(Rmat).t()*kron(Lambda, ZKtZt*solve(ZKZt+tolparinv*eye(Z.n_rows, Z.n_rows),ZKtZt))*vectorise(Rmat));
  ////Rcpp::Rcout  << funcmin << std::endl;
  return funcmin;
}







// [[Rcpp::depends("RcppArmadillo")]]


List simplex_minfuncforcormv(arma::vec start,const double tolparconv,const  int  maxiter, const arma::mat data, std::string funname, arma::mat Omegainv, arma::mat Z,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambda, const double tolparinv)
{
  const  double scale=1;
  int npar=start.n_elem;
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    
    // minfuncforcor(const arma::vec  params,const arma::mat data, std::string funname, arma::mat Omegainv, arma::mat Z,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambda)
    //minfuncforcor(as_scalar(vvec),,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda);
    f[j] =  minfuncforcor(vvec,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    
    fr = minfuncforcor(vr,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      
      
      fe = minfuncforcor(ve,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        
        
        fc= minfuncforcor(vc,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        
        fc= minfuncforcor(vc,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        
        f[vg] =minfuncforcor(vvec,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
        ksimplex++;
        
        
        vvec=(v.row(vh)).t();
        f[vh] =minfuncforcor(vvec,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  
  minval =minfuncforcor(vvec,data,funname,Omegainv,Z,ZKtZt,Rmat,Lambda, tolparinv);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}


/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//Omegainv,Ztemp,ZKZt,Rmat, as<arma::mat>(sigmahatlist(i))
double minfuncforsig(const arma::vec  params,const arma::mat data, std::string funname, arma::mat Omegainv,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambdat, const double tolparinv) {
  
  arma::vec paramstemp=params;
  arma::mat Lambda=callViaStringSigma(paramstemp,data, funname);
  
  ////Rcpp::Rcout  << det(Omegatemp) << std::endl;
  double funcmin;
  
  funcmin=trace(Omegainv*kron(Lambda, ZKtZt));
  funcmin=funcmin+as_scalar(vectorise(Rmat).t()*kron(Lambdat*solve(Lambda+tolparinv*eye(Lambda.n_rows,Lambda.n_rows),Lambdat), ZKtZt)*vectorise(Rmat));
  
  ////Rcpp::Rcout  << funcmin << std::endl;
  return funcmin;
}







// [[Rcpp::depends("RcppArmadillo")]]


List simplex_minfuncforsigmv(arma::vec start,const double   tolparconv,const  int  maxiter, const arma::mat data, std::string funname, arma::mat Omegainv,arma::mat ZKtZt,arma::mat Rmat, arma::mat Lambdat, const double tolparinv)
{
  const  double scale=1;
  int npar=start.n_elem;
  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  unsigned int i,j,m,row;
  int ksimplex;   	      /* track the number of function evaluations */
  unsigned int itr;	      /* track the number of iterations */
  
  //v;     /* holds vertices of simplex */
  double pn,qn;   /* values used to create initial simplex */
  //f;      /* value of function at each vertex */
  double fr;      /* value of function at reflection point */
  double fe;      /* value of function at expansion point */
  double fc;      /* value of function at contraction point */
  //vr;     /* reflection - coordinates */
  //ve;     /* expansion - coordinates */
  //vc;     /* contraction - coordinates */
  //vm;     /* centroid - coordinates */
  double minval;
  
  double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  arma::mat v(npar+1,npar);
  arma::vec f(npar+1);
  arma::vec vr(npar);
  arma::vec ve(npar);
  arma::vec vc(npar);
  arma::vec vm(npar);
  arma::vec vvec(npar);
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(static_cast<double>(npar+1))-1+npar)/(npar*sqrt(2.0));
  qn = scale*(sqrt(static_cast<double>(npar+1))-1)/(npar*sqrt(2.0));
  
  for (i=0; i<npar; i++) {
    v(0,i) = start[i];
  }
  
  for (i=1; i<=npar; i++) {
    for (j=0; j<npar; j++) {
      if (i-1 == j) {
        v(i,j) = pn + start[j];
      }
      else {
        v(i,j) = qn + start[j];
      }
    }
  }
  //Use the constraints, just scale so that 1>v[j]>0, sum(v[j])=1
  
  /* find the initial function values */
  for (j=0; j<=npar; j++) {
    vvec=(v.row(j)).t();
    
    
    // Omegainv,ZKtZt,Rmat,Lambdat)
    f[j] =  minfuncforsig(vvec,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
  }
  
  ksimplex = npar+1;
  
  /* begin the main loop of the minimization */
  for (itr=1; itr<=maxiter; itr++) {
    /* find the index of the largest value */
    vg=0;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vg]) {
        vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0; j<=npar; j++) {
      if (f[j] < f[vs]) {
        vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0; j<=npar; j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
        vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0; j<=npar-1; j++) {
      cent=0.0;
      for (m=0; m<=npar; m++) {
        if (m!=vg) {
          cent += v(m,j);
        }
      }
      vm[j] = cent/npar;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0; j<=npar-1; j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j]+ALPHA*(vm[j]-v(vg,j));
    }
    
    fr =minfuncforsig(vr,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
    
    ksimplex++;
    
    if (fr < f[vh] && fr >= f[vs]) {
      for (j=0; j<=npar-1; j++) {
        v(vg,j) = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0; j<=npar-1; j++) {
        /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
        ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
      }
      
      
      fe =minfuncforsig(ve,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
      ksimplex++;
      
      /* by making fe < fr as opposed to fe < f[vs],
                                              Rosenbrocks function takes 63 iterations as opposed
      to 64 when using double variables. */
      
      if (fe < fr) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j)= ve[j];
        }
        f[vg] = fe;
      }
      else {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vr[j];
        }
        f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
        }
        
        
        fc= minfuncforsig(vc,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
        ksimplex++;
      }
      else {
        /* perform inside contraction */
        for (j=0; j<=npar-1; j++) {
          /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
          vc[j] = vm[j]-BETA*(vm[j]-v(vg,j));
        }
        
        fc= minfuncforsig(vc,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
        ksimplex++;
      }
      
      
      if (fc < f[vg]) {
        for (j=0; j<=npar-1; j++) {
          v(vg,j) = vc[j];
        }
        f[vg] = fc;
      }
      /* at this point the contraction is not successful,
      we must halve the distance from vs to all the
      vertices of the simplex and then continue.
      10/31/97 - modified to account for ALL vertices.
      */
      else {
        for (row=0; row<=npar; row++) {
          if (row != vs) {
            for (j=0; j<=npar-1; j++) {
              v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0;
            }
          }
        }
        
        f[vg] =minfuncforsig(vvec,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
        ksimplex++;
        
        
        vvec=(v.row(vh)).t();
        f[vh] =minfuncforsig(vvec,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
        ksimplex++;
        /* print out the value at each iteration */
        
        
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0; j<=npar; j++) {
      fsum += f[j];
    }
    favg = fsum/(npar+1);
    s = 0.0;
    for (j=0; j<=npar; j++) {
      s += pow((f[j]-favg),2.0)/(npar);
    }
    s = sqrt(s);
    if (s < tolparconv) break;
    
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0; j<=npar; j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  
  for (j=0; j<npar; j++) {
    
    start[j] = v(vs,j);
  }
  vvec=(v.row(vs)).t();
  
  
  minval =minfuncforsig(vvec,data,funname,Omegainv,ZKtZt,Rmat,Lambdat,tolparinv);
  
  List outaa(2);
  outaa(0)=minval;
  outaa(1)=vvec;
  return outaa;
  
}








// [[Rcpp::depends("RcppArmadillo")]]


List mmmkmvcorfuncsigmafuncmvopt(const arma::mat & Y,const arma::mat & X, const arma::uvec & corfunc,const arma::uvec & corfuncfixed, const arma::uvec & sigfunc,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const Rcpp::List & R,const Rcpp::List & Siglist,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool Hinv=false) {
  /////////////
  //Rcpp::Rcout  << 1 << std::endl;
  
  List brentout;
  int k=Zlist.size();
  List corfuncparamslist(k+1);
  List sigfuncparamslist(k+1);
  arma::vec dimvec(k);
  int dimen=0;
  int q = X.n_cols;
  int d = Y.n_cols;
  arma::colvec betahat(q*d);
  arma::mat Betahat(q, d);
  int n = X.n_rows;
  arma::vec ehat(n*d);
  arma::mat Ehat(n,d);
  arma::mat Rmat(n,d);
  arma::mat bigX=kron(eye(d,d), X);
  
  arma::vec Kparamstemp;
  std::string Kfunctiontemp;
  
  arma::mat Kdatatemp;
  arma::mat  Ktemp;
  
  
  arma::vec Sparamstemp;
  std::string Sfunctiontemp;
  
  arma::mat Sdatatemp;
  arma::mat Stemp;
  
  
  arma::mat  Ztemp;
  arma::mat ZKZt(n,n);
  arma::mat Omega(n*d,n*d);
  arma::mat Omegainv(n*d,n*d);
  arma::mat LLt(d,d);
  arma::mat Lt(d,d);
  arma::mat Ltinv(d,d);
  List sigmahatlist(k+1);
  arma::mat sigma(d,d);
  arma::mat sigma0(d,d);
  if (d==1) {
    sigma=1;
  }
  else {
    sigma=cov(Y)/(k+1);
  }
  arma::mat Onesn=ones(n,1);
  arma::mat IdOnesn=kron(eye(d,d),Onesn);
  
  Rcpp::List ZKZtlist(k+1);
  
  for(unsigned int i=0; i<k; i++) {
    if (corfunc(i)) {
      //Rcpp::Rcout  << 1999999 << std::endl;
      Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
      Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(Klist(i))(1));
      corfuncparamslist(i)=Kparamstemp;
      Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
      //Rcpp::Rcout  << Ktemp.n_rows << std::endl << Ktemp.n_cols << std::endl;
    } else {
      Ktemp=as<arma::mat>(Klist(i));
    }
    dimvec(i) = Ktemp.n_rows ;
    dimen = dimen+dimvec(i);
    Ztemp=as<arma::mat>(Zlist(i));
    ZKZtlist(i)=Ztemp*Ktemp*Ztemp.t();
    //Rcpp::Rcout  << 1999999 << std::endl;
    
  }
  //Rcpp::Rcout  <<2<< std::endl;
  
  if (corfunc(k)) {
    //Rcpp::Rcout  << 1999999 << std::endl;
    Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(R(0))(0));
    Kparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(R(0))(1));
    corfuncparamslist(k)=Kparamstemp;
    Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(R(0))(2));
    Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
  } else {
    Ktemp=as<arma::mat>(R(0));
  }
  
  ZKZtlist(k)=W*Ktemp*W.t();
  
  
  
  //Rcpp::Rcout  << 3 << std::endl;
  
  
  for(unsigned int i=0; i<(k+1); i++) {
    if (sigfunc(i)) {
      //Rcpp::Rcout  << 3.1 << std::endl;
      // Rcpp::Rcout  << 1999999 << std::endl;
      Sfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Siglist(i))(0));
      // Rcpp::Rcout  << 1999999 << std::endl;
      Sparamstemp= Rcpp::as<arma::vec>(Rcpp::as<Rcpp::List>(Siglist(i))(1));
      sigfuncparamslist(i)=Sparamstemp;
      Sdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Siglist(i))(2));
      Stemp=callViaStringSigma(Sparamstemp,Sdatatemp,Sfunctiontemp);
      sigmahatlist(i)=Stemp;
      //Rcpp::Rcout  << Stemp.n_rows << std::endl << Stemp.n_cols << std::endl;
    } else {
      //Rcpp::Rcout  << 3.2 << std::endl;
      
      sigmahatlist(i)=sigma;
      //Rcpp::Rcout  << sigma << std::endl;
      
    }
    
  }
  //Rcpp::Rcout  << 4 << std::endl;
  
  unsigned int itercountouter=0;
  do {
    itercountouter=itercountouter+1;
    
    sigma0=as<arma::mat>(sigmahatlist(k));////this is for convergence check
    
    Omega=zeros(n*d,n*d);
    for(unsigned int i=0; i<(k+1); i++) {
      Omega=Omega+kron(as<arma::mat>(sigmahatlist(i)),as<arma::mat>(ZKZtlist(i)));
    }
    //Rcpp::Rcout  << itercountouter << std::endl;
    Omegainv=inv_sympdsamm(Omega,tolparinv);
    betahat=solve(bigX.t()*Omegainv*bigX, bigX.t()*Omegainv*reshape(Y,n*d,1));
    
    
    ehat=reshape(Y, n*d,1)-bigX*betahat;
    
    Ehat=reshape(ehat,n,d);
    
    Rmat=reshape(Omegainv*ehat,n,d);
    
    for(unsigned int i=0; i<(k+1); i++) {
      //Rcpp::Rcout  << itercountouter << std::endl;
      //Rcpp::Rcout  << i << std::endl;
      
      ZKZt=as<arma::mat>(ZKZtlist(i));
      LLt=IdOnesn.t()*(kron(ones(d,d),ZKZt)%Omegainv)*IdOnesn;
      Lt=diagmat(pow(diagvec(LLt)+tolparinv,.5));
      Lt=cholsammupper(Lt,tolparinv);
      Ltinv=inv(trimatu(Lt+tolparinv*eye(LLt.n_rows,LLt.n_rows)));
      if (!(sigfunc(i))) {
        sigmahatlist(i)=Ltinv*symMroot(Lt*(as<arma::mat>(sigmahatlist(i))*Rmat.t()*ZKZt*Rmat*as<arma::mat>(sigmahatlist(i)))*Lt.t())*Ltinv.t();
      } else {
        
        Sfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Siglist(i))(0));
        Sparamstemp= Rcpp::as<arma::vec>(sigfuncparamslist(i));
        Sdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Siglist(i))(2));
        
        
        if (i<k) {
          Ztemp=Rcpp::as<arma::mat>(Zlist(i));
          ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
        } else {
          Ztemp=W;
          ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
        }
        
        //Rcpp::Rcout  << 99999999999<< std::endl;
        
        for(unsigned int iii=1; iii<(Sparamstemp.n_elem+1); iii++) {
          //Rcpp::Rcout  << itercountouter << std::endl;
          
          brentout=simplex_minfuncforsigmv(Sparamstemp,tolparconv,100, Sdatatemp, Sfunctiontemp, Omegainv,ZKZt,Rmat, as<arma::mat>(sigmahatlist(i)), tolparinv);
          //Rcpp::Rcout  << brentout << std::endl;
          Sparamstemp=Rcpp::as<arma::vec>(brentout(1));
          sigfuncparamslist(i)=Sparamstemp;
          
          
          sigmahatlist(i)=callViaStringSigma(Sparamstemp,Sdatatemp,Sfunctiontemp);
          
        }
        
      }
      
      
      if (itercountouter % 10==0) {
        if (corfunc(i)) {
          if (!corfuncfixed(i)) {
            if (i<k) {
              Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
              Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
              Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
            } else {
              Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(R(0))(0));
              Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
              Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(R(0))(2));
            }
            
            
            if (i<k) {
              Ztemp=Rcpp::as<arma::mat>(Zlist(i));
              ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
            } else {
              Ztemp=W;
              ZKZt=Rcpp::as<arma::mat>(ZKZtlist(i));
            }
            
            
            brentout=simplex_minfuncforcormv(Kparamstemp,tolparconv,100, Kdatatemp, Kfunctiontemp, Omegainv,Ztemp,ZKZt,Rmat, as<arma::mat>(sigmahatlist(i)), tolparinv);
            Kparamstemp=Rcpp::as<arma::vec>(brentout(1));
            corfuncparamslist(i)=Kparamstemp;
            
            ZKZtlist(i)= Ztemp*callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp)*Ztemp.t();
            
          }
        }
      }
    }
    
  }      while (itercountouter<maxiter && fabs(sum(diagvec(as<arma::mat>(sigmahatlist(k))-sigma0)))/(fabs(1+sum(diagvec(sigma0))))>tolparconv);
  
  
  
  Betahat=reshape(betahat,q, d);
  arma::mat HobsInv;
  arma::vec eehat;
  arma::mat Omegainve;
  arma::mat varvecG;
  arma::mat gpred;
  List Gpredlist(k);
  List PEVlist(k);
  arma::mat Vet=as<arma::mat>(sigmahatlist(k));
  Omega=zeros(n*d,n*d);
  
  for(unsigned int i=0; i<(k+1); i++) {
    Omega=Omega+kron(as<arma::mat>(ZKZtlist(i)),as<arma::mat>(sigmahatlist(i)));
  }
  Omegainv=inv_sympdsamm(Omega,tolparinv);
  
  arma::mat Xforvec=kron(X,eye(d,d));
  arma::mat varBhat =inv(Xforvec.t()*Omegainv*Xforvec);
  List Klistnew(k);
  arma::mat Rnew;
  int nk;
  for(unsigned int i=0; i<(k); i++) {
    if (corfunc(i)) {
      Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(Klist(i))(0));
      Kparamstemp= Rcpp::as<arma::vec>(corfuncparamslist(i));
      Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(Klist(i))(2));
      Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
      Klistnew(i)=Ktemp;
    } else {
      Ktemp=as<arma::mat>(Klist(i));
      Klistnew(i)=Ktemp;
    }
    nk=  Ktemp.n_cols;
    arma::vec eehat =arma::vectorise((Y-X*Betahat).t());
    Omegainve = Omegainv * eehat;
    varvecG = kron(Ktemp,as<arma::mat>(sigmahatlist(i)));
    gpred =  varvecG*(kron(as<arma::mat>(Zlist(i)).t(),eye(d,d)))*Omegainve;
    Gpredlist(i)= (reshape(gpred, d, nk)).t();
    
  }
  
  if (corfunc(k)) {
    Kfunctiontemp= Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(R(0))(0));
    Kparamstemp=Rcpp::as<arma::vec>(corfuncparamslist(k));
    Kdatatemp= Rcpp::as<arma::mat>(Rcpp::as<Rcpp::List>(R(0))(2));
    Ktemp=callViaString(Kparamstemp,Kdatatemp,Kfunctiontemp);
    Rnew=Ktemp;
  } else {
    Ktemp=as<arma::mat>(R(0));
    Rnew=Ktemp;
  }
  
  
  
  if (geterrors==true) {
    for(unsigned int i=0; i<(k); i++) {
      arma::mat P= Omegainv - Omegainv *Xforvec* varBhat*Xforvec.t()*Omegainv;
      arma::mat Zforvec=kron(as<arma::mat>(Zlist(i)), eye(d,d));
      arma::mat ZKforvec=Zforvec * varvecG;
      arma::vec varGhat = diagvec(ZKforvec.t()* P* ZKforvec);
      PEVlist(i)= diagvec(varvecG) - varGhat;
    }
    
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klistnew,sigmahatlist,Betahat, W, Rnew);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("corfuncparamslist") =corfuncparamslist,Rcpp::Named("sigfuncparamslist") =sigfuncparamslist,Rcpp::Named("Hinv")= HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist,Rcpp::Named("PEVlist") =PEVlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("corfuncparamslist") =corfuncparamslist,Rcpp::Named("sigfuncparamslist") =sigfuncparamslist);
      return outaa;
    }
  } else {
    arma::vec loglik=loglikfuncmmmkmv(Y,X,Zlist, Klistnew,sigmahatlist,Betahat, W, Rnew);
    List outaa;
    if (Hinv) {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("corfuncparamslist") =corfuncparamslist,Rcpp::Named("sigfuncparamslist") =sigfuncparamslist,Rcpp::Named("Hinv")= HobsInv);
      return outaa;
    }
    else {
      outaa=List::create(Rcpp::Named("sigmahatlist") =sigmahatlist,Rcpp::Named("betahat") =betahat, Rcpp::Named("Gpredlist") =Gpredlist, Rcpp::Named("loglik") =loglik, Rcpp::Named("corfuncparamslist") =corfuncparamslist,Rcpp::Named("sigfuncparamslist") =sigfuncparamslist);
      return outaa;
    }
  }
  
}

////////////////////////////////////
using namespace arma;
using namespace Rcpp;



///////////////////////////////###############################EMM

// [[Rcpp::depends("RcppArmadillo")]]

bool isidentity(arma::mat a) {
  bool flag=true;
  for (unsigned int i = 0; i < a.n_rows; i++)
  {
    for (unsigned int j = 0; j <a.n_cols; j++)
    {
      if (a(i,i) != 1 && a(i,j) != 0)
      {
        flag = false;
        break;
      }
    }
  }
  return flag;
}

/////////////////////################################################










// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List  SAMM(const arma::mat & Y,const arma::mat & X, const Rcpp::List & Zlist, const Rcpp::List & Klist, double lambda, const arma::mat & W,const Rcpp::List & R,const Rcpp::List & Siglist, const arma::uvec & corfunc, const arma::uvec & corfuncfixed, const arma::uvec & sigfunc,const std::string mmalg, const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool Hinv=false) {
  List outSAM;
  ////univariate/multivariate
  bool univariate=false;
  int ncolY=Y.n_cols;
  if (!(ncolY>1)) {
    univariate=true;
  }
  //Rcpp::Rcout  << "nkernel:" << std::endl << Zlist.length() << std::endl;
  
  ///one kernel / multi kernel
  bool onekernel=true;
  int nKernel=Zlist.length();
  if (nKernel>1) {
    onekernel=false;
  }
  ///corfunc
  bool corfunctions=false;
  for (unsigned int i=0; i<(nKernel+1); i++) {
    if (corfunc(i)) {
      corfunctions=true;
    }
  }
  
  ///shrinkage
  bool shrink=false;
  if (lambda>0) {
    shrink=true;
  }
  
  ///R diag or not
  bool Rident=false;
  if (!corfunc(nKernel)) {
    Rident=isidentity(Rcpp::as<arma::mat>(R(0)));
    // Rcpp::Rcout  << "Rident:" << std::endl << Rident << std::endl;
  }
  
  bool Wident=false;
  if (!corfunc(nKernel)) {
    Wident=isidentity(W);
    // Rcpp::Rcout  << "Rident:" << std::endl << Rident << std::endl;
  }
  
  //Rcpp::Rcout  << univariate << std::endl;
  //Rcpp::Rcout  << onekernel << std::endl;
  //Rcpp::Rcout  << corfunctions << std::endl;
  
  ///####################1
  if (univariate) {
    if (onekernel) {
      if (!corfunctions) {
        if ((Rident && Wident)) {
          //Rcpp::Rcout  << 1 << std::endl;
          if (mmalg=="emm_reml") {
            //emm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const arma::mat & K, const bool & REML=true, const double & tolparconv=1e-7, const int maxiter=1000, bool geterrors=false)
            outSAM=emm(vectorise(Y), X,Zlist, Klist, true,tolparconv, maxiter,geterrors,Hinv);
          }
          else if (mmalg=="emm_ml") {
            //emm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const arma::mat & K, const bool & REML=true, const double & tolparconv=1e-7, const int maxiter=1000, bool geterrors=false)
            outSAM=emm(vectorise(Y), X,Zlist, Klist, false,tolparconv, maxiter,geterrors,Hinv);
          }
          
          else if (mmalg=="dmm_ml") {
            // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
            outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,false,Hinv);
          }
          else if (mmalg=="dmm_reml") {
            // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
            outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,true,Hinv);
          }
          else if (mmalg=="dermm_reml1") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
          }
          else if (mmalg=="dermm_reml2") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
          }
          
          else if (mmalg=="mm_ml") {
            List Zlist0(2);
            List Klist0(2);
            Zlist0(0)=Rcpp::as<arma::mat>(Zlist(0));
            Zlist0(1)=W;
            Klist0(0)=Rcpp::as<arma::mat>(Klist(0));
            Klist0(1)=Rcpp::as<arma::mat>(R(0));
            
            //
            outSAM=mm(vectorise(Y), X,Zlist0, Klist0,tolparconv, tolparinv,maxiter, geterrors,Hinv);
          }
          
          else if (mmalg=="emmmk_reml") {
            //List emmmk(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool REML=true, bool geterrors=false){
            outSAM=emmmk(vectorise(Y), X,Zlist, Klist, tolparconv,maxiter, true,geterrors,Hinv);
          }
          else if (mmalg=="emmmk_ml") {
            //List emmmk(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool REML=true, bool geterrors=false){
            outSAM=emmmk(vectorise(Y), X,Zlist, Klist, tolparconv,maxiter, false,geterrors,Hinv);
          }
          
          else if (mmalg=="mmmk_ml") {
            //List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, arma::mat & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmk(vectorise(Y), X,Zlist, Klist, W,Rmat, tolparconv,tolparinv,maxiter, geterrors, lambda,Hinv);
          }
          
          
          else if (mmalg=="emmmv_ml") {
            //List emmmv(arma::mat & Y,arma::mat & X,const arma::mat & Z, const arma::mat & K,const double & tolparconv=1e-8, double tolparinv=1e-8,const int & maxiter=500, bool geterrors=false){
            arma::mat Zmat=Rcpp::as<arma::mat>(Zlist(0));
            arma::mat Kmat=Rcpp::as<arma::mat>(Klist(0));
            outSAM=emmmv(Y,X,Zmat, Kmat, tolparconv,tolparinv,maxiter, geterrors,Hinv);
          }
          else if (mmalg=="mmmv_ml") {
            //List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200,const  bool & geterrors=false){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmv(Y,X,Zlist, Klist, W, Rmat,tolparconv,tolparinv,maxiter, geterrors,Hinv);
          }
          
          
        }
      }
    }
  }
  
  
  ///####################2
  if (univariate) {
    if (onekernel) {
      if (!corfunctions) {
        if (!(Rident && Wident)) {
          //Rcpp::Rcout  << 2<< std::endl;
          if (mmalg=="dmm_ml") {
            // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
            outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,false,Hinv);
          }
          else if (mmalg=="dmm_reml") {
            // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
            outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,true,Hinv);
          }
          else if (mmalg=="mm_ml") {
            //mm(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false)
            List Zlist0(2);
            List Klist0(2);
            Zlist0(0)=Rcpp::as<arma::mat>(Zlist(0));
            Zlist0(1)=W;
            Klist0(0)=Rcpp::as<arma::mat>(Klist(0));
            Klist0(1)=Rcpp::as<arma::mat>(R(0));
            
            //
            outSAM=mm(vectorise(Y), X,Zlist0, Klist0,tolparconv, tolparinv,maxiter, geterrors,Hinv);
          }
          else if (mmalg=="dermm_reml1") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
          }
          else if (mmalg=="dermm_reml2") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
          }
          else if (mmalg=="mmmk_ml") {
            //List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, arma::mat & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmk(vectorise(Y), X,Zlist, Klist, W,Rmat, tolparconv,tolparinv,maxiter, geterrors, lambda,Hinv);
          }
          else if (mmalg=="mmmv_ml") {
            //List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200,const  bool & geterrors=false){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmv(Y,X,Zlist, Klist, W, Rmat,tolparconv,tolparinv,maxiter, geterrors,Hinv);
          }
          
          //List mmmkmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200, bool geterrors=false){
          else if (mmalg=="mmmkmv_ml") {
            //List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200,const  bool & geterrors=false){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmkmv(Y,X,Zlist, Klist, W, Rmat,tolparconv,tolparinv,maxiter, geterrors,Hinv);
          }
        }
      }
    }
  }
  
  
  ///####################3
  if (univariate) {
    if (onekernel) {
      if ((corfunc(0) && !corfunc(1))) {
        //Rcpp::Rcout  << 3 << std::endl;
        if (mmalg=="dmm_ml") {
          // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
          outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,false,Hinv);
        }
        else if (mmalg=="dmm_reml") {
          // dmm(const arma::colvec & y,const arma::mat & X,const arma::mat & Z, const Rcpp::List & Kfunc, const arma::mat & W, const arma::mat & R,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, bool geterrors=false, bool REML=false){
          outSAM=dmm(vectorise(Y), X,Rcpp::as<arma::mat>(Zlist(0)), as<Rcpp::List>(Klist(0)), W,Rcpp::as<arma::mat>(R(0)),tolparconv, tolparinv,maxiter, geterrors,true,Hinv);
        }
        
        else if (mmalg=="dermm_reml1") {
          //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
          outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
        }
        else if (mmalg=="dermm_reml2") {
          //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
          outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
        }
      }
    }
  }
  
  
  
  ////################4
  
  if (univariate) {
    if (!onekernel) {
      if (!corfunctions) {
        if (!shrink) {
          if ((Rident && Wident)) {
            //Rcpp::Rcout  << 4 << std::endl;
            if (mmalg=="emmmk_reml") {
              //List emmmk(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool REML=true, bool geterrors=false){
              outSAM=emmmk(vectorise(Y), X,Zlist, Klist, tolparconv,maxiter, true,geterrors,Hinv);
            }
            else if (mmalg=="emmmk_ml") {
              //List emmmk(const arma::vec  y,const arma::mat  X,const Rcpp::List  Zlist, const Rcpp::List  Klist, const double  tolparconv=1e-9, const int  maxiter=1000, bool REML=true, bool geterrors=false){
              outSAM=emmmk(vectorise(Y), X,Zlist, Klist, tolparconv,maxiter, false,geterrors,Hinv);
            }
            else if (mmalg=="dermm_reml1") {
              //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
              outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
            }
            else if (mmalg=="dermm_reml2") {
              //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
              outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
            }
            else if (mmalg=="mmmk_ml") {
              //List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, arma::mat & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
              arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
              outSAM=mmmk(vectorise(Y), X,Zlist, Klist, W,Rmat, tolparconv,tolparinv,maxiter, geterrors, lambda,Hinv);
            }
            else if (mmalg=="mmmkmv_ml") {
              //List mmmkmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200, bool geterrors=false){
              arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
              outSAM=mmmkmv(Y, X,Zlist, Klist, W,Rmat, tolparconv,tolparinv, maxiter, geterrors,Hinv);
            }
            
          }
        }
      }
    }
  }
  ////################5
  
  if (univariate) {
    if (!onekernel) {
      if (!shrink) {
        if (!corfunctions) {
          if (!(Rident && Wident)) {
            //Rcpp::Rcout  << 5 << std::endl;
            if (mmalg=="mmmk_ml") {
              //List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, arma::mat & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
              arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
              outSAM=mmmk(vectorise(Y), X,Zlist, Klist, W,Rmat, tolparconv,tolparinv,maxiter, geterrors, lambda,Hinv);
            }
            else if (mmalg=="dermm_reml1") {
              //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
              outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
            }
            else if (mmalg=="dermm_reml2") {
              //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
              outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
            }
            else if (mmalg=="mmmkmv_ml") {
              //List mmmkmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200, bool geterrors=false){
              arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
              outSAM=mmmkmv(Y, X,Zlist, Klist, W,Rmat, tolparconv,tolparinv, maxiter, geterrors,Hinv);
            }
          }
        }
      }
    }
  }
  
  
  ////################4.2
  
  
  ////################5.2
  
  if (univariate) {
    if (onekernel) {
      if (!shrink) {
        if (corfunctions) {
          //Rcpp::Rcout  << 6 << std::endl;
          if (mmalg=="mmmk_ml") {
            //List mmmkcorfuncmvopt(const arma::colvec & y,const arma::mat & X, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, Rcpp::List & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            outSAM=mmmkcorfuncmvopt(vectorise(Y), X,corfunc, corfuncfixed, Zlist, Klist, W,R,tolparconv,tolparinv,maxiter, geterrors,0,Hinv);
          }
          else if (mmalg=="dermm_reml1") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
          }
          else if (mmalg=="dermm_reml2") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
          }
          
        }
      }
    }
  }
  
  
  
  ////################5.2
  
  if (univariate) {
    if (!onekernel) {
      if (!shrink) {
        if (corfunctions) {
          //Rcpp::Rcout  << 6 << std::endl;
          if (mmalg=="mmmk_ml") {
            //List mmmkcorfuncmvopt(const arma::colvec & y,const arma::mat & X, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, Rcpp::List & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            outSAM=mmmkcorfuncmvopt(vectorise(Y), X,corfunc, corfuncfixed, Zlist, Klist, W,R,tolparconv,tolparinv,maxiter, geterrors,0,Hinv);
          }
          else if (mmalg=="dermm_reml1") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 1,Hinv);
          }
          else if (mmalg=="dermm_reml2") {
            //dermm(const arma::vec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W, const Rcpp::List & R, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const double & tolparconv=1e-10, const double & tolparinv=1e-10,const int & maxiter=1000, const bool & geterrors=false,const int  & methodopt=0)
            outSAM=dermm(vectorise(Y), X,Zlist, Klist, W,R,corfunc, corfuncfixed, tolparconv, tolparinv,maxiter, geterrors, 2,Hinv);
          }
          
        }
      }
    }
  }
  
  
  /////////////////////////////
  ////################4
  
  if (univariate) {
    if (!onekernel) {
      if (!corfunctions) {
        if (shrink) {
          //Rcpp::Rcout  << 7 << std::endl;
          if (mmalg=="mmmk_ml") {
            //List mmmk(const arma::colvec & y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, arma::mat & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            //Rcpp::Rcout  <<lambda << std::endl;
            //  List mmmk(y,X,Zlist,Klist,W,R,1e-10 ,const double & tolparinv=1e-10,const int & maxiter=200, bool geterrors=false, const double & lambda=0){
            
            outSAM=mmmk(vectorise(Y), X,Zlist, Klist, W,Rmat, tolparconv,tolparinv,maxiter, geterrors, lambda,Hinv);
          }
          
        }
      }
    }
  }
  
  
  
  
  ////################5.2
  
  if (univariate) {
    if (!onekernel) {
      if (shrink) {
        if (corfunctions) {
          //Rcpp::Rcout  <<8 << std::endl;
          if (mmalg=="mmmk_ml") {
            //List mmmkcorfuncmvopt(const arma::colvec & y,const arma::mat & X, const arma::uvec & corfunc,const arma::uvec & corfuncfixed,const Rcpp::List & Zlist, const Rcpp::List & Klist, arma::mat & W, Rcpp::List & R,const double & tolparconv=1e-10, const int & maxiter=200, bool geterrors=false, double lambda=0){
            outSAM=mmmkcorfuncmvopt(vectorise(Y), X,corfunc, corfuncfixed, Zlist, Klist, W,R,tolparconv,tolparinv,maxiter, geterrors,lambda,Hinv);
          }
          
        }
      }
    }
  }
  
  
  
  
  
  
  
  
  
  
  ////################6
  
  if (!univariate) {
    if (onekernel) {
      if (!corfunctions) {
        if ((Rident && Wident)) {
          // Rcpp::Rcout  << 9 << std::endl;
          if (mmalg=="emmmv_ml") {
            arma::mat Zmat=Rcpp::as<arma::mat>(Zlist(0));
            arma::mat Kmat=Rcpp::as<arma::mat>(Klist(0));
            outSAM=emmmv(Y,X,Zmat, Kmat, tolparconv,tolparinv,maxiter, geterrors,Hinv);
          }
          
          else if (mmalg=="mmmv_ml") {
            //List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200,const  bool & geterrors=false){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmv(Y,X,Zlist, Klist, W, Rmat,tolparconv,tolparinv,maxiter, geterrors,Hinv);
            
          }
        }
      }
    }
  }
  //////#########
  
  ////################7
  
  if (!univariate) {
    if (onekernel) {
      if (!corfunctions) {
        if (!(Rident && Wident)) {
          //Rcpp::Rcout  << 10 << std::endl;
          if (mmalg=="mmmv_ml") {
            //List mmmv(const arma::mat & Y,const arma::mat & X,const Rcpp::List & Zlist, const Rcpp::List & Klist, const arma::mat & W,const arma::mat & R,const double & tolparconv=1e-9, const double & tolparinv=1e-9,const int & maxiter=200,const  bool & geterrors=false){
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmv(Y,X,Zlist, Klist, W, Rmat,tolparconv,tolparinv,maxiter, geterrors,Hinv);
            
          }
        }
      }
    }
  }
  //////#########
  
  
  ////################7
  
  if (!univariate) {
    if (!onekernel) {
      if (!corfunctions) {
        if (!(Rident && Wident)) {
          //Rcpp::Rcout  << 11 << std::endl;
          if (mmalg=="mmmkmv_ml") {
            arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
            outSAM=mmmkmv(Y, X,Zlist, Klist, W,Rmat, tolparconv,tolparinv, maxiter, geterrors,Hinv);
          }
        }
      }
    }
  }
  
  if (!univariate) {
    if (!corfunctions) {
      if (mmalg=="mmmkmv_ml") {
        arma::mat Rmat=Rcpp::as<arma::mat>(R(0));
        outSAM=mmmkmv(Y, X,Zlist, Klist, W,Rmat, tolparconv,tolparinv, maxiter, geterrors,Hinv);
      }
    }
  }
  //////#########
  //////#########
  
  
  
  
  return outSAM;
}


///////////////////////////////////////////////////////////////////


//////////////////////////////////////////////



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat ar1cov_cppforR(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho = (2/M_PI)*atan(params(0));
  arma::mat H = disteucarma(times, times);
  arma::mat V(data(0,0),data(0,0));
  for (unsigned int i = 0; i < data(0,0); i++) {
    for (unsigned int j = 0; j < data(0,0); j++) {
      V(i,j)=powf(rho,H(i,j));
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat ar1hetcov_cppforR(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho = (2/M_PI)*atan(params(0));
  arma::mat H = disteucarma(times, times);
  arma::mat V(data(0,0),data(0,0));
  for (unsigned int i = 0; i < data(0,0); i++) {
    for (unsigned int j = 0; j < data(0,0); j++) {
      V(i,j)=powf(rho,H(i,j));
    }
  }
  arma::vec vvec(data(0,0));
  vvec(0)=1;
  for (unsigned int i = 1; i < data(0,0); i++) {
    vvec(i)=exp(params(i));
  }
  return diagmat(vvec)*V*diagmat(vvec);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat arma11cov_cppforR(const arma::vec & params, const arma::mat  & data) {
  arma::vec times =linspace(1,data(0,0),data(0,0));
  double rho =(2/M_PI)*atan(params(0));
  double lambda =(2/M_PI)*atan(params(1));
  
  arma::mat H = disteucarma(times, times);
  
  arma::mat V(data(0,0),data(0,0));
  for (int i = 0; i < data(0,0); i++) {
    for (int j = 0; j < data(0,0); j++) {
      if (abs(i-j)>1) {
        V(i,j)=lambda*powf(rho,H(i,j));
      } else {
        if (i==j) {
          V(i,j)=1;
        } else {
          V(i,j)=lambda;
        }
      }
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmcov_cppforR(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =data(0,0);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < data(0,0); i++) {
    V(i,i)=1;
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmhetcov_cppforR(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =data(0,0);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < data(0,0); i++) {
    V(i,i)=1;
  }
  arma::vec vvec(data(0,0));
  vvec(0)=1;
  for (unsigned int i = 1; i < data(0,0); i++) {
    vvec(i)=exp(params(i));
  }
  return diagmat(vvec)*V*diagmat(vvec);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat lincombcov_cppforR(const arma::vec & params, const arma::mat  & data) {
  int k=params.n_elem;
  arma::vec weights=params-min(params);
  weights=weights/sum(weights);
  arma::mat V(data.n_rows,data.n_rows);
  for (unsigned int i = 0; i < k; i++) {
    V=V+weights(i)*data.submat(0,i*data.n_rows,data.n_rows-1,(i+1)*data.n_rows-1);
  }
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat unstrcov_cppforR(const arma::vec & params, const arma::mat  & data) {
  int d1=data(0,0);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat V=diagmat(d2vec)*D1*diagmat(d2vec);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat diagcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int k=params.n_elem;
  arma::vec v=ones(k+1);
  for (unsigned int i = 1; i < (k+1); ++i) {
    v(i)=exp(params(i-1));
  }
  arma::mat V=diagmat(v);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat unstrKronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(diagmat(d2vec)*D1*diagmat(d2vec),D2);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat diagKronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1=zeros(d1,d1);
  D1.diag()=join_cols(ones(1,1),exp(params));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat ar1KronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat ar1hetKronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1hetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat arma11KronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=arma11cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmKronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmhetKronKcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmhetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D1,D2);
  return V;
}






// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat UnstrKronUnstrcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d1vec(d1);
  d1vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d1vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2(d2,d2);
  
  for (unsigned int i = 0; i < d2; i++) {
    for (unsigned int j = i; j < d2; j++) {
      if (i==j) {
        D2(i,j)=1;
      }
      else {
        D2(i,j)=(2/M_PI)*atan(params(ii));
        D2(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d2);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d2; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  
  arma::mat V=kron(diagmat(d1vec)*D1*diagmat(d1vec),diagmat(d2vec)*D2*diagmat(d2vec));
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat rbfcov_cppforR(const arma::vec & params, const arma::mat & data) {
  arma::mat C=disteucarma(data,data);
  return exp(-exp(params(0))*pow(C, 2));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat expcov_cppforR(const arma::vec & params, const arma::mat & data) {
  arma::mat C=disteucarma(data,data);
  return exp(-(exp(params(0)))*C);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat rbfdistcov_cppforR(const arma::vec & params, const arma::mat & data) {
  return exp(-exp(params(0))*pow(data, 2));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat expdistcov_cppforR(const arma::vec & params, const arma::mat & data) {
  return exp(-(exp(params(0)))*data);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat relmatcov_cppforR(const arma::vec & params, const arma::mat & data) {
  arma::vec pks(data.n_cols);
  double c=0;
  for (unsigned int iter = 0; iter < data.n_cols; ++iter) {
    pks(iter)=sum(data.col(iter)+ones(data.n_rows))/(2*data.n_rows);
    c=c+2*pks(iter)*(1-pks(iter));
  }
  arma::mat W=data+1-2*ones(data.n_rows,1)*pks.t();
  arma::mat Amat=(1/c)*W*W.t();
  return Amat+params(0)*eye(Amat.n_cols,Amat.n_cols);
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat ConstMatcov_cppforR(const arma::vec & params, const arma::mat & data) {
  return data;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKronunstrcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d1);
  d2vec(0)=1;
  for (unsigned int i = 1; i < d1; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,diagmat(d2vec)*D1*diagmat(d2vec));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKrondiagcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1=zeros(d1,d1);
  D1.diag()=join_cols(ones(1,1),exp(params));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKronar1cov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKronar1hetcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=ar1hetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKronarma11cov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=arma11cov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKroncompsymmcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat KKroncompsymmhetcov_cppforR(const arma::vec & params, const arma::mat & data) {
  int d2=data(0,1);
  arma::mat D1=compsymmhetcov_cpp(params,data.submat(0,0,0,0));
  arma::mat D2=data.submat(1,2,d2,d2+1);
  arma::mat V=kron(D2,D1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat sppowcov_cppforR(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      V(i,j)=pow(rho, data(i,j));
    }
  }
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat splincov_cppforR(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      if (rho*data(i,j)<=1) {
        V(i,j)=(1-rho*data(i,j));
      } else {
        V(i,j)=0;
      }
    }
  }
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat splinlogcov_cppforR(const arma::vec & params, const arma::mat & data) {
  double rho=(2/M_PI)*atan(params(0));
  arma::mat V(data.n_cols,data.n_cols);
  for (unsigned int i = 0; i < data.n_cols; i++) {
    for (unsigned int j = 0; j < data.n_cols; j++) {
      if (std::log(rho*data(i,j))<=1) {
        V(i,j)=(1-rho*std::log(data(i,j)));
      } else {
        V(i,j)=0;
      }
    }
  }
  return V;
}




//////////////////////////////SIGMA Functions

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat diagSig_cppforR(const arma::vec & params, const arma::mat & data) {
  arma::mat V=diagmat(exp(params));
  return V;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat IdentSig_cppforR(const arma::vec & params, const arma::mat & data) {
  arma::mat V=exp(params(0))*eye(as_scalar(data),as_scalar(data));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat UnstrKronIdentSig_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D1(d1,d1);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d1; i++) {
    for (unsigned int j = i; j < d1; j++) {
      if (i==j) {
        D1(i,j)=1;
      }
      else {
        D1(i,j)=(2/M_PI)*atan(params(ii));
        D1(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d1vec(d1);
  for (unsigned int i = 0; i < d1; i++) {
    d1vec(i)=exp(params(ii));
    ii=ii+1;
  }
  arma::mat V=kron(diagmat(d1vec)*D1*diagmat(d1vec),eye(d2,d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat IdentKronUnstrSig_cppforR(const arma::vec & params, const arma::mat & data) {
  int d1=data(0,0);
  int d2=data(0,1);
  arma::mat D2(d2,d2);
  
  unsigned int ii=0;
  for (unsigned int i = 0; i < d2; i++) {
    for (unsigned int j = i; j < d2; j++) {
      if (i==j) {
        D2(i,j)=1;
      }
      else {
        D2(i,j)=(2/M_PI)*atan(params(ii));
        D2(j,i)=(2/M_PI)*atan(params(ii));
        ii=ii+1;
      }
    }
  }
  
  arma::vec d2vec(d2);
  for (unsigned int i = 0; i < d2; i++) {
    d2vec(i)=exp(params(ii));
    ii=ii+1;
    
  }
  arma::mat V=kron(eye(d1,d1),diagmat(d2vec)*D2*diagmat(d2vec));
  return V;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat FA1hetSig_cppforR(const arma::vec & params, const arma::mat & data) {
  int dim1=as_scalar(data);
  arma::vec d1=params.subvec(0,(dim1-1));
  arma::vec d2=params.subvec((dim1),(2*dim1-1));
  arma::mat V=d1*d1.t()+diagmat(exp(d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat FA1homSig_cppforR(const arma::vec & params, const arma::mat & data) {
  int dim1=as_scalar(data);
  arma::vec d1=params.subvec(0,(dim1-1));
  double d2=params(dim1);
  arma::mat V=d1*d1.t()+exp(d2)*eye(dim1,dim1);
  return V;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmhomSig_cppforR(const arma::vec & params, const arma::mat  & data) {
  double sigma = exp(params(0));
  double rho = (2/M_PI)*atan(params(1));
  int dim =as_scalar(data);
  arma::mat V=rho*ones(dim,dim);
  for (unsigned int i = 0; i < as_scalar(data); i++) {
    V(i,i)=1;
  }
  return sigma*V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat compsymmhetSig_cppforR(const arma::vec & params, const arma::mat  & data) {
  double rho = (2/M_PI)*atan(params(0));
  int dim =as_scalar(data);
  arma::mat V=rho*ones(dim,dim);
  arma::vec d1=params.subvec(1,(dim));
  return diagmat(d1)*V*diagmat(d1);
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat FAhetSig_cppforR(const arma::vec & params, const arma::mat  & data) {
  arma::mat Lambda=zeros(data(0,0),data(0,1));
  unsigned int ii=0;
  for (unsigned int i = 0; i < data(0,1); i++) {
    for (unsigned int j = i; j < data(0,0); j++) {
      Lambda(j,i)=params(ii);
      ii=ii+1;
    }
  }
  arma::vec d2(data(0,0));
  for (unsigned int j = 0; j < data(0,0); j++) {
    d2(j)=params(ii);
    ii=ii+1;
  }
  arma::mat V=Lambda*Lambda.t()+diagmat(exp(d2));
  return V;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat FAhomSig_cppforR(const arma::vec & params, const arma::mat  & data) {
  arma::mat Lambda=zeros(data(0,0),data(0,1));
  unsigned int ii=0;
  for (unsigned int i = 0; i < data(0,1); i++) {
    for (unsigned int j = i; j < data(0,0); j++) {
      Lambda(j,i)=params(ii);
      ii=ii+1;
    }
  }
  double d2=params(ii);
  arma::mat V=Lambda*Lambda.t()+exp(d2)*eye(data(0,0),data(0,0));
  return V;
}


