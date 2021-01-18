// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R

// [[Rcpp::export]]
Eigen::MatrixXd csqei(const Eigen::Map<Eigen::VectorXd> w_v, const Eigen::MatrixXd & mx, const Eigen::Map<Eigen::VectorXd> rs_rs, 
                      const Eigen::Map<Eigen::VectorXd> rs_cs,const Eigen::MatrixXi & ind, const Eigen::Map<Eigen::VectorXd> av) {
  
  int nr = mx.rows();
  int nc = mx.cols();
  
  Eigen::MatrixXd w_v_x = mx.array().colwise() * w_v.array();
  Eigen::MatrixXd temp(nr,nc);
  for(int j=0; j<nr; j++)
    temp.row(j) = w_v_x.row(ind(j,0));
  temp = temp.colwise().reverse().eval();
  
  for(int j=1; j<nr; ++j)
    temp.row(j) += temp.row(j-1);
  for(int j=0; j<nr; j++)
    w_v_x.row(j) = temp.row(rs_rs(j));
  w_v_x = w_v_x.array().colwise() * av.array();
  for(int j=1; j<nr; ++j)
    w_v_x.row(j) += w_v_x.row(j-1);
  for(int j=0; j<nr; j++)
    temp.row(j) = w_v_x.row(rs_cs(j));	
  for(int j=0; j<nr; j++)
    w_v_x.row(j) = temp.row(ind(j,1));
  w_v_x = w_v_x.array().colwise() * w_v.array();
  return w_v_x;
  
}



// [[Rcpp::export]]
Eigen::MatrixXd wma_cp(const Eigen::Map<Eigen::VectorXd> w, const Eigen::Map<Eigen::VectorXd> cs_p, const Eigen::MatrixXi & ind,
                       const Eigen::Map<Eigen::VectorXd> a) {
  
  int n = w.size();
  
  Eigen::ArrayXd a2 = a.array()*a.array();
  int as = a2.size();
  Eigen::ArrayXd a2csum(as);
  double cst = 0;
  for(int j = 0; j<as;j++)
  {
    cst += a2(j);
    a2csum(j) = cst;
  }
  
  Eigen::MatrixXd q(n,n);
  for(int j=0;j<n;j++)
  {
    for(int k=j;k<n;k++)
    {
      int min = cs_p(ind(j,1));
      int ind_2 = ind(k,1);
      if(min>cs_p(ind_2))
        min = cs_p(ind_2);
      
      double q_t = w(j)*w(k)*a2csum(min);
      q(j,k) = q_t;
      q(k,j) = q_t;
    }
  }
  
  return q;
}


//
// [[Rcpp::export]]
Eigen::VectorXd score_test(const Eigen::Map<Eigen::VectorXd> deriv, const Eigen::Map<Eigen::VectorXd> bw_v,
               const Eigen::Map<Eigen::VectorXd> w, const Eigen::Map<Eigen::VectorXd> rs_rs,const Eigen::Map<Eigen::VectorXd> rs_cs,
               const Eigen::Map<Eigen::VectorXd> cs_p, const Eigen::MatrixXi & ind,
               const Eigen::Map<Eigen::VectorXd> a, const Eigen::Map<Eigen::VectorXd> a2, const Eigen::VectorXd & tau, 
               const Eigen::Map<Eigen::MatrixXd> v,const Eigen::MatrixXd & cov, const Eigen::MatrixXd & x) {

    int n = bw_v.size();
    int n_c = cov.cols();
    if(cov.rows()==0)
    {
      n_c = 0;
    }
    
    int ns = x.cols();  
    Eigen::VectorXd tst(ns);
    for(int i=0; i<ns; i++)
    {
      Eigen::RowVectorXd xh = bw_v.array()*x.col(i).array() - csqei(w,x.col(i),rs_rs,rs_cs,ind,a2).col(0).array(); 
      Eigen::RowVectorXd v21(n+n_c);
      if(n_c>0)
      {
        v21.head(n_c) = xh*cov;
      }
      v21.tail(n) = xh;
      double sc_v = (xh*x.col(i) - v21*v*v21.transpose()).value();
      double deriv_x = x.col(i).transpose()*deriv;
      tst(i) = deriv_x*deriv_x/sc_v;
    }
    
    return tst;
}

// [[Rcpp::export]]
double logdet_ch(const Eigen::MatrixXd & X_m, const Eigen::MatrixXd & rad_m, const Eigen::VectorXd & bma_d,
                       const Eigen::VectorXd & bpa_d, const Eigen::VectorXd & cj_v) {
  
  int t = rad_m.cols();
  int q = cj_v.size();
  
  double app = 0;
  double ratio = bpa_d[0]/bma_d[0];
  
  Eigen::MatrixXd u = cj_v[0]*rad_m;
  Eigen::MatrixXd w0 = rad_m;
  Eigen::MatrixXd w1 = 2/bma_d[0]*(X_m*rad_m) - ratio*rad_m;
  Eigen::MatrixXd w2 = w1;
  u = u+cj_v[1]*w1;
  for(int j=2; j<q; j++)
  {
    w2 = 4/bma_d[0]*(X_m*w1) - w0 - 2*ratio*w1;
    u = u + cj_v[j]*w2;
    w0 = w1;
    w1 = w2;
  }
  u = rad_m.array()*u.array();
  app = u.sum()/t;
  
  return app;
}

// [[Rcpp::export]]
double logdet_lanczos(const Eigen::Map<Eigen::MatrixXd> X_m, const Eigen::Map<Eigen::MatrixXd> rad_m, const Eigen::VectorXi & m_d) {
  
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::ArrayXXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m-1,t);
  Eigen::ArrayXXd temp(n,t);
  Eigen::MatrixXd T;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  temp = w*v1.array();
  alpha.row(0) = temp.colwise().sum();
  w = w - v1.array().rowwise()*alpha.row(0);
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    temp = w*w;
    beta.row(j-1) = sqrt(temp.colwise().sum());
    
    for(int k=0;k<t;k++)
    {
      if(beta(j-1,k)<1e-10)
      {
        if(iter(k)==m)
          iter(k) = j;
      }
    }
    
    v2 = w.rowwise()*beta.row(j-1).inverse();
    
    w = X_m*v2;
    temp = w*v2.array();
    alpha.row(j) = temp.colwise().sum();
    w = w - v2.array().rowwise()*alpha.row(j) - v1.array().rowwise()*beta.row(j-1);
    v1 = v2;
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    T = Eigen::MatrixXd::Zero(it,it);
    T.diagonal() = alpha.col(i).head(it);
    if(it>1)
    {
      T.diagonal(1) = beta.col(i).head(it-1);
      T.diagonal(-1) = beta.col(i).head(it-1);
    }
    es.compute(T);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    eivl = eivt*eivl;
    logdet = logdet + eivl.sum();
  }
  
  return logdet*n/t;
}

// [[Rcpp::export]]
double logdet_lanczos_sp(const Eigen::MappedSparseMatrix<double> X_m, const Eigen::Map<Eigen::MatrixXd> rad_m, const Eigen::VectorXi & m_d) {
  
  int t = rad_m.cols();
  int n = rad_m.rows();
  int m = m_d[0];
  
  double logdet = 0;
  
  Eigen::MatrixXd v1;
  Eigen::MatrixXd v2;
  Eigen::ArrayXXd w(n,t);
  Eigen::ArrayXXd alpha(m,t);
  Eigen::ArrayXXd beta(m-1,t);
  Eigen::ArrayXXd temp(n,t);
  Eigen::MatrixXd T;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  v1 = rad_m;
  w = X_m*v1;
  temp = w*v1.array();
  alpha.row(0) = temp.colwise().sum();
  w = w - v1.array().rowwise()*alpha.row(0);
  
  Eigen::VectorXd iter(t);
  iter.fill(m);
  for(int j=1;j<m;j++)
  {
    temp = w*w;
    beta.row(j-1) = sqrt(temp.colwise().sum());
    
    for(int k=0;k<t;k++)
    {
      if(beta(j-1,k)<1e-10)
      {
        if(iter(k)==m)
          iter(k) = j;
      }
    }
    
    v2 = w.rowwise()*beta.row(j-1).inverse();
    
    w = X_m*v2;
    temp = w*v2.array();
    alpha.row(j) = temp.colwise().sum();
    w = w - v2.array().rowwise()*alpha.row(j) - v1.array().rowwise()*beta.row(j-1);
    v1 = v2;
  }
  
  for(int i = 0; i<t; i++)
  { 
    int it = iter(i);
    T = Eigen::MatrixXd::Zero(it,it);
    T.diagonal() = alpha.col(i).head(it);
    if(it>1)
    {
      T.diagonal(1) = beta.col(i).head(it-1);
      T.diagonal(-1) = beta.col(i).head(it-1);
    }
    es.compute(T);
    Eigen::ArrayXd eivt = Eigen::square(es.eigenvectors().row(0).array());
    Eigen::ArrayXd eivl = Eigen::log(es.eigenvalues().array());
    eivl = eivt*eivl;
    logdet = logdet + eivl.sum();
  }
  
  return logdet*n/t;
}
