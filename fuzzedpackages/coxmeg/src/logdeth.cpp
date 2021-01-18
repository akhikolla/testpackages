// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SimplicialLDLT<SpMat> SpChol;
typedef Eigen::LDLT<Eigen::MatrixXd> Chol;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
double logdeth(Eigen::SparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> dv, const Eigen::Map<Eigen::VectorXd> bw_v,
                  const Eigen::Map<Eigen::VectorXd> w, const Eigen::Map<Eigen::VectorXd> cs_p, const Eigen::MatrixXi & v4,
                  const Eigen::Map<Eigen::VectorXd> a, const Eigen::VectorXd & tau, const Eigen::VectorXi & inv, const Eigen::VectorXi & detap) {
  
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
  
  if(detap(0)==1)
  {
    Eigen::ArrayXd qd(n);
    
    for(int j=0;j<n;j++)
    {
      qd(j) = w(j)*w(j)*a2csum(cs_p(v4(j,1)));
    }
   
    if(inv(0)==1)
    {
      for (int k=0; k<A.outerSize(); ++k)
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
          if(it.row()==it.col())
          {
            it.valueRef() = dv(it.row()) + (bw_v(it.row()) - qd(it.row()))*tau(0);
            break;
          }
        }
        
      SpChol solver(A);
      Eigen::ArrayXd logdet = solver.vectorD();
      logdet = logdet.log();
      return logdet.sum();
    }else{
      for (int k=0; k<A.outerSize(); ++k)
        for (SpMat::InnerIterator it(A,k); it; ++it)
        {
          if(it.row()==it.col())
          {
            it.valueRef() = dv(it.row()) + 1/(bw_v(it.row()) - qd(it.row()))/tau(0);
            break;
          }
        }
        
        SpChol solver(A);
      Eigen::ArrayXd logdet = solver.vectorD();
      logdet = logdet.log();
      Eigen::ArrayXd dd = bw_v.array() - qd;
      dd = dd.log();
      return logdet.sum() + dd.sum();
    }
    
  }else{
    Eigen::MatrixXd q(n,n);
    
    for(int j=0;j<n;j++)
    {
      for(int k=j;k<n;k++)
      {
        int min = cs_p(v4(j,1));
        int ind_2 = v4(k,1);
        if(min>cs_p(ind_2))
          min = cs_p(ind_2);
        double q_t = w(j)*w(k)*a2csum(min);
        q(j,k) = q_t;
        q(k,j) = q_t;
      }
    }
    
    q = Eigen::MatrixXd(A)/tau(0) - q;
    q.diagonal() += bw_v; 
    Chol solver(q);
    Eigen::ArrayXd logdet = solver.vectorD();
    logdet = logdet.log();
    return logdet.sum();
  }
  
}

