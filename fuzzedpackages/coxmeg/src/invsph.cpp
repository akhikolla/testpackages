// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SimplicialLDLT<SpMat> SpChol;

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List invsph(Eigen::SparseMatrix<double> & A, const Eigen::Map<Eigen::VectorXd> der, const Eigen::Map<Eigen::VectorXd> dv,
                  const Eigen::Map<Eigen::VectorXd> v1, const Eigen::MatrixXd & mx, const Eigen::Map<Eigen::VectorXd> v2, 
                  const Eigen::Map<Eigen::VectorXd> v3,const Eigen::MatrixXd & v4,
                  const Eigen::Map<Eigen::VectorXd> av, const Eigen::Map<Eigen::VectorXd> bw, 
                  const Eigen::VectorXd & f,const Eigen::VectorXd & inv,const Eigen::VectorXd & tau,
                  const int sol) {
  
  int nr = v1.size();
  int nc = mx.cols();
  double tol = 1.0e-8;

  Eigen::ArrayXd av2 = av.array()*av.array();
  Eigen::MatrixXd w_v_x;
  Eigen::MatrixXd hx;
  
  if(mx.rows()==0)
  {
    nc = 0;
  }
  
  int nc_1 = nc+1;
  Eigen::MatrixXd B(nr, nc_1);
  B.col(0) = der;
  
  if(nc>0)
  {
    w_v_x = mx.array().colwise() * v1.array();
    Eigen::MatrixXd temp(nr,nc);
    for(int j=0; j<nr; j++)
      temp.row(j) = w_v_x.row(v4(j,0));
    temp = temp.colwise().reverse().eval();
    
    for(int j=1; j<nr; ++j)
      temp.row(j) += temp.row(j-1);
    for(int j=0; j<nr; j++)
      w_v_x.row(j) = temp.row(v2(j));
    w_v_x = w_v_x.array().colwise() * av2;
    for(int j=1; j<nr; ++j)
      w_v_x.row(j) += w_v_x.row(j-1);
    for(int j=0; j<nr; j++)
      temp.row(j) = w_v_x.row(v3(j));	
    for(int j=0; j<nr; j++)
      w_v_x.row(j) = temp.row(v4(j,1));
    w_v_x = w_v_x.array().colwise() * v1.array();
    
    hx = mx.array().colwise() * bw.array();
    hx = hx.array() - w_v_x.array();
    for(int i = 0; i < nc; i++)
    {
      B.col(i+1) = hx.col(i);
    }
  }
  
  if(inv(0)==0)
  {
    B = A*B;
  }
  
  double tau_i = tau(0);
  Eigen::ArrayXd sd_v = bw.array();
  if(inv(0)==0)
  {
    sd_v = 1/sd_v;
    tau_i = 1/tau_i;
  }
  
  SpMat A_bk;
  A_bk = A;
    
  for (int k=0; k<A_bk.outerSize(); ++k)
    for (SpMat::InnerIterator it(A_bk,k); it; ++it)
    {
      if(it.row()==it.col())
      {
        it.valueRef() = dv(it.row()) + sd_v(it.row())*tau_i;
        break;
      }
    }
  
  SpChol solver_ch;
  Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper > solver_cg;
  Eigen::MatrixXd wb_sig_i_hx_der;
  if(sol<2)
  {
    solver_ch.compute(A_bk);
    wb_sig_i_hx_der = solver_ch.solve(B);
  }else{
    solver_cg.setTolerance(tol);
    solver_cg.compute(A_bk);
    wb_sig_i_hx_der = solver_cg.solve(B);
  }
  
  if(inv(0)==0)
  {
    for(int k=0; k<nc_1; k++)
      wb_sig_i_hx_der.col(k) = sd_v*wb_sig_i_hx_der.col(k).array();
  }else{
    wb_sig_i_hx_der = tau_i*wb_sig_i_hx_der;
  }
  
  if(f(0)==0)
  {
    
    return Rcpp::List::create(Rcpp::Named("xh") = hx.transpose(),
                              Rcpp::Named("wb_sig_i_hx_der") = wb_sig_i_hx_der);
  }else{
    Eigen::MatrixXd wb_sig_i_hx_der2 = wb_sig_i_hx_der;
    int order_i = f(0);
    
    for(int i = 0; i<order_i; i++)
    {
      w_v_x = wb_sig_i_hx_der2.array().colwise() * v1.array();
      for(int j=0; j<nr; j++)
        B.row(j) = w_v_x.row(v4(j,0));
      B = B.colwise().reverse().eval();
      
      for(int j=1; j<nr; ++j)
        B.row(j) += B.row(j-1);
      for(int j=0; j<nr; j++)
        w_v_x.row(j) = B.row(v2(j));
      w_v_x = w_v_x.array().colwise() * av2;
      for(int j=1; j<nr; ++j)
        w_v_x.row(j) += w_v_x.row(j-1);
      for(int j=0; j<nr; j++)
        B.row(j) = w_v_x.row(v3(j));	
      for(int j=0; j<nr; j++)
        w_v_x.row(j) = B.row(v4(j,1));
      w_v_x = w_v_x.array().colwise() * v1.array();
      
      if(inv(0)==0)
      {
        w_v_x = A*w_v_x;
      }
      if(sol<2)
      {
        wb_sig_i_hx_der2 = solver_ch.solve(w_v_x);
      }else{
        wb_sig_i_hx_der2 = solver_cg.solve(w_v_x);
      }
      
      if(inv(0)==0)
      {
        for(int k=0; k<nc_1; k++)
          wb_sig_i_hx_der2.col(k) = sd_v*wb_sig_i_hx_der2.col(k).array();
      }else{
        wb_sig_i_hx_der2 = tau_i*wb_sig_i_hx_der2;
      }
      wb_sig_i_hx_der = wb_sig_i_hx_der+wb_sig_i_hx_der2;
    }
    
    return Rcpp::List::create(Rcpp::Named("xh") = hx.transpose(),
                              Rcpp::Named("wb_sig_i_hx_der") = wb_sig_i_hx_der);
    
  }
  
}

