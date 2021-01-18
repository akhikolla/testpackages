#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat Ac(arma::mat X, int n);
arma::mat Acs(arma::vec z, int n);
arma::vec Pinv(arma::vec z, int n);
arma::mat projAXB(arma::mat X0, double alpha, int n);
arma::mat projToSDC(arma::mat M);


List getEigenValues(arma::mat M) {
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, M);
 
  
  return List::create(
    _["val"]=eigval, 
    _["vec"]=eigvec
  );
}

// [[Rcpp::export]]
List sdp1_admm(arma::mat As, int K, List opts) {
  
  double rho = (opts.containsElementNamed("rho") ?  opts["rho"] : .1);
  int    T   = (opts.containsElementNamed("T") ?  opts["T"] : 10000);
  double tol = (opts.containsElementNamed("tol") ?  opts["tol"] : 1e-5);
  int report_interval = (opts.containsElementNamed("report_interval") ?  opts["report_interval"] : 100);
  
  int    n = As.n_rows;
  arma::vec delta = arma::zeros(T);
  
  arma::mat As_rescaled = (1./rho)*As, 
            U = arma::zeros(n,n),
            V = arma::zeros(n,n),
            X = arma::zeros(n,n),
            Xold = arma::zeros(n,n),
            Y = arma::zeros(n,n),
            Z = arma::zeros(n,n);
  
  double alpha = (n*1.)/K;
  

  int t = 0;
  bool CONVERGED = false;
  while (!CONVERGED && t<T) {
    Xold = X;
    X = projAXB( 0.5*(Z-U+Y-V+As_rescaled), alpha, n);
    Z = max(X+U, arma::zeros(n,n));
    Y = projToSDC(X+V);
    U = U+X-Z;
    V = V+X-Y;
   
    delta(t) = norm(X-Xold,2);
    CONVERGED = delta(t) < tol;
    
    if ((t+1) % report_interval == 0) {
      Rprintf("%4d | %15e\n", t+1, delta(t));  
    }
    
    t++;
  }
  
  return List::create(
      _["X"]=X,
      _["delta"]=delta,
      _["T_term"]=t
  );
}


arma::mat projToSDC(arma::mat M) {
  int n = M.n_rows;

  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, M);
  
  for (int i=0; i < eigval.n_elem; i++){
    if ( eigval(i) < 0 ){ 
      eigval(i) = 0;
    }
  }
  
  M = eigvec * arma::diagmat(eigval) * eigvec.t();
  // arma::vec x = arma::eig_sym(M);
  // std::cout << x(3);
  return M;
}


arma::mat projAXB(arma::mat X0, double alpha, int n) {
//   arma::vec b (2*n);
//   b.ones();
  arma::vec b = arma::ones(2*n);

  b(arma::span(0,n-1)) = 2*(alpha-1)*arma::ones(n);
  return X0 - Acs( Pinv( Ac(X0, n)-b,n ), n);
}


arma::mat Acs(arma::vec z, int n) {
  arma::vec mu = z.head(n);
  arma::vec nu = z.tail(n);
  arma::mat Z(n,n);
  
  for (int i=0; i < n; i++) {
    Z(i,i) = nu(i);
  }
  
  for (int i=0; i < n; i++) {
    for (int j=i+1; j < n; j++) {
        Z(i,j) = mu(i) + mu(j);
        Z(j,i) = Z(i,j);
    }
  }
  
  return Z;
}


arma::vec Pinv(arma::vec z, int n) {
  arma::vec mu = z.head(n);
  arma::vec nu = z.tail(n);
  
  return arma::join_vert( (1./(2*(n-2)))*(mu - arma::ones(n)*arma::sum(mu)/(2*n-2)), nu);
}


arma::mat Ac(arma::mat X, int n) {
  return arma::join_vert( 2*(X - arma::diagmat( X ))*arma::ones(n), arma::diagvec(X) );
}

