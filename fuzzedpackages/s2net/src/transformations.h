// Tranformations

// #include "RcppArmadillo.h"

#define TYPE_TRANSFORM_JT 0
#define TYPE_TRANSFORM_ExtJT 1


arma::mat transform_ExtJT(const arma::mat & X, double gamma2, double gamma3){
    int n = X.n_rows;
    int p = X.n_cols;
    
    // center X
    arma::rowvec u = arma::mean(X, 0);
    arma::mat X_u(n, p);
    for(int col = 0; col < p; col++){
        X_u.col(col) = X.col(col) - u(col);
    }
    
    arma::mat U;
    arma::vec s;
    arma::mat V;
    
    arma::svd(U, s, V, X_u);
    U = U.cols(0, std::min(n, p) - 1);

    return std::sqrt(gamma2)*U*arma::diagmat(1.0/arma::sqrt(s % s + gamma2))*U.t()*X_u + gamma3*arma::ones(n)*u;
}


arma::mat transform_JT( const arma::mat & X, double gamma2){
    int n = X.n_rows;
    int p = X.n_cols;
    arma::mat U;
    arma::vec s;
    arma::mat V;
    
    arma::svd(U, s, V, X);
    U = U.cols(0, std::min(n, p)-1);
    return std::sqrt(gamma2)*U*arma::diagmat(1.0/arma::sqrt(s % s + gamma2))*U.t()*X;
    //return arma::diagmat(1.0/arma::sqrt(s % s + gamma))*U.t()*X;
}
