#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

const double log2pi = std::log(2.0 * M_PI);

//' Function to calculate the similarity matrix based on the
//' cluster membership indicator of each iteration.
//' @param mat A matrix of cluster membership indicators.
//' @return returns a similarity matrix.
//' @export
//' @examples
//' n = 90 ## number of subjects
//' iters = 200 ## number of iterations
//' ## matrix of cluster membership indicators
//' ## perfect clustering with three clusters
//' mat = matrix(rep(1:3,each=n/3),nrow=n,ncol=iters)
//' sim = calSim(t(mat))
//' ## plot similarity matrix
//' x <- rep(1:n,times=n)
//' y <- rep(1:n,each=n)
//' z <- as.vector(sim)
//' levelplot(z~x*y,col.regions=rev(gray.colors(n^2)), xlab = "Subject ID",ylab = "Subject ID")
// [[Rcpp::export]]
arma::mat calSim(arma::mat mat) {
    arma::mat ans(mat.n_cols,mat.n_cols);
    for(int i=0;i < mat.n_cols-1;i++){
        for(int j =i+1; j < mat.n_cols; j++){
            ans(i,j) = sum(mat.col(i) == mat.col(j));
            ans(j,i) = ans(i,j);
        }
    }
    ans = ans/mat.n_rows;
    ans.diag().ones();
    return ans;
}

//' Internal function to sample from multivariate normal distribution.
//' @keywords internal
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//' Internal function to sample Lambda
//' @keywords internal
// [[Rcpp::export]]
arma::mat samLamV3Cpp(arma::mat A, arma::mat eta, arma::vec sig, double sigl,
                      arma::mat lam){
    for(int i=1; i < lam.n_rows; i++){

        if(i < lam.n_cols){
            arma::mat I;
            I.eye(i,i);
            arma::mat tmpSig = arma::inv( 1.0/sig[i]*eta.rows(0,i-1)*arma::trans(eta.rows(0,i-1))
                                              + 1.0/sigl*I);
            arma::vec tmpMu = 1.0/sig[i]*tmpSig*eta.rows(0,i-1)*(A.col(i)- arma::trans(eta.row(i)) ) ;
            lam(arma::span(i,i),arma::span(0,i-1)) = mvrnormArma(1, tmpMu, tmpSig);
        }else{
            arma::mat I;
            I.eye(lam.n_cols,lam.n_cols);
            arma::mat tmpSig = arma::inv( 1.0/sig[i]*eta*arma::trans(eta)
                                              + 1.0/sigl*I);
            arma::vec tmpMu = 1.0/sig[i]*tmpSig*eta*A.col(i);
            lam.row(i) = mvrnormArma(1, tmpMu, tmpSig);
        }
    }
    return(lam);
}

//' Internal function to calculate the density of multivariate normal distribution.
//' @keywords internal
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) {
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
    }

    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

//' Internal function to find matched index.
//' @keywords internal
// [[Rcpp::export]]
arma::uvec myfind(IntegerVector evec, int e){
    arma::uvec u(sum(evec==e));
    int i,j=0;
    for(i=0;i < evec.size() ; i++ ){
        if(evec(i)==e){u(j)=i;j++;}
    }
    return(u);
}

//' Internal function to sample cluster membership indicator
//' @keywords internal
// [[Rcpp::export]]
IntegerVector polyurncpp(IntegerVector e,arma::vec muA0,arma::mat sigma0,arma::mat A,
                arma::mat sigma,arma::mat sigmaInv,arma::mat sigma0Inv, double c) {

    int i,j;
    arma::mat tmpSig = sigma0;
    arma::vec tmpMu = muA0;

    for(i = 0;i < e.size(); i++){
        e[i] = max(e) + 1;
        IntegerVector eset = unique(e);
        NumericVector pp(eset.size());
        for(j = 0; j < eset.size(); j++){
            if(eset[j] == max(e)){
                tmpSig = sigma0;
                tmpMu = muA0;
                pp[j] = dmvnrm_arma(A.row(i),arma::trans(tmpMu),tmpSig+sigma,true)[0] + log(c);
                continue;
            }
            arma::uvec ind = myfind(e,eset[j]);
            arma::vec unit(ind.size());
            unit.ones();
            tmpSig = arma::inv(sigma0Inv + ind.size()*sigmaInv);
            tmpMu = tmpSig*(sigmaInv*arma::trans(A.rows(ind))*unit + sigma0Inv*muA0);
            pp[j] = dmvnrm_arma(A.row(i),arma::trans(tmpMu),tmpSig+sigma,true)[0] +
                log(ind.size()*1.0);
        }
        //Rf_PrintValue(wrap(eset));
        //Rf_PrintValue(wrap(pp));
        pp = pp - max(pp);
        pp = exp(pp - log(sum(exp(pp))));
        pp = pp/sum(pp);
        e[i] = RcppArmadillo::sample(eset,1,true,pp)[0];
    }

    return e;
}

