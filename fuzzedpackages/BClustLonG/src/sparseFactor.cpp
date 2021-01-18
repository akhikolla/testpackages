#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]

const double log2pi = std::log(2.0 * M_PI);

//' Internal function for sorting
//' @noRd
// [[Rcpp::export]]
IntegerVector stlSort(IntegerVector x) {
    IntegerVector y = clone(x);
    std::sort(y.begin(), y.end());
    return y;
}

//' Internal function for concatenating
//' @noRd
// [[Rcpp::export]]
IntegerVector myc(IntegerVector x, IntegerVector y) {
    IntegerVector z(x.size()+y.size());
    if(x.size()>0 and y.size() > 0){
        z[Range(0,x.size()-1)] = x;
        z[Range(x.size(),y.size()+x.size()-1)] = y;
    }else if(x.size()==0 and y.size() >0 ){
        z=y;
    }else if(y.size()==0 and x.size() >0){
        z=x;
    }else{}

    return stlSort(z);
}

//' Internal function for density function of MVN
//' @noRd
// [[Rcpp::export]]
arma::vec dmvnrmArma(arma::mat x,
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

//' Internal function for finding matched index
//' @noRd
// [[Rcpp::export]]
arma::uvec myfind(IntegerVector evec, int e){
    arma::uvec u(sum(evec==e));
    int i,j=0;
    for(i=0;i < evec.size() ; i++ ){
        if(evec(i)==e){u(j)=i;j++;}
    }
    return(u);
}

//' Internal function for generating MVN random variables
//' @noRd
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//' Internal function for sampling lambda in the factor analysis
//' @noRd
// [[Rcpp::export]]
arma::mat samLamV2Cpp(arma::mat A, arma::mat eta, arma::vec sig,
                      arma::mat lam, arma::mat phi, arma::rowvec tau){
    arma::mat D(lam.n_cols,lam.n_cols);
    D.eye();
    for(int i=0; i < lam.n_rows; i++){
        D.diag() = phi.row(i)%tau;
        arma::mat tmpSig = arma::inv( 1.0/sig[i]*eta*arma::trans(eta) + D);
        arma::vec tmpMu = 1.0/sig[i]*tmpSig*eta*A.col(i);
        lam.row(i) = mvrnormArma(1, tmpMu, tmpSig);
    }
    return(lam);
}

//' Internal function for sampling memberhsip indicators in sparse both
//' @noRd
// [[Rcpp::export]]
IntegerVector polyurncppBoth(IntegerVector e,arma::mat A,arma::vec muA0,
                             arma::mat sigmaA, arma::mat sigmaAInv,arma::mat sigmaA0,
                             arma::mat sigmaA0Inv, arma::mat B,arma::vec muB0,
                             arma::mat sigmaB,arma::mat sigmaBInv,arma::mat sigmaB0,
                             arma::mat sigmaB0Inv, double c) {

    int i,j;
    arma::mat tmpSigA = sigmaA0;
    arma::vec tmpMuA = muA0;
    arma::mat tmpSigB = sigmaB0;
    arma::vec tmpMuB = muB0;

    for(i = 0;i < e.size(); i++){
        e[i] = max(e) + 1;
        IntegerVector eset = unique(e);
        NumericVector pp(eset.size());
        for(j = 0; j < eset.size(); j++){
            if(eset[j] == max(e)){
                tmpSigA = sigmaA0;
                tmpMuA = muA0;
                tmpSigB = sigmaB0;
                tmpMuB = muB0;
                pp[j] = dmvnrmArma(A.row(i),arma::trans(tmpMuA),tmpSigA+sigmaA,true)[0]+
                    dmvnrmArma(B.row(i),arma::trans(tmpMuB),tmpSigB+sigmaB,true)[0] + log(c);
                continue;
            }
            arma::uvec ind = myfind(e,eset[j]);
            arma::vec unit(ind.size());
            unit.ones();
            tmpSigA = arma::inv_sympd(sigmaA0Inv + ind.size()*sigmaAInv);
            tmpMuA = tmpSigA*(sigmaAInv*arma::trans(A.rows(ind))*unit + sigmaA0Inv*muA0);
            tmpSigB = arma::inv_sympd(sigmaB0Inv + ind.size()*sigmaBInv);
            tmpMuB = tmpSigB*(sigmaBInv*arma::trans(B.rows(ind))*unit + sigmaB0Inv*muB0);

            pp[j] = dmvnrmArma(A.row(i),arma::trans(tmpMuA),tmpSigA+sigmaA,true)[0]+
                dmvnrmArma(B.row(i),arma::trans(tmpMuB),tmpSigB+sigmaB,true)[0] +
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


//' Internal function for sampling memberhsip indicators in sparse int
//' @noRd
// [[Rcpp::export]]
IntegerVector polyurncppInt(IntegerVector e,arma::vec muA0,arma::mat sigma0,
                            arma::mat A, arma::mat sigma,arma::mat sigmaInv,
                            arma::mat sigma0Inv, double c) {

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
                pp[j] = dmvnrmArma(A.row(i),arma::trans(tmpMu),tmpSig+sigma,true)[0] + log(c);
                continue;
            }
            arma::uvec ind = myfind(e,eset[j]);
            arma::vec unit(ind.size());
            unit.ones();
            tmpSig = arma::inv_sympd(sigma0Inv + ind.size()*sigmaInv);
            tmpMu = tmpSig*(sigmaInv*arma::trans(A.rows(ind))*unit + sigma0Inv*muA0);
            pp[j] = dmvnrmArma(A.row(i),arma::trans(tmpMu),tmpSig+sigma,true)[0] +
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


//' Function to calculate the similarity matrix based on the
//' cluster membership indicator of each iteration.
//' @param mat Matrix of cluster membership indicator from all iterations
//' @examples
//' n = 90 ##number of subjects
//' iters = 200 ##number of iterations
//' ## matrix of cluster membership indicators
//' ## perfect clustering with three clusters
//' mat = matrix(rep(1:3,each=n/3),nrow=n,ncol=iters)
//' sim = calSim(t(mat))
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

