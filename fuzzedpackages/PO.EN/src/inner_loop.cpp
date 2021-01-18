//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>

//'@name sum_compute_single_rcpp
//'@title Softthresholding computation function
//‘@usage sum_compute_single_rcpp(y, start_beta, ytx, XtX_reduce, lambda, alpha, xtx,  i_iter)
//'Compute the maximization in the EM-algorithm
//'@return updated coefficients vector

using namespace Rcpp ;
double softholding(const double & z,const double & g ){
   // std::cout<<"abs(z): "<<abs(z)<<std::endl;
    if( fabs(z)>g){
        if(z>0){
            return(z-g);
        }else{
            return(z+g);
        }
    }else{
        return(0);
    }
}

// [[Rcpp::export]]
Rcpp::List sum_compute_single_rcpp(const arma::vec &y, const arma::vec & start_beta,const arma::vec & ytx,const arma::mat & XtX_reduce,const double &lambda,const double & alpha,const arma::vec & xtx, unsigned int & i_iter){
     double N=y.n_rows;
//  std::cout<<N<<std::endl;
    unsigned int p=XtX_reduce.n_cols;
    arma::vec xx(1);
    arma::vec reduce_beta(p-1);
     arma::vec working_beta=start_beta;
  //  std::cout<<p<<std::endl;
  //   std::cout<<p<<std::endl;
   //w std::cout<<"lambda: "<<lambda<<std::endl;
//  std::cout<<XtX_reduce.submat(10,10,p-2,10)<<std::endl;
    for(unsigned int i=0;i<i_iter;i++){
    working_beta[0]=mean(y);
    for(unsigned int j=1; j<p-1;j++){
        reduce_beta.subvec(0,j-1)=working_beta.subvec(0,j-1 );
        reduce_beta.subvec(j,p-2)=working_beta.subvec(j+1,p-1);
        xx=trans(XtX_reduce.submat(0,j,p-2,j))*reduce_beta;
        working_beta[j]=softholding ((ytx[j]-xx[0])/N,4*lambda*alpha )/(xtx[j]/N+4*lambda*(1.00-alpha));
        
      
    }
    reduce_beta=working_beta.subvec(0,p-2);
    xx=trans(XtX_reduce.submat(0,p-1,p-2,p-1))*reduce_beta;
    working_beta[p-1]=softholding((ytx[p-1]-xx[0])/N,4*lambda*alpha)/(xtx[p-1]/N+4*lambda*(1.00-alpha));
    
    }
    
    return(Rcpp::List::create(Rcpp::Named("beta")=working_beta));
   

}
