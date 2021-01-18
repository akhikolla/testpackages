
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List IM4ECpp(Function oneIM4E,NumericMatrix train_xx,NumericVector train_yy,double epsilon=0.01,
    double sig=1,double lambda=1,int max_iter=10,bool removesmall=false){
    double p = train_xx.ncol();
    NumericVector w0 = NumericVector(p);
    for ( int i = 0; i < p;i++){
    w0[i] = 1/p;
    }
    double c0 = 0;
    double c_before = c0;
    double c_after = c0+1;
    NumericVector w_before = w0;
    NumericVector w_after = w0;
    NumericVector w;
    int iter = 0;
    if (removesmall == false){
    while ( (fabs(c_before - c_after)>epsilon)&&(iter< max_iter)  ){
    w_before = w_after;
    c_before = c_after;
    List res = oneIM4E(train_xx,train_yy,w_before,sig,lambda);
    w_after = res["w"];
    c_after = res["C"];
    iter++;
    }
    }else{
    while ( (fabs(c_before - c_after)>epsilon)&&(iter< max_iter)  ){
    w_before = w_after;
    c_before = c_after;
    List res = oneIM4E(train_xx,train_yy,w_before,sig,lambda);
    w_after = res["w"];
    for ( int j = 0 ; j < w_after.size(); j++ ){
    if(w_after[j] < 1/p){
    w_after[j] = 0;
    }
    }
    c_after = res["C"];
    iter++;
    }
    }
    List newList;
    newList["w"]=w_after;
    newList["c"]=c_after;
    newList["iter_num"]=iter;
    return newList;
    }
