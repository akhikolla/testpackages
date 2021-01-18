//Author: Manos Papadakis

// This file was generated by compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <vector>
#include "Rfast.h"

using namespace Rcpp;
using namespace std;

static NumericVector negative_and_positive_min(NumericVector &x){
    double mnp=INT_MAX,mnn=-1,v;
    for(auto xx=x.begin();xx!=x.end();++xx){
        v=*xx;
        if(v<0){
            if(v<mnn){
                mnn=v;
            }
        }else{
            if(v<mnp){
                mnp=v;
            }
        }
    }
    return NumericVector::create(mnn,mnp);
}

static NumericVector negative_and_positive_max(NumericVector &x){
    double mxp=0,mxn=INT_MIN,v;
    for(auto xx=x.begin();xx!=x.end();++xx){
        v=*xx;
        if(v<0){
            if(v>mxn){
                mxn=v;
            }
        }else{
            if(v>mxp){
                mxp=v;
            }
        }
    }
    return NumericVector::create(mxn,mxp);
}

static NumericVector negative_and_positive_min_max(NumericVector &x){
    double mnp=INT_MAX,mnn=-1,mxp=0,mxn=INT_MIN,v;
    for(auto xx=x.begin();xx!=x.end();++xx){
        v=*xx;
        if(v<0){
            if(v<mnn)
                mnn=v;
            else if(v>mxn)
                mxn=v;
        }else{
            if(v>mxp)
                mxp=v;
            else if(v<mnp)
                mnp=v;
        }
    }
    return NumericVector::create(mnn,mxn,mnp,mxp);
}

NumericVector negative(NumericVector x,string method = "min"){
    NumericVector s=0;
    if(method == "min"){
        s = negative_or_positive<mless<double> ,mless<double>>(x); 
    }else if(method == "max"){
        s = negative_or_positive<mless<double> ,mgreater<double>>(x);
    }else if(method == "min.max"){
        s = negative_or_positive_min_max<mless<double>>(x); 
    }else{
        stop("Error: Unsupported method.");
    }
    return s;
}

NumericVector positive(NumericVector x,string method = "min"){
    NumericVector s=0;
    if(method == "min"){
        s = negative_or_positive<mgreater<double> ,mless<double>>(x); 
    }else if(method == "max"){
        s = negative_or_positive<mgreater<double> ,mgreater<double>>(x);
    }else if(method == "min.max"){
        s = negative_or_positive_min_max<mgreater<double>>(x);
    }else{
        stop("Error: Unsupported method.");
    }
    return s;
}

NumericVector positive_negative(NumericVector x,string method = "min"){
    NumericVector s;
    if(method == "min"){
        s = negative_and_positive_min(x); 
    }else if(method == "max"){
        s = negative_and_positive_max(x);
    }else if(method == "min.max"){
        s = negative_and_positive_min_max(x);
    }else{
        stop("Error: Unsupported method.");
    }
    return s;
}

RcppExport SEXP Rfast_negative(SEXP xSEXP,SEXP methodSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    __result = negative(x,method);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_positive(SEXP xSEXP,SEXP methodSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    __result = positive(x,method);
    return __result;
END_RCPP
}

RcppExport SEXP Rfast_positive_negative(SEXP xSEXP,SEXP methodSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< string >::type method(methodSEXP);
    __result = positive_negative(x,method);
    return __result;
END_RCPP
}
