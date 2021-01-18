#include <Rcpp.h>

using namespace Rcpp;

// vectorSum
double VectorSum(NumericVector x);
RcppExport SEXP GiniDistance_VectorSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(VectorSum(x));
    return __result;
END_RCPP
}

// order a vector
IntegerVector orderc(NumericVector x);
RcppExport SEXP GiniDistance_orderc(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(orderc(x));
    return __result;
END_RCPP
}

// Rcpp distance matrix
NumericMatrix rcpp_Eu_distance(NumericMatrix mat);
RcppExport SEXP GiniDistance_rcpp_Eu_distance(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    __result = Rcpp::wrap(rcpp_Eu_distance(mat));
    return __result;
END_RCPP
}


// Rcpp vector GiniDistance 
double rcpp_covg(NumericVector x, NumericVector y);
RcppExport SEXP GiniDistance_rcpp_covg(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(rcpp_covg(x,y));
    return __result;
END_RCPP
}

// Rcpp vector GiniDistance 
double rcpp_covg_alpha(NumericVector x, NumericVector y, double alpha);
RcppExport SEXP GiniDistance_rcpp_covg_alpha(SEXP xSEXP, SEXP ySEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(rcpp_covg_alpha(x,y, alpha));
    return __result;
END_RCPP
}


// Rcpp Kernel Distance matrix
NumericMatrix rcpp_Kernel_Distance(NumericVector vec, double sigma);
RcppExport SEXP GiniDistance_rcpp_Kernel_Distance(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(rcpp_Kernel_Distance(x,sigma));
    return __result;
END_RCPP
}

//RcppKernelGiniDistance for vector
double rcpp_Kcovg(NumericVector x, NumericVector y, double sigma);
RcppExport SEXP GiniDistance_rcpp_Kcovg(SEXP xSEXP, SEXP ySEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(rcpp_Kcovg(x,y, sigma));
    return __result;
END_RCPP
}

//RcppKernelGiniDistance for matrix
double Rcpp_KCovg(NumericMatrix x, NumericVector y, double sigma); 
RcppExport SEXP GiniDistance_Rcpp_KCovg(SEXP xSEXP, SEXP ySEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_KCovg(x,y, sigma));
    return __result;
END_RCPP
}

//RcppKernelGiniDistance for matrix raise to power alpha
double Rcpp_KCovg_Alpha(NumericMatrix x, NumericVector y, double sigma, double alpha); 
RcppExport SEXP GiniDistance_Rcpp_KCovg_Alpha(SEXP xSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_KCovg_Alpha(x,y, sigma, alpha));
    return __result;
END_RCPP
}


// Rcpp matrix GiniDistance  
double Rcpp_Covg(NumericMatrix x, NumericVector y);
RcppExport SEXP GiniDistance_Rcpp_Covg(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(Rcpp_Covg(x,y));
    return __result;
END_RCPP
}

// Rcpp matrix GiniDistance raise to power alpha 
double Rcpp_Covg_Alpha(NumericMatrix x, NumericVector y, double alpha);
RcppExport SEXP GiniDistance_Rcpp_Covg_Alpha(SEXP xSEXP, SEXP ySEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_Covg_Alpha(x,y, alpha));
    return __result;
END_RCPP
}


// Rcpp GiniCovariance Multiclasses  
double Rcpp_gCov(NumericMatrix x);
RcppExport SEXP GiniDistance_Rcpp_gCov(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(Rcpp_gCov(x));
    return __result;
END_RCPP
}

// Rcpp GiniCovariance Multiclasses to power alpha  
double Rcpp_gCov_Alpha(NumericMatrix x, double alpha);
RcppExport SEXP GiniDistance_Rcpp_gCov_Alpha(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_gCov_Alpha(x, alpha));
    return __result;
END_RCPP
}


// Rcpp Kernel GiniCovariance Multiclasses to power alpha  
double Rcpp_KgCov_Alpha(NumericMatrix x, double sigma, double alpha);
RcppExport SEXP GiniDistance_Rcpp_KgCov_Alpha(SEXP xSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_KgCov_Alpha(x, sigma, alpha));
    return __result;
END_RCPP
}

// Rcpp Kernel GiniCorrelation Multiclasses to power alpha  
double Rcpp_KgCor_Alpha(NumericMatrix x, double sigma, double alpha);
RcppExport SEXP GiniDistance_Rcpp_KgCor_Alpha(SEXP xSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_KgCor_Alpha(x, sigma, alpha));
    return __result;
END_RCPP
}

// Rcpp GiniCovariance Multiclasses to power alpha  
double Rcpp_gCor_Alpha(NumericMatrix x, double alpha);
RcppExport SEXP GiniDistance_Rcpp_gCor_Alpha(SEXP xSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_gCor_Alpha(x, alpha));
    return __result;
END_RCPP
}


// Rcpp GiniCorrelation Multiclasses  
double Rcpp_gCor(NumericMatrix x);
RcppExport SEXP GiniDistance_Rcpp_gCor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(Rcpp_gCor(x));
    return __result;
END_RCPP
}

//RcppKernelGiniCovariance
double Rcpp_KgCov(NumericMatrix x, double sigma);
RcppExport SEXP GiniDistance_Rcpp_KgCov(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_KgCov(x, sigma));
    return __result;
END_RCPP
}

//RcppKernelGiniCorrelation
double Rcpp_KgCor(NumericMatrix x, double sigma);
RcppExport SEXP GiniDistance_Rcpp_KgCor(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_KgCor(x, sigma));
    return __result;
END_RCPP
}

//Rcpp Jacknife variance
double Rcpp_HatV_gCov(NumericMatrix x);
RcppExport SEXP GiniDistance_Rcpp_HatV_gCov(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_gCov(x));
    return __result;
END_RCPP
}

//Rcpp Jacknife variance to power alpha
double Rcpp_HatV_gCov_Alpha(NumericMatrix x, double alpha);
RcppExport SEXP GiniDistance_Rcpp_HatV_gCov_Alpha(SEXP xSEXP, SEXP alphaSEXP ) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_gCov_Alpha(x, alpha));
    return __result;
END_RCPP
}

//Rcpp Jacknife variance of correlation to power alpha
double Rcpp_HatV_gCor_Alpha(NumericMatrix x, double alpha);
RcppExport SEXP GiniDistance_Rcpp_HatV_gCor_Alpha(SEXP xSEXP, SEXP alphaSEXP ) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_gCor_Alpha(x, alpha));
    return __result;
END_RCPP
}


//Rcpp Jacknife variance of Gini correlation
double Rcpp_HatV_gCor(NumericMatrix x);
RcppExport SEXP GiniDistance_Rcpp_HatV_gCor(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_gCor(x));
    return __result;
END_RCPP
}

//jackknife variance of Kernel Gini covariance  
double Rcpp_HatV_KgCov(NumericMatrix x, double sigma);
RcppExport SEXP GiniDistance_Rcpp_HatV_KgCov(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_KgCov(x, sigma));
    return __result;
END_RCPP
}

//jackknife variance of Kernel Gini covariance to power alpha  
double Rcpp_HatV_KgCov_Alpha(NumericMatrix x, double sigma, double alpha);
RcppExport SEXP GiniDistance_Rcpp_HatV_KgCov_Alpha(SEXP xSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_KgCov_Alpha(x, sigma, alpha));
    return __result;
END_RCPP
}

//jackknife variance of Kernel Gini correlation to power alpha  
double Rcpp_HatV_KgCor_Alpha(NumericMatrix x, double sigma, double alpha);
RcppExport SEXP GiniDistance_Rcpp_HatV_KgCor_Alpha(SEXP xSEXP, SEXP sigmaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_KgCor_Alpha(x, sigma, alpha));
    return __result;
END_RCPP
}


//jackknife variance of Kernel Gini correlation  
double Rcpp_HatV_KgCor(NumericMatrix x, double sigma);
RcppExport SEXP GiniDistance_Rcpp_HatV_KgCor(SEXP xSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    __result = Rcpp::wrap(Rcpp_HatV_KgCor(x, sigma));
    return __result;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"GiniDistance_VectorSum", (DL_FUNC) &GiniDistance_VectorSum, 1},
    {"GiniDistance_orderc", (DL_FUNC) &GiniDistance_orderc, 1},
    {"GiniDistance_rcpp_Eu_distance", (DL_FUNC) &GiniDistance_rcpp_Eu_distance, 1},
    {"GiniDistance_rcpp_Kernel_Distance", (DL_FUNC) &GiniDistance_rcpp_Kernel_Distance, 2},
    {"GiniDistance_rcpp_Kcovg", (DL_FUNC) &GiniDistance_rcpp_Kcovg, 3},
    {"GiniDistance_Rcpp_KCovg", (DL_FUNC) &GiniDistance_Rcpp_KCovg, 3},
    {"GiniDistance_Rcpp_KCovg_Alpha", (DL_FUNC) &GiniDistance_Rcpp_KCovg_Alpha, 4},
    {"GiniDistance_Rcpp_KgCov", (DL_FUNC) &GiniDistance_Rcpp_KgCov, 2},
    {"GiniDistance_Rcpp_KgCor", (DL_FUNC) &GiniDistance_Rcpp_KgCor, 2},
    {"GiniDistance_rcpp_covg", (DL_FUNC) &GiniDistance_rcpp_covg, 2},
    {"GiniDistance_rcpp_covg_alpha", (DL_FUNC) &GiniDistance_rcpp_covg_alpha, 3},
    {"GiniDistance_Rcpp_Covg", (DL_FUNC) &GiniDistance_Rcpp_Covg, 2},
    {"GiniDistance_Rcpp_Covg_Alpha", (DL_FUNC) &GiniDistance_Rcpp_Covg_Alpha, 3},
    {"GiniDistance_Rcpp_gCov_Alpha", (DL_FUNC) &GiniDistance_Rcpp_gCov_Alpha, 2},
    {"GiniDistance_Rcpp_KgCov_Alpha", (DL_FUNC) &GiniDistance_Rcpp_KgCov_Alpha, 3},
    {"GiniDistance_Rcpp_KgCor_Alpha", (DL_FUNC) &GiniDistance_Rcpp_KgCor_Alpha, 3},
    {"GiniDistance_Rcpp_gCor_Alpha", (DL_FUNC) &GiniDistance_Rcpp_gCor_Alpha, 2},
    {"GiniDistance_Rcpp_gCov", (DL_FUNC) &GiniDistance_Rcpp_gCov, 1},
    {"GiniDistance_Rcpp_gCor", (DL_FUNC) &GiniDistance_Rcpp_gCor, 1},
    {"GiniDistance_Rcpp_HatV_gCov", (DL_FUNC) &GiniDistance_Rcpp_HatV_gCov, 1},
    {"GiniDistance_Rcpp_HatV_gCov_Alpha", (DL_FUNC) &GiniDistance_Rcpp_HatV_gCov_Alpha, 2},
    {"GiniDistance_Rcpp_HatV_gCor_Alpha", (DL_FUNC) &GiniDistance_Rcpp_HatV_gCor_Alpha, 2},
    {"GiniDistance_Rcpp_HatV_gCor", (DL_FUNC) &GiniDistance_Rcpp_HatV_gCor, 1},
    {"GiniDistance_Rcpp_HatV_KgCov", (DL_FUNC) &GiniDistance_Rcpp_HatV_KgCov, 2},
    {"GiniDistance_Rcpp_HatV_KgCov_Alpha", (DL_FUNC) &GiniDistance_Rcpp_HatV_KgCov_Alpha, 3},
    {"GiniDistance_Rcpp_HatV_KgCor", (DL_FUNC) &GiniDistance_Rcpp_HatV_KgCor, 2},
    {"GiniDistance_Rcpp_HatV_KgCor_Alpha", (DL_FUNC) &GiniDistance_Rcpp_HatV_KgCor_Alpha, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_GiniDistance(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
