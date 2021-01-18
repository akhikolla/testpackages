// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvn_internal
Rcpp::List mvn_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd covM, bool useLog2, int N);
RcppExport SEXP _tlrmvnmvt_mvn_internal(SEXP aSEXP, SEXP bSEXP, SEXP covMSEXP, SEXP useLog2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type covM(covMSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(mvn_internal(a, b, covM, useLog2, N));
    return rcpp_result_gen;
END_RCPP
}
// mvn_internal2
Rcpp::List mvn_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, int N);
RcppExport SEXP _tlrmvnmvt_mvn_internal2(SEXP aSEXP, SEXP bSEXP, SEXP geomSEXP, SEXP kernelTypeSEXP, SEXP paraSEXP, SEXP nuggetSEXP, SEXP useLog2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type kernelType(kernelTypeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type para(paraSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(mvn_internal2(a, b, geom, kernelType, para, nugget, useLog2, N));
    return rcpp_result_gen;
END_RCPP
}
// tlrmvn_internal
Rcpp::List tlrmvn_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd covM, bool useLog2, int m, double epsl, int N);
RcppExport SEXP _tlrmvnmvt_tlrmvn_internal(SEXP aSEXP, SEXP bSEXP, SEXP covMSEXP, SEXP useLog2SEXP, SEXP mSEXP, SEXP epslSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type covM(covMSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsl(epslSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(tlrmvn_internal(a, b, covM, useLog2, m, epsl, N));
    return rcpp_result_gen;
END_RCPP
}
// tlrmvn_internal2
Rcpp::List tlrmvn_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, int m, double epsl, int N);
RcppExport SEXP _tlrmvnmvt_tlrmvn_internal2(SEXP aSEXP, SEXP bSEXP, SEXP geomSEXP, SEXP kernelTypeSEXP, SEXP paraSEXP, SEXP nuggetSEXP, SEXP useLog2SEXP, SEXP mSEXP, SEXP epslSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type kernelType(kernelTypeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type para(paraSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsl(epslSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(tlrmvn_internal2(a, b, geom, kernelType, para, nugget, useLog2, m, epsl, N));
    return rcpp_result_gen;
END_RCPP
}
// mvt_internal
Rcpp::List mvt_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd mu, double nu, Eigen::MatrixXd covM, bool useLog2, int N);
RcppExport SEXP _tlrmvnmvt_mvt_internal(SEXP aSEXP, SEXP bSEXP, SEXP muSEXP, SEXP nuSEXP, SEXP covMSEXP, SEXP useLog2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type covM(covMSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(mvt_internal(a, b, mu, nu, covM, useLog2, N));
    return rcpp_result_gen;
END_RCPP
}
// mvt_internal2
Rcpp::List mvt_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd mu, double nu, Eigen::MatrixXd geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, int N);
RcppExport SEXP _tlrmvnmvt_mvt_internal2(SEXP aSEXP, SEXP bSEXP, SEXP muSEXP, SEXP nuSEXP, SEXP geomSEXP, SEXP kernelTypeSEXP, SEXP paraSEXP, SEXP nuggetSEXP, SEXP useLog2SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type kernelType(kernelTypeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type para(paraSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(mvt_internal2(a, b, mu, nu, geom, kernelType, para, nugget, useLog2, N));
    return rcpp_result_gen;
END_RCPP
}
// tlrmvt_internal
Rcpp::List tlrmvt_internal(Eigen::VectorXd a, Eigen::VectorXd b, double nu, Eigen::VectorXd mu, Eigen::MatrixXd covM, bool useLog2, int m, double epsl, int N);
RcppExport SEXP _tlrmvnmvt_tlrmvt_internal(SEXP aSEXP, SEXP bSEXP, SEXP nuSEXP, SEXP muSEXP, SEXP covMSEXP, SEXP useLog2SEXP, SEXP mSEXP, SEXP epslSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type covM(covMSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsl(epslSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(tlrmvt_internal(a, b, nu, mu, covM, useLog2, m, epsl, N));
    return rcpp_result_gen;
END_RCPP
}
// tlrmvt_internal2
Rcpp::List tlrmvt_internal2(Eigen::VectorXd a, Eigen::VectorXd b, double nu, Eigen::VectorXd mu, Eigen::MatrixXd geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, int m, double epsl, int N);
RcppExport SEXP _tlrmvnmvt_tlrmvt_internal2(SEXP aSEXP, SEXP bSEXP, SEXP nuSEXP, SEXP muSEXP, SEXP geomSEXP, SEXP kernelTypeSEXP, SEXP paraSEXP, SEXP nuggetSEXP, SEXP useLog2SEXP, SEXP mSEXP, SEXP epslSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type mu(muSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type kernelType(kernelTypeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type para(paraSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    Rcpp::traits::input_parameter< bool >::type useLog2(useLog2SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type epsl(epslSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(tlrmvt_internal2(a, b, nu, mu, geom, kernelType, para, nugget, useLog2, m, epsl, N));
    return rcpp_result_gen;
END_RCPP
}
// zorder
Eigen::VectorXi zorder(const Eigen::MatrixXd& geom);
RcppExport SEXP _tlrmvnmvt_zorder(SEXP geomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type geom(geomSEXP);
    rcpp_result_gen = Rcpp::wrap(zorder(geom));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tlrmvnmvt_mvn_internal", (DL_FUNC) &_tlrmvnmvt_mvn_internal, 5},
    {"_tlrmvnmvt_mvn_internal2", (DL_FUNC) &_tlrmvnmvt_mvn_internal2, 8},
    {"_tlrmvnmvt_tlrmvn_internal", (DL_FUNC) &_tlrmvnmvt_tlrmvn_internal, 7},
    {"_tlrmvnmvt_tlrmvn_internal2", (DL_FUNC) &_tlrmvnmvt_tlrmvn_internal2, 10},
    {"_tlrmvnmvt_mvt_internal", (DL_FUNC) &_tlrmvnmvt_mvt_internal, 7},
    {"_tlrmvnmvt_mvt_internal2", (DL_FUNC) &_tlrmvnmvt_mvt_internal2, 10},
    {"_tlrmvnmvt_tlrmvt_internal", (DL_FUNC) &_tlrmvnmvt_tlrmvt_internal, 9},
    {"_tlrmvnmvt_tlrmvt_internal2", (DL_FUNC) &_tlrmvnmvt_tlrmvt_internal2, 12},
    {"_tlrmvnmvt_zorder", (DL_FUNC) &_tlrmvnmvt_zorder, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_tlrmvnmvt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
