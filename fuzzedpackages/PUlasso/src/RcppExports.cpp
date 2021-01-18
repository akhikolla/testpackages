// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// LU_dense_cpp
Rcpp::List LU_dense_cpp(Eigen::Map<Eigen::MatrixXd> X_, Eigen::VectorXd& z_, Eigen::VectorXd& icoef_, Eigen::ArrayXd& gsize_, Eigen::ArrayXd& pen_, Eigen::ArrayXd& lambdaseq_, bool user_lambdaseq_, int pathLength_, double lambdaMinRatio_, double pi_, int max_nUpdates_, int maxit_, Eigen::VectorXd& wei_, bool weiOption_, double tol_, double inner_tol_, bool useStrongSet_, bool verbose_, double stepSize_, double stepSizeAdj_, int batchSize_, int updateFreq_, std::vector<double> samplingProbabilities_, bool useLipschitz_, std::string method_, int trace_, bool skipFitting_);
RcppExport SEXP _PUlasso_LU_dense_cpp(SEXP X_SEXP, SEXP z_SEXP, SEXP icoef_SEXP, SEXP gsize_SEXP, SEXP pen_SEXP, SEXP lambdaseq_SEXP, SEXP user_lambdaseq_SEXP, SEXP pathLength_SEXP, SEXP lambdaMinRatio_SEXP, SEXP pi_SEXP, SEXP max_nUpdates_SEXP, SEXP maxit_SEXP, SEXP wei_SEXP, SEXP weiOption_SEXP, SEXP tol_SEXP, SEXP inner_tol_SEXP, SEXP useStrongSet_SEXP, SEXP verbose_SEXP, SEXP stepSize_SEXP, SEXP stepSizeAdj_SEXP, SEXP batchSize_SEXP, SEXP updateFreq_SEXP, SEXP samplingProbabilities_SEXP, SEXP useLipschitz_SEXP, SEXP method_SEXP, SEXP trace_SEXP, SEXP skipFitting_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type icoef_(icoef_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type gsize_(gsize_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type pen_(pen_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type lambdaseq_(lambdaseq_SEXP);
    Rcpp::traits::input_parameter< bool >::type user_lambdaseq_(user_lambdaseq_SEXP);
    Rcpp::traits::input_parameter< int >::type pathLength_(pathLength_SEXP);
    Rcpp::traits::input_parameter< double >::type lambdaMinRatio_(lambdaMinRatio_SEXP);
    Rcpp::traits::input_parameter< double >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< int >::type max_nUpdates_(max_nUpdates_SEXP);
    Rcpp::traits::input_parameter< int >::type maxit_(maxit_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type wei_(wei_SEXP);
    Rcpp::traits::input_parameter< bool >::type weiOption_(weiOption_SEXP);
    Rcpp::traits::input_parameter< double >::type tol_(tol_SEXP);
    Rcpp::traits::input_parameter< double >::type inner_tol_(inner_tol_SEXP);
    Rcpp::traits::input_parameter< bool >::type useStrongSet_(useStrongSet_SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_(verbose_SEXP);
    Rcpp::traits::input_parameter< double >::type stepSize_(stepSize_SEXP);
    Rcpp::traits::input_parameter< double >::type stepSizeAdj_(stepSizeAdj_SEXP);
    Rcpp::traits::input_parameter< int >::type batchSize_(batchSize_SEXP);
    Rcpp::traits::input_parameter< int >::type updateFreq_(updateFreq_SEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type samplingProbabilities_(samplingProbabilities_SEXP);
    Rcpp::traits::input_parameter< bool >::type useLipschitz_(useLipschitz_SEXP);
    Rcpp::traits::input_parameter< std::string >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< int >::type trace_(trace_SEXP);
    Rcpp::traits::input_parameter< bool >::type skipFitting_(skipFitting_SEXP);
    rcpp_result_gen = Rcpp::wrap(LU_dense_cpp(X_, z_, icoef_, gsize_, pen_, lambdaseq_, user_lambdaseq_, pathLength_, lambdaMinRatio_, pi_, max_nUpdates_, maxit_, wei_, weiOption_, tol_, inner_tol_, useStrongSet_, verbose_, stepSize_, stepSizeAdj_, batchSize_, updateFreq_, samplingProbabilities_, useLipschitz_, method_, trace_, skipFitting_));
    return rcpp_result_gen;
END_RCPP
}
// LU_sparse_cpp
Rcpp::List LU_sparse_cpp(Eigen::SparseMatrix<double>& X_, Eigen::VectorXd& z_, Eigen::VectorXd& icoef_, Eigen::ArrayXd& gsize_, Eigen::ArrayXd& pen_, Eigen::ArrayXd& lambdaseq_, bool user_lambdaseq_, int pathLength_, double lambdaMinRatio_, double pi_, int max_nUpdates_, int maxit_, Eigen::VectorXd& wei_, bool weiOption_, double tol_, double inner_tol_, bool useStrongSet_, bool verbose_, double stepSize_, double stepSizeAdj_, int batchSize_, int updateFreq_, std::vector<double> samplingProbabilities_, bool useLipschitz_, std::string method_, int trace_, bool skipFitting_);
RcppExport SEXP _PUlasso_LU_sparse_cpp(SEXP X_SEXP, SEXP z_SEXP, SEXP icoef_SEXP, SEXP gsize_SEXP, SEXP pen_SEXP, SEXP lambdaseq_SEXP, SEXP user_lambdaseq_SEXP, SEXP pathLength_SEXP, SEXP lambdaMinRatio_SEXP, SEXP pi_SEXP, SEXP max_nUpdates_SEXP, SEXP maxit_SEXP, SEXP wei_SEXP, SEXP weiOption_SEXP, SEXP tol_SEXP, SEXP inner_tol_SEXP, SEXP useStrongSet_SEXP, SEXP verbose_SEXP, SEXP stepSize_SEXP, SEXP stepSizeAdj_SEXP, SEXP batchSize_SEXP, SEXP updateFreq_SEXP, SEXP samplingProbabilities_SEXP, SEXP useLipschitz_SEXP, SEXP method_SEXP, SEXP trace_SEXP, SEXP skipFitting_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type icoef_(icoef_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type gsize_(gsize_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type pen_(pen_SEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd& >::type lambdaseq_(lambdaseq_SEXP);
    Rcpp::traits::input_parameter< bool >::type user_lambdaseq_(user_lambdaseq_SEXP);
    Rcpp::traits::input_parameter< int >::type pathLength_(pathLength_SEXP);
    Rcpp::traits::input_parameter< double >::type lambdaMinRatio_(lambdaMinRatio_SEXP);
    Rcpp::traits::input_parameter< double >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< int >::type max_nUpdates_(max_nUpdates_SEXP);
    Rcpp::traits::input_parameter< int >::type maxit_(maxit_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type wei_(wei_SEXP);
    Rcpp::traits::input_parameter< bool >::type weiOption_(weiOption_SEXP);
    Rcpp::traits::input_parameter< double >::type tol_(tol_SEXP);
    Rcpp::traits::input_parameter< double >::type inner_tol_(inner_tol_SEXP);
    Rcpp::traits::input_parameter< bool >::type useStrongSet_(useStrongSet_SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose_(verbose_SEXP);
    Rcpp::traits::input_parameter< double >::type stepSize_(stepSize_SEXP);
    Rcpp::traits::input_parameter< double >::type stepSizeAdj_(stepSizeAdj_SEXP);
    Rcpp::traits::input_parameter< int >::type batchSize_(batchSize_SEXP);
    Rcpp::traits::input_parameter< int >::type updateFreq_(updateFreq_SEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type samplingProbabilities_(samplingProbabilities_SEXP);
    Rcpp::traits::input_parameter< bool >::type useLipschitz_(useLipschitz_SEXP);
    Rcpp::traits::input_parameter< std::string >::type method_(method_SEXP);
    Rcpp::traits::input_parameter< int >::type trace_(trace_SEXP);
    Rcpp::traits::input_parameter< bool >::type skipFitting_(skipFitting_SEXP);
    rcpp_result_gen = Rcpp::wrap(LU_sparse_cpp(X_, z_, icoef_, gsize_, pen_, lambdaseq_, user_lambdaseq_, pathLength_, lambdaMinRatio_, pi_, max_nUpdates_, maxit_, wei_, weiOption_, tol_, inner_tol_, useStrongSet_, verbose_, stepSize_, stepSizeAdj_, batchSize_, updateFreq_, samplingProbabilities_, useLipschitz_, method_, trace_, skipFitting_));
    return rcpp_result_gen;
END_RCPP
}
// deviances_dense_cpp
Eigen::MatrixXd deviances_dense_cpp(Eigen::MatrixXd& coefMat_, Eigen::Map<Eigen::MatrixXd>& X_, Eigen::VectorXd& z_, double pi_, const Eigen::VectorXd& wei_, bool weiOption_);
RcppExport SEXP _PUlasso_deviances_dense_cpp(SEXP coefMat_SEXP, SEXP X_SEXP, SEXP z_SEXP, SEXP pi_SEXP, SEXP wei_SEXP, SEXP weiOption_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type coefMat_(coefMat_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< double >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wei_(wei_SEXP);
    Rcpp::traits::input_parameter< bool >::type weiOption_(weiOption_SEXP);
    rcpp_result_gen = Rcpp::wrap(deviances_dense_cpp(coefMat_, X_, z_, pi_, wei_, weiOption_));
    return rcpp_result_gen;
END_RCPP
}
// deviances_sparse_cpp
Eigen::MatrixXd deviances_sparse_cpp(Eigen::MatrixXd& coefMat_, Eigen::SparseMatrix<double>& X_, Eigen::VectorXd& z_, double pi_, const Eigen::VectorXd& wei_, bool weiOption_);
RcppExport SEXP _PUlasso_deviances_sparse_cpp(SEXP coefMat_SEXP, SEXP X_SEXP, SEXP z_SEXP, SEXP pi_SEXP, SEXP wei_SEXP, SEXP weiOption_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type coefMat_(coefMat_SEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double>& >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< double >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wei_(wei_SEXP);
    Rcpp::traits::input_parameter< bool >::type weiOption_(weiOption_SEXP);
    rcpp_result_gen = Rcpp::wrap(deviances_sparse_cpp(coefMat_, X_, z_, pi_, wei_, weiOption_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PUlasso_LU_dense_cpp", (DL_FUNC) &_PUlasso_LU_dense_cpp, 27},
    {"_PUlasso_LU_sparse_cpp", (DL_FUNC) &_PUlasso_LU_sparse_cpp, 27},
    {"_PUlasso_deviances_dense_cpp", (DL_FUNC) &_PUlasso_deviances_dense_cpp, 6},
    {"_PUlasso_deviances_sparse_cpp", (DL_FUNC) &_PUlasso_deviances_sparse_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_PUlasso(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
