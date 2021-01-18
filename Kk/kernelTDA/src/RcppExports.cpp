// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// svmtrain_R
List svmtrain_R(Eigen::MatrixXd x, int r, int c, Eigen::VectorXd y, int rowindex, int colindex, int svm_type, int kernel_type, int degree, double gamma, double coef0, double cost, double nu, Eigen::VectorXi weightlabels, Eigen::VectorXd weights, int nweights, double cache, double tolerance, double epsilon, int shrinking, int cross, int sparse, int probability, int nclasses, int nr, Eigen::VectorXi index, Eigen::VectorXi labels, Eigen::VectorXi nSV, Eigen::VectorXd rho, Eigen::VectorXd coefs, double sigma, Eigen::VectorXd probA, Eigen::VectorXd probB, Eigen::VectorXd cresults, double ctotal1, double ctotal2);
RcppExport SEXP _kernelTDA_svmtrain_R(SEXP xSEXP, SEXP rSEXP, SEXP cSEXP, SEXP ySEXP, SEXP rowindexSEXP, SEXP colindexSEXP, SEXP svm_typeSEXP, SEXP kernel_typeSEXP, SEXP degreeSEXP, SEXP gammaSEXP, SEXP coef0SEXP, SEXP costSEXP, SEXP nuSEXP, SEXP weightlabelsSEXP, SEXP weightsSEXP, SEXP nweightsSEXP, SEXP cacheSEXP, SEXP toleranceSEXP, SEXP epsilonSEXP, SEXP shrinkingSEXP, SEXP crossSEXP, SEXP sparseSEXP, SEXP probabilitySEXP, SEXP nclassesSEXP, SEXP nrSEXP, SEXP indexSEXP, SEXP labelsSEXP, SEXP nSVSEXP, SEXP rhoSEXP, SEXP coefsSEXP, SEXP sigmaSEXP, SEXP probASEXP, SEXP probBSEXP, SEXP cresultsSEXP, SEXP ctotal1SEXP, SEXP ctotal2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type rowindex(rowindexSEXP);
    Rcpp::traits::input_parameter< int >::type colindex(colindexSEXP);
    Rcpp::traits::input_parameter< int >::type svm_type(svm_typeSEXP);
    Rcpp::traits::input_parameter< int >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< double >::type cost(costSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type weightlabels(weightlabelsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type nweights(nweightsSEXP);
    Rcpp::traits::input_parameter< double >::type cache(cacheSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type shrinking(shrinkingSEXP);
    Rcpp::traits::input_parameter< int >::type cross(crossSEXP);
    Rcpp::traits::input_parameter< int >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< int >::type probability(probabilitySEXP);
    Rcpp::traits::input_parameter< int >::type nclasses(nclassesSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type nSV(nSVSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type probA(probASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type probB(probBSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type cresults(cresultsSEXP);
    Rcpp::traits::input_parameter< double >::type ctotal1(ctotal1SEXP);
    Rcpp::traits::input_parameter< double >::type ctotal2(ctotal2SEXP);
    rcpp_result_gen = Rcpp::wrap(svmtrain_R(x, r, c, y, rowindex, colindex, svm_type, kernel_type, degree, gamma, coef0, cost, nu, weightlabels, weights, nweights, cache, tolerance, epsilon, shrinking, cross, sparse, probability, nclasses, nr, index, labels, nSV, rho, coefs, sigma, probA, probB, cresults, ctotal1, ctotal2));
    return rcpp_result_gen;
END_RCPP
}
// svmpredict_R
List svmpredict_R(int decisionvalues, int probability, Eigen::MatrixXd v, int r, int c, int rowindex, int colindex, Eigen::VectorXd coefs, Eigen::VectorXd rho, int compprob, Eigen::VectorXd probA, Eigen::VectorXd probB, int nclasses, int totnSV, Eigen::VectorXi labels, Eigen::VectorXi nSV, int sparsemodel, int svm_type, int kernel_type, int degree, double gamma, double coef0, Eigen::MatrixXd x, int xr, Eigen::VectorXi xrowindex, Eigen::VectorXi xcolindex, int sparsex, Eigen::VectorXd ret, Eigen::VectorXd dec, Eigen::VectorXd prob);
RcppExport SEXP _kernelTDA_svmpredict_R(SEXP decisionvaluesSEXP, SEXP probabilitySEXP, SEXP vSEXP, SEXP rSEXP, SEXP cSEXP, SEXP rowindexSEXP, SEXP colindexSEXP, SEXP coefsSEXP, SEXP rhoSEXP, SEXP compprobSEXP, SEXP probASEXP, SEXP probBSEXP, SEXP nclassesSEXP, SEXP totnSVSEXP, SEXP labelsSEXP, SEXP nSVSEXP, SEXP sparsemodelSEXP, SEXP svm_typeSEXP, SEXP kernel_typeSEXP, SEXP degreeSEXP, SEXP gammaSEXP, SEXP coef0SEXP, SEXP xSEXP, SEXP xrSEXP, SEXP xrowindexSEXP, SEXP xcolindexSEXP, SEXP sparsexSEXP, SEXP retSEXP, SEXP decSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type decisionvalues(decisionvaluesSEXP);
    Rcpp::traits::input_parameter< int >::type probability(probabilitySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type rowindex(rowindexSEXP);
    Rcpp::traits::input_parameter< int >::type colindex(colindexSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type compprob(compprobSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type probA(probASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type probB(probBSEXP);
    Rcpp::traits::input_parameter< int >::type nclasses(nclassesSEXP);
    Rcpp::traits::input_parameter< int >::type totnSV(totnSVSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type nSV(nSVSEXP);
    Rcpp::traits::input_parameter< int >::type sparsemodel(sparsemodelSEXP);
    Rcpp::traits::input_parameter< int >::type svm_type(svm_typeSEXP);
    Rcpp::traits::input_parameter< int >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< int >::type degree(degreeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type coef0(coef0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type xrowindex(xrowindexSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXi >::type xcolindex(xcolindexSEXP);
    Rcpp::traits::input_parameter< int >::type sparsex(sparsexSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ret(retSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type dec(decSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(svmpredict_R(decisionvalues, probability, v, r, c, rowindex, colindex, coefs, rho, compprob, probA, probB, nclasses, totnSV, labels, nSV, sparsemodel, svm_type, kernel_type, degree, gamma, coef0, x, xr, xrowindex, xcolindex, sparsex, ret, dec, prob));
    return rcpp_result_gen;
END_RCPP
}
// wasserstein_distance
double wasserstein_distance(const Rcpp::NumericMatrix& Diag1, const Rcpp::NumericMatrix& Diag2, int q, double internal_p, double delta, double initial_eps, double eps_factor);
RcppExport SEXP _kernelTDA_wasserstein_distance(SEXP Diag1SEXP, SEXP Diag2SEXP, SEXP qSEXP, SEXP internal_pSEXP, SEXP deltaSEXP, SEXP initial_epsSEXP, SEXP eps_factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Diag1(Diag1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Diag2(Diag2SEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type internal_p(internal_pSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type initial_eps(initial_epsSEXP);
    Rcpp::traits::input_parameter< double >::type eps_factor(eps_factorSEXP);
    rcpp_result_gen = Rcpp::wrap(wasserstein_distance(Diag1, Diag2, q, internal_p, delta, initial_eps, eps_factor));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kernelTDA_svmtrain_R", (DL_FUNC) &_kernelTDA_svmtrain_R, 36},
    {"_kernelTDA_svmpredict_R", (DL_FUNC) &_kernelTDA_svmpredict_R, 30},
    {"_kernelTDA_wasserstein_distance", (DL_FUNC) &_kernelTDA_wasserstein_distance, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_kernelTDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
