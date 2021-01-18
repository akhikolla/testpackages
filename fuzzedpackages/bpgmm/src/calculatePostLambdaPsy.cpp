// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include "calculatePostLambdaPsy.h"
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List Calculate_PostLambdaPsy(int m,
                                   int p,
                                   Rcpp::S4 hparam,
                                   Rcpp::List CxyList,
                                   Rcpp::S4 thetaYList,
                                   arma::vec qVec,
                                   arma::vec constraint) {

  // double alpha1 = hparam.slot("alpha1");
  double alpha2 = hparam.slot("alpha2");
  double bbeta  = hparam.slot("bbeta");
  double delta  = hparam.slot("delta");

  List Cxxk      = CxyList["Cxxk"];
  List Cxyk      = CxyList["Cxyk"];
  List Cyyk      = CxyList["Cyyk"];
  List Cytytk    = CxyList["Cytytk"];
  List Cxtytk    = CxyList["Cxtytk"];
  List CxL1k     = CxyList["CxL1k"];
  List Cxmyk     = CxyList["Cxmyk"];

  arma::mat sumCxmyk  = CxyList["sumCxmyk"];
  arma::mat sumCyyk   = CxyList["sumCyyk"];

  List A         = CxyList["A"];
  arma::vec nVec = CxyList["nVec"];

  List M      = thetaYList.slot("M");
  List psy    = thetaYList.slot("psy");

  List lambda(m);

  // std::cout << constraint << std::endl;
  // std::cout << constraint[0] << std::endl;
  // std::cout << constraint[1] << std::endl;
  // std::cout << constraint[2] << std::endl;


  // Obtaining namespace of rmvnorm function
  Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
  // Picking up rmvnorm() function from rmvnorm package
  Rcpp::Function rmvnorm = mvtnorm["rmvnorm"];


  Rcpp::Environment base = Rcpp::Environment::namespace_env("base");
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function kronecker = base["kronecker"];
  Rcpp::Function c = base["c"];
  Rcpp::Function matrix = base["matrix"];
  Rcpp::Function t = base["t"];
  Rcpp::Function diag = base["diag"];
  // Rcpp::Function rgamma = stats["rgamma"];


  // Rcpp::NumericVector test_norm = c(sumCxmyk * sumCyyk.i());
  // test_norm = c(sumCxmyk * sumCyyk.i());
  // std::cout  << test_norm << std::endl;


  // Rcpp::NumericVector test_norm1 = kronecker(sumCyyk.i(), psy[1]);
  // std::cout  << test_norm << std::endl;

  // std::cout << sumCyyk << sumCyyk.i() << std::endl;
  // arma::vec test_norm;
  // test_norm = rmvnorm(_["n"] = 1, _["mean"]=arma::vec(1), _["sigma"]=arma::vec(1));

  if ((constraint[0] == 1) & (constraint[1] == 1) & (constraint[2] == 1)) {
    // std::cout << "Model 1" << std::endl;
    // Model 1
    sumCxmyk = sumCxmyk.zeros();
    sumCyyk = sumCyyk.zeros();

    for (int k=0; k<m; ++k) {
      arma::mat Cxmyk_k = Cxmyk[k];
      arma::mat Cyyk_k = Cyyk[k];
      sumCxmyk = sumCxmyk + Cxmyk_k;
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);
      alpha2_eye = alpha2 / m * alpha2_eye;  
      sumCyyk = sumCyyk + Cyyk_k + alpha2_eye;
      // std::cout << sumCxmyk << std::endl;
      // std::cout << sumCyyk << std::endl;
    }


    
    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumCxmyk * sumCyyk.i());
        // std::cout << sumCxmyk << std::endl;
        // std::cout << sumCyyk.i() << std::endl;
        // std::cout << mean_vec << std::endl;

        Rcpp::NumericMatrix sigma_mat = kronecker(sumCyyk.i(), psy[k]);
        // std::cout << sigma_mat << std::endl;


        lambda[k] = rmvnorm(Named("n", 1),
                            Named("mean", mean_vec),
                            Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        // std::cout << lambdak << std::endl;
        lambda[k] = matrix(lambdak, p, qVec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];

      shapePara += 0.5 * p * (nVec[k] + qVec[k] / m + (2 * delta - 2) / (m * p) + 1);

      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      // std::cout << "bbeta" << std::endl <<  2*bbeta << std::endl;
      // std::cout << "bbeta/(m * p)" << std::endl <<  2*bbeta/(m * p) << std::endl;

      bbeta_eye = 2 * bbeta/(m * p) * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "A_ka/m" << std::endl <<  A_ka/m << std::endl;

      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka / m) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) * 0.5 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;


    }
    shapePara += 1;
    double ratePara = arma::sum(ratePara_vec);
    double scalePara = 1 / ratePara;
    // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
    // std::cout << "ratePara: " << std::endl << ratePara << std::endl;
    // std::cout << "scalePara: " << std::endl << scalePara << std::endl;
    double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    // double invpsy = rgamma(Named("n", 1), Named("shape", shapePara), Named("rate", ratePara));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    for (int k=0; k<m; ++k) {
      arma::mat post_psy_eye(p, p, arma::fill::eye);
      // std::cout << "post_psy_eye: " << std::endl << post_psy_eye << std::endl;
      // std::cout << "post_psy_eye: " << std::endl << 1/invpsy * post_psy_eye << std::endl;

      post_psy(k) = 1/invpsy * post_psy_eye;
    }
    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if ((constraint[0] == 1) & (constraint[1] == 1) & (constraint[2] == 0)) {
    // std::cout << "Model 2" << std::endl;

    // Model 2

    sumCxmyk = sumCxmyk.zeros();
    sumCyyk = sumCyyk.zeros();

    for (int k=0; k<m; ++k) {
      arma::mat Cxmyk_k = Cxmyk[k];
      arma::mat Cyyk_k = Cyyk[k];
      sumCxmyk = sumCxmyk + Cxmyk_k;
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);
      alpha2_eye = alpha2 / m * alpha2_eye;  
      sumCyyk = sumCyyk + Cyyk_k + alpha2_eye;
      // std::cout << sumCxmyk << std::endl;
      // std::cout << sumCyyk << std::endl;
    }

    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumCxmyk * sumCyyk.i());
        // std::cout << sumCxmyk << std::endl;
        // std::cout << sumCyyk.i() << std::endl;
        // std::cout << mean_vec << std::endl;

        Rcpp::NumericMatrix sigma_mat = kronecker(sumCyyk.i(), psy[k]);
        // std::cout << sigma_mat << std::endl;


        lambda[k] = rmvnorm(Named("n", 1),
                            Named("mean", mean_vec),
                            Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        // std::cout << lambdak << std::endl;
        lambda[k] = matrix(lambdak, p, qVec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "  << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: " << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara += 0.5 * (nVec[k] + qVec[k] / m + (2 * delta - 2)/m + 1);
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = (2 * bbeta / m) * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka/m) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) * 0.5 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;


    }
    shapePara += 1;
    arma::vec ratePara = ratePara_vec;
    arma::vec scalePara = 1 / ratePara;
    // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
    // std::cout << "ratePara: " << std::endl << ratePara << std::endl;
    // std::cout << "scalePara: " << std::endl << scalePara << std::endl;

    arma::vec invpsy(p);
    for (int j=0; j<p; ++j) {
      invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));
      // std::cout << "invpsy: " << std::endl << invpsy[j] << std::endl;

    }
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;


    for (int k=0; k<m; ++k) {
      // arma::mat post_psy_eye(p, p, arma::fill::eye);
      // std::cout << "post_psy_eye: " << std::endl << post_psy_eye << std::endl;
      // std::cout << "post_psy_eye: " << std::endl << 1/invpsy * post_psy_eye << std::endl;

      // post_psy(k) = trans(invpsy) * post_psy_eye;
      post_psy(k) = diagmat(1/invpsy);
      // std::cout << "post_psy(k): " << std::endl << trans(invpsy) * post_psy_eye << std::endl;
      // std::cout << "1/invpsy: " << std::endl << diagmat(1/invpsy) << std::endl;

    }
    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);
  } else if ((constraint[0] == 1) & (constraint[1] == 0) & (constraint[2] == 1)) {
    // std::cout << "Model 3" << std::endl;
    arma::mat Cxmyk_0 = Cxmyk[0];
    arma::mat Cyyk_0 = Cyyk[0];

    // std::cout << "Cxmyk_0: " << std::endl << Cxmyk_0.n_rows  << std::endl << Cxmyk_0.n_cols << std::endl;
    arma::mat sumPhiCxy(Cxmyk_0.n_rows, Cxmyk_0.n_cols, fill::zeros);
    arma::mat sumPhiCyy(Cyyk_0.n_rows, Cyyk_0.n_cols, fill::zeros);

    // std::cout << "sumPhiCxy: " << std::endl << sumPhiCxy << std::endl;
    // std::cout << "sumPhiCyy: " << std::endl << sumPhiCyy << std::endl;

    for (int k=0; k<m; ++k) {
      arma::mat psy_k = psy[k];
      // std::cout << "psy[k]: " << std::endl << psy_k(1, 1)     << std::endl;
      // std::cout << "1 / psy_k(1, 1): " << std::endl << 1 / psy_k(1, 1) << std::endl;
      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);


      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);
      alpha2_eye = (alpha2 / m) * alpha2_eye;
      // std::cout << "alpha2: " << std::endl << alpha2 << std::endl;
      // std::cout << "alpha2_eye: " << std::endl << alpha2_eye << std::endl;

      // std::cout << "Cxmyk_ka: " << std::endl << Cxmyk_ka << std::endl;
      // std::cout << "1 / psy_k(1, 1): " << std::endl << 1 / psy_k(1, 1) << std::endl;
      // std::cout << "1 / psy_k(1, 1) * Cxmyk_ka: " << std::endl << Cxmyk_ka / psy_k(1, 1) << std::endl;

      sumPhiCxy += Cxmyk_ka / psy_k(1, 1);
      sumPhiCyy += ( Cyyk_ka + alpha2_eye ) / psy_k(1, 1);

      // std::cout << "sumPhiCxy: " << std::endl << sumPhiCxy << std::endl;
      // std::cout << "sumPhiCyy: " << std::endl << sumPhiCyy << std::endl;
    }

    for (int k=0; k<m; ++k) {
      if (k == 0) {

        Rcpp::NumericVector mean_vec = c(sumPhiCxy * sumPhiCyy.i());
        // std::cout << sumPhiCxy << std::endl;
        // std::cout << sumPhiCyy.i() << std::endl;
        // std::cout << mean_vec << std::endl;

        arma::mat p_eye(p, p, arma::fill::eye);
        Rcpp::NumericMatrix sigma_mat = kronecker(sumPhiCyy.i(), p_eye);
        // std::cout << sigma_mat << std::endl;


        lambda[k] = rmvnorm(Named("n", 1),
                            Named("mean", mean_vec),
                            Named("sigma", sigma_mat));
        Rcpp::NumericMatrix lambdak = lambda[k];
        // std::cout << lambdak << std::endl;
        lambda[k] = matrix(lambdak, p, qVec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara = 0.5 * p * (nVec[k] + qVec[k] / m + (2 * delta - 2) / p + 1) + 1;
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / p) * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      ratePara_vec = arma::diagvec(ratePara_k) * 0.5;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) * 0.5 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;
      double ratePara = arma::sum(ratePara_vec);
      double scalePara = 1 / ratePara;
      double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));
      arma::mat post_psy_eye(p, p, arma::fill::eye);
      post_psy(k) = 1/invpsy * post_psy_eye;

    }
    // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
    // std::cout << "ratePara: " << std::endl << ratePara << std::endl;
    // std::cout << "scalePara: " << std::endl << scalePara << std::endl;
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    // double invpsy = rgamma(Named("n", 1), Named("shape", shapePara), Named("rate", ratePara));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;


    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if ((constraint[0] == 1) & (constraint[1] == 0) & (constraint[2] == 0)) {
    // std::cout << "Model 4" << std::endl;
    arma::mat sumVar;
    arma::mat B;

    for (int k=0; k<m; ++k) {
      arma::mat psy_k = psy[k];
      // double qVec_k = qVec[k];
      // std::cout << "psy_k: " << std::endl << psy_k << std::endl;
      // std::cout << "qVec: " << std::endl << qVec << std::endl;

      // std::cout << "qVec_k: " << std::endl << qVec_k << std::endl;
      // std::cout << "psy_k.i(): " << std::endl << psy_k.i() << std::endl;
      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);

      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);
      alpha2_eye = (alpha2 / m) * alpha2_eye;
      // std::cout << "alpha2_eye: " << std::endl << alpha2_eye << std::endl;
      // std::cout << "Cyyk_ka: " << std::endl << Cyyk_ka << std::endl;

      // std::cout << "Cyyk_k + alpha2_eye: " << std::endl << Cyyk_ka + alpha2_eye << std::endl;
      // std::cout << "psy_k.i(): " << std::endl << psy_k.i() << std::endl;

      Rcpp::NumericMatrix sumVar_plus = kronecker(Cyyk_ka + alpha2_eye, psy_k.i());
      arma::mat sumVar_plusa = Rcpp::as<arma::mat>(sumVar_plus);
      // std::cout << "sumVar_plusa: " << std::endl << sumVar_plusa << std::endl;
      if (k == 0) {
        sumVar = sumVar_plusa;
        B = psy_k.i() * Cxmyk_ka;
      } else {
        sumVar += sumVar_plusa;
        B += psy_k.i() * Cxmyk_ka;

      }
      // std::cout << "sumVar: " << std::endl << sumVar << std::endl;
      // std::cout << "B: " << std::endl << B << std::endl;

    }

    arma::mat lambdaVar = sumVar.i();
    Rcpp::NumericVector c_B = c(B);
    arma::vec c_Ba = Rcpp::as<arma::vec>(c_B);

    // std::cout << "B: " << std::endl << B << std::endl;
    // std::cout << "c_Ba: " << std::endl << c_Ba << std::endl;
    // std::cout << "lambdaVar: " << std::endl << lambdaVar.n_rows << lambdaVar.n_cols << std::endl;

    arma::mat lambdaMean = trans(c_Ba) * lambdaVar;
    // std::cout << "lambdaMean: " << std::endl << lambdaMean << std::endl;
    // std::cout << "lambdaMean: " << std::endl << lambdaMean.n_rows << lambdaMean.n_cols << std::endl;

    for (int k=0; k<m; ++k) {
      if (k == 0) {
        lambda[k] = rmvnorm(Named("n", 1),
                            Named("mean", lambdaMean),
                            Named("sigma", lambdaVar));
        Rcpp::NumericMatrix lambdak = lambda[k];
        // std::cout << lambdak << std::endl;
        lambda[k] = matrix(lambdak, p, qVec[k]);

      } else{
        lambda[k] = lambda[0];
      }
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara = 0.5 * (nVec[k] + qVec[k] / m + 2 * delta - 1) + 1;
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * bbeta * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka / m) * trans(tildaLambda_ka) + bbeta_eye;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      arma::vec scalePara = 1 / ratePara_k;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;
      // std::cout << "scalePara: " << std::endl << scalePara << std::endl;

      arma::vec invpsy(p);
      for (int j=0; j<p; ++j) {

        // std::cout << "arma::vec: " << std::endl << arma::randg( 1, distr_param(shapePara, scalePara[j]) ) << std::endl;

        invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));;
        // std::cout << "invpsy[j]: " << std::endl << invpsy[j] << std::endl;

      }

      //

      // std::cout << "diagmat(1/invpsy): " << std::endl << diagmat(1/invpsy) << std::endl;


      post_psy(k) = diagmat(1/invpsy);



    }

    // for (int k=0; k<m; ++k) {
    //   arma::mat post_psyk = post_psy(k);
    //
    //   std::cout << "post_psy(k): " << std::endl << post_psyk << std::endl;
    //
    // }

    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if ((constraint[0] == 0) & (constraint[1] == 1) & (constraint[2] == 1)) {
    // std::cout << "Model 5" << std::endl;
    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);

      // std::cout << alpha2_eye << std::endl;
      alpha2_eye *= alpha2;
      // std::cout << alpha2_eye << std::endl;

      // std::cout << Cyyk_ka    << std::endl;
      Cyyk_ka += alpha2_eye;
      // std::cout << Cyyk_ka    << std::endl;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());

      // std::cout << mean_vec   << std::endl;

      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);
      // std::cout << sigma_mat << std::endl;


      lambda[k] = rmvnorm(Named("n", 1),
                          Named("mean", mean_vec),
                          Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      // std::cout << lambdak << std::endl;
      lambda[k] = matrix(lambdak, p, qVec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara += 0.5 * p * (nVec[k] + qVec[k] + (2 * delta - 2)/(m * p) + 1);
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / (m * p)) * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) * 0.5 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;


    }
    shapePara += 1;
    double ratePara = arma::sum(ratePara_vec);
    double scalePara = 1 / ratePara;
    // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
    // std::cout << "ratePara: " << std::endl << ratePara << std::endl;
    // std::cout << "scalePara: " << std::endl << scalePara << std::endl;
    double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    // double invpsy = rgamma(Named("n", 1), Named("shape", shapePara), Named("rate", ratePara));
    // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

    for (int k=0; k<m; ++k) {
      arma::mat post_psy_eye(p, p, arma::fill::eye);
      // std::cout << "post_psy_eye: " << std::endl << post_psy_eye << std::endl;
      // std::cout << "post_psy_eye: " << std::endl << 1/invpsy * post_psy_eye << std::endl;

      post_psy(k) = 1/invpsy * post_psy_eye;
    }
    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);


  } else if ((constraint[0] == 0) & (constraint[1] == 1) & (constraint[2] == 0)) {
    // std::cout << "Model 6" << std::endl;

    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);

      // std::cout << alpha2_eye << std::endl;
      alpha2_eye *= alpha2;
      // std::cout << alpha2_eye << std::endl;

      // std::cout << Cyyk_ka    << std::endl;
      Cyyk_ka += alpha2_eye;
      // std::cout << Cyyk_ka    << std::endl;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());

      // std::cout << mean_vec   << std::endl;

      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);
      // std::cout << sigma_mat << std::endl;


      lambda[k] = rmvnorm(Named("n", 1),
                          Named("mean", mean_vec),
                          Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      // std::cout << lambdak << std::endl;
      lambda[k] = matrix(lambdak, p, qVec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      // double nVec_k = nVec[k];
      // double qVec_k = qVec[k];


      // std::cout << "nVec[k]: "  << std::endl << nVec[k] << std::endl;
      // std::cout << "qVec[k]: "  << std::endl << qVec[k] << std::endl;
      // std::cout << "delta: "  << std::endl << delta << std::endl;


      shapePara += 0.5 * (nVec[k] + qVec[k] + (2 * delta - 2)/m + 1);
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;
      // ratePara_vec += 1；

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      // std::cout << "Cxxk[k]: " << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_k: " << std::endl << Cxtytk_k << std::endl;
      // std::cout << "tildaLambda_k: " << std::endl << tildaLambda_k << std::endl;

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / m) * bbeta_eye;
      // std::cout << "Cxxk[k]: (p * p)" << std::endl << Cxxk_k << std::endl;
      // std::cout << "Cxtytk_ka: (p * m)" << std::endl << Cxtytk_ka << std::endl;
      // std::cout << "tildaLambda_ka: (m * p)" << std::endl << trans(tildaLambda_ka) << std::endl;
      // std::cout << "A_ka" << std::endl <<  A_ka << std::endl;
      // std::cout << "bbeta_eye" << std::endl <<  bbeta_eye << std::endl;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      // ratePara_k = arma::diagvec(ratePara_k) * 0.5;
      ratePara_vec += arma::diagvec(ratePara_k) * 0.5;
      // std::cout << "ratePara_k: " << std::endl << ratePara_k << std::endl;

      // std::cout << "ratePara_k: " << std::endl << arma::diagvec(ratePara_k) * 0.5 << std::endl;
      // std::cout << "ratePara_vec: " << std::endl << ratePara_vec << std::endl;

      // std::cout << "test: (p * p)" << std::endl << arma::diagvec(test) << std::endl;


    }
    shapePara += 1;

    arma::vec scalePara = 1 / ratePara_vec;
    arma::vec invpsy(p);
    for (int j=0; j<p; ++j) {

      invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara[j]) ));;
      // std::cout << "invpsy[j]: " << std::endl << invpsy[j] << std::endl;

    }

    for (int k=0; k<m; ++k) {

      post_psy(k) = diagmat(1/invpsy);
      // std::cout << "diagmat(1/invpsy): " << std::endl << diagmat(1/invpsy) << std::endl;
    }

    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if ((constraint[0] == 0) & (constraint[1] == 0) & (constraint[2] == 1)) {
    // std::cout << "Model 7" << std::endl;


    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);

      // std::cout << alpha2_eye << std::endl;
      alpha2_eye *= alpha2;
      // std::cout << alpha2_eye << std::endl;

      // std::cout << Cyyk_ka    << std::endl;
      Cyyk_ka += alpha2_eye;
      // std::cout << Cyyk_ka    << std::endl;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());

      // std::cout << mean_vec   << std::endl;

      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);
      // std::cout << sigma_mat << std::endl;


      lambda[k] = rmvnorm(Named("n", 1),
                          Named("mean", mean_vec),
                          Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      // std::cout << lambdak << std::endl;
      lambda[k] = matrix(lambdak, p, qVec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      shapePara = 0.5 * p * (nVec[k] + qVec[k] +  (2 * delta - 2)/p + 1) + 1;
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * (bbeta / p) * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      ratePara_vec = arma::diagvec(ratePara_k) * 0.5;
      double ratePara = arma::sum(ratePara_vec);
      // std::cout << "ratePara: "  << std::endl << ratePara << std::endl;

      double scalePara = 1 / ratePara;
      double invpsy = sum(arma::randg( 1, distr_param(shapePara, scalePara) ));

      arma::mat post_psy_eye(p, p, arma::fill::eye);

      // std::cout << "post_psy_eye: " << std::endl << post_psy_eye << std::endl;
      // std::cout << "post_psy_eye: " << std::endl << 1/invpsy * post_psy_eye << std::endl;

      post_psy(k) = 1/invpsy * post_psy_eye;

    }

    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);

    return(res);

  } else if ((constraint[0] == 0) & (constraint[1] == 0) & (constraint[2] == 0)) {
    // std::cout << "Model 8" << std::endl;

    for (int k=0; k<m; ++k) {


      Rcpp::NumericMatrix Cxmyk_k = Cxmyk[k];
      Rcpp::NumericMatrix Cyyk_k = Cyyk[k];

      arma::mat Cxmyk_ka = Rcpp::as<arma::mat>(Cxmyk_k);
      arma::mat Cyyk_ka = Rcpp::as<arma::mat>(Cyyk_k);
      arma::mat alpha2_eye(qVec[k], qVec[k], arma::fill::eye);

      // std::cout << alpha2_eye << std::endl;
      alpha2_eye *= alpha2;
      // std::cout << alpha2_eye << std::endl;

      // std::cout << Cyyk_ka    << std::endl;
      Cyyk_ka += alpha2_eye;
      // std::cout << Cyyk_ka    << std::endl;

      Rcpp::NumericVector mean_vec = c(Cxmyk_ka * Cyyk_ka.i());

      // std::cout << mean_vec   << std::endl;

      Rcpp::NumericMatrix sigma_mat = kronecker(Cyyk_ka.i(), psy[k]);
      // std::cout << sigma_mat << std::endl;


      lambda[k] = rmvnorm(Named("n", 1),
                          Named("mean", mean_vec),
                          Named("sigma", sigma_mat));
      Rcpp::NumericMatrix lambdak = lambda[k];
      // std::cout << lambdak << std::endl;
      lambda[k] = matrix(lambdak, p, qVec[k]);
    }

    // post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k

    List tildaLambda(m);
    for (int k=0; k<m; ++k) {
      Rcpp::NumericMatrix lambda_k = lambda[k];
      Rcpp::NumericVector m_k = M[k];
      // std::cout << "lambda_k" << lambda_k << std::endl;
      // std::cout << "m_k: " << m_k << std::endl;

      arma::vec m_ka = m_k;
      // std::cout << "m_ka: "       << std::endl << m_ka << std::endl;
      arma::mat lambda_ka = Rcpp::as<arma::mat>(lambda_k);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      lambda_ka.insert_cols(0, m_ka);
      // std::cout << "lambda_ka: "  << std::endl << lambda_ka << std::endl;
      tildaLambda[k] = lambda_ka;
    }

    // Post psy

    List post_psy(m);

    double shapePara = 0;
    arma::vec ratePara_vec(p, arma::fill::zeros);
    // std::cout << "ratePara_vec: "  << std::endl << ratePara_vec << std::endl;

    for (int k=0; k<m; ++k) {

      shapePara = 0.5 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1;
      // std::cout << "shapePara: "  << std::endl << shapePara << std::endl;

      Rcpp::NumericMatrix Cxxk_k = Cxxk[k];
      Rcpp::NumericMatrix Cxtytk_k = Cxtytk[k];
      Rcpp::NumericMatrix Cytytk_k = Cytytk[k];

      Rcpp::NumericMatrix tildaLambda_k = tildaLambda[k];
      Rcpp::NumericMatrix A_k = A[k];

      arma::mat Cxxk_ka = Rcpp::as<arma::mat>(Cxxk_k);
      arma::mat Cxtytk_ka = Rcpp::as<arma::mat>(Cxtytk_k);
      arma::mat Cytytk_ka = Rcpp::as<arma::mat>(Cytytk_k);

      arma::mat tildaLambda_ka = Rcpp::as<arma::mat>(tildaLambda_k);
      arma::mat A_ka = Rcpp::as<arma::mat>(A_k);
      arma::mat bbeta_eye(p, p, arma::fill::eye);

      bbeta_eye = 2 * bbeta * bbeta_eye;

      arma::mat ratePara_k = Cxxk_ka - 2 * Cxtytk_ka * trans(tildaLambda_ka) + tildaLambda_ka * (Cytytk_ka + A_ka) * trans(tildaLambda_ka) + bbeta_eye;
      Rcpp::NumericVector ratePara = c(arma::diagvec(ratePara_k) * 0.5);
      // std::cout << "ratePara: "  << std::endl << ratePara << std::endl;

      // arma:vec scalePara_vec = 1.0 / ratePara_vec;
      // std::cout << "scalePara_vec: "  << std::endl << scalePara_vec << std::endl;

      // arma::vec invpsy(p);
      // for (int j=0; j<p; ++j) {

        // std::cout << "arma::vec: " << std::endl << arma::randg( 1, distr_param(shapePara, scalePara[j]) ) << std::endl;

        // invpsy[j] = sum(arma::randg( 1, distr_param(shapePara, scalePara_vec[j]) ));;
        // std::cout << "invpsy[j]: " << std::endl << invpsy[j] << std::endl;

      // }

      // std::cout << "diagmat(1/invpsy): " << std::endl << diagmat(1/invpsy) << std::endl;
      // Rcpp::NumericVector invpsy = rgamma(Name("n", p),
      //                                     Name("shape", shapePara),
      //                                     Name("rate", ratePara));


      // std::cout << "p: " << std::endl << p << std::endl;
      // std::cout << "shapePara: " << std::endl << shapePara << std::endl;
      // std::cout << "ratePara: " << std::endl << ratePara << std::endl;

      arma::vec invpsy(p);
      for (int j=0; j<p; ++j) {
        // std::cout << "ratePara[j]: " << std::endl << ratePara[j] << std::endl;
        double ratePara_j = ratePara[j];
        double scalePara_j = 1.0/ratePara_j;
        invpsy[j] = sum(Rcpp::rgamma(1, shapePara, scalePara_j));

      }
      // std::cout << "invpsy: " << std::endl << invpsy << std::endl;

      post_psy(k) = diagmat(1.0/invpsy);
      // std::cout << "diagmat(1/invpsy): " << std::endl << diagmat(1/invpsy) << std::endl;

    }

    List res = Rcpp::List::create(Named("lambda") = lambda,
                                  Named("psy")    = post_psy);
    return(res);
  }

  List res = Rcpp::List::create(Named("lambda") = lambda,
                                Named("psy")    = psy);
  return(res);
}
