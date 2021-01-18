#include "distribution.h"
#include "cumulativeR.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;
#include <algorithm>    // std::sort

// [[Rcpp::depends(RcppEigen)]]

CumulativeR::CumulativeR(void) {
  distribution dist;
}

Eigen::VectorXd CumulativeR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_logit( eta(0) );
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = Logistic::cdf_logit( eta(j) ) -  Logistic::cdf_logit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_normal( eta(0) );
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_normal( eta(j) ) - cdf_normal( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_cauchit( eta(0) );
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_cauchit( eta(j) ) - cdf_cauchit( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_gompertz( eta(0) );
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_gompertz( eta(j) ) - cdf_gompertz( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_student( eta(0) , freedom_degrees);
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_student( eta(j) , freedom_degrees ) - cdf_student( eta(j-1) , freedom_degrees ); }
  return in_open_corner(ordered_pi);
}

Eigen::VectorXd CumulativeR::inverse_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd ordered_pi( eta.size() );
  ordered_pi[0] = cdf_gumbel( eta(0) );
  for(int j=1; j<eta.size(); ++j)
  { ordered_pi[j] = cdf_gumbel( eta(j) ) - cdf_gumbel( eta(j-1) ); }
  return in_open_corner(ordered_pi);
}


Eigen::MatrixXd CumulativeR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  {F_1(j,j) =  pdf_logit(eta(j));}
  return (F_1 * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  { F_1(j,j) = pdf_normal( eta(j) ); }
  return (F_1 * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  { F_1(j,j) = pdf_cauchit( eta(j) ); }
  return (F_1 * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_student(const Eigen::VectorXd& eta,const double& freedom_degrees) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  { F_1(j,j) = pdf_student( eta(j) , freedom_degrees); }
  return (F_1 * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  { F_1(j,j) = pdf_gompertz( eta(j) ); }
  return (F_1 * R);
}

Eigen::MatrixXd CumulativeR::inverse_derivative_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(eta.rows(), eta.rows());
  R.block(0, 1, eta.rows()-1, eta.rows()-1) -= Eigen::MatrixXd::Identity(eta.rows() -1, eta.rows()-1);
  Eigen::MatrixXd F_1 = Eigen::MatrixXd::Zero(eta.rows(),eta.rows());
  for(int j=0; j<eta.rows(); ++j)
  { F_1(j,j) = pdf_gumbel( eta(j) ); }
  return (F_1 * R);
}

// // [[Rcpp::export("GLMcum")]]
// List GLMcum(Formula formula,
//             CharacterVector categories_order,
//             CharacterVector proportional,
//             DataFrame data,
//             std::string distribution,
//             double freedom_degrees,
//             Eigen::VectorXd beta_init,
//             std::string threshold){
//
//   if (!(threshold == "standard" || threshold == "equidistant" )){
//     Rcpp::stop("Unrecognized threshold restriction; options are: standard and equidistant");
//   }
//
//   class distribution dist_cum;
//
//   std::string ratio = "cumulative";
//
//   const int N = data.nrows() ; // Number of observations
//   List Full_M = dist_cum.All_pre_data_or(formula, data,
//                                          categories_order,
//                                          proportional,
//                                          threshold, ratio);
//
//   Eigen::MatrixXd Y_init = Full_M["Response_EXT"];
//   Eigen::MatrixXd X_EXT = Full_M["Design_Matrix"];
//   CharacterVector levs1 = Full_M["Levels"];
//   categories_order = Full_M["categories_order"];
//   CharacterVector explanatory_complete = Full_M["Complete_effects"];
//   int N_cats = Full_M["N_cats"];
//
//   int P_c = explanatory_complete.length();
//   int P_p = 0;
//   if(proportional[0] != "NA"){P_p = proportional.length();}
//   // if(threshold == "equidistant"){P_p = P_p + 2;}
//
//   int Q = Y_init.cols();
//   int K = Q + 1;
//   // Rcout << "X_EXT" << std::endl;
//   // Rcout << X_EXT << std::endl;
//
//   // // // Beta initialization with zeros
//   Eigen::MatrixXd BETA2;
//   BETA2 = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
//   Eigen::VectorXd BETA3 = BETA2;
//   Eigen::MatrixXd BETA = BETA2;
//   if(beta_init.size() >= 2 ){
//     BETA = beta_init;
//   }
//   int iteration = 0;
//   // double check_tutz = 1.0;
//   double Stop_criteria = 1.0;
//   Eigen::MatrixXd X_M_i ;
//   Eigen::VectorXd Y_M_i ;
//   Eigen::VectorXd eta ;
//   Eigen::VectorXd pi ;
//   Eigen::MatrixXd D ;
//   Eigen::MatrixXd Cov_i ;
//   Eigen::MatrixXd W_in ;
//   Eigen::MatrixXd Score_i_2 ;
//   Eigen::MatrixXd F_i_2 ;
//   Eigen::VectorXd LogLikIter;
//   LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
//   Eigen::MatrixXd cov_beta;
//   Eigen::VectorXd Std_Error;
//   double LogLik;
//   Eigen::MatrixXd pi_ma(N, K);
//   Eigen::MatrixXd F_i_final = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
//
//   // for (int iteration=1; iteration < 18; iteration++){
//   // while (check_tutz > 0.0001){
//   double epsilon = 0.0001 ;
//   while (Stop_criteria >( epsilon / N)){
//
//     Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
//     Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
//     LogLik = 0.;
//     // Loop by subject
//     for (int i=0; i < N; i++){
//       // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
//       X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//       Y_M_i = Y_init.row(i);
//       eta = X_M_i * BETA;
//
//       CumulativeR cum;
//
//       // Vector pi depends on selected distribution
//       if(distribution == "logistic"){
//         pi = cum.inverse_logistic(eta);
//         D = cum.inverse_derivative_logistic(eta);
//       }else if(distribution == "normal"){
//         pi = cum.inverse_normal(eta);
//         D = cum.inverse_derivative_normal(eta);
//       }else if(distribution == "cauchit"){
//         pi = cum.inverse_cauchit(eta);
//         D = cum.inverse_derivative_cauchit(eta);
//       }else if(distribution == "gompertz"){
//         pi = cum.inverse_gompertz(eta);
//         D = cum.inverse_derivative_gompertz(eta);
//       }else if(distribution == "student"){
//         pi = cum.inverse_student(eta,freedom_degrees);
//         D = cum.inverse_derivative_student(eta,freedom_degrees);
//       }else if(distribution == "gumbel"){
//         pi = cum.inverse_gumbel(eta);
//         D = cum.inverse_derivative_gumbel(eta);
//       }else{
//         Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
//       }
//
//       Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
//       W_in = D * Cov_i.inverse();
//       Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
//       Score_i = Score_i + Score_i_2;
//       F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
//       F_i = F_i + F_i_2;
//       LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );
//
//       pi_ma.row(i) = pi.transpose();
//
//     }
//
//     Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(pi_ma.rows());
//     pi_ma.col(Q) = Ones1 - pi_ma.rowwise().sum() ;
//
//     LogLikIter.conservativeResize(iteration+2, 1);
//     LogLikIter(iteration+1) = LogLik;
//     Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
//     Eigen::VectorXd beta_old = BETA;
//
//     // MatrixXd inverse;
//     FullPivLU<MatrixXd> lu(F_i);
//     bool invertible = lu.isInvertible();
//
//     if(!invertible) {
//       Rcpp::stop("Fisher matrix is not invertible");
//     }
//
//
//     BETA = BETA + (F_i.inverse() * Score_i);
//     // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
//     iteration = iteration + 1;
//     F_i_final = F_i;
//
//   }
//
//   cov_beta = F_i_final.inverse();
//   Std_Error = cov_beta.diagonal();
//   Std_Error = Std_Error.array().sqrt() ;
//
//   // Rcout << "cum1" << std::endl;
//
//   std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
//   std::vector<std::string> level_text=as<std::vector<std::string>>(categories_order);
//
//   StringVector names;
//
//   if(threshold == "equidistant"){
//     // if(P_c>1){ // hay alguna complete
//     StringVector names1(2*P_c + P_p);
//     int ind_name = 0;
//     for(int var = 0 ; var < explanatory_complete.size() ; var++){
//       names1[ind_name] = dist_cum.concatenate(text[var], level_text[0]);
//       names1[ind_name+1] = dist_cum.concatenate(text[var], "distance");
//       ind_name = ind_name + 2;
//     }
//     if(P_p>0){ // hay alguna proporcional
//       for(int var = 0 ; var < P_p ; var++){
//         names1[ind_name] = proportional[var];
//         ind_name = ind_name + 1;
//       }
//     }
//     names = names1;
//   }else{
//     StringVector names1(Q*P_c + P_p);
//     if(P_c > 0){
//       for(int var = 0 ; var < explanatory_complete.size() ; var++){
//         for(int cat = 0 ; cat < Q ; cat++){
//           names1[(Q*var) + cat] = dist_cum.concatenate(text[var], level_text[cat]);
//         }
//       }
//     }
//     if(P_p > 0){
//       for(int var_p = 0 ; var_p < proportional.size() ; var_p++){
//         names1[(Q*P_c) + var_p] = proportional[var_p];
//       }
//     }
//     names = names1;
//   }
//   // if (threshold == "equidistant"){
//   //   names[Q*P_c + P_p-2] = dist_cum.concatenate("(Intercept)", level_text[0]);
//   //   names[Q*P_c + P_p-1] = "(Intercept) distance";
//   // }
//
//   // Rcout << names << std::endl;
//
//   // // TO NAMED THE RESULT BETAS
//   NumericMatrix coef = wrap(BETA);
//   rownames(coef) = names; // this is the problem
//
//   // Rcout << coef << std::endl;
//   // AIC
//   // double AIC = (-2*LogLik) + (2 *coef.length());
//
//   // AIC
//   // double BIC = (-2*LogLik) + (coef.length() * log(N) );
//
//   int df = (N*Q) - coef.length();
//
//   Eigen::MatrixXd predict_glmcated = X_EXT * BETA;
//
//   Eigen::VectorXd Ones2 = Eigen::VectorXd::Ones(Y_init.rows());
//   Eigen::VectorXd vex1 = (Y_init.rowwise().sum()) ;
//   Y_init.conservativeResize( Y_init.rows(), K);
//   Y_init.col(Q) = (vex1 - Ones2).array().abs() ;
//
//   Eigen::MatrixXd residuals = Y_init - pi_ma;
//
//   Eigen::VectorXd pi_ma_vec(Eigen::Map<Eigen::VectorXd>(pi_ma.data(), pi_ma.cols()*pi_ma.rows()));
//   Eigen::VectorXd Y_init_vec(Eigen::Map<Eigen::VectorXd>(Y_init.data(), Y_init.cols()*Y_init.rows()));
//
//   Eigen::VectorXd div_arr = Y_init_vec.array() / pi_ma_vec.array();
//
//   Eigen::VectorXd dev_r(Y_init.rows());
//
//   int el_1 = 0;
//
//   for (int element = 0 ; element < div_arr.size() ;  element++){
//     if (div_arr[element] != 0){
//       dev_r[el_1] = div_arr[element];
//       el_1 = el_1 +1 ;
//     }
//   }
//
//   Eigen::ArrayXd dev_log = dev_r.array().log();
//
//   double deviance = dev_log.sum();
//   deviance = -2*deviance;
//
//   List output_list = List::create(
//     Named("coefficients") = coef,
//     Named("stderr") = Std_Error,
//     Named("iteration") = iteration,
//     // Named("ratio") = ratio,
//     // Named("AIC") = AIC,
//     // Named("pinv") = pinv,
//     Named("cov_beta") = cov_beta,
//     Rcpp::Named("df of the model") = df,
//     // Rcpp::Named("predict_glmcated") = predict_glmcated,
//     // Rcpp::Named("fitted") = pi_ma,
//     // Rcpp::Named("pi_ma_vec") = pi_ma_vec,
//     // Rcpp::Named("Y_init_vec") = Y_init_vec,
//     // Rcpp::Named("dev_log") = dev_log,
//     Rcpp::Named("deviance") = deviance,
//     // Rcpp::Named("residuals") = residuals,
//     Named("LogLikelihood") = LogLik,
//     // Named("freedom_degrees") = freedom_degrees,
//     // Named("Y_init") = Y_init,
//     // Named("LogLikIter") = LogLikIter,
//     Named("formula") = formula,
//     Named("categories_order") = categories_order,
//     Named("proportional") = proportional,
//     Named("N_cats") = N_cats,
//     Named("nobs_glmcat") = N,
//     Named("distribution") = distribution,
//     Named("freedom_degrees") = freedom_degrees
//   );
//
//   output_list.attr("class") = "glmcat";
//
//   return output_list;
// }
// RCPP_MODULE(cumulativemodule){
//   Rcpp::function("GLMcum", &GLMcum,
//                  List::create(_["formula"],
//                               _["categories_order"] = CharacterVector::create(NA_STRING),
//                               _["proportional"] = CharacterVector::create(NA_STRING),
//                               _["data"],
//                                _["distribution"] = "logistic",
//                                _["freedom_degrees"] = 1.0,
//                                _["beta_init"] = NumericVector::create(1),
//                                _["threshold"] = "standard"
//                  ),
//                  "Family of cumulative models");
//   Rcpp::class_<CumulativeR>("CumulativeR")
//     .constructor()
//   // .method( "GLMcum", &CumulativeR::GLMcum )
//   ;
// }







