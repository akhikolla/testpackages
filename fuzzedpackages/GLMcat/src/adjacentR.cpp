#include "distribution.h"
#include "adjacentR.h"
using namespace std;
using namespace Rcpp ;

// [[Rcpp::depends(RcppEigen)]]

AdjacentR::AdjacentR(void) {
  distribution dist;
}

Eigen::VectorXd AdjacentR::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_logit( eta(eta.size()-1) ) / ( 1-cdf_logit( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_logit( eta(j-1) ) / ( 1-cdf_logit( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_logistic(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_logit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_logit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_logit(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_normal( eta(eta.size()-1) ) / ( 1-cdf_normal( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_normal( eta(j-1) ) / ( 1-cdf_normal( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_normal(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_normal( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_normal(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_cauchit( eta(eta.size()-1) ) / ( 1-cdf_cauchit( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_cauchit( eta(j-1) ) / ( 1-cdf_cauchit( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_cauchit(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_cauchit( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_cauchit(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_cauchit(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_gompertz( eta(eta.size()-1) ) / ( 1-cdf_gompertz( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_gompertz( eta(j-1) ) / ( 1-cdf_gompertz( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_gompertz(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_gompertz(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_gompertz( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_gompertz(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_gompertz(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_gumbel( eta(eta.size()-1) ) / ( 1-cdf_gumbel( eta(eta.size()-1) ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_gumbel( eta(j-1) ) / ( 1-cdf_gumbel( eta(j-1) ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_gumbel(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_gumbel(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_gumbel( eta(j) ) /( std::max(1e-10, std::min(1-1e-6, cdf_gumbel(eta(j)))) * std::max(1e-10, std::min(1-1e-6, 1-cdf_gumbel(eta(j)))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}

Eigen::VectorXd AdjacentR::inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const
{
  Eigen::VectorXd pi( eta.size() );
  pi[eta.size()-1] = cdf_student( eta(eta.size()-1),freedom_degrees ) / ( 1-cdf_student( eta(eta.size()-1),freedom_degrees ) );
  double norm = 1 + pi[eta.size()-1];
  for(int j=(eta.size()-1); j>0; --j)
  {
    pi[j-1] = pi[j] * cdf_student( eta(j-1) ,freedom_degrees) / ( 1-cdf_student( eta(j-1),freedom_degrees ) );
    norm += pi[j-1];
  }
  return in_open_corner(pi/norm);
}

Eigen::MatrixXd AdjacentR::inverse_derivative_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const
{
  Eigen::VectorXd pi = AdjacentR::inverse_student(eta, freedom_degrees);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(pi.rows(),pi.rows());
  for(int j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_student( eta(j) ,freedom_degrees) /( std::max(1e-10, std::min(1-1e-6, pdf_student(eta(j),freedom_degrees))) * std::max(1e-10, std::min(1-1e-6, 1-pdf_student(eta(j),freedom_degrees))) ); }

  return D * Eigen::TriangularView<Eigen::MatrixXd, Eigen::UpLoType::Lower>(Ones) * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose() );

}



// distribution dist_adj;

// // [[Rcpp::export(".GLMadj")]]
// List GLMadj(Formula formula,
//             CharacterVector categories_order,
//             CharacterVector proportional,
//             DataFrame data,
//             std::string distribution,
//             double freedom_degrees){
//
//   const int N = data.nrows() ; // Number of observations
//
//   List Full_M = dist_adj.All_pre_data_or(formula, data,
//                                       categories_order, proportional);
//
//   Eigen::MatrixXd Y_init = Full_M["Response_EXT"];
//   Eigen::MatrixXd X_EXT = Full_M["Design_Matrix"];
//   CharacterVector levs1 = Full_M["Levels"];
//   CharacterVector explanatory_complete = Full_M["Complete_effects"];
//
//   int P_c = explanatory_complete.length();
//   int P_p = 0;
//   if(proportional[0] != "NA"){P_p = proportional.length();}
//   // int P =  P_c +  P_p ; // Number of explanatory variables without intercept
//
//   int Q = Y_init.cols();
//   int K = Q + 1;
//
//   // // // Beta initialization with zeros
//   Eigen::MatrixXd BETA;
//   BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
//   //
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
//
//     // Loop by subject
//     for (int i=0; i < N; i++){
//       // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
//       X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
//       Y_M_i = Y_init.row(i);
//       eta = X_M_i * BETA;
//
//       AdjacentR adj;
//
//       // Vector pi depends on selected distribution
//       if(distribution == "logistic"){
//         pi = adj.inverse_logistic(eta);
//         D = adj.inverse_derivative_logistic(eta);
//       }else if(distribution == "normal"){
//         pi = adj.inverse_normal(eta);
//         D = adj.inverse_derivative_normal(eta);
//       }else if(distribution == "cauchit"){
//         pi = adj.inverse_cauchit(eta);
//         D = adj.inverse_derivative_cauchit(eta);
//       }else if(distribution == "gompertz"){
//         pi = adj.inverse_gompertz(eta);
//         D = adj.inverse_derivative_gompertz(eta);
//       }else if(distribution == "gumbel"){
//         pi = adj.inverse_gumbel(eta);
//         D = adj.inverse_derivative_gumbel(eta);
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
//     BETA = BETA + (F_i.inverse() * Score_i);
//     // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
//     iteration = iteration + 1;
//     F_i_final = F_i;
//   }
//
//   // cov_beta = (((X_EXT.transpose() * F_i_final) * X_EXT).inverse());
//   cov_beta = F_i_final.inverse();
//   Std_Error = cov_beta.diagonal();
//   Std_Error = Std_Error.array().sqrt() ;
//
//   std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
//   std::vector<std::string> level_text=as<std::vector<std::string>>(categories_order);
//   StringVector names(Q*P_c + P_p);
//   if(P_c > 0){
//     for(int var = 0 ; var < explanatory_complete.size() ; var++){
//       for(int cat = 0 ; cat < Q ; cat++){
//         names[(Q*var) + cat] = dist_adj.concatenate(text[var], level_text[cat]);
//       }
//     }
//   }
//   if(P_p > 0){
//     for(int var_p = 0 ; var_p < proportional.size() ; var_p++){
//       names[(Q*P_c) + var_p] = proportional[var_p];
//     }
//   }
//   // TO NAMED THE RESULT BETAS
//   NumericMatrix coef = wrap(BETA);
//   rownames(coef) = names;
//
//   // AIC
//   double AIC = (-2*LogLik) + (2 *coef.length());
//
//   // AIC
//   double BIC = (-2*LogLik) + (coef.length() * log(N) );
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
//   return List::create(
//     // Named("Nb. iterations") = iteration-1 ,
//     Named("coefficients") = coef,
//     Named("AIC") = AIC,
//     Named("BIC") = BIC,
//     // Named("cov_beta") = cov_beta,
//     Named("stderr") = Std_Error,
//     Rcpp::Named("df") = df,
//     Rcpp::Named("predict_glmcated") = predict_glmcated,
//     Rcpp::Named("fitted") = pi_ma,
//     Rcpp::Named("pi_ma_vec") = pi_ma_vec,
//     Rcpp::Named("Y_init_vec") = Y_init_vec,
//     Rcpp::Named("dev_log") = dev_log,
//     Rcpp::Named("deviance") = deviance,
//     Rcpp::Named("residuals") = residuals,
//     Named("LogLikelihood") = LogLik
//   );
// }



// RCPP_MODULE(adjacentmodule){
//   Rcpp::function("GLMadj", &GLMadj,
//                  List::create(_["formula"] = R_NaN,
//                               _["categories_order"] = CharacterVector::create( "A", NA_STRING),
//                               _["proportional"] = CharacterVector::create(NA_STRING),
//                               _["data"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf),
//                               _["distribution"] = "a",
//                               _["freedom_degrees"] = 1.0),
//                               "Adjacent model");
//   // Rcpp::class_<AdjacentR>("AdjacentR")
//   // .constructor()
//   // // .method( "GLMadj", &AdjacentR::GLMadj)
//   // ;
// }

