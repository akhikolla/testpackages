#include "distribution.h"
#include "reference.h"
#include "adjacentR.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

//' Family of models for Discrete Choice
//' @description Discrete choice model: Requires data in long form.
//' For each individual (or decision maker), there are multiple observations (rows),
//' one for each of the alternatives the individual could have chosen.
//' We call the group of observations for an individual a “case”.
//' Each case represents a single statistical observation although it comprises
//' multiple observations.
//' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model. A particularity for the formula is that for the case-specific variables, the user can define a specific effect for a category.
//' @param case_id a string with the name of the column that identifies each case.
//' @param alternatives a string with the name of the column that identifies the vector of alternatives the individual could have chosen.
//' @param reference a string indicating the reference category
//' @param alternative_specific a character vector with the name of the explanatory variables that are different for each case, these are the alternative specific variables. By default, the case specific variables are the explanatory variables that are not identify in here, but that are part of the formula.
//' @param data a dataframe object in R, with the dependent variable as factor.
//' @param distribution a string indicating the F distribution, options are: logistic, normal, cauchit, student (any df), gompertz, gumbel.
//' @param freedom_degrees an optional scalar to indicate the degrees of freedom for the Student distribution.
//' @examples
//' library(GLMcat)
//' data(TravelChoice)
//' Discrete_CM(formula = choice ~ hinc + gc + invt,
//' case_id = "indv",alternatives = "mode",reference = "air",
//' data = TravelChoice,  alternative_specific = c("gc", "invt"),
//' distribution = "logistic")
//' @note For these models it is not allowed to exclude the intercept.
//' @export
// [[Rcpp::export("Discrete_CM")]]
List Discrete_CM(Formula formula,
                 String case_id,
                 String alternatives,
                 CharacterVector reference,
                 CharacterVector alternative_specific,
                 DataFrame data,
                 std::string distribution,
                 double freedom_degrees
){

  class distribution dist1;

  List Full_M = dist1.select_data_nested(formula,
                                         case_id,
                                         alternatives,
                                         reference,
                                         alternative_specific,
                                         data
                                           //   ,
                                           // ratio
  );

  Eigen::MatrixXd Y_init = Full_M["Response_M"];
  Eigen::MatrixXd X_EXT = Full_M["Design_Matrix"];

  CharacterVector Names_design = Full_M["Names_design"];

  int Q = Y_init.cols();
  int K = Q + 1;
  int N = K * Y_init.rows();

  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1); // Beta initialization with zeros
  int iteration = 0;
  double Stop_criteria = 1.0;
  Eigen::MatrixXd X_M_i ;
  Eigen::VectorXd Y_M_i ;
  Eigen::VectorXd eta ;
  Eigen::VectorXd pi ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd Cov_i ;
  Eigen::MatrixXd W_in ;
  Eigen::MatrixXd Score_i_2 ;
  Eigen::MatrixXd F_i_2 ;
  Eigen::VectorXd LogLikIter;
  LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
  double LogLik;

  Eigen::MatrixXd F_i_final = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
  Eigen::MatrixXd cov_beta;
  Eigen::VectorXd Std_Error;

  double epsilon = 0.0001 ;
  // for (int iteration=1; iteration < 18; iteration++){
  while (Stop_criteria >( epsilon / N)){
    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;


    for (int i=0; i < N/K; i++){
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      // if(ratio == "reference"){
      ReferenceF ref;
      if(distribution == "logistic"){
        pi = ref.inverse_logistic(eta);
        D = ref.inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = ref.inverse_normal(eta);
        D = ref.inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = ref.inverse_cauchit(eta);
        D = ref.inverse_derivative_cauchit(eta);
      }else if(distribution == "student"){
        pi = ref.inverse_student(eta, freedom_degrees);
        D = ref.inverse_derivative_student(eta, freedom_degrees);
      }
      // }else{
      //   AdjacentR adj;
      //   if(distribution == "logistic"){
      //     pi = adj.inverse_logistic(eta);
      //     D = adj.inverse_derivative_logistic(eta);
      //   }else if(distribution == "normal"){
      //     pi = adj.inverse_normal(eta);
      //     D = adj.inverse_derivative_normal(eta);
      //   }else if(distribution == "cauchit"){
      //     pi = adj.inverse_cauchit(eta);
      //     D = adj.inverse_derivative_cauchit(eta);
      //   }else if(distribution == "student"){
      //     pi = adj.inverse_student(eta, freedom_degrees);
      //     D = adj.inverse_derivative_student(eta, freedom_degrees);
      //   }
      // }

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      W_in = D * Cov_i.inverse();
      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );
    }
    // To stop when LogLik is smaller than the previous
    if(iteration>1){
      if (LogLikIter[iteration] > LogLik)
        break;
    }

    LogLikIter.conservativeResize( LogLikIter.rows() +1 , 1);
    LogLikIter(LogLikIter.rows() - 1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);

    iteration = iteration + 1;
    // Rcout << "BETA" << std::endl;
    // Rcout << BETA << std::endl;

    // Rcout << "iteration" << std::endl;
    // Rcout << iteration << std::endl;
    // Rcout << "LogLik" << std::endl;
    // Rcout << LogLik << std::endl;
    F_i_final = F_i;

    // Rcout << "LogLikIter" << std::endl;
    // Rcout << LogLikIter << std::endl;

  }

  cov_beta = F_i_final.inverse();
  Std_Error = cov_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  // Eigen::MatrixXd X_M_i_1 = X_EXT.block(0*Q , 0 , Q , X_EXT.cols());
  // Eigen::VectorXd Y_M_i_1 = Y_init.row(0);

  NumericMatrix BETA_2 = wrap(BETA);
  rownames(BETA_2) = Names_design;

  List output_list_dis = List::create(
    Named("Nb. iterations") = iteration-1 ,
    Named("coefficients") = BETA_2,
    Named("LogLikelihood") = LogLikIter(LogLikIter.rows() - 1),
    Named("LogLikIter") =  LogLikIter,
    Named("X_M_i") =  X_M_i,
    Named("stderr") =  Std_Error
  );

  output_list_dis.attr("class") = "glmcat";
  return output_list_dis;
}

RCPP_MODULE(discretemodule){
  Rcpp::function("Discrete_CM", &Discrete_CM,
                 List::create(_["formula"] = R_NaN,
                              _["case_id"] = "a",
                              _["alternatives"] = "a",
                              _["reference"] = R_NaN,
                              _["alternative_specific"] = CharacterVector::create( NA_STRING),
                              _["data"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf),
                              _["distribution"] = "a",
                              _["freedom_degrees"] = 1.0),
                              "Discrete Choice Model");

}
