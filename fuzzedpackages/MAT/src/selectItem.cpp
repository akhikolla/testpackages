#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
arma::mat selectItem_cpp(arma::mat ipar, 
                    arma::uvec available,
                    arma::uvec given, 
                    arma::colvec th,  
                    const int p,
                    arma::mat sigma, 
                    const double D = 1.7,
                    String method = "D",
                    String selectionType = "FISHER",
                    arma::colvec c_weights = 0,  
                    bool content_balancing=false,
                    int topN = 1) {
    
// input for selectItem

int i;
int ni = ipar.n_rows;
arma::mat info_i = arma::zeros<arma::mat>(p,p);
arma::vec available_num = arma::conv_to<arma::vec>::from(available);
//arma::mat available_num_matrix = arma::mat(ni,1);
//available_num_matrix = available_num ;
arma::vec given_num = arma::conv_to<arma::vec>::from(given);
arma::mat items_given; 
arma::mat info_given = arma::zeros<arma::mat>(p,p);
arma::mat info_matrix = arma::zeros<arma::mat>(p,p);
arma::mat sigma_inv = inv(sigma);
arma::mat info = arma::zeros<arma::mat>(ni,1);


if (sum(given_num)>0) {
  
  items_given = ipar.rows(find(given_num==1));
  
  //input for makeFI
  
  int ig = items_given.n_rows;
  arma::mat FI = arma::zeros<arma::mat>(p,p);
  arma::mat FI_temp = arma::zeros<arma::mat>(p,p);
  arma::mat a = items_given.cols(0,p-1);
  arma::colvec d = items_given.col(p);
  arma::colvec c = items_given.col(p+1);
  arma::mat P;
  arma::mat cf;
  double num;
  bool addsigma = false;
  
  //makeFI from previous items given
  
    for (i=0; i<ig; i++) {
    
    P = c.row(i) + (1-c.row(i))/(1+exp(-D*(a.row(i)*th + d.row(i))));
    cf = (1-P)*pow(P-c.row(i),2.0)/(P*pow(1-c.row(i),2.0));
    FI_temp = trans(a.row(i))*a.row(i);
    num = arma::as_scalar(pow(D,2.0)*cf);
    FI_temp = num*FI_temp;
    FI = FI + FI_temp;
    
    }
    
    if (addsigma == true)
    {
    FI = FI + inv(sigma);
    }
    
    info_given = FI;
                      }
 
//input for calcFI
  
  arma::mat FI_calc = arma::zeros<arma::mat>(p,p);
  arma::mat a_calc = ipar.cols(0,p-1);
  arma::colvec d_calc = ipar.col(p);
  arma::colvec c_calc = ipar.col(p+1);
  arma::mat P_calc;
  arma::mat cf_calc;
  double num_calc;
  
//calculate information for current item ;

for (i=0; i<ni; i++) {
  if (available_num(i) == 1 ) {
    if (method == "R") { 
      info.row(i) = arma::as_scalar(arma::randu<arma::vec>(1));
    }
    else {
        
      info_i = arma::zeros<arma::mat>(p,p);
      info_matrix = arma::zeros<arma::mat>(p,p);
      
      P_calc = c_calc.row(i) + (1-c_calc.row(i))/(1+exp(-D*(a_calc.row(i)*th + d_calc.row(i))));
      cf_calc = (1-P_calc)*pow(P_calc-c_calc.row(i),2.0)/(P_calc*pow(1-c_calc.row(i),2.0));
      FI_calc = trans(a_calc.row(i))*a_calc.row(i);
      num_calc = arma::as_scalar(pow(D,2.0)*cf_calc);
      FI_calc = num_calc*FI_calc;
    
      info_i = FI_calc;
      info_matrix = info_given + info_i;
      if (selectionType == "BAYESIAN") {
        info_matrix = info_matrix + sigma_inv;
      }
      if (method=="D") { 
        info.row(i) = arma::as_scalar(det(info_matrix));
      }
      else if (method=="A") {
        arma::mat a_opt_temp = pinv(info_matrix); 
        info.row(i) = arma::as_scalar(1/sum(a_opt_temp.diag()));
      }
      else if (method=="C") {
        info.row(i) = arma::as_scalar(-trans(c_weights)*pinv(info_matrix)*c_weights);
      }
    }
  }    
}
  
    return (info);
//  return List::create(Named("info")=info);
 
}

    
    
 

