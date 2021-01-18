#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#include <iostream>
#include <iomanip>
using namespace std;

// function prototypes
arma::vec grad_loss_alpha_beta(arma::vec & alpha,
				  arma::vec & beta,
				  double delta,
				  arma::mat & SampCov,
				  arma::mat & S,
				  arma::mat & Sinv,
				  int index_group);

double grad_loss_delta(arma::vec & alpha,
                       arma::vec & beta,
                       double delta,
                       arma::mat & SampCov,
                       arma::mat & S,
                       arma::mat & Sinv);


double grad_loss_vx(arma::vec & alpha,
                        arma::vec & beta,
                        double delta,
                        arma::mat & SampCov,
                        arma::mat & Sinv);

double grad_loss_vy(arma::vec & alpha,
                        arma::vec & beta,
                        double delta,
                        arma::mat & SampCov,
                        arma::mat & Sinv);


arma::mat compute_B(arma::vec & alpha, arma::vec & beta, double delta);

arma::mat inverse_ImpCov(arma::vec & alpha, arma::vec & beta, double delta, arma::mat & Sinv);

void update_ImpCov(arma::mat & ImpCov, arma::vec & alpha, arma::vec & beta, double delta);

double loss_unpenalized(arma::mat & ImpCov, arma::mat & SampCov, 
			arma::vec & alpha, arma::vec & beta, 
			double delta, arma::mat & Sinv, double logdetMedCov);

double soft_threshold(double z, double lambda); 

arma::vec update_theta(arma::vec & gradient, arma::vec & theta_ctr, double step, 
                       double lambda, double fracLasso);

double update_delta(double grad_delta, double delta_old, double step, 
		    double lambda,  double wt_delta);

double penalty(arma::vec & alpha,
               arma::vec & beta,
               double delta,
               double lambda,
               double fracLasso,
	       double wt_delta);

void print_vec(arma::vec x);
void print_mat(arma::mat x);

double logdet_medcov(arma::mat MedCov);

// end of function prototypes

arma::vec grad_loss_alpha_beta(arma::vec & alpha,
				  arma::vec & beta,
				  double delta,
				  arma::mat & SampCov,
				  arma::mat & S,
				  arma::mat & Sinv,
				  int index_group){

  // compute gradient of loss function for specified index_group for alpha and beta,
  // assuming same index_group for alpha & beta in a mediation pathway
  
  arma::vec grad(2, fill::zeros);
  
  arma::mat deriv15;
  
  arma::mat A2(size(S), fill::zeros);
  
  arma::mat B = compute_B(alpha, beta, delta);
  
  arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, Sinv);
  
  arma::mat Ident  = eye(size(SampCov));
  
  arma::mat C = Ident  - ImpCov_inv * SampCov;
  
  arma::mat E = B * S * B.t();
  
  // alpha's is in first  col, beginning at row 2, ending at next to last row
  // beta is in last row, beginning in col 2, eding in next to last col
  
  // deriv for alpha
  A2(index_group+1,0) = 1.0;
  deriv15 = B * A2 * E;
  deriv15 = deriv15 + deriv15.t();
  grad(0)  = trace(ImpCov_inv * deriv15 * C);
  
  // deriv for beta
  A2(index_group+1,0) = 0.0;
  A2(A2.n_cols-1, index_group+1) = 1.0;
  deriv15 = B * A2 * E; 
  deriv15 = deriv15 + deriv15.t();
  grad(1)  = trace(ImpCov_inv * deriv15 * C);
  
  return grad;
  
}

double grad_loss_delta(arma::vec & alpha,
                       arma::vec & beta,
                       double delta,
                       arma::mat & SampCov,
                       arma::mat & S,
                       arma::mat & Sinv){
  
  // compute gradient of loss function for delta, in 
  // last row of first col of A mat
  
  double grad = 0.0;
  
  arma::mat deriv15;
  
  arma::mat A2(size(S), fill::zeros);
  
  arma::mat B = compute_B(alpha, beta, delta);
  
  arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, Sinv);
  
  arma::mat Ident  = eye(size(SampCov));
  
  arma::mat C = Ident  - ImpCov_inv * SampCov;
  
  arma::mat E = B * S * B.t();
  
  // delta  first  col, last row
    
  // deriv for delta

  A2( (A2.n_rows-1) ,0) = 1.0;

  deriv15 = B * A2 * E;
  deriv15 = deriv15 + deriv15.t();
  grad  = trace(ImpCov_inv * deriv15 * C);
  
  return grad;
}

double grad_loss_vx(arma::vec & alpha,
                        arma::vec & beta,
                        double delta,
                        arma::mat & SampCov,
                        arma::mat & Sinv){
  
  // compute gradient of loss function for vx
  
  arma::vec grad(2, fill::zeros);
  
  arma::mat deriv15;
  
  arma::mat S2(size(Sinv), fill::zeros);
  
  arma::mat B = compute_B(alpha, beta, delta);
 
  arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, Sinv);

  arma::mat Ident  = eye(size(SampCov));
  
  arma::mat C = Ident  - ImpCov_inv * SampCov;
    
  // gradient for var(x)

  S2(0,0) = 1.0;
  deriv15 =  B * S2 * B.t(); 
  grad(0)  = trace(ImpCov_inv * deriv15 * C);

  return grad(0);
}

double grad_loss_vy(arma::vec & alpha,
                        arma::vec & beta,
                        double delta,
                        arma::mat & SampCov,
                        arma::mat & Sinv){
  
  // compute gradient of loss function for vy
  
  arma::vec grad(2, fill::zeros);
  
  arma::mat deriv15;
  
  arma::mat S2(size(Sinv), fill::zeros);
  
  arma::mat B = compute_B(alpha, beta, delta);
 
  arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, Sinv);

  arma::mat Ident  = eye(size(SampCov));
  
  arma::mat C = Ident  - ImpCov_inv * SampCov;
    
  
  S2(S2.n_rows-1, S2.n_cols-1) = 1.0;
  deriv15 =  B * S2 * B.t(); 
  grad(0)  = trace(ImpCov_inv * deriv15 * C);
  
  return grad(0);
  
}

arma::mat compute_B(arma::vec & alpha, arma::vec & beta, double delta){
  
  // compute  B = inv(I-A), but compute B with analytic solution
  
  // adding small value to diag(1). Assume that
  // diag(1) is replaced with k, where k > 1

  double k = 1.1;

  int nvar = alpha.size() + 2;
  
  arma::mat B(nvar, nvar);
  B.eye();

  for(int i = 0; i < B.n_rows; i++){
    B(i,i) = 1.0/k;
  }

  
  double k2 = k*k;
  double k3 = k2*k;

  double sum = 0.0;

  for(int i = 0; i < alpha.size(); i++){
    B(i+1,0) = alpha(i)/k2;
    B( (nvar-1), i+1 ) = beta(i)/k2;
    sum += alpha(i) * beta(i);
  }

  
  B( (nvar-1), 0 )  = delta/k2 + sum/k3;
  
  return B;
}



arma::mat inverse_ImpCov(arma::vec & alpha, arma::vec & beta, double delta, arma::mat & Sinv){
  

  int nvar = alpha.size() + 2;

  arma::mat ImpCov_inv(size(Sinv), fill::zeros);
  
  // (I-A) matrix
  arma::mat ImA(size(Sinv));
  ImA.eye();
  
  for(int i = 0; i < alpha.size(); i++){
    ImA(i+1, 0) = -alpha(i);
    ImA((nvar-1), i+1) = -beta(i);
  }
  ImA((nvar-1), 0) = -delta;
  
 
  // the step below takes most time to compute

  ImpCov_inv = ImA.t() * Sinv * ImA;
  

  return ImpCov_inv;
}

double loss_unpenalized(arma::mat & ImpCov, arma::mat &SampCov, 
                        arma::vec & alpha, arma::vec & beta, 
                        double delta, arma::mat & Sinv, double logdetMedCov){
  
  arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, Sinv);
  int nrow = ImpCov.n_rows;
  double vx = ImpCov(0,0);
  double vy = ImpCov(nrow-1, nrow-1);
  double logdet_ImpCov = logdetMedCov + log(vx) + log(vy);
  double loss = logdet_ImpCov  + trace(SampCov * ImpCov_inv);

  return loss;
}


double soft_threshold(double z, double lambda){
  
  double result;
  
  double sgn = 1.0;
  if(z < 0.0) sgn = -1.0;
  double diff = fabs(z) - lambda;
  if(diff < 0.0) {
    result = 0.0;
  } else { 
    result = sgn * diff;
  }
  return result;
}

arma::vec update_theta(arma::vec & gradient, arma::vec & theta, double step, 
                       double lambda, double fracLasso){
  
  arma::vec diff(theta.size(), fill::zeros);
  arma::vec soft_diff(theta.size(), fill::zeros);
  arma::vec theta_new(theta.size(), fill::zeros);
  
  double soft2 = 0.0;
  
  diff = theta - step*gradient;

  for(int i = 0; i < diff.size(); i++){ 
    soft_diff(i) = soft_threshold(diff(i), step*lambda*fracLasso);
    soft2 += soft_diff(i) * soft_diff(i);
  }
  
  double denom = sqrt(soft2);
  double mult = (1.0 - step*(1.0-fracLasso)*lambda/denom);
  
  if( (mult <= 0.0)  || (denom < step*(1.0-fracLasso)*lambda) ){
    for(int i = 0; i <  diff.size(); i++){
      theta_new(i) = 0.0;
    }
  } else {
    for(int i = 0; i <  diff.size(); i++){
      theta_new(i) = mult * soft_diff(i);
    }
  }
  
  return theta_new;
}

double update_delta(double grad_delta, double delta_old, double step, 
		    double lambda, double wt_delta){
  double diff  = delta_old - step*grad_delta;
  
  double delta_new = soft_threshold(diff, step*lambda*wt_delta);
  return delta_new;
}

double penalty(arma::vec & alpha,
               arma::vec & beta,
               double delta,
               double lambda,
               double fracLasso,
	       double wt_delta){
  
  // sparse group lasso 
  
  // group = alpha and beta for x -> (alpha) m -> (beta) y
  
   
  double pen_lasso = 0.0;
  double pen_group = 0.0;
  
  for(int i = 0; i < alpha.size(); i++){
    pen_lasso += abs(alpha(i)) + abs(beta(i));
    pen_group  += sqrt(alpha(i)*alpha(i) + beta(i)*beta(i));
  }
  
  double penalty =  lambda*wt_delta*abs(delta) + (lambda*(1.0 - fracLasso))*pen_group + 
    (lambda*fracLasso) * pen_lasso;
  
  return penalty;
}
/***** for testing only ****
void print_vec(arma::vec x){

  for(int i = 0; i < x.size(); i++){
    cout << x(i)  << ", ";
  }
  return;
}
void print_mat(arma::mat x){

  for(int i = 0; i < x.n_rows; i++){
    for(int j = 0; j < x.n_cols; j++){
      cout << setprecision(3) << x(i,j)  << ", ";
    }
    cout << endl;
  }
  cout << endl;

  return;
}
*****/


double logdet_medcov(arma::mat MedCov){

  // compute the contribution of MedCov to logdet of implied cov matrix
 
  arma::mat lower = chol(MedCov, "lower");
  int nmed = MedCov.n_rows;
  double logdet = 0.0;
  for(int i = 0; i < nmed; i++){
    logdet += log(lower(i,i));
  }
  logdet = 2.0 * logdet;
  return logdet;
}

		  
  
//[[Rcpp::export]]
List  rcpp_regmed(arma::vec alpha,
                  arma::vec beta,
                  double delta,
		  double vary,
		  double varx,
		  arma::mat SampCov,
		  arma::mat MedCov,
		  double sample_size,
		  double fracLasso,
		  double lambda,
		  double wt_delta,
		  int max_iter,
		  int max_iter_inner,
		  double tol=1e-5,
		  double vary_step_size = .05, // DJS added
		  double step_multiplier =.5,  // DJS changed default to .5
		  bool verbose=false){
  
  // fit regularized mediation model via structural equation model
  
  bool converge = false;
 
  arma::mat ImpCov(size(SampCov), fill::zeros);
  
  //  S matrix arranged as (vx, 0..0,   0
  //                        0, MedCov, 0
  //                        0, 0..0,  vy)
  
  // S has dim nmed + 2
  
  arma::mat S(size(SampCov), fill::zeros);
  int nmed  = MedCov.n_rows;
  S(0,0) = varx;
  S(nmed+1, nmed+1) = vary;
  for(int i = 0; i < nmed; i++){
    for(int j = 0; j < nmed; j++){
      S(i+1, j+1) = MedCov(i,j);
    }
  }
  
  arma::mat MedCov_inv = pinv(MedCov);
  
  // setup Sinv arranged as S
  arma::mat Sinv(size(SampCov), fill::zeros);
  Sinv(0,0) = 1.0/varx;
  Sinv(nmed+1, nmed+1) = 1.0/vary;
  for(int i = 0; i < nmed; i++){
    for(int j = 0; j < nmed; j++){
      Sinv(i+1, j+1) = MedCov_inv(i,j);
    }
  } 

  // compute the contribution of MedCov to the logdet of implied cov
  // matrix. This part is constant because MedCov is assumed fixed and
  // mediators independent from x and y (corresponding to arrangement of 
  // S matrix)

  double logdetMedCov = logdet_medcov(MedCov);

  // work arrays for iterations
  // B = inv(I-A)
  arma::mat B(size(SampCov), fill::zeros); 
  
  arma::vec gradient(2, fill::zeros);
  arma::vec theta_old(2, fill::zeros);
  arma::vec theta_new(2, fill::zeros);
  arma::vec theta_ctr(2, fill::zeros);
  arma::vec theta_diff(2, fill::zeros);
  
  double diff2 = 0.0;
  double gdiff = 0.0;
  double grad_delta = 0.0;
  double delta_old = 0.0;
  double delta_new = 0.0;
  double delta_diff = 0.0;
  double loss_majorize = 0.0;
 
  arma::vec alpha_temp(size(alpha), fill::zeros);
  arma::vec beta_temp(size(beta), fill::zeros);

  int iter_inner = 0;
  
  double varx_old = 0.0;
  double varx_new = 0.0;
  double varx_diff = 0.0;
  
  double vary_old = 0.0;
  double vary_new = 0.0;
  double vary_diff = 0.0;
  
  double loss_old = 0.0;
  double loss_new = 0.0;
  double step = 1.0;
  double pen_loss_begin = 0.0;
  double pen_loss_old = 0.0;
  double pen_loss_end = 0.0;

  // compute implied cov
  B = compute_B(alpha, beta, delta);
  ImpCov = B * S * B.t();
 
  loss_old =  loss_unpenalized(ImpCov, SampCov, alpha, beta, delta, Sinv, logdetMedCov);
  
  pen_loss_old =  loss_old + penalty(alpha, beta, delta, lambda, fracLasso, wt_delta);
  pen_loss_begin = pen_loss_old;
  
  double pen_loss_new = 0.0;
  double grad_varx = 0.0;
  double grad_vary = 0.0;

  bool converge_inner = false;
   
  // outer loop for all parameters
  int iter=1;
  converge = false;
 
  while( (iter <= max_iter) & !converge ){
    
    if(verbose & ((iter % 100) == 0)) Rcout << "iter = " << iter << endl;

    //============================== update delta ======================================//

    iter_inner = 0;
    converge_inner = false;
    while( (iter_inner < max_iter_inner) & !converge_inner){

      // update gradient

      delta_old = delta;
      grad_delta =  grad_loss_delta(alpha, beta, delta, SampCov, S, Sinv);
      loss_new =  loss_unpenalized(ImpCov, SampCov, alpha, beta, delta, Sinv, logdetMedCov);
    
      // optimize step size
  
      step = 0.5; // DJS changed from 1.0 to 0.5

      for(int i=0; i<10; i++){

        delta_new = update_delta(grad_delta, delta_old, step, lambda, wt_delta);
	
	loss_new =  loss_unpenalized(ImpCov, SampCov, alpha, beta, delta_new, Sinv, logdetMedCov);

	delta_diff = delta_new - delta_old;
	gdiff = grad_delta * delta_diff;
	diff2 = delta_diff * delta_diff;
	loss_majorize = loss_old + gdiff + diff2 / (2.0*step);

	if(loss_new <= loss_majorize){
	  break;
	}
	step = step_multiplier * step;
      } // end of step optimization
      
      // Nesterov step
      delta = delta_old + delta_diff * double(iter_inner + 1)/double(iter_inner + 4) ;

    
      // since delta updated, update loss 
      B = compute_B(alpha,beta, delta);
      ImpCov = B * S * B.t();
      loss_new = loss_unpenalized(ImpCov, SampCov, alpha, beta, delta, Sinv,logdetMedCov);
     
      // check convergence
      pen_loss_new = loss_new + penalty(alpha, beta, delta, lambda, fracLasso, wt_delta);
     
      if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0) ) {
	converge_inner = true;
      }

      pen_loss_old = pen_loss_new;
      loss_old = loss_new;
      iter_inner ++;

    } // end inner loop for delta


    //========= update  mediator groups (1, ..., nmed) ===================================//

    
    for(int index_group = 0; index_group < nmed; index_group ++){
      
      alpha_temp = alpha;
      beta_temp = beta;
      
      // inner loop to optimize over a single group
      iter_inner = 0;
      converge_inner = false;
      
      while( (iter_inner < max_iter_inner) & !converge_inner ){

	theta_old(0) = alpha(index_group);
	theta_old(1) = beta(index_group);
	theta_ctr = theta_old;
      
        // update vector of gradients for alpha, beta
	gradient =  grad_loss_alpha_beta(alpha, beta, delta, SampCov, S, Sinv, index_group);

        // optimize step size
        
        step = 0.5; // DJS changed from 1.0 to 0.5

        for(int i=0; i<10; i++){
      
          theta_new = update_theta(gradient, theta_ctr, step, lambda, fracLasso);
	  alpha_temp(index_group) = theta_new(0);
          beta_temp(index_group)  = theta_new(1);	
          loss_new =  loss_unpenalized(ImpCov, SampCov, alpha_temp, beta_temp, delta, Sinv, logdetMedCov);
	  
	  theta_diff = theta_new - theta_old;
          gdiff = as_scalar( gradient.t() * theta_diff );
          diff2 = as_scalar( theta_diff.t() * theta_diff );
          loss_majorize = loss_old + gdiff + diff2 / (2.0*step);
         
          if(loss_new <= loss_majorize){
            break;
          }
          
          step = step_multiplier * step;

	} // end step optimize
        
        // Nesterov step
	theta_new = theta_old + theta_diff*double(iter_inner + 1)/double(iter_inner + 4);
        alpha(index_group) = theta_new(0);
        beta(index_group) = theta_new(1);

	// since alpha, beta updated, update loss
        B = compute_B(alpha,beta, delta);
        ImpCov = B * S * B.t();
        loss_new = loss_unpenalized(ImpCov, SampCov, alpha, beta, delta, Sinv, logdetMedCov);
	
	// check convergence
        pen_loss_new =  loss_new + penalty(alpha, beta, delta, lambda, fracLasso, wt_delta);
        
	if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0) ) {
          converge_inner = true;
        }


        pen_loss_old = pen_loss_new;
	loss_old = loss_new;

	iter_inner ++;

      } // end inner loop to optimize for a single group
      
  

    } // end loop over mediator groups
   
    // skip update of var(x) since doesn't change
 

    //================================= update var(y) =======================================//

  
    iter_inner = 0;
    converge_inner = false;

    while( (iter_inner < max_iter_inner) & !converge_inner ){
       
      // update gradient

      vary_old = vary;
      grad_vary =  grad_loss_vy(alpha, beta, delta, SampCov, Sinv);
 
      // optimize step size
      
      step = vary_step_size; // DJS changed from 1.0 to vary_step_size
      
      for(int i=0; i<10; i++){
	
	vary_new = vary_old - step*grad_vary;
	
	// bound var away from 0
	if(vary_new < 0.0001){
	  step = step_multiplier * step;
	  continue;
	}
	
	S(S.n_rows-1,S.n_cols-1) = vary_new;
	Sinv(Sinv.n_rows-1,Sinv.n_cols-1) = 1.0/vary_new;
	ImpCov = B * S * B.t();

	loss_new =  loss_unpenalized(ImpCov, SampCov, alpha, beta, delta_new, Sinv, logdetMedCov);

	// DJS changed to loss_old below
	if(loss_new <= loss_old){
	  break;
	}
	step = step_multiplier * step;
      }
          
      // if after step opt, loss_new > loss_old, revert to loss_old and parm_old
       if(loss_old < loss_new){
	vary = vary_old;
	S(S.n_rows-1,S.n_cols-1) = vary_old;
	Sinv(Sinv.n_rows-1,Sinv.n_cols-1) = 1.0/vary_old;
	ImpCov = B * S * B.t();
	loss_new =  loss_unpenalized(ImpCov, SampCov, alpha, beta, delta_new, Sinv, logdetMedCov);
       }  else{
           vary = vary_new;
       }
      // check convergence
        
 
      pen_loss_new =  loss_new + penalty(alpha, beta, delta, lambda, fracLasso, wt_delta);
      
      if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0) ) {
	converge_inner = true;
      }

      pen_loss_old = pen_loss_new;
      loss_old = loss_new;
      iter_inner ++;

    }
   

    // this completes sequentially optimizing over each parameter
    // check for global convergence by comparing the penalized loss function
    // at the beginning of a global iter to the end of a global iter
         
    if (fabs(pen_loss_new - pen_loss_begin) < tol * (fabs(pen_loss_new) + 1.0) ) {
	converge = true;
    }

    iter++;
    pen_loss_begin = pen_loss_new;

    
  } // end outer loop
  


  // compute BIC using threshold of eps to determine if 
  // a parameter is non-zero

  double eps = 0.001;
  double df =  0.0;
  for(int i = 0; i < alpha.size(); i++){
    if(fabs(alpha(i))  >  eps) {df++; }
    if(fabs(beta(i))   >  eps) {df++; }
  }
  
  if(fabs(delta) > eps) {df++; }
  if(fabs(vary) > eps)  {df++; }

  // DJS: should not have multiplied loss by 2
  double bic =  loss_new*sample_size + log(sample_size)*df;
 
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("delta") = delta,
			    Rcpp::Named("var.x") = varx,
			    Rcpp::Named("var.y") = vary,
                            Rcpp::Named("converge") = converge,
			    Rcpp::Named("iter") = iter,
			    Rcpp::Named("loss") = loss_new,
			    Rcpp::Named("penloss") = pen_loss_new,
			    Rcpp::Named("bic") = bic,
			    Rcpp::Named("df") = df);

}

