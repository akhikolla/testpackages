#include "EMVS.h"

// [[Rcpp::export(".conj_EMVS")]]
SEXP conj_EMVS(SEXP Y_,
	  SEXP X_,
	  SEXP v0s_,
	  SEXP v1_,
	  SEXP type_,
	  SEXP beta_init_,
	  SEXP sigma_init_,
	  SEXP epsilon_,
	  SEXP temperature_,
	  SEXP theta_,
	  SEXP a_,
	  SEXP b_,
	  SEXP v1_g_,
	  SEXP direction_){

  Rcout << "Loading Data " << endl;
 
  vec Y = as<vec>(Y_);  
  mat X = as<mat>(X_);  
  vec v0s = as<vec>(v0s_);
  double v1 = as<double>(v1_);
  const string type = as<string>(type_);
  vec beta_init = as<vec>(beta_init_);
  double sigma_init = as<double>(sigma_init_);
  double epsilon = as<double>(epsilon_);
  double temperature = as<double>(temperature_);
  double theta = as<double>(theta_);
  double a = as<double>(a_);
  double b = as<double>(b_);
  double v1_g = as<double>(v1_g_);
  const string direction = as<string>(direction_);

  const int L = v0s.n_elem;
  const int dim = X.n_cols;
 
  vec intersects = zeros<vec>(L);
  vec sigmas = zeros<vec>(L);
  mat betas = zeros<mat>(L,dim);
  mat posts = zeros<mat>(L,dim);  
  vec log_post = zeros<vec>(L);
  vec niters = zeros<vec>(L);
  vec thetas = zeros<vec>(L);

  mat E_step;
  mat XtY=X.t()*Y;
  mat XtX=X.t()*X;
  mat Xt = X.t(); 
  uvec index;
  vec beta_k, beta_new,theta_ks;
  vec Y_copy=Y;
  double sigma_k, c, w, theta_k;
  vec inv_var;
  vec post;
  double v0,  eps;
  int niter;
 
  beta_new = beta_init;

  // Begin EMVS regularization

  Rcout << "EMVS Begins" << endl;

  for(int i = 0; i < L; ++i){
  
    if ((direction.compare("forward")==0)|(direction.compare("null")==0)) {
    v0 = v0s[i];
    } else if (direction.compare("backward")==0){
      v0=v0s[L-1-i];  
    }
        Rcout << "v0 = " << v0 << endl;

    // Initialization of parameters

    if ((direction.compare("forward")==0)| (direction.compare("backward")==0)) {
	   beta_k = beta_new;
	  } else {
      beta_k = beta_init;
    }

    beta_new = beta_k;
    sigma_k = sigma_init;
    
    
     if(type.compare("betabinomial") == 0){
       theta_k = 0.5;
    } else if (type.compare("fixed") == 0){
       theta_k = theta;
    }

    eps = epsilon + 1;
    niter = 1;

    while(eps>epsilon){

  
      // ******** E-step ************ //

      // Update inclusion probabilities

    
      E_step = conj_E_beta_binom(beta_k, sigma_k, v0, v1, theta_k, temperature);

      inv_var = E_step.col(0);
  
      post = E_step.col(1);
  
  
      
      // ******** M-step ************ //

      beta_k = conj_M_beta(XtY, X, Y, XtX, inv_var);
      
    
      sigma_k = conj_M_sigma(Y,X,beta_k, inv_var, 1, 1);
    
    
      if(type.compare("betabinomial") == 0){
	     theta_k = M_p(post, a, b);
      } 

      //eps = max(abs(beta_new - beta_k));
      eps=accu((beta_new - beta_k)%(beta_new - beta_k));
      beta_new = beta_k;
      niter++;
      Rcout << "Epsilon " << eps << endl;

    }
    
    index = find(post > 0.5);
    
    // Store values:
    if ((direction.compare("forward")==0)|(direction.compare("null")==0)){
      posts.row(i) = post.t();
      betas.row(i) = beta_new.t();
      sigmas[i] = sigma_k;
      niters[i] = niter;
      
      if (index.n_elem < 1000){
        log_post[i] = log_g(index,X,Y,1,1,0,v1_g,type,a,b);
      } else {
        log_post[i]=datum::nan;
      } 
      
      if(( type.compare("betabinomial") == 0) || (type.compare("fixed") == 0)){
        thetas[i]=theta_k;}
      
      c = sqrt(v1/v0s[i]);
      if(( type.compare("betabinomial") == 0) || (type.compare("fixed") == 0)){
        w = (1-theta_k)/theta_k;
        intersects[i] = sigmas[i];
        intersects[i] *= sqrt(v0s[i]);
        intersects[i] *= sqrt(2*log(w*c)*c*c/(c*c-1));
      } 
      
    } else if (direction.compare("backward")==0){
      posts.row(L-1-i) = post.t();
      betas.row(L-1-i) = beta_new.t();
      sigmas[L-1-i] = sigma_k;
      niters[L-1-i] = niter;
      
      if (index.n_elem < 1000){
        log_post[L-1-i] = log_g(index,X,Y,1,1,0,v1_g,type,a,b);
      } else {
        log_post[L-1-i]=datum::nan;
      }
      
      if(( type.compare("betabinomial") == 0) || (type.compare("fixed") == 0)){
        thetas[L-1-i]=theta_k;}
      
      c = sqrt(v1/v0s[L-1-i]);
      if(( type.compare("betabinomial") == 0) || (type.compare("fixed") == 0)){
        w = (1-theta_k)/theta_k;
        intersects[L-1-i] = sigmas[L-1-i];
        intersects[L-1-i] *= sqrt(v0s[L-1-i]);
        intersects[L-1-i] *= sqrt(2*log(w*c)*c*c/(c*c-1));
      } 
      
    }
  
  }
  

  // Wrap the results into a list and return.
  List list;
  list["betas"] = betas;
  list["log_g_function"] = log_post;
  list["intersects"] = intersects;
  list["sigmas"] = sigmas;
  list["v1"] = v1;
  list["v0"] = v0s;
  list["niters"] = niters;
  list["prob_inclusion"] = posts;
  list["type"] = type;
  list["direction"]=direction;
  list["theta"] = thetas;
  
  Rcout << "Done! " << endl;
  return list;
}

