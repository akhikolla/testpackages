functions{
  
  vector kronecker_simplified_J(int I, int J, vector x){ 
    matrix[I,J] X;
    vector[I*J] x_vec;
    for(i in 1:I){
      X[i] = rep_row_vector(x[i],J);}
    x_vec = to_vector(X');
    return x_vec;
  }
  
  vector kronecker_simplified_I(int I, int J, vector x){ 
    matrix[I,J] X;
    vector[I*J] x_vec;
    for(i in 1:I){
      X[i] = to_row_vector(x);}
    x_vec = to_vector(X');
    return x_vec;
  }
  
  vector rescale_data(vector x, real min_x, real max_x){
    int I = rows(x);
    real max_Vx = max(x);
    real min_Vx = min(x);
    return rep_vector(min_x,I) + ( (x-min_Vx) ./ rep_vector(max_Vx-min_Vx,I)) .* rep_vector(max_x-min_x,I);
  }
  
}

data{
  int I; // number of individuals
  int N; // length of Y-trajectories
  int J; // number of trials (the same for each individual)
  int KK; // total number of categorical levels minus one (number of columns of the model matrix Z)
  vector[I*J] Y[N]; // NxIJ matrix of Y-trajectories
  vector[I*J] DY[N]; //NxIJ matrix for delta-values
  vector<lower=0>[I] sigmaz; // sigmax parameter for the latent dynamics
  matrix<lower=0,upper=pi()>[I*J,3] bnds; // matrix of bounds (lb,ub,ub-lb) for sampling the Y-trajectories
  matrix[I*J,KK] D; // IJxKK matrix of the exp design 
  vector<lower=0>[I*J] a; // it does not work in this version
  vector<lower=0>[I*J] lambda_vec; // lambda parameter to compute the kappa of the von-mises distribution
  matrix[I*J,I] Am; // working matrix to be used in the Kalman filter loop
  real<lower=0> kappa_lb; //lower boud (fixed by the user) for the kappa parameter
  real<lower=0> kappa_ub; //upper boud (fixed by the user) for the kappa parameter
  matrix[KK,5] priors_matrix; //prior distributions and associated parameters
}

parameters{
  vector[KK] gamma;
}

transformed parameters{
  vector[I] z_pred[N]; // latent states (mean: E[X]) in the prediction stage of the Kalman filter loop
  vector[I] z_upd[N]; // latent states (mean: E[X]) in the update stage of the Kalman filter loop
  vector<lower=0>[I] lambda_pred[N]; // latent states (variance: VAR[X]) in the prediction stage of the Kalman filter loop
  vector<lower=0>[I] lambda_upd[N]; // latent states (variance: VAR[X]) in the prediction stage of the Kalman filter loop
  vector<lower=0>[I*J] sigma_kf[N]; // working variables for the Kalman filter loop
  vector[I*J] kappa_vec; // working variables for the Kalman filter loop
  vector[I*J] y_star[N]; // predicted Y-trajectories in the Kalman filter loop
  vector[I*J] G; // working variables for the Kalman filter loop
  vector[I*J] z_vec; // working variables for the Kalman filter loop
  vector[I*J] lambda_pred_vec; // working variables for the Kalman filter loop
  vector[I*J] b; // working variables for the Kalman filter loop
  
  b = D*gamma; 
  
  // ************************************** adapted-KALMAN FILTER loop ************************************** //
  z_pred[1] = rep_vector(1e-04,I);
  lambda_pred[1] = rep_vector(1,I);
  kappa_vec = sqrt(1 ./ rep_vector(kappa_lb,I*J));
  
  for(n in 1:N){
    if(n>1){
      // Prediction
      z_pred[n] = z_upd[n-1];
      lambda_pred[n] = lambda_upd[n-1] + sigmaz;
      
      kappa_vec = sqrt(1 ./ rescale_data((exp(lambda_vec .* DY[n])),kappa_lb,kappa_ub)); 
    }
    // Working variables
    z_vec = kronecker_simplified_J(I,J,z_pred[n]);
    lambda_pred_vec = kronecker_simplified_J(I,J,lambda_pred[n]);
    y_star[n] = bnds[,1] + (bnds[,3] ./ ( 1 + exp(b - z_vec)));
    sigma_kf[n] = lambda_pred_vec + kappa_vec;
    G  = lambda_pred_vec ./ sigma_kf[n];

    // Updating
    z_upd[n] = z_pred[n] + (((Y[n]-y_star[n]) .* G)' * Am)';
    lambda_upd[n] = lambda_pred[n] - ((G .* sigma_kf[n] .* G)' * Am)';
  }
  // ************************************** ************************** ************************************** //

}

model{
  // prior for the model parameters
  
  /// betas of the stimuli equation
  for(k in 1:KK){
    if(priors_matrix[k,1]==1) gamma[k] ~ lognormal(priors_matrix[k,2],priors_matrix[k,3]);
    else if(priors_matrix[k,1]==2) gamma[k] ~ normal(priors_matrix[k,2],priors_matrix[k,3]);
    else if(priors_matrix[k,1]==201) gamma[k] ~ normal(priors_matrix[k,2],priors_matrix[k,3])T[priors_matrix[k,4],priors_matrix[k,5]];
    else if(priors_matrix[k,1]==3) gamma[k] ~ chi_square(priors_matrix[k,2]);
    else if(priors_matrix[k,1]==4) gamma[k] ~ inv_chi_square(priors_matrix[k,2]);
    else if(priors_matrix[k,1]==5) gamma[k] ~ gamma(priors_matrix[k,2],priors_matrix[k,3]);
    else if(priors_matrix[k,1]==6) gamma[k] ~ pareto(priors_matrix[k,2],priors_matrix[k,3]);
    else if(priors_matrix[k,1]==7) gamma[k] ~ uniform(priors_matrix[k,2],priors_matrix[k,3]);
  }
  
  
  // marginal likelihood of the model
  for(n in 1:N){
    Y[n] ~ multi_normal(y_star[n],diag_matrix(sigma_kf[n]));
  }
}


generated quantities{
  // Fixed-Interval backward smoother 
  vector[I] z_s_pred[N];
  vector[I] z_s_upd[N]; 
  vector[I] lambda_s_pred[N];
  vector<lower=0>[I] lambda_s_upd[N]; 
  vector[I] G_s;
  
  z_s_upd[N] = z_upd[N];
  lambda_s_upd[N] = lambda_upd[N];
  
  z_s_pred[1] = rep_vector(0,I);
  lambda_s_pred[1] = rep_vector(1,I);
  
  for(n in 1:(N-1)){
    int nn = N-n;
    
    z_s_pred[nn+1] = z_upd[nn];
    lambda_s_pred[nn+1] = lambda_upd[nn] + sigmaz;
    
    G_s = lambda_upd[nn] ./ lambda_s_pred[nn+1];
    z_s_upd[nn] = z_upd[nn] + G_s .* (z_s_upd[nn+1]-z_s_pred[nn+1]);
    lambda_s_upd[nn] = lambda_upd[nn] + G_s .* (lambda_s_upd[nn+1]-lambda_s_pred[nn+1]) .* G_s;
  }
}

