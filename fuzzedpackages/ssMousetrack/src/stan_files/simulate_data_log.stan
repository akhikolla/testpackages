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
  
  vector compute_deltay(vector x, int I, real pT, real pD, real pC){
    vector[I] dx;
    for(i in 1:I){
      if(x[i]>=pC) dx[i] = fabs(x[i]-pD);
      else dx[i] = fabs(x[i]-pT);}
    return dx;
  }
}

data{
  int I; // number of individuals
  int N; // length of Y-trajectories
  int J; // number of trials (the same for each individual)
  int KK; // total number of categorical levels minus one (number of columns of the model matrix Z)
  vector[I*J] Y[N]; // NxIJ matrix of Y-trajectories
  vector<lower=0>[I] sigmaz; // sigmax parameter for the latent dynamics
  matrix<lower=0,upper=pi()>[I*J,3] bnds; // matrix of bounds (lb,ub,ub-lb) for sampling the Y-trajectories
  matrix[I*J,KK] D; // IJxKK matrix for delta-values
  vector<lower=0>[I*J] lambda_vec; // lambda parameter to compute the kappa of the von-mises distribution
  matrix[KK,5] priors_matrix; //prior distributions and associated parameters
  real pT; // position of TARGET
  real pD; // position of DISTRACTOR
  real pC; // (pT+pD)/2
  real<lower=0> kappa_lb; // lower bound (fixed by the user) for the kappa parameter
  real<lower=0> kappa_ub; // upper bound (fixed by the user) for the kappa parameter
}

parameters{
  vector[KK] gamma;
}

transformed parameters{
  vector[I*J] b;

  b = D*gamma;
}

model{
  
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
}

generated quantities{
  vector[I] z[N]; // NxI matrix of latent dynamics
  vector[I*J] mu[N]; // NxIJ matrix of von-mises means
  vector[I*J] y_sim[N]; // NxIJ matrix of simulated Y-trajectories
  vector[I*J] dy_sim[N]; // NxIJ matrix of delta values
  vector[I*J] z_vec; // working variable
  vector<lower=kappa_lb,upper=kappa_ub>[I*J] kappas[N]; // working variable

  // Sampling the latent dynamics
  z[1] = rep_vector(1e-04,I);
  for(n in 2:N){
    for(i in 1:I){
      z[n,i] = normal_rng(z[n-1,i],sigmaz[i]);
    }
  }
  
  // Sampling the Y-trajectories
  for(n in 1:N){
    if(n==1){
      mu[n] = rep_vector(pi()/2,I*J);
      dy_sim[n] = compute_deltay(mu[n],I*J,pT,pD,pC); //in this case, dy=(pi/2) as mu[1] == pC
      kappas[n] = rep_vector(kappa_ub,I*J);
    }else{
      z_vec = kronecker_simplified_J(I,J,z[n]);
      mu[n] = bnds[,1] + (bnds[,3] ./ ( 1 + exp(b - z_vec))); // generalized logistic function  
      dy_sim[n] = compute_deltay(mu[n],I*J,pT,pD,pC) + multi_normal_rng(rep_vector(0,I*J),diag_matrix(rep_vector(2.5e-2,I*J)));
      kappas[n] = rescale_data((exp(lambda_vec .* dy_sim[n])),kappa_lb,kappa_ub);
    }
    for(q in 1:(I*J))
      y_sim[n,q] = von_mises_rng(mu[n,q],kappas[n,q]);
  }

}

