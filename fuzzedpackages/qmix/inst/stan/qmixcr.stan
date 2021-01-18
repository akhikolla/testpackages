functions{

  real pald2(real mu, real p){
    real prob;
    if (mu <0){
      prob = p*exp((mu)*(1-p));
    } else {
      prob = 1-(1-p)*exp(-(mu)*(p));
    }
    return(prob);
}

 real cald_lpdf(real x, real theta, real p, real sigma){
  real prob;
  real out;
  if (x - theta < 0){
    prob = p*(1-p)/sigma*exp(-(x - theta)*(p - 1)/sigma);
  } else {
    prob = p*(1-p)/sigma*exp(- (x - theta)*p/sigma);
  }
  out = log(prob);
  return(out);
 }

}

data {
  int N; 
  int D; 
  vector[N] Y; 
  matrix[N,D] X; 
  int k;
}


parameters {
  vector[D] beta[k];
  real<lower=0> sigma[k];
  //real<lower=0> sigma;
  simplex[k] theta;
  vector<lower = 0, upper = 1>[k] p_tmp;
}

transformed parameters{
  vector[k] p;
  p = sort_asc(p_tmp);
}


model{
  vector[k] prob_tmp;
  for (i in 1:k){
    beta[i] ~ normal(0,10);
  }
  sigma ~ cauchy(0, 1);
  theta ~ dirichlet(rep_vector(1.0, k));
  p_tmp ~ uniform(0,1);

  for (i in 1:N){
    for (j in 1:k){
      prob_tmp[j] = log(theta[j]) + cald_lpdf(Y[i] | dot_product(X[i,],beta[j]), p[j], sigma[j]);
      //prob_tmp[j] = log(theta[j]) + cald_lpdf(Y[i] | dot_product(X[i,],beta[j]), p[j], sigma);
    }
    target += log_sum_exp(prob_tmp);
  }
  
  
}



