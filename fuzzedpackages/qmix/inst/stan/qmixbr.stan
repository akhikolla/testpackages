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
  
  
}

data {
  
  int N;
  int D;
  vector<lower=0,upper=1>[N] Y;
  matrix[N,D] X;
  int k;
  real offset;
  

  
}

transformed data{
}

parameters {
  vector[D] beta[k];
  vector<lower = 0, upper = 1>[k] p_tmp; 
  simplex[k] theta;
  
}

transformed parameters{
  vector[k] p;
  p = sort_asc(p_tmp);
  
}

model{
  real lik;
  for (i in 1:k){
    beta[i] ~ normal(0,10);  
  }
  
  theta ~ dirichlet(rep_vector(1.0, k));
  p_tmp ~ uniform(0,1);
  
  
  
  for (i in 1:N){
    if (Y[i] == 1){
      lik = 0;
      for (j in 1:k){
        lik = lik + theta[j] * pald2((dot_product(X[i,],beta[j])),p[j]);
      }
      lik =  lik + offset;  
    }
    if (Y[i] == 0){
      lik = 0;
      for (j in 1:k){
        lik = lik + theta[j] * (1 - pald2((dot_product(X[i,],beta[j])),p[j]));
      }
      lik = lik + offset;  
    }
    
    target += log(lik);
  }
}



