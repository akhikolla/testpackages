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
  vector[k] p;
  real offset;
  
}

transformed data{
}

parameters {
  vector[D] beta[k];
  simplex[k] theta;
  
}

transformed parameters{
  
}

model{
  real lik;
  for (i in 1:k){
    beta[i] ~ normal(0,10);  
  }
  
  
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



