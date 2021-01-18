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
  vector<lower=-1,upper=1>[N] Y;
  matrix[N,D] X; 
  real offset;
  real<lower = 0, upper = 1> q;
}


parameters {
  vector[D] beta;

}

model{
  real lik;
  // alpha ~ normal(0,10);
  beta ~ normal(0,10);

  for (i in 1:N){
    if (Y[i] == 1){
      lik = pald2((dot_product(X[i,],beta)),q) + offset;  
    }
    if (Y[i] == 0){
      lik = 1 - pald2((dot_product(X[i,],beta)),q) + offset;  
    }
    target += log(lik);
  }
}



