//
// Stan model to evaluated the cases of Covid-19 - Poisson model
// model: generalized static logistics (d=0;sazonal_effect)

//functions{

  // mean of the number of new cases
//  real mu_new_cases(real time, real a, real b, real c, real d, real f){
//    return(log(f)+log(a-d)+log(c)-(c*time)-(f+1)*log( b+exp(-c*time) ) );
//  }

//}


data {

  //-----------------------------
  // observed data
  int<lower=1> n; // number of observations
  int<lower=0> y[n]; // counts of new cases
  real<lower=0> pop;
  real<lower=0,upper=1> p;
  int<lower=0, upper=7> w1;
  int<lower=0, upper=7> w2;
  int<lower=0, upper=7> w3;
  //-----------------------------
}


parameters {
  real<lower=0> c;
  real<lower=1> f;
  real b01;

  real<lower=0, upper=p*pop*exp(f*b01)> a;
  real<lower=0> d_1;
  real<lower=0> d_2;
  real<lower=0> d_3;
}

transformed parameters{
  real<lower=0> b;
  real<lower=0, upper=pop> mu[n];

  b = exp(b01);
  for(t in 1:n){
    mu[t]= exp( log(f)+log(a)+log(c)-(c*t)-(f+1)*log( b+exp(-c*t) ))*d_1^( ((t-w1) % 7 == 0) * (w1>0) )*d_2^( ((t-w2) % 7 == 0) * (w2>0) )*d_3^( ((t-w3) % 7 == 0) * (w3>0) );
}
}


model {
  //----------------------------
  // likelihood function
    y ~ poisson(mu); // observed model
  //----------------------
   // prior distributions
   a ~ gamma(0.1, 0.1);
   //b ~ gamma(0.1, 0.1);  // sqrt(1/0.1)
   b01 ~ normal(0, sqrt(20));
   c ~ gamma(2,9);
   f ~ gamma(0.01,0.01);   // shape, scale
   //beta_1 ~ beta(2, 1);
   //beta_2 ~ beta(2, 1);
   d_1 ~ gamma(2,1);
   d_2 ~ gamma(2,1);
   d_3 ~ gamma(2,1);
}
