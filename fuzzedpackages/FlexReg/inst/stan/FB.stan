data {
	int<lower=1> N; // number of observations
	vector[N]    y; // outcomes
	int<lower=1> K; // number of covariates+intercept for mu
	matrix[N,K]    X; // covariate
	int<lower=1> link_code_mu;
	int<lower=1> link_prior_beta;
	real hyperprior_beta;  // Prior standard deviation
	int<lower=1> prior_code_phi;
	real<lower=0> hyper_phi;
}

parameters {
	vector[K] beta;
	real<lower=0> phi;
	real<lower=0,upper=1> w;
	real<lower=0,upper=1> p;
}


transformed parameters {
	vector<lower=0,upper=1>[N]  mu;
	vector<lower=0,upper=1>[N]  lambda1;
	vector<lower=0,upper=1>[N]  lambda2;
	vector<lower=0>[N] b1;
	vector<lower=0>[N] b2;
	vector<lower=0>[N] a1;
	vector<lower=0>[N] a2;
	vector<lower=0>[N] wtilde;

	if(link_code_mu == 1)
		mu = inv_logit(X * beta);
	else if(link_code_mu == 2)
		mu = Phi(X * beta);
	else if(link_code_mu == 3)
		mu = inv_cloglog(X * beta);
	else if(link_code_mu == 4)
		mu = exp(-exp(X * beta));

	for (i in 1:N){
	  wtilde[i] = w*fmin(mu[i]/p,(1-mu[i])/(1-p));
	  lambda1[i] = mu[i]+(1-p)*wtilde[i];
	  lambda2[i] = mu[i]-p*wtilde[i];
	//con probabilità p
	b2[i] = (1-lambda1[i])*phi;
	a2[i] = lambda1[i]*phi;
	//con probabilità 1-p
	b1[i] = (1-lambda2[i])*phi;
	a1[i] = lambda2[i]*phi;
	}
}

model {
	//priors
	if(prior_code_phi == 1)
	phi ~ gamma(hyper_phi,hyper_phi);
	else if(prior_code_phi == 2)
	phi ~ uniform(0.0,hyper_phi);

	for (l in 1:K) {
	if(link_prior_beta == 1)
		beta[l] ~ normal(0, hyperprior_beta);
	else if(link_prior_beta == 2)
		beta[l] ~ cauchy(0, hyperprior_beta);
	}
	// likelihood of log(y)

	for(i in 1:N)
	target += log_mix(p, beta_lpdf(y[i] | a2[i],b2[i]),
                         beta_lpdf(y[i] | a1[i], b1[i]));
}

generated quantities{
	vector[N] log_lik;
	for(i in 1:N){
		log_lik[i] = log_mix(p, beta_lpdf(y[i] | a2[i],b2[i]),
                         beta_lpdf(y[i] | a1[i], b1[i]));
	}
}   
