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

//## function to return the number of observations in a group
  int group_size(int[] ref, int value) {
    int count;
    count = 0;
    for (ii in 1:size(ref))
      if (ref[ii]==value)
        count = count + 1;
      return count;
  }
  
  //## function to subset an integer array (return just those observations in a given group)
  int[] subset_intarray(int[] y, int[] ref, int value) {
    int jj;
    int res[group_size(ref, value)];
    if (size(ref) != size(y))
      reject("illegal input: non-matching dimensions")
    jj = 1;
    for(ii in 1:size(ref)) {
      if (ref[ii] == value) {
        res[jj] = y[ii];
        jj = jj+1;
      }
    }
    return res;
  }

}

/////////////////////////////////////////////////////////////////////////////////////
// Model
/////////////////////////////////////////////////////////////////////////////////////

data {
	int N; //# number of observations
	int D_common; //# number of common covariates
	vector<lower=-1,upper=1>[N] Y; //# data of choices y = 0 -> -1, y = 1 -> 1
	matrix[N,D_common] X_common; //# common covariates
	int N_indx; //# number of groups
	int ind[N]; //# group index
  int N_person; // number of individuals
  int person[N]; // person index
  int N_wave; // number of waves
  int wave[N];// wave index
  real q; // quantile
  real offset;
}

transformed data{
  int n_group[N_indx]; //# number of observations in the group
  for (ii in 1:N_indx) {
    n_group[ii] = group_size(ind, ii);
  }
}

parameters {
	vector[D_common] beta;
  vector[N_person] beta_ind;
  vector[N_wave] beta_wave;
  real<lower=0> sigma_beta_ind;

}

transformed parameters{
	
}

model{
	int pos;
	vector[N] xb_common;
  
  sigma_beta_ind ~ cauchy(0,1);

	beta ~ normal(0,10);
  beta_ind ~ normal(0,sigma_beta_ind);
  beta_wave ~ normal(0,10);

  

  	for (i in 1:N){
  		xb_common[i] = X_common[i,]*beta + beta_ind[person[i]] + beta_wave[wave[i]];
  	}
  	
	pos = 1;

  	for (i in 1:N_indx){
      real lik0;
  		real lik;
      real lik2;
      real lik3;
      real lik4;
  		vector[n_group[i]] y_g;
  		vector[n_group[i]] xb_common_g;
  		vector[n_group[i]] ystar_g;

  		y_g = segment(Y, pos, n_group[i]);
  		xb_common_g = segment(xb_common,pos,n_group[i]);

  		lik = 1;
  		for (j in 1:(n_group[i]-1) ){ 
              //needs to order the data with the y==1 at the end of each choice set
  			lik =	lik * (1 - pald2(-(xb_common_g[n_group[i]] - xb_common_g[j] ),q));
  		}

      lik3=1;
      for (j in 1:(n_group[i]-1)){
        lik0 = 1;
        for (k in 1:(n_group[i])){
          if (j != k){
            lik0 = lik0 * (1 - pald2( - (xb_common_g[j] - xb_common_g[k] ),q)) ;
          }
        }
        lik3 = lik3*(1-lik0);
      }

      lik2 = (lik+offset) * lik3;

     	target += log(lik2);
     	pos = pos + n_group[i];
  }

}



