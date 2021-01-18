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

  int group_size(int[] ref, int value) {
    int count;
    count = 0;
    for (ii in 1:size(ref))
      if (ref[ii]==value)
        count = count + 1;
      return count;
  }

  int[] subset_intarray(int[] y, int[] ref, int value) {
    int jj;
    int res[group_size(ref, value)];
    if (size(ref) != size(y))
      reject("illegal input")
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


data {
	int N;
	int D_common;
	vector<lower=-1,upper=1>[N] Y;
	matrix[N,D_common] X_common;
	int N_indx; 
	int ind[N]; 
  real offset;
  real<lower = 0, upper = 1> q;
}

transformed data{
  int n_group[N_indx]; 
  for (ii in 1:N_indx) {
    n_group[ii] = group_size(ind, ii);
  }
}

parameters {
	vector[D_common] beta;
}


model{
	int pos;
	vector[N] xb_common;
	beta ~ normal(0,10);

  for (i in 1:N){
  	xb_common[i] = X_common[i,]*beta;
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

  		y_g = segment(Y, pos, n_group[i]);
  		xb_common_g = segment(xb_common,pos,n_group[i]);

  		lik = 1;
  		for (j in 1:(n_group[i]-1) ){ 
  			lik =	lik * pald2((xb_common_g[n_group[i]] - xb_common_g[j] ),q);
  		}

      lik3=1;
      for (j in 1:(n_group[i]-1)){
        lik0 = 1;
        for (k in 1:(n_group[i])){
          if (j != k){
            lik0 = lik0 * pald2((xb_common_g[j] - xb_common_g[k] ),q) ;
          }
        }
        lik3 = lik3*(1-lik0);
      }

     	lik2 = rising_factorial(n_group[i],0) * (lik+offset) * lik3;
     	target += log(lik2);
     	pos = pos + n_group[i];
  }

}



