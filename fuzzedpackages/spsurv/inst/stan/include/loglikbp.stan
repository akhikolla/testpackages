// Functions block (optional).
// Here we implement a likelihood function and other useful functions.
functions{

 vector loglik_null(vector beta, vector gamma, vector status, matrix X,
   matrix b, matrix B, int M, int dist, int [] id, vector z){

    int n = num_elements(status);
    vector[n] log_lik;
    vector[n] bp;
    vector[n] BP;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      bp = b * gamma;
      BP = B * gamma;
    }
    else if(dist == 1){ //gamma
      bp = b * gamma .* z[id];
      BP = B * gamma .* z[id];
    }
    else{//additive
      bp = b * gamma .* exp(z[id]);
      BP = B * gamma .* exp(z[id]);
    }
    if(M == 0){
        log_lik = (log(bp) - log(1 + BP)) .* status - log(1 + BP); // prop odds
    }
    else{
        log_lik = log(bp) .* status - BP; // aft or prop hazards
    }
    return log_lik;
  }

 vector loglik_ph(vector beta, vector gamma, vector status,
 matrix X, matrix b, matrix B, int dist, int [] id, vector z){

    int n = num_elements(status);
    vector[n] log_lik;
    vector[n] h0;
    vector[n] H0;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      h0 = b * gamma;
      H0 = B * gamma;
    }
    else if(dist == 1){ //gamma
      h0 = b * gamma .* z[id];
      H0 = B * gamma .* z[id];
    }
    else{//additive
      h0 = b * gamma .* exp(z[id]);
      H0 = B * gamma .* exp(z[id]);
    }
    log_lik = (log(h0) + log(theta)) .* status - (H0 .* theta);
 return log_lik;
}

 vector loglik_po(vector beta, vector gamma, vector status, matrix X,
  matrix b, matrix B, int dist, int [] id, vector z){

    int n = num_elements(status);
    vector[n] log_lik;
    vector[n] r0;
    vector[n] R0;
    vector[n] theta = exp(X * beta);

    if(dist == 0){
      r0 = b * gamma;
      R0 = B * gamma;
    }
    else if(dist == 1){
      r0 = b * gamma .* z[id];
      R0 = B * gamma .* z[id];
    }
    else{//additive
      r0 = b * gamma .* exp(z[id]);
      R0 = B * gamma .* exp(z[id]);
    }
	  log_lik = (log((r0 .* theta)) - log(1 + R0 .* theta)) .* status - log(1 + R0 .* theta);
 return log_lik;
 }

 vector loglik_aft(vector time,  vector beta, vector gamma, vector status, matrix X,
    matrix b, matrix B, int dist, int [] id, vector z){

  int n = num_elements(status);
  int m = num_elements(gamma);
  vector[n] log_lik;
  vector[n] h0;
  vector[n] H0;
  vector[n] theta = exp(X * beta);
  real tau_aft = max(time ./ theta) + 0.000001;
  vector[n] y = time ./theta;
  vector[n] y_alt = y ./tau_aft;
  matrix[n, m] b2;
  matrix[n, m] B2;
  int j = 1;

    while(j < m +1){
      for(i in 1:n){
        b2[i, j] = beta_lpdf(y_alt[i] | j, (m - j + 1));
        B2[i, j] = beta_lcdf(y_alt[i]| j, (m - j + 1));
      }
      j = j+1;
    }
    b2 = exp(b2) ./ tau_aft;
    B2 = exp(B2);

    if(dist == 0){
      h0 = b2 * gamma;
      H0 = B2 * gamma;
    }
    else if(dist == 1){ //gamma multiplicative
      h0 = b2 * gamma .* z[id];
      H0 = B2 * gamma .* z[id];
    }
    else{ //additive
      h0 = b2 * gamma .* exp(z[id]);
      H0 = B2 * gamma .* exp(z[id]);
    }
	  log_lik =  ((log(h0)-log(theta)) .* status) - (H0);
 return log_lik;
 }
}
// Final line empty to avoid warnings.
