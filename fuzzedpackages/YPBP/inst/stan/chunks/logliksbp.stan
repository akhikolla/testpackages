functions{

  real dbeta(real x, real shape1, real shape2){
    return(exp(beta_lpdf(x|shape1, shape2)));
  }

  real pbeta(real x, real shape1, real shape2){
    return(exp(beta_lcdf(x|shape1, shape2)));
  }



  //-------------------------------------------------------------
  // Likelihood function for model M1:
  vector loglik1_bp(vector status, matrix Z, matrix g, matrix G, real tau,
                    vector gamma, vector psi, vector phi){

    int q = cols(Z);
    int n = rows(Z);
    vector[n] lht0;
    vector[n] Ht0;
    //vector[n] Ft0;
    vector[n] St0;
    vector[n] lp_short;
    vector[n] lp_long;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] loglik;
    vector[n] theta;
    vector[n] logmix;
    vector[n] aux;

    lht0 = log(g*gamma) - log(tau);
    Ht0 = G*gamma;
    //Ft0 = -expm1(-Ht0);
    St0 = exp(-Ht0);
    lp_short = Z*psi;
    lp_long = Z*phi;
    theta = exp(lp_long);
    aux = lp_long - Ht0;

    for(i in 1:n){
      //logmix[i] = log_mix(Ft0[i], lp_short[i], lp_long[i]);
      logmix[i] = log_mix(St0[i], lp_long[i], lp_short[i]);
      log_ht[i] = lp_short[i] + lp_long[i] - logmix[i] + lht0[i];
      log_St[i] = -theta[i]*(logmix[i] - aux[i]);
      loglik[i] = status[i]*log_ht[i] + log_St[i];
    }
    return loglik;
  }


  //-------------------------------------------------------------
  // Likelihood function for model M2:
  vector loglik2_bp(vector status, matrix Z, matrix g, matrix G, real tau,
                    vector gamma, vector psi, vector phi){

    int q = cols(Z);
    int n = rows(Z);
    vector[n] Rt0;
    vector[n] lrt0;
    vector[n] lp_short;
    vector[n] theta_L;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] ratio;
    vector[n] aux;
    vector[n] loglik;

    lrt0 = log(g*gamma) - log(tau);
    Rt0 = G*gamma;
    lp_short = Z*psi;
    theta_L = exp( Z*phi);
    ratio = exp( Z*(psi-phi));

    for(i in 1:n){
      aux[i] = ratio[i]*Rt0[i];
      log_ht[i] = lp_short[i] - log1p(aux[i] ) + lrt0[i];
      log_St[i] = - theta_L[i]*log1p(aux[i]);
      loglik[i] = status[i]*log_ht[i] + log_St[i];
    }

    return loglik;
  }


  //-------------------------------------------------------------
  // Likelihood function for model M3:
  vector loglik3_bp(vector status, matrix Z, matrix X, matrix g, matrix G,  real tau,
                    vector gamma, vector psi, vector phi, vector beta){

    int q = cols(Z);
    int p = cols(X);
    int n = rows(Z);

    vector[n] lht0;
    vector[n] Ht0;
    //vector[n] Ft0;
    vector[n] St0;
    vector[n] lp_short;
    vector[n] lp_long;
    vector[n] lp_const;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] ratio;
    vector[n] loglik;
    vector[n] theta;
    vector[n] logmix;
    vector[n] aux;

    lht0 = log(g*gamma) - log(tau);
    Ht0 = G*gamma;
    //Ft0 = -expm1(-Ht0);
    St0 = exp(-Ht0);
    lp_short = Z*psi;
    lp_long = Z*phi;
    lp_const = X*beta;
    //theta = exp(lp_long);
    aux = lp_long - Ht0;

    for(i in 1:n){
      //logmix[i] = log_mix(Ft0[i], lp_short[i], lp_long[i]);
      logmix[i] = log_mix(St0[i], lp_long[i], lp_short[i]);
      log_ht[i] = lp_short[i] + lp_long[i] + lp_const[i] - logmix[i] + lht0[i];
      //log_St[i] = -theta[i]*exp(lp_const[i])*(logmix[i] - aux[i]);
      log_St[i] = -exp(lp_long[i] + lp_const[i])*(logmix[i] - aux[i]);
      loglik[i] = status[i]*log_ht[i] + log_St[i];
    }

    return loglik;
  }


  //-------------------------------------------------------------
  // Likelihood function for model M4:
  vector loglik4_bp(vector status, matrix Z, matrix X, matrix g, matrix G,  real tau,
                    vector gamma, vector psi, vector phi, vector beta){


    int q = cols(Z);
    int p = cols(X);
    int n = rows(Z);

    vector[n] Rt0;
    vector[n] lrt0;
    vector[n] lp_short;
    vector[n] lp_const;
    vector[n] theta_L;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] ratio;
    vector[n] aux;
     vector[n] loglik;

    lrt0 = log(g*gamma) - log(tau);
    Rt0 = G*gamma;
    lp_short = Z*psi;
    lp_const = X*beta;
    theta_L = exp( Z*phi);
    ratio = exp( Z*(psi-phi));

    for(i in 1:n){
      aux[i] = ratio[i]*Rt0[i];
      log_ht[i] = lp_short[i] + lp_const[i] - log1p(aux[i]) + lrt0[i];
      log_St[i] = - theta_L[i]*exp(lp_const[i])*log1p(aux[i]);
      loglik[i] = status[i]*log_ht[i] + log_St[i];
    }

    return loglik;
  }

}
