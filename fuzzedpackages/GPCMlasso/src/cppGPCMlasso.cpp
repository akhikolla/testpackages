// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <math.h>    
#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

/////////////////////////////////////////////////////////
// Definition of some help functions
/////////////////////////////////////////////////////////


// response function for acat
vec responseFun(vec eta){
  int q = eta.n_rows;
  eta(find(eta>10)) = ones(size(find(eta>10)))*10;
  eta(find(eta<-10)) = -ones(size(find(eta<-10)))*10;
  mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
  eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
  vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
  for(int k=1; k<q ;k++){
    pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
  }
  pi_vec = (pi_vec-0.5)*0.99999+0.5;
  return pi_vec;
}

// create inverse sigma for acat
mat createSigmaInv(vec mu){
  mat Sigma = diagmat(mu) - mu * trans(mu);
  mat SigmaInv;
    try{
    SigmaInv = inv(Sigma);
    }
    catch(...){
     SigmaInv = pinv(Sigma);
    }
  return SigmaInv;
}

// create derivative matrix for acat
mat createD(vec mu){
  int q = mu.n_rows;

  mat D2 = zeros(q,q) - diagmat(1/mu);

  if(q==2){
    D2(1,0) = 1/mu(1);
  }else{
    D2(span(1,q-1),span(0,q-2)) = D2(span(1,q-1),span(0,q-2)) + diagmat(1/mu(span(1,q-1)));
  }

  D2(span::all,q-1) = -ones(q,1)/(1-as_scalar(sum(mu)));
  D2(q-1,q-1) = as_scalar(-(1-sum(mu(span(0,q-2))))/((1-sum(mu))*mu(q-1)));

  mat D;
  try{
    D = inv(D2);
  }
  catch(...){
    D = pinv(D2);
  }

  return D;
}

// inverse of L1 penalty
mat L1(mat xi, double cvalue){
  mat ret = 1/sqrt(xi%xi+cvalue);
  return ret;
}

// general function to create derivative of penalties
mat A(mat beta,
      mat acoefs,
      double lambda,
      vec weight,
      double cvalue){
  vec nbs = trans(acoefs)*beta;
  uvec fd2 = find(nbs!=0);
  mat fd = zeros(nbs.n_elem,1);
  fd.rows(fd2) = ones(fd2.n_elem);
  
  mat appro = fd%weight%L1(nbs,cvalue)*lambda;
  mat Aret = trans(acoefs)%repmat(sqrt(appro),1,acoefs.n_rows);
  
  return(trans(Aret)*Aret);
}

/////////////////////////////////////////////////////////
// Definition of PCMlasso functions
/////////////////////////////////////////////////////////

// [[Rcpp::export]]
double loglikPCMlasso(arma::vec alpha,
                arma::vec Y,
                arma::mat X,
                arma::mat Z,
                int Q,
                arma::vec q,
                int n,
                int I,
                int px,
                arma::vec GHweights,
                arma::vec GHnodes,
                arma::mat acoefs,
                double lambda,
                double lambda2,
                double cvalue,
                int cores,
                arma::vec weight,
                int n_sigma,
                double scale_fac
){

  // initalize loglikelihood
  vec f = zeros(n);

  // current value of L1 penalty
  double P1 = accu(weight%(1/(L1(trans(acoefs)*alpha,cvalue))))*lambda;

  int sumq = accu(q);
  
  // get standard deviations and matrix of different sigmas
  vec sigma = alpha(span(px-n_sigma,px-1));
  
  // create long vector of different sigmas;
  vec sigma_long = ones(sumq)*sigma(0);
  if(n_sigma>1){
    int pos_sigma = q(0);
    for(double k=1; k<I ;k++){
      sigma_long(span(pos_sigma,pos_sigma+q(k)-1)) =  ones(q(k))*sigma(k);
      pos_sigma = pos_sigma + q(k);
    }
  }

  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda2;

    // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));
  
  
  // initialize number of different threads
#ifdef _OPENMP
  // if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  vec yi; mat Xi; mat eta_i2; vec eta_i;
  mat Xij; vec eta_ij; 
  vec mu_i; vec yi_k; vec eta_ijk; vec mu_k; vec help_pow;
  vec prods_i; int pos; int pos_person; 
  double j; double k; int qk; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, qk, j, k,  pos, yi, Xi, eta_i2, eta_i, Xij, eta_ij,  mu_i, yi_k, eta_ijk, mu_k, help_pow, prods_i, pos_person) shared(f)
#endif
    for(i=0; i<n ;i++){
    pos = i*sumq;
    
    // initialize response for person i
    yi = Y(span(pos,pos+sumq-1));
    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+sumq-1),span::all));
    }
    Xi = Xi%repmat(sigma_long,1,Xi.n_cols);

    // current values of eta without random effect
    eta_i = Xi*alpha;

    // initialize likelihood contributions of person i for all Q knots
    prods_i = ones(Q);

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      eta_ij = eta_i+sigma_long*GHnodes(j);

      pos_person = 0;
      // iterate through items
      for(k=0; k<I ;k++){
        qk = q(k);
        yi_k = yi(span(pos_person,qk+pos_person-1));
        yi_k = join_cols(yi_k,1-sum(yi_k,0));
        
        // mu and eta for item k, knot j, person i
        eta_ijk = eta_ij(span(pos_person,qk+pos_person-1));
        mu_k = responseFun(eta_ijk);
        
        // update current cij by probability of current response
        mu_k = join_cols(mu_k,1-sum(mu_k,0));
        
        // likelihood of response for person i for item k
        help_pow = mu_k % yi_k - (yi_k-1);
        prods_i(j) = prods_i(j)*prod(help_pow);
        pos_person = pos_person + qk;
      }
    }

    // update likelihood by weighted contributions of person i
    f(i) = -log(sum(prods_i%GHweights));
  }
    double fx = sum(f)/scale_fac + P1 + P2;
  return fx;
}

// [[Rcpp::export]]
arma::vec scorePCMlasso(arma::vec alpha,
                        arma::vec Y,
                        arma::mat X,
                        arma::mat Z,
                        int Q,
                        arma::vec q,
                        int n,
                        int I,
                        int px,
                        arma::vec GHweights,
                        arma::vec GHnodes,
                        arma::mat acoefs,
                        double lambda,
                        double lambda2,
                        double cvalue,
                        int cores,
                        arma::vec weight,
                        int n_sigma,
                        double scale_fac
                        ) {
  

  // initialize score vector
  vec s = zeros(px);

  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  mat cij_mat = zeros(n,Q);

  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(px,n);

    // current deriv of L1 penalty
  vec P1 = A(alpha,
             acoefs,
             lambda,
             weight,
             cvalue)*alpha;
  
  // length of response vector per person
  int sumq = accu(q);

  // extract discrimination parameters(s)
  vec sigma = alpha(span(px-n_sigma,px-1));
  
  // create long vector of different sigmas; also create design matrix for sigma
  vec sigma_long = ones(sumq)*sigma(0);
  mat design_sigma;
  if(n_sigma==1){
    design_sigma = ones(sumq,1);
  }else{
    design_sigma = join_cols(ones(q(0)),zeros(sumq-q(0)));
    int pos_sigma = q(0);
    for(double k=1; k<I ;k++){
      sigma_long(span(pos_sigma,pos_sigma+q(k)-1)) =  ones(q(k))*sigma(k);
      design_sigma = join_rows(design_sigma,join_cols(join_cols(zeros(pos_sigma),ones(q(k))),zeros(sumq-pos_sigma-q(k))));
      pos_sigma = pos_sigma + q(k);
    }
  }

  // current deriv of L2 penalty
  vec P2 = 2*alpha*lambda2;

  // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));

  // initialize number of different threads
#ifdef _OPENMP
  // if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  
  vec yi; mat Xi; mat eta_i2; vec eta_i;
  mat Xij; vec eta_ij; mat D_i; mat SigmaInv_i; 
  vec mu_i; vec yi_k; vec eta_ijk; vec mu_k; vec help_pow;
  vec prods_i; double cij_help; int pos; int pos_person; 
  double help_me; int j; int k; int qk; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, qk, j, k, help_me, pos, yi, Xi, eta_i2, eta_i, Xij, eta_ij, D_i, SigmaInv_i, mu_i, yi_k, eta_ijk, mu_k, help_pow, prods_i, cij_help, pos_person) shared(help_mat, cij_mat)
#endif  
  for(i=0; i<n ;i++){
    pos = i*sumq;
    
    // extract response for person i
    yi = Y(span(pos,pos+sumq-1));
    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+sumq-1),span::all));
    }
    eta_i2 = repmat(Xi*alpha,1,n_sigma);
    Xi = Xi%repmat(sigma_long,1,Xi.n_cols);
 
    // current values of eta without random effect
    eta_i = Xi*alpha;

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      Xij = join_rows(Xi,design_sigma*GHnodes(j)+design_sigma%eta_i2);

      eta_ij = eta_i+sigma_long*GHnodes(j);

      // derivative and inverse sigma
      D_i = zeros(sumq,sumq);
      SigmaInv_i = zeros(sumq,sumq);
      // initalize mu for person i and knot j
      mu_i = zeros(sumq);
      
      // initialize cij weights
      cij_help = GHweights(j);

      pos_person = 0;
      // iterate through items
      for(k=0; k<I ;k++){
        qk = q(k);
        yi_k = yi(span(pos_person,qk+pos_person-1));
        yi_k = join_cols(yi_k,1-sum(yi_k,0));

        // mu and eta for item k, knot j, person i
        eta_ijk = eta_ij(span(pos_person,qk+pos_person-1));
        mu_k = responseFun(eta_ijk);

        // D and SigmaInv for item k, knot j, person i
        if(qk==1){
          mu_i(pos_person) = mu_k(0);
          help_me = exp(eta_ijk(0))/pow(1+exp(eta_ijk(0)),2);
          D_i(pos_person,pos_person)= -help_me;
          SigmaInv_i(pos_person,pos_person) = 1/help_me;
        }else{
          mu_i(span(pos_person,qk+pos_person-1)) = mu_k;
          D_i(span(pos_person,qk+pos_person-1),span(pos_person,qk+pos_person-1)) = createD(mu_k);
          SigmaInv_i(span(pos_person,qk+pos_person-1),span(pos_person,qk+pos_person-1)) = createSigmaInv(mu_k);
        }

        // update current cij by probability of current response
        mu_k = join_cols(mu_k,1-sum(mu_k,0));
        
        cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        pos_person = pos_person + qk;
      }

      // after looping through all items, update cij_mat
      cij_mat(i,j) = cij_help;

      // update all contribution of current knot-combination and person i to column of person i
      help_mat(span::all,i) = help_mat(span::all,i) + (trans(Xij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;

    }
  }
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  // normalize row of each parameter by cij_norm
  help_mat = help_mat%(ones(help_mat.n_rows,1)*trans(1/cij_norm));

  // sum score contributions over all persons per covariate
  s = -sum(help_mat,1)/scale_fac + P1 + P2;

    return s;
}

//////////////////////////////



// [[Rcpp::export]]
double loglikDIFlasso(arma::vec alpha,
                      arma::vec Y,
                      arma::mat X,
                      arma::mat Z,
                      int Q,
                      arma::vec q,
                      int n,
                      int I,
                      int px,
                      arma::vec GHweights,
                      arma::vec GHnodes,
                      arma::mat acoefs,
                      double lambda,
                      double lambda2,
                      double cvalue,
                      int cores,
                      arma::vec weight,
                      int n_sigma,
                      double scale_fac
){

  // initalize loglikelihood
  vec f = zeros(n);

  // current value of L1 penalty
  double P1 = accu(weight%(1/(L1(trans(acoefs)*alpha,cvalue))))*lambda;

  // get standard deviations and matrix of different sigmas
  vec sigma = alpha(span(px-n_sigma,px-1));
  vec sigma_mat;
  if(n_sigma==1){
    sigma_mat = repmat(sigma,I,1);
  }else{
    sigma_mat = sigma;
  }

  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda2;

    // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));

  // initialize number of different threads
#ifdef _OPENMP
//   if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  
  vec yi; mat Xi; vec eta_i; vec prods_i;
  vec eta_ij; vec mu_ij; 
  int pos; int j; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, j, pos, yi, Xi, eta_i, eta_ij, mu_ij,prods_i) shared(f)
#endif  
  for(i=0; i<n ;i++){
    pos = i*I;

    // initialize response for person i
    yi = Y(span(pos,pos+I-1));
    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+I-1),span::all));
    }
    Xi = Xi%repmat(sigma_mat,1,Xi.n_cols);
    // current values of eta without random effect
    eta_i = Xi*alpha;

    // initialize likelihood contributions of person i for all Q knots
    prods_i = ones(Q);

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      eta_ij = eta_i+sigma_mat*GHnodes(j);
      mu_ij = exp(eta_ij)/(1+exp(eta_ij));
      prods_i(j) = prod(mu_ij%yi+(1-mu_ij)%(1-yi));
    }
    // update likelihood by weighted contributions of person i
    f(i) = -log(sum(prods_i%GHweights));
  }

  double fx = sum(f)/scale_fac + P1 + P2;
  return fx;
}


// [[Rcpp::export]]
arma::vec scoreDIFlasso(arma::vec alpha,
                        arma::vec Y,
                        arma::mat X,
                        arma::mat Z,
                        int Q,
                        arma::vec q,
                        int n,
                        int I,
                        int px,
                        arma::vec GHweights,
                        arma::vec GHnodes,
                        arma::mat acoefs,
                        double lambda,
                        double lambda2,
                        double cvalue,
                        int cores,
                        arma::vec weight,
                        int n_sigma,
                        double scale_fac) {

  // initialize score vector
  vec s = zeros(px);

  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  mat cij_mat = zeros(n,Q);

  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(px,n);

  // current deriv of L1 penalty
  vec P1 = A(alpha,
             acoefs,
             lambda,
             weight,
             cvalue)*alpha;

  // get standard deviations and matrix of different sigmas
  vec sigma = alpha(span(px-n_sigma,px-1));
  vec sigma_mat;
  if(n_sigma==1){
    sigma_mat = repmat(sigma,I,1);
  }else{
    sigma_mat = sigma;
  }

  // current deriv of L2 penalty
  vec P2 = 2*alpha*lambda2;

  // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));

  // create design matrix for sigma
  mat design_sigma;
  if(n_sigma==1){
    design_sigma = ones(I,1);
  }else{
    design_sigma = diagmat(ones(I));
  }

  // initialize number of different threads
#ifdef _OPENMP
//   if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  
  vec yi; mat Xi; mat eta_i2; vec eta_i; vec prods_i;
  mat Xij; vec eta_ij; vec mu_ij; 
  double cij_help; int pos; int j; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, j, pos, yi, Xi, eta_i2, eta_i, Xij, eta_ij, mu_ij, prods_i, cij_help) shared(help_mat, cij_mat)
#endif
  for(i=0; i<n ;i++){
    pos = i*I;

    // extract response for person i
    yi = Y(span(pos,pos+I-1));

    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+I-1),span::all));
    }
    eta_i2 = Xi*alpha;
    Xi = Xi%repmat(sigma_mat,1,Xi.n_cols);

    // current values of eta without random effect
    eta_i = Xi*alpha;

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      Xij = join_rows(Xi,design_sigma*GHnodes(j)+design_sigma%repmat(eta_i2,1,n_sigma));
      eta_ij = eta_i+sigma_mat*GHnodes(j);

      // get mu
      mu_ij = exp(eta_ij)/(1+exp(eta_ij));

      // initialize cij weights
      cij_help = GHweights(j)*prod(mu_ij%yi+(1-mu_ij)%(1-yi));

      // after looping through all items, update cij_mat
      cij_mat(i,j) = cij_help;

      // update all contribution of current knot-combination and person i to column of person i
      help_mat(span::all,i) = help_mat(span::all,i) + (trans(Xij)*(yi-mu_ij))*cij_help;

    }
  }
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  // normalize row of each parameter by cij_norm
  help_mat = help_mat%(ones(help_mat.n_rows,1)*trans(1/cij_norm));

  // sum score contributions over all persons per covariate
  s = -sum(help_mat,1)/scale_fac + P1 + P2;

  return s;
}




// [[Rcpp::export]]
Rcpp::List loglikscorePCMlasso(arma::vec alpha,
                      arma::vec Y,
                      arma::mat X,
                      arma::mat Z,
                      int Q,
                      arma::vec q,
                      int n,
                      int I,
                      int px,
                      arma::vec GHweights,
                      arma::vec GHnodes,
                      arma::mat acoefs,
                      double lambda,
                      double lambda2,
                      double cvalue,
                      int cores,
                      arma::vec weight,
                      int n_sigma,
                      double scale_fac
){

  // initalize loglikelihood
  vec f = zeros(n);

  // current value of L1 penalty
  double P1 = accu(weight%(1/(L1(trans(acoefs)*alpha,cvalue))))*lambda;

  // current deriv of L1 penalty
  vec P1der = A(alpha,
             acoefs,
             lambda,
             weight,
             cvalue)*alpha;

  // length of response vector per person
  int sumq = accu(q);

  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  mat cij_mat = zeros(n,Q);

  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(px,n);
  
  // extract discrimination parameters(s)
  vec sigma = alpha(span(px-n_sigma,px-1));
  
  // create long vector of different sigmas; also create design matrix for sigma
  vec sigma_long = ones(sumq)*sigma(0);
  mat design_sigma;
  if(n_sigma==1){
    design_sigma = ones(sumq,1);
  }else{
    design_sigma = join_cols(ones(q(0)),zeros(sumq-q(0)));
    int pos_sigma = q(0);
    for(double k=1; k<I ;k++){
      sigma_long(span(pos_sigma,pos_sigma+q(k)-1)) =  ones(q(k))*sigma(k);
      design_sigma = join_rows(design_sigma,join_cols(join_cols(zeros(pos_sigma),ones(q(k))),zeros(sumq-pos_sigma-q(k))));
      pos_sigma = pos_sigma + q(k);
    }
  }

    // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda2;

  // current deriv of L2 penalty
  vec P2der = 2*alpha*lambda2;

  // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));

  // initialize number of different threads
#ifdef _OPENMP
//   if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  vec yi; mat Xi; mat eta_i2; vec eta_i;
  mat Xij; vec eta_ij; mat D_i; mat SigmaInv_i; 
  vec mu_i; vec yi_k; vec eta_ijk; vec mu_k; vec help_pow;
  vec prods_i; double cij_help; int pos; int pos_person; 
  double help_me; int j; int k; int qk; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, qk, j, k, help_me, pos, yi, Xi, eta_i2, eta_i, Xij, eta_ij, D_i, SigmaInv_i, mu_i, yi_k, eta_ijk, mu_k, help_pow, prods_i, cij_help, pos_person) shared(help_mat, cij_mat, f)
#endif
  for(i=0; i<n ;i++){
    pos = i*sumq;
    
    // initialize response for person i
    yi = Y(span(pos,pos+sumq-1));

    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+sumq-1),span::all));
    }

    eta_i2 = repmat(Xi*alpha,1,n_sigma);
    Xi = Xi%repmat(sigma_long,1,Xi.n_cols);

    // current values of eta without random effect
    eta_i = Xi*alpha;


    // initialize likelihood contributions of person i for all Q knots
    prods_i = ones(Q);

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      Xij = join_rows(Xi,design_sigma*GHnodes(j)+design_sigma%eta_i2);
      eta_ij = eta_i+sigma_long*GHnodes(j);

      // initialize derivative and inverse sigma
      D_i = zeros(sumq,sumq);
      SigmaInv_i = zeros(sumq,sumq);

      // initalize mu for person i and knot j
      mu_i = zeros(sumq);

      // initialize cij weights
      cij_help = GHweights(j);
      
      pos_person = 0;
      // iterate through items
      for(k=0; k<I ;k++){
        qk = q(k);
        yi_k = yi(span(pos_person,qk+pos_person-1));
        yi_k = join_cols(yi_k,1-sum(yi_k,0));

        // mu and eta for item k, knot j, person i
        eta_ijk = eta_ij(span(pos_person,qk+pos_person-1));
        vec mu_k = responseFun(eta_ijk);

        // D and SigmaInv for item k, knot j, person i
        if(qk==1){
          mu_i(pos_person) = mu_k(0);
          help_me = exp(eta_ijk(0))/pow(1+exp(eta_ijk(0)),2);
          D_i(pos_person,pos_person)= -help_me;
          SigmaInv_i(pos_person,pos_person) = 1/help_me;
        }else{
          mu_i(span(pos_person,qk+pos_person-1)) = mu_k;
          D_i(span(pos_person,qk+pos_person-1),span(pos_person,qk+pos_person-1)) = createD(mu_k);
          SigmaInv_i(span(pos_person,qk+pos_person-1),span(pos_person,qk+pos_person-1)) = createSigmaInv(mu_k);
        }

        // update mu with prob of last category
        mu_k = join_cols(mu_k,1-sum(mu_k,0));

        // likelihood of response for person i for item k
        help_pow = mu_k % yi_k - (yi_k-1);
        prods_i(j) = prods_i(j)*prod(help_pow);

        // update current cij by probability of current response
        cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
        pos_person = pos_person + qk;
      }
      // after looping through all items, update cij_mat
      cij_mat(i,j) = cij_help;

      // update all contribution of current knot-combination and person i to column of person i
      help_mat(span::all,i) = help_mat(span::all,i) +
        (trans(Xij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
    }

    // update likelihood by weighted contributions of person i
    f(i) = -log(sum(prods_i%GHweights));
  }

  
  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  // normalize row of each parameter by cij_norm
  help_mat = help_mat%(ones(help_mat.n_rows,1)*trans(1/cij_norm));

  // sum score contributions over all persons per covariate
  vec s = -sum(help_mat,1)/scale_fac + P1der + P2der;

  double fx = sum(f)/scale_fac + P1 + P2;

  return List::create(Named("objective") = fx,
                      Named("gradient") = s);
}


// [[Rcpp::export]]
Rcpp::List loglikscoreDIFlasso(arma::vec alpha,
                      arma::vec Y,
                      arma::mat X,
                      arma::mat Z,
                      int Q,
                      arma::vec q,
                      int n,
                      int I,
                      int px,
                      arma::vec GHweights,
                      arma::vec GHnodes,
                      arma::mat acoefs,
                      double lambda,
                      double lambda2,
                      double cvalue,
                      int cores,
                      arma::vec weight,
                      int n_sigma,
                      double scale_fac
){

  // initalize loglikelihood
  vec f = zeros(n);
  
  // initialize score vector
  vec s = zeros(px);

  // initialize cij matrix for all persons and all knots,
  // will be needed for normalization per person afterwards
  mat cij_mat = zeros(n,Q);

  // initialize matrix containing score contributions per person and per parameter
  // has to be kept because one has to use cij_mat before one can sum over persons
  mat help_mat = zeros(px,n);

  // current value of L1 penalty
  double P1 = accu(weight%(1/(L1(trans(acoefs)*alpha,cvalue))))*lambda;

  // current deriv of L1 penalty
  vec P1der = A(alpha,
                acoefs,
                lambda,
                weight,
                cvalue)*alpha;


  // get standard deviations and matrix of different sigmas
  vec sigma = alpha(span(px-n_sigma,px-1));
  vec sigma_mat;
  if(n_sigma==1){
    sigma_mat = repmat(sigma,I,1);
  }else{
    sigma_mat = sigma;
  }

  // current value of L2 penalty
  double P2 = accu(alpha%alpha)*lambda2;

  // current deriv of L2 penalty
  vec P2der = 2*alpha*lambda2;

  // shorten alpha
  alpha = alpha(span(0,px-n_sigma-1));

  // create design matrix for sigma
  mat design_sigma;
  if(n_sigma==1){
    design_sigma = ones(I,1);
  }else{
    design_sigma = diagmat(ones(I));
  }

  // initialize number of different threads
#ifdef _OPENMP
//   if(cores > 1)
    omp_set_num_threads(cores);
#endif
  
  
  vec yi; mat Xi; mat eta_i2; vec eta_i; vec prods_i;
  mat Xij; vec eta_ij; vec mu_ij; 
  double cij_help; int pos; int j; int i;
  
  // iterate through persons
#ifdef _OPENMP
#pragma omp parallel for private(i, j, pos, yi, Xi, eta_i2, eta_i, Xij, eta_ij, mu_ij, prods_i, cij_help) shared(help_mat, cij_mat, f)
#endif
  for(i=0; i<n ;i++){
    pos = i*I;

    // initialize response for person i
    yi = Y(span(pos,pos+I-1));
    // design for person i
    Xi = X;
    if(Z.n_rows>1){
      Xi = join_rows(X, Z(span(pos,pos+I-1),span::all));
    }
    eta_i2 = Xi*alpha;
    Xi = Xi%repmat(sigma_mat,1,Xi.n_cols);
    // current values of eta without random effect
    eta_i = Xi*alpha;

    // initialize likelihood contributions of person i for all Q knots
    prods_i = ones(Q);

    // iterate through knots
    for(j=0; j<Q ;j++){
      // create linear predictor for person i and knot j
      Xij = join_rows(Xi,design_sigma*GHnodes(j)+design_sigma%repmat(eta_i2,1,n_sigma));
      eta_ij = eta_i+sigma_mat*GHnodes(j);

      // get mu
      mu_ij = exp(eta_ij)/(1+exp(eta_ij));
      prods_i(j) = prod(mu_ij%yi+(1-mu_ij)%(1-yi));

    // initialize cij weights
    cij_help = GHweights(j)*prod(mu_ij%yi+(1-mu_ij)%(1-yi));

    // after looping through all items, update cij_mat
    cij_mat(i,j) = cij_help;

    // update all contribution of current knot-combination and person i to column of person i
    help_mat(span::all,i) = help_mat(span::all,i) + (trans(Xij)*(yi-mu_ij))*cij_help;

    }
    // update likelihood by weighted contributions of person i
    f(i) = - log(sum(prods_i%GHweights));
  }

  // normalization value per person
  vec cij_norm = sum(cij_mat,1);
  // normalize row of each parameter by cij_norm
  help_mat = help_mat%(ones(help_mat.n_rows,1)*trans(1/cij_norm));

  // sum score contributions over all persons per covariate
  s = -sum(help_mat,1)/scale_fac + P1der + P2der;


  double fx = sum(f)/scale_fac + P1 + P2;

  return List::create(Named("objective") = fx,
                      Named("gradient") = s);
}


