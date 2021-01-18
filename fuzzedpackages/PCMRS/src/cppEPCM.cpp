
// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <vector>
#include <algorithm>
#include <assert.h>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::plugins(openmp)]]

// square root of a matrix
mat matsqrt2(mat A){
  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, A);
  
  colvec d = (eigval+abs(eigval))*0.5;
  colvec d2 = sqrt(d);
  mat B = eigvec*diagmat(d2)*trans(eigvec);
  
  return B;
}

// response function for acat
vec responseFun(vec eta){
  int q = eta.n_rows;
  mat eta_help1 = ones(q+1)*trans(join_cols(zeros(1),eta));
  eta_help1 = eta_help1 % trimatl(ones(eta_help1.n_rows,eta_help1.n_cols));
  vec pi_vec = ones(q)/as_scalar(sum(prod(exp(eta_help1),1)));
  for(int k=1; k<q ;k++){
    pi_vec(k) =pi_vec(k-1)*exp(eta(k-1));
  }
  return pi_vec;
}

// create inverse sigma for acat
mat createSigmaInv(vec mu){
  mat Sigma = diagmat(mu) - mu * trans(mu);
  mat SigmaInv = inv(Sigma);
  return SigmaInv;
}

// create derivative matrix for acat
mat createD(vec mu){
  int q = mu.n_rows;
  mat mu_k = join_cols(mu,1-sum(mu,1));
  mat D2 = zeros(q,q) - diagmat(1/mu);
  
  if(q==2){
    D2(1,0) = 1/mu(1);
  }else{
    D2(span(1,q-1),span(0,q-2)) = D2(span(1,q-1),span(0,q-2)) + diagmat(1/mu(span(1,q-1)));
  }

  D2(span::all,q-1) = -ones(q,1)/(1-as_scalar(sum(mu)));
  D2(q-1,q-1) = as_scalar(-(1-sum(mu(span(0,q-2))))/((1-sum(mu))*mu(q-1)));
  
  mat D = inv(D2);
  return D;
}

  
  // [[Rcpp::export]]
  arma::vec scoreEPCM(arma::vec alpha,
           arma::vec Y,
           arma::mat X,
           int Q,
           int q,
           int n,
           int I,
           int pall,
           int pX,
           arma::mat GHprobs,
           arma::mat GHweights,
           arma::vec GHnodes,
           int scaled,
           double cores) { 
             
      // initialize score vector
      vec s = zeros(pall);

      // initialize cij matrix for all persons and all knots, 
      // will be needed for normalization per person afterwards
      // in contrast, for implementation of PCM this was calculated BEFORE looping through persons
      // should be faster this way
      mat cij_mat = zeros(n,Q*Q);
      
      // initialize matrix containing score contributions per person and per parameter
      // has to be kept because one has to use cij_mat before one can sum over persons
      mat help_mat = zeros(pall,n);

      // initialize design for one person, without random effects
      mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,3+2*pX));
      mat etai = Z*alpha;
      // current linear predictor without random effects, per item and category
      etai = reshape(etai,q,I);

      // create design matrix for both random effects
      vec resp_style_help = zeros(q);

      for(int r=0; r<q ;r++){
        resp_style_help(r) = double(q-1)/2 - r;
      }
      
      if(scaled==0){
        resp_style_help = sign(resp_style_help);
      }
      mat design_rnd = join_rows(ones(q), resp_style_help);

      // create sigma matrix from current parameters
      
      double co_var = alpha(pall-2*pX-2)*sqrt(alpha(pall-2*pX-1))*sqrt(alpha(pall-2*pX-3));
      mat sigma = zeros(2,2);        
      sigma(0,0) = alpha(pall-2*pX-3);
      sigma(1,0) = co_var;
      sigma(0,1) = co_var;
      sigma(1,1) = alpha(pall-2*pX-1);
      
      // square root of current sigma matrix
      mat sigma12 = trans(chol(sigma));


      mat etaBoth = zeros(n,2);
      if(pX>0){
        vec betaTheta = alpha(span(pall-2*pX,pall-pX-1));
        vec betaGamma = alpha(span(pall-pX,pall-1));
        vec etaTheta = X * betaTheta;
        vec etaGamma = X * betaGamma;
        etaBoth = join_rows(etaTheta,etaGamma);
      }

      // initialize number of different threads
      #ifdef _OPENMP
      if(cores > 1)
            omp_set_num_threads(cores);
      #endif
      
  
    // loop through persons
    #pragma omp parallel for
    for(int i=0; i<n ;i++){ 

      // initialize a running vector over all Q*Q knots
      int pos_knot = 0;
      
      int pos = i*q*I;
      // get response of person i
      vec yi = Y(span(pos,pos+q*I-1));

      mat Xihelp;
      if(pX>0){
      Xihelp = X(i,span::all);
      if(q>1){
        for(int f=1; f<q ;f++){ 
          Xihelp = join_cols(Xihelp,X(i,span::all));
      }
      }

      mat Xihelp2 = Xihelp;
      
        for(int f=0; f<pX ;f++){ 
          Xihelp2(span::all,f) = Xihelp2(span::all,f)%resp_style_help;
        }
      
      
      Xihelp = join_rows(Xihelp,Xihelp2);
      }

      // initialize design for person i
      mat Zij = Z;

      
      // loop through all knots, for both random effects
      for(int j=0; j<Q ;j++){  
        for(int jj=0; jj<Q ;jj++){
          // initialize derivatives, inverse sigma and mu for all items and current knots
          mat D_i = zeros(q*I,q*I);
          mat SigmaInv_i = zeros(q*I,q*I);
          vec mu_i = zeros(q*I);

          // initialize current random effects or knots
          vec rnd_act = zeros(2);
          rnd_act(0) = GHnodes(j);
          rnd_act(1) = GHnodes(jj);
          
          rnd_act = sigma12*rnd_act + trans(etaBoth(i,span::all));

          // initialize current cij value with weight and prob of current knot-combination
          double cij_help = GHweights(j,jj)*GHprobs(j,jj);
            
          // loop through all items  
          for(int k=0; k<I ;k++){  
            
            // response of person i for current item k
            vec yi_k = yi(span(k*q,q+k*q-1));
            yi_k = join_cols(yi_k,1-sum(yi_k,0));
            
            // create current eta and mu
            vec eta_k = etai(span::all,k) + design_rnd*rnd_act;
            vec mu_k = responseFun(eta_k);
            mu_i(span(k*q,q+k*q-1)) = mu_k;

            // create colmuns for random effects in design matrix, 
            // one column per parameter, therefore three columns
            mat design_act2 = zeros(q,3);
            
            design_act2(span::all,0) = (design_rnd(span::all,0)*GHnodes(j))/(2*sqrt(sigma(0,0)));
            design_act2(span::all,1) = GHnodes(j)*sqrt(sigma(1,1))*design_rnd(span::all,1) - (design_rnd(span::all,1)*GHnodes(jj)*sigma(1,1)*alpha(pall-2*pX-2))/sqrt(sigma(1,1)-sigma(1,1)*pow(alpha(pall-2*pX-2),2));
            design_act2(span::all,2) = design_rnd(span::all,1)*GHnodes(j)*alpha(pall-2*pX-2)/(2*sqrt(sigma(1,1))) + (GHnodes(jj)*(1-pow(alpha(pall-2*pX-2),2))*design_rnd(span::all,1))/(2*sqrt(sigma(1,1)-sigma(1,1)*pow(alpha(pall-2*pX-2),2)));
            
            
            // update the respective part of design matrix
            if(pX>0){
              Zij(span(k*q,q+k*q-1),span(pall-2*pX-3,pall-1)) = join_rows(design_act2,Xihelp);
            }else{
              Zij(span(k*q,q+k*q-1),span(pall-2*pX-3,pall-1)) = design_act2;
            }

  
            
            // update derivative and inverse sigma for current item           
            D_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createD(mu_k);
            SigmaInv_i(span(k*q,q+k*q-1),span(k*q,q+k*q-1)) = createSigmaInv(mu_k);

            // update current cij by probability of current response
            mu_k = join_cols(mu_k,1-sum(mu_k,0));

            cij_help = cij_help * prod(mu_k % yi_k - (yi_k-1));
          }
          // after looping through all items, update cij_mat
          cij_mat(i,pos_knot) = cij_help;
       
          // update all contribution of current knot-combination and person i to column of person i
          help_mat(span::all,i) = help_mat(span::all,i) + (trans(Zij)*D_i*SigmaInv_i*(yi-mu_i))*cij_help;
          
          pos_knot = pos_knot + 1;
        }}
    }
    
      // normalization value per person
      vec cij_norm = sum(cij_mat,1);
      
      // normalize row of each parameter by cij_norm
      for(int e=0; e<pall ;e++){  
        help_mat(e,span::all) = help_mat(e,span::all)%trans(1/cij_norm);
      }

      // sum score contributions over all persons per covariate
      s = -sum(help_mat,1);

      return s;
           }
           
 // [[Rcpp::export]]
  double loglikEPCM(arma::vec alpha,
           arma::vec Y,
           arma::mat X,
           int Q,
           int q,
           int n,
           int I,
           int pall,
           int pX,
           arma::mat GHprobs,
           arma::mat GHweights,
           arma::vec GHnodes,
           int scaled,
           int cores) { 
    
      
      // initialize loglikelihood       
      double f = 0;   

      // initialize design for one person, without random effects
      mat Z = -join_rows(diagmat(ones(q*I)),zeros(q*I,3+2*pX));
      mat etai = Z*alpha;
      // current linear predictor without random effects, per item and category
      etai = reshape(etai,q,I);

     
      // create design matrix for both random effects
      vec resp_style_help = zeros(q);

        for(int r=0; r<q ;r++){
          resp_style_help(r) = double(q-1)/2 - r;
        }

      if(scaled==0){
        resp_style_help = sign(resp_style_help);
      }
      mat design_rnd = join_rows(ones(q), resp_style_help);

      // create sigma matrix from current parameters
      double co_var = alpha(pall-2*pX-2)*sqrt(alpha(pall-2*pX-1))*sqrt(alpha(pall-2*pX-3));
      mat sigma = zeros(2,2);        
      sigma(0,0) = alpha(pall-2*pX-3);
      sigma(1,0) = co_var;
      sigma(0,1) = co_var;
      sigma(1,1) = alpha(pall-2*pX-1);

      // square root of current sigma matrix
      mat sigma12 = trans(chol(sigma));

      
      mat etaBoth = zeros(n,2);
      if(pX>0){
            vec betaTheta = alpha(span(pall-2*pX,pall-pX-1));
            vec betaGamma = alpha(span(pall-pX,pall-1));
            vec etaTheta = X * betaTheta;
            vec etaGamma = X * betaGamma;
            etaBoth = join_rows(etaTheta,etaGamma);
      }
      
      // initialize number of different threads
      #ifdef _OPENMP
        if(cores > 1)
          omp_set_num_threads(cores);
      #endif

      // loop through all persons
      #pragma omp parallel for reduction(-:f)
      for(int i=0; i<n ;i++){
      int pos = i*q*I;
        // get response of person i
        vec yi = Y(span(pos,pos+q*I-1));
 
        mat prods_i = ones(Q,Q);
          // loop through all knots, for both random effects
          for(int j=0; j<Q ;j++){ 
            for(int jj=0; jj<Q ;jj++){ 
              // initialize current random effects or knots
              vec rnd_act = zeros(2);
              rnd_act(0) = GHnodes(j);
              rnd_act(1) = GHnodes(jj);

              rnd_act = sigma12*rnd_act + trans(etaBoth(i,span::all));

              
              // loop through all items  
              for(int k=0; k<I ;k++){  
                // response of person i for current item k
                vec yi_k = yi(span(k*q,q+k*q-1));
                yi_k = join_cols(yi_k,1-sum(yi_k,0));
                
                // get eta and mu of person i for current item k
                vec eta_k = etai(span::all,k) + design_rnd*rnd_act;
                vec mu_k = responseFun(eta_k);
                mu_k = join_cols(mu_k,1-sum(mu_k,0));
                
                // create prob for current item, update prob matrix for respective knots
                vec help_pow = mu_k % yi_k - (yi_k-1);
                prods_i(j,jj) = prods_i(j,jj)*prod(help_pow);
              }
            }
          }
                      f  -= log(accu(prods_i%GHweights%GHprobs));
         // accumulate all likelihood contributions, weights by respective weights and probs of knots
      }

  return f;
}
