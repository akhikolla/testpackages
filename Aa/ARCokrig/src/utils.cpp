#include <Rmath.h>
#include <R_ext/Utils.h>   // interrupt the Gibbs sampler from R

#include <RcppArmadillo.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(cpp11)]]

// #include <math.h>

using namespace Rcpp;





// function definition

// void expcor(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       cormat(i,j) = 1.0;
//       for(k=0; k<dist.n_slices; k++){
//         cormat(i,j) *= exp(-dist(i,j,k)*beta(k)); 
//       }
//     }
//   }
  
// }


void expcor_deriv(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat, int l){
  // derivative wrt inverse range parameter 

  int i=0, j=0, k=0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        if(k==l){
          cormat(i,j) *= exp(-beta(k)*dist(i,j,k)) * (-dist(i,j,k)); 
        }else{
          cormat(i,j) *= exp(-dist(i,j,k)*beta(k));
        }
      }
    }
  }

}


// void expcor_nonsep(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double dtemp = 0.0;
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       dtemp = 0.0;
//       for(k=0; k<dist.n_slices; k++){
//         dtemp += (dist(i,j,k)*beta(k))*(dist(i,j,k)*beta(k)); 
//       }
//       cormat(i,j) = exp(-sqrt(dtemp));
//     }
//   }
  
// }


// void matern_3_2_cor(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double constant, sqrt_3=sqrt(3.0);
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       cormat(i,j) = 1.0;
//       for(k=0; k<dist.n_slices; k++){
//         constant = sqrt_3*dist(i,j,k)*beta(k);
//         cormat(i,j) *= (1.0 + constant)*exp(-constant); 
//       }

//     }
//   }

// }


void matern_3_2_cor_deriv(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat, int l)
{

  // derivative wrt inverse range parameter 
  // Eigen::VectorXd beta = phi.array().inverse();

  int i=0, j=0, k=0;
  double constant, sqrt_3=sqrt(3.0);
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        constant = sqrt_3*dist(i,j,k);
        if(k==l){
          cormat(i,j) *= exp(-constant*beta(k))*(constant - 
                        (1.0+constant*beta(k))*constant);
        }else{
          cormat(i,j) *= (1.0 + constant*beta(k))*exp(-constant*beta(k));
        }
         
      }

    }
  }

}

// void matern_3_2_cor_nonsep(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double constant, sqrt_3=sqrt(3.0);
//   double dtemp=0.0;
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       dtemp = 0.0;
//       for(k=0; k<dist.n_slices; k++){
//         dtemp += (dist(i,j,k)*beta(k))*(dist(i,j,k)*beta(k));
//       }
//       constant = sqrt_3*sqrt(dtemp);
//       cormat(i,j) = (1.0 + constant)*exp(-constant); 
//     }
//   }

// }

// void matern_5_2_cor(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double constant, sqrt_5 = sqrt(5.0);
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       cormat(i,j) = 1.0;
//       for(k=0; k<dist.n_slices; k++){
//         constant = sqrt_5*dist(i,j,k)*beta(k);
//         cormat(i,j) *= (1.0 + constant + constant*constant/3.0)*exp(-constant); 
//       }

//     }
//   }

// }


void matern_5_2_cor_deriv(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat, int l)
{

  // derivative wrt inverse range parameter 
  // Eigen::VectorXd beta = phi.array().inverse();

  int i=0, j=0, k=0;
  double constant, sqrt_5 = sqrt(5.0), delta;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        constant = sqrt_5*dist(i,j,k)*beta(k);
        if(k==l){
          delta = constant + constant*constant/3.0;
          cormat(i,j) *= exp(-constant)*(delta/beta(k) + constant*constant/(beta(k)*3.0)
                         -(1.0+delta)*sqrt_5*dist(i,j,k)); 
        }else{
          cormat(i,j) *= (1.0 + constant + constant*constant/3.0)*exp(-constant); 
        }
      }

    }
  }

}


// void matern_5_2_cor_nonsep(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double constant, sqrt_5 = sqrt(5.0);
//   double dtemp = 0.0;
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       dtemp = 0.0;
//       for(k=0; k<dist.n_slices; k++){
//         dtemp += (dist(i,j,k)*beta(k))*(dist(i,j,k)*beta(k));
//       }
//       constant = sqrt_5*sqrt(dtemp);
//       cormat(i,j) = (1.0 + constant + constant*constant/3.0)*exp(-constant); 
//     }
//   }

// }


// void sqexpcor(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
// {
  
//   int i=0, j=0, k=0;
//   double temp = 0.0;
//   for(i=0; i<dist.n_rows; i++){
//     for(j=0; j<dist.n_cols; j++){
//       temp = 0.0;
//       for(k=0; k<dist.n_slices; k++){
//         temp += dist(i,j,k)*dist(i,j,k)*(beta(k)*beta(k)); 
//       }
//       cormat(i,j) = exp(-temp);
//     }
//   }



// }


void sqexpcor_deriv(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat, int l){
  // derivative wrt inverse range parameter 

  int i=0, j=0, k=0;
  double alpha = 2.0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        if(k==l){
          cormat(i,j) *= exp(-beta(k)*dist(i,j,k)) * (-pow(dist(i,j,k), alpha)) * alpha * pow(beta(k), alpha-1.0); 
        }else{
          cormat(i,j) *= exp(-pow(dist(i,j,k)*beta(k), alpha));
        }
      }
    }
  }

}




void powercor_deriv(Eigen::Ref<Eigen::VectorXd> beta, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat, int l){
  // derivative wrt inverse range parameter 

  int i=0, j=0, k=0;
  double alpha = 1.9;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        if(k==l){
          cormat(i,j) *= exp(-beta(k)*dist(i,j,k)) * (-pow(dist(i,j,k), alpha)) * alpha * pow(beta(k), alpha-1.0); 
        }else{
          cormat(i,j) *= exp(-pow(dist(i,j,k)*beta(k), alpha));
        }
      }
    }
  }

}







// [[Rcpp::export]]
Eigen::MatrixXd buildcov_deriv(Eigen::VectorXd& beta, arma::cube& dist, int l,
       const String& covmodel, const bool& nugget){

  Eigen::MatrixXd cormat_deriv(dist.n_rows, dist.n_cols);


  if(covmodel=="exp"){
    expcor_deriv(beta, dist, cormat_deriv, l-1);
  }else if(covmodel=="matern_3_2"){
    matern_3_2_cor_deriv(beta, dist, cormat_deriv, l-1);
  }else if(covmodel=="matern_5_2"){
    matern_5_2_cor_deriv(beta, dist, cormat_deriv, l-1);
  }else if(covmodel=="Gaussian"){
    sqexpcor_deriv(beta, dist, cormat_deriv, l-1);
  }else if(covmodel=="powexp"){
    powercor_deriv(beta, dist, cormat_deriv, l-1);
  }else{
    Rcpp::Rcout<<"\nNot Implemented yet!\n";
  }



  return cormat_deriv;
}





double log_reference_prior(Eigen::VectorXd& beta, arma::cube& dist, Eigen::MatrixXd& RInv,
        Eigen::MatrixXd& X, const String& covmodel, const bool& nugget){

  int n = dist.n_rows;
  int p = beta.size();
  int q = X.cols();

  // compute Q 
  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  ldlt.compute(X.transpose()*RInv*X);
  Eigen::MatrixXd Q = RInv * (Eigen::MatrixXd::Identity(n, n) - X*ldlt.solve(X.transpose())*RInv);

  Eigen::MatrixXd Rl_deriv(n, n);

  List W(p);

  if(nugget){ // if nugget exists
    if(covmodel=="exp"){
      for(int l=0; l<p-1; l++){
        expcor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="matern_3_2"){
      for(int l=0; l<p-1; l++){
        matern_3_2_cor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="matern_5_2"){
      for(int l=0; l<p-1; l++){
        matern_5_2_cor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="Gaussian"){
      for(int l=0; l<p-1; l++){
        sqexpcor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="powexp"){
      for(int l=0; l<p-1; l++){
        powercor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else{
      Rcpp::Rcout<<"\nThe correlation function is not implemented yet!\n";
    }

    W[p-1] = Q;

  }else{ // no nugget

    if(covmodel=="exp"){
      for(int l=0; l<p; l++){
        expcor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="matern_3_2"){
      for(int l=0; l<p; l++){
        matern_3_2_cor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="matern_5_2"){
      for(int l=0; l<p; l++){
        matern_5_2_cor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="Gaussian"){
      for(int l=0; l<p; l++){
        sqexpcor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else if(covmodel=="powexp"){
      for(int l=0; l<p; l++){
        powercor_deriv(beta, dist, Rl_deriv, l);
        W[l] = Rl_deriv * Q;
      }
    }else{
      Rcpp::Rcout<<"\nThe correlation function is not implemented yet!\n";
    }

  }



  Eigen::MatrixXd FisherIR(p+1, p+1);
  Eigen::MatrixXd W_l(n, n), W_k(n, n);
  FisherIR(0,0) = n - q;
  for(int l=0; l<p; l++){
    W_l = as<Eigen::MatrixXd>(W[l]);
    FisherIR(0, l+1) = W_l.trace();
    FisherIR(l+1, 0) = W_l.trace();

    for(int k=0; k<p; k++){
      W_k = as<Eigen::MatrixXd>(W[k]);
      FisherIR(l+1, k+1) = (W_l*W_k).trace();
      FisherIR(k+1, l+1) = (W_l*W_k).trace();
    }
  }

  //Eigen::LDLT<Eigen::MatrixXd> ldlt;
  ldlt.compute(FisherIR);

  double detIR = ldlt.vectorD().array().log().sum();

  return 0.5*detIR;

}





Eigen::VectorXd log_Jeffreys_prior(Eigen::VectorXd& beta, arma::cube& dist, Eigen::MatrixXd& RInv,
        Eigen::MatrixXd& X, const String& covmodel, const bool& nugget){

  int n = dist.n_rows;
  int p = beta.size();
  // int q = X.cols();
  Eigen::MatrixXd Rl_deriv(n, n);

  List U(p);

  if(nugget){ // if nugget exists
    if(covmodel=="exp"){
      for(int l=0; l<p-1; l++){
        expcor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="matern_3_2"){
      for(int l=0; l<p-1; l++){
        matern_3_2_cor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="matern_5_2"){
      for(int l=0; l<p-1; l++){
        matern_5_2_cor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="Gaussian"){
      for(int l=0; l<p-1; l++){
        sqexpcor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="powexp"){
      for(int l=0; l<p-1; l++){
        powercor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else{
      Rcpp::Rcout<<"\nThe correlation function is not implemented yet!\n";
    }

    U[p-1] = RInv;

  }else{ // no nugget

    if(covmodel=="exp"){
      for(int l=0; l<p; l++){
        expcor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="matern_3_2"){
      for(int l=0; l<p; l++){
        matern_3_2_cor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="matern_5_2"){
      for(int l=0; l<p; l++){
        matern_5_2_cor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="Gaussian"){
      for(int l=0; l<p; l++){
        sqexpcor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else if(covmodel=="powexp"){
      for(int l=0; l<p; l++){
        powercor_deriv(beta, dist, Rl_deriv, l);
        U[l] = Rl_deriv * RInv;
      }
    }else{
      Rcpp::Rcout<<"\nThe correlation function is not implemented yet!\n";
    }

  }



  Eigen::MatrixXd FisherIJ(p+1, p+1);
  Eigen::MatrixXd U_l(n, n), U_k(n, n);
  FisherIJ(0,0) = n;
  for(int l=0; l<p; l++){
    U_l = as<Eigen::MatrixXd>(U[l]);
    FisherIJ(0, l+1) = U_l.trace();
    FisherIJ(l+1, 0) = U_l.trace();

    for(int k=0; k<p; k++){
      U_k = as<Eigen::MatrixXd>(U[k]);
      FisherIJ(l+1, k+1) = (U_l*U_k).trace();
      FisherIJ(k+1, l+1) = (U_l*U_k).trace();
    }
  }

  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  ldlt.compute(FisherIJ);

  double detIJ = ldlt.vectorD().array().log().sum();

  Eigen::VectorXd log_IJ(2);

  log_IJ(0) = 0.5*detIJ;
  
  
  ldlt.compute(X.transpose()*RInv*X);
  log_IJ(1) = log_IJ(0) + 0.5*ldlt.vectorD().array().log().sum();

  return log_IJ;

}


// [[Rcpp::export]]
double log_objective_prior(Eigen::VectorXd& beta, arma::cube& dist, Eigen::MatrixXd& RInv,
        Eigen::MatrixXd& X, const String& covmodel, const bool& nugget, const String& prior){

  double log_prior=0.0;
  if(prior=="Reference"){
    log_prior = log_reference_prior(beta, dist, RInv, X, covmodel, nugget);
  }else if(prior=="Jeffreys"){
    Eigen::VectorXd tmp = log_Jeffreys_prior(beta, dist, RInv, X, covmodel, nugget);
    log_prior = tmp(0);
  }else if(prior=="Ind_Jeffreys"){
    Eigen::VectorXd tmp = log_Jeffreys_prior(beta, dist, RInv, X, covmodel, nugget);
    log_prior = tmp(1);
  }

  return log_prior;
}







// function definition
void expcor(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        cormat(i,j) *= exp(-dist(i,j,k)/phi(k)); 
      }
    }
  }
  

}


void expcor_nonsep(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double dtemp = 0.0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      dtemp = 0.0;
      for(k=0; k<dist.n_slices; k++){
        dtemp += (dist(i,j,k)/phi(k))*(dist(i,j,k)/phi(k)); 
      }
      cormat(i,j) = exp(-sqrt(dtemp));
    }
  }
  
}


void matern_3_2_cor(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double constant, sqrt_3=sqrt(3.0);
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        constant = sqrt_3*dist(i,j,k)/phi(k);
        cormat(i,j) *= (1.0 + constant)*exp(-constant); 
      }

    }
  }

}


void matern_3_2_cor_nonsep(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double constant, sqrt_3=sqrt(3.0);
  double dtemp=0.0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      dtemp = 0.0;
      for(k=0; k<dist.n_slices; k++){
        dtemp += (dist(i,j,k)/phi(k))*(dist(i,j,k)/phi(k));
      }
      constant = sqrt_3*sqrt(dtemp);
      cormat(i,j) = (1.0 + constant)*exp(-constant); 
    }
  }

}

void matern_5_2_cor(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double constant, sqrt_5 = sqrt(5.0);
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      cormat(i,j) = 1.0;
      for(k=0; k<dist.n_slices; k++){
        constant = sqrt_5*dist(i,j,k)/phi(k);
        cormat(i,j) *= (1.0 + constant + constant*constant/3.0)*exp(-constant); 
      }

    }
  }

}


void matern_5_2_cor_nonsep(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double constant, sqrt_5 = sqrt(5.0);
  double dtemp = 0.0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      dtemp = 0.0;
      for(k=0; k<dist.n_slices; k++){
        dtemp += (dist(i,j,k)/phi(k))*(dist(i,j,k)/phi(k));
      }
      constant = sqrt_5*sqrt(dtemp);
      cormat(i,j) = (1.0 + constant + constant*constant/3.0)*exp(-constant); 
    }
  }

}


void sqexpcor(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double temp = 0.0;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      temp = 0.0;
      for(k=0; k<dist.n_slices; k++){
        temp += dist(i,j,k)*dist(i,j,k)/(phi(k)*phi(k)); 
      }
      cormat(i,j) = exp(-temp);
    }
  }


}


void powercor(Eigen::Ref<Eigen::VectorXd> phi, arma::cube& dist, Eigen::Ref<Eigen::MatrixXd> cormat)
{
  
  int i=0, j=0, k=0;
  double temp = 0.0;
  double alpha = 1.9;
  for(i=0; i<dist.n_rows; i++){
    for(j=0; j<dist.n_cols; j++){
      temp = 0.0;
      for(k=0; k<dist.n_slices; k++){
        temp += pow(dist(i,j,k)/phi(k), alpha); 
      }
      cormat(i,j) = exp(-temp);
    }
  }


}


// [[Rcpp::export]]
Eigen::MatrixXd buildcov(Eigen::VectorXd& phi, arma::cube& dist,  
              const String& covmodel, const bool& nugget){

  Eigen::MatrixXd cormat(dist.n_rows, dist.n_cols);
  if(covmodel=="exp"){
    expcor(phi, dist, cormat);
  }else if(covmodel=="matern_3_2"){
    matern_3_2_cor(phi, dist, cormat);
  }else if(covmodel=="matern_5_2"){
    matern_5_2_cor(phi, dist, cormat);
  }else if(covmodel=="Gaussian"){
    sqexpcor(phi, dist, cormat);
  }else if(covmodel=="powexp"){
    powercor(phi, dist, cormat);
  }else if(covmodel=="aniso_exp"){
    expcor_nonsep(phi, dist, cormat);
  }else if(covmodel=="aniso_matern_3_2"){
    matern_3_2_cor_nonsep(phi, dist, cormat);
  }else if(covmodel=="aniso_matern_5_2"){
    matern_5_2_cor_nonsep(phi, dist, cormat);
  }else{
    Rcpp::Rcout<<"\nThe correlation function is not implemented yet!\n";
  }

  if(nugget){
    // cormat.diagonal().array() += phi[dist.n_slices];
    double dtemp = 0.0;
    for(int i=0; i<dist.n_rows; i++){
      for(int j=0; j<dist.n_cols; j++){
        dtemp = 0.0;
        for(int k=0; k<dist.n_slices; k++){
          dtemp += dist(i,j,k);
        }
        if(dtemp==0){
          cormat(i,j) += phi[dist.n_slices];
        }
      }
    }
  }



  return cormat;

}




// [[Rcpp::export]]
arma::cube compute_distance(arma::mat& input1, arma::mat& input2){
  
  int m1 = input1.n_rows;
  int m2 = input2.n_rows;
  int Dim_x = input1.n_cols;

  arma::cube dist = arma::zeros(m1, m2, Dim_x);

  // #pragma omp parallel for collapse(2) 
  for(int k=0; k<Dim_x; k++){
    // #pragma omp parallel for collapse(2) 
    for(int i=0; i<m1; i++){
      for(int j=0; j<m2; j++){
        dist(i,j,k) = abs(input1(i,k) - input2(j,k));
      }
    }
  }


  return dist;
}




// This routine generate samples from multivariate t distribution.
// [[Rcpp::export]]
arma::cube sample_mvt(const arma::mat& mu,
          const arma::mat& L, 
          const arma::vec& sigma, double df,
          int nsample 
                 ){


  int N = mu.n_cols;
  int n = mu.n_rows;
  // Eigen::LDLT<Eigen::MatrixXd> ldlt;
  // ldlt.compute(c_star);
  // Eigen::MatrixXd L = ldlt.matrixL();
  // arma::mat L = chol(c_star);

  arma::cube y(nsample, n, N, arma::fill::zeros);

#ifdef USE_R
    GetRNGstate();
#endif

    for(int k=0; k<nsample; k++){

    for(int i=0; i<N; i++){
      y(arma::span(k,k), arma::span(0, n-1), arma::span(i,i)) = mu.col(i) + 
             sqrt(sigma(i)) * L *  as<arma::vec>(Rcpp::rnorm(n, 0.0, 1.0)) 
             / sqrt(as<double>(Rcpp::rchisq(1, df)) / df);
      // y.col(i) = mu.col(i) + sqrt(sigma(i)) * L * as<Eigen::VectorXd>(Rcpp::rt(n, df));
    }
  }

#ifdef USE_R
  PutRNGstate();
#endif


  return y;
}






// [[Rcpp::export]]
double compute_S(const Eigen::Map<Eigen::MatrixXd>& output,
                 const Eigen::Map<Eigen::MatrixXd>& Q) {
  
  int k_outputs = output.cols();
  Eigen::VectorXd S(k_outputs);

  // omp_set_num_threads(2);
  // #pragma omp parallel for 
  for(int i=0; i<k_outputs; i++){
    S(i) = output.col(i).transpose()*Q*output.col(i);
  }
  
  return S.array().log().sum();
  
}







// [[Rcpp::export]]
Eigen::VectorXd compute_Svec(const Eigen::Map<Eigen::MatrixXd>& output,
                 const Eigen::Map<Eigen::MatrixXd>& Q) {
  
  int k_outputs = output.cols();
  Eigen::VectorXd S(k_outputs);
  
  //omp_set_num_threads(2);
  //#pragma omp parallel for 
  for(int i=0; i<k_outputs; i++){
    S(i) = output.col(i).transpose()*Q*output.col(i);
  }
  
  return S;
  
}



// [[Rcpp::export]]
double compute_S_sum(const Eigen::Map<Eigen::MatrixXd>& y_t,
                     const Eigen::Map<Eigen::MatrixXd>& H_t,
                     const Eigen::Map<Eigen::MatrixXd>& y_t1,
                     const Eigen::Map<Eigen::MatrixXd>& RInv,
                     const Eigen::Map<Eigen::MatrixXd>& K
                     ) {
   
  int k_outputs = y_t.cols();
  int n_input = y_t.rows();
  int q = H_t.cols();


  Eigen::VectorXd S(k_outputs);
  Eigen::VectorXd S_quad(k_outputs);

  Eigen::MatrixXd Q = RInv - K;

  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  Eigen::MatrixXd HRH(q,q), HRInv(q, n_input), HRInvW(q,1), XRX(q,q);
  Eigen::VectorXd logdetX(k_outputs);

  HRInv = H_t.transpose() * RInv;
  HRH = HRInv * H_t;

  double WRW, ytQyt, ytKyt1, yKOKy, ytRyt1, yROKy, yRORy, S_sum;

  // omp_set_num_threads(2);
  // #pragma omp parallel for 
  for(int i=0; i<k_outputs; i++){
    WRW = y_t1.col(i).transpose() * RInv * y_t1.col(i);
    HRInvW = HRInv * y_t1.col(i);

    ldlt.compute(HRH*WRW - HRInvW*HRInvW.transpose());
    logdetX(i) = ldlt.vectorD().array().log().sum();

    // compute t(y_t1)%*%QH%*%y_t1
    S(i) = y_t1.col(i).transpose() * Q * y_t1.col(i);

    // compute S^2
    ytQyt = y_t.col(i).transpose() * Q * y_t.col(i);
    // ytKW =  y_t.col(i).transpose() * K * y_t1.col(i);
    ytKyt1 = y_t.col(i).transpose() * K * y_t1.col(i);
    yKOKy = ytKyt1 * ytKyt1 / S(i);
    ytRyt1 = y_t.col(i).transpose()*RInv*y_t1.col(i);
    yROKy = ytRyt1 * ytKyt1 / S(i);
    yRORy = ytRyt1 * ytRyt1 / S(i);

    S_quad(i) = ytQyt - yKOKy + 2.0*yROKy - yRORy;

  }

  double S_2_log = S_quad.array().log().sum();
  S_sum = logdetX.sum() + (n_input-q-1)*S_2_log;
  
  return S_sum;
  
}




// This routine computes the predictive mean nad predictive variance 
// when the design is hierarchically nested in multivariate autogressive
// cokriging models.
// [[Rcpp::export]]
List compute_prediction(const Eigen::Map<Eigen::MatrixXd>& y_t,
                     const Eigen::Map<Eigen::MatrixXd>& Ht,
                     const Eigen::Map<Eigen::MatrixXd>& y_t1,
                     const Eigen::Map<Eigen::MatrixXd>& yhat_t1,
                     const Eigen::Map<Eigen::MatrixXd>& vhat_t1,
                     const Eigen::Map<Eigen::MatrixXd>& RInv,
                     const Eigen::Map<Eigen::MatrixXd>& Hnew,
                     const Eigen::Map<Eigen::MatrixXd>& Wnew_t1,
                     const Eigen::Map<Eigen::MatrixXd>& Rmo,
                     const Eigen::Map<Eigen::MatrixXd>& R_sk
                     ) {
   
  int k_outputs = y_t.cols(); // number of spatial locations
  int n = y_t.rows(); // number of model runs
  int q = Ht.cols(); // number of fixed basis functions
  int m = Hnew.rows(); 

  Eigen::MatrixXd A(q,q), RY(n,k_outputs), HY(q, k_outputs), 
                  K(n,n), Q(n,n), RmoR(m,n), RH(n,q);

  // compute intermediate terms
  Eigen::LDLT<Eigen::MatrixXd> ldlt;

  RH = RInv * Ht;

  ldlt.compute(Ht.transpose()*RH);
  A = ldlt.solve(Eigen::MatrixXd::Identity(q,q));
  RY = RInv*y_t;
  // AHR = A*RH.transpose();
  HY = Ht.transpose()*RY; // q-by-k matrix
  K = RH * A * RH.transpose();
  Q = RInv - K;
  RmoR = Rmo*RInv;


  Eigen::VectorXd sigma(k_outputs);
  // create temp variables
  // double omega; 
  Eigen::VectorXd omega(k_outputs);
  Eigen::MatrixXd X(n, q+1), kappa(m, k_outputs), temp(m, q+1), Fnew(m,q+1);
  double ytQyt, ytKyt1, yKOKy, ytRyt1, yROKy, yRORy, constant;
  constant = (n-q-1.0) / (n-q-3.0);

  X.block(0,0,n,q) = Ht;
  Fnew.block(0,0,m,q) = Hnew;


  Eigen::MatrixXd b = Eigen::MatrixXd::Identity(q+1, k_outputs);

  b.block(0,0,q,k_outputs) = A*HY;


  Eigen::MatrixXd eye(q+1,q+1), XRXInv(q+1,q+1);
  eye.setIdentity();

  Eigen::VectorXd var_betahat(m);

  for(int i=0; i<k_outputs; i++){
    // compute 1 / t(W_t1)%*%QH%*%W_t1
    omega(i) = 1.0 / (y_t1.col(i).dot(Q * y_t1.col(i)));

    // compute S^2
    ytQyt = y_t.col(i).transpose() * (Q * y_t.col(i));
    ytKyt1 = y_t.col(i).transpose() * (K * y_t1.col(i));
    yKOKy = ytKyt1 * ytKyt1 * omega(i);
    ytRyt1 = y_t.col(i).transpose()*(RInv*y_t1.col(i));
    yROKy = ytRyt1 * ytKyt1 * omega(i);
    yRORy = ytRyt1 * ytRyt1 * omega(i);

    // sigma is k_output by 1 vector
    sigma(i) = (ytQyt - yKOKy + 2.0*yROKy - yRORy)/(n-q-1); 

    // compute XRXInv
    // X.block(0,q,n,1) = y_t1;
    X.col(q) = y_t1.col(i);

    // XRX = X.transpose() * RInv * X;
    ldlt.compute(X.transpose() * RInv * X);
    XRXInv = ldlt.solve(eye);  // for trace calcuation later

    // compute beta 
    b.col(i) = XRXInv * (X.transpose() * (RY.col(i)));


    // Fnew.block(0,q,m,1) = Wnew_t1.col(i);
    Fnew.col(q) = Wnew_t1.col(i);
    temp = Fnew - RmoR*X; // m-by-(q+1)
    var_betahat = (temp*XRXInv*temp.transpose()).diagonal();
    kappa.col(i) = R_sk + var_betahat + 
                   vhat_t1.col(i)*omega(i);

  }


  // define outputs
  Eigen::MatrixXd krige=Eigen::MatrixXd::Zero(m, k_outputs);
  Eigen::MatrixXd krige_var=Eigen::MatrixXd::Zero(m, k_outputs);

  // get kriging mean
  Eigen::VectorXd gamma = b.row(q); //b.block(q,0,1,k_outputs);
  krige = (Hnew - RmoR*Ht) * b.block(0,0,q,k_outputs)  
          + (yhat_t1 - RmoR*y_t1) * gamma.asDiagonal() + RmoR * y_t;

  // // get kriging variance 
  krige_var = vhat_t1 * (gamma.array()*gamma.array()).matrix().asDiagonal() + 
              constant * kappa * sigma.asDiagonal();
                     
  return List::create(_["krige"] = krige,
                      _["krige.var"] = krige_var
                      );
  
  // return List::create(_["b"] = b,
  //                     _["sigma"] = sigma,
  //                     _["kappa"] = kappa
  //                     );
}





// conditional simulation for multivaraite autoregressive cokriging models
// [[Rcpp::export]]
Eigen::MatrixXd conditional_simulation(const Eigen::Map<Eigen::MatrixXd>& y_t,
                     const Eigen::Map<Eigen::MatrixXd>& Ht,
                     const Eigen::Map<Eigen::MatrixXd>& y_t1,
                     const Eigen::Map<Eigen::MatrixXd>& RInv,
                     const Eigen::Map<Eigen::MatrixXd>& Hnew,
                     const Eigen::Map<Eigen::MatrixXd>& Wnew_t1,
                     const Eigen::Map<Eigen::MatrixXd>& Rmo,
                     const Eigen::Map<Eigen::MatrixXd>& R_sk
                     ) {
   
   
  int k_outputs = y_t.cols(); // number of spatial locations
  int n = y_t.rows(); // number of model runs 60
  int q = Ht.cols(); // number of fixed basis functions  
  int m = Hnew.rows(); 

  Eigen::MatrixXd A(q,q), RY(n,k_outputs), AHR(q, n), HY(q, k_outputs), 
                  K(n,n), Q(n,n), RmoR(m,n);

  // compute intermediate terms
  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  ldlt.compute(Ht.transpose()*RInv*Ht);
  A = ldlt.solve(Eigen::MatrixXd::Identity(q,q));
  RY = RInv*y_t;
  AHR = A*Ht.transpose()*RInv;
  HY = Ht.transpose()*RY; // q-by-k matrix
  K = RInv*Ht*AHR;
  Q = RInv - K;
  RmoR = Rmo*RInv;
  // RmoRY = RmoR*y_t;

  Eigen::VectorXd sigma(k_outputs);
  // create temp variables
  double omega; 
  Eigen::MatrixXd X(n, q+1), c_star(m, m), temp(m, q+1), Fnew(m,q+1);
  double ytQyt, ytKyt1, yKOKy, ytRyt1, yROKy, yRORy;
  // constant = (n-q-1) / (n-q-3);

  Eigen::MatrixXd Sig_ymy(m, m), L(m,m), ym(m, k_outputs);
  Eigen::MatrixXd mu_ymy(m, k_outputs);
  // Eigen::VectorXd mu_ymy(m);

  X.setZero();
  Fnew.setZero();
  c_star.setZero();
  Sig_ymy.setZero();
  temp.setZero();
  mu_ymy.setZero();
  ym.setZero();


  X.block(0,0,n,q) = Ht;
  Fnew.block(0,0,m,q) = Hnew;


  Eigen::MatrixXd b = Eigen::MatrixXd::Identity(q+1, k_outputs);

  b.block(0,0,q,k_outputs) = A*HY;


  Eigen::MatrixXd eye(q+1,q+1), XRXInv(q+1,q+1);
  eye.setIdentity();

#ifdef USE_R
    GetRNGstate();
#endif

  for(int i=0; i<k_outputs; i++){
    // compute t(W_t1)%*%QH%*%W_t1
    omega = 1.0 / (y_t1.col(i).transpose() * Q * y_t1.col(i));

    // compute S^2
    ytQyt = y_t.col(i).transpose() * Q * y_t.col(i);
    ytKyt1 = y_t.col(i).transpose() * K * y_t1.col(i);
    yKOKy = ytKyt1 * ytKyt1 * omega;
    ytRyt1 = y_t.col(i).transpose()*RInv*y_t1.col(i);
    yROKy = ytRyt1 * ytKyt1 * omega;
    yRORy = ytRyt1 * ytRyt1 * omega;

    // sigma is k_output by 1 vector
    sigma(i) = (ytQyt - yKOKy + 2.0*yROKy - yRORy)/(n-q-1); 

    // compute XRXInv
    // X.block(0,q,n,1) = y_t1.col(i);
    X.col(q) = y_t1.col(i);

    // XRX = X.transpose() * RInv * X;
    ldlt.compute(X.transpose() * RInv * X);
    XRXInv = ldlt.solve(eye);  

    // compute beta 
    b.col(i) = XRXInv * (X.transpose() * (RY.col(i))); // it is the same as in R

    // compute conditional mean
    Fnew.col(q) = Wnew_t1.col(i);
    mu_ymy.col(i) = Fnew*b.col(i) + RmoR*(y_t.col(i)-X*b.col(i));


    // compute conditional variance
    // Fnew.block(0,q,m,1) = Wnew_t1.col(i);
    temp = Fnew - RmoR*X; // m-by-(q+1)
    c_star = R_sk + temp*XRXInv*temp.transpose();

    Sig_ymy = c_star;
    ldlt.compute(Sig_ymy);
    L = ldlt.matrixL();

    // generate random numbers from a multivariate t distribution

    ym.col(i) = mu_ymy.col(i) + sqrt(sigma(i)) * L * as<Eigen::VectorXd>(Rcpp::rnorm(m, 0.0, 1.0)) 
               / sqrt(as<double>(Rcpp::rchisq(1, n-q-1)) / (n-q-1));

    // ym.col(i) = mu_ymy.col(i) + L * as<Eigen::VectorXd>(Rcpp::rt(m, n-q-1.0));

    //interface with R
    // #ifdef USE_R
    // if (i%20==0) R_CheckUserInterrupt();
    // #endif

  }

#ifdef USE_R
  PutRNGstate();
#endif

                
  return ym;
  
  
  // return List::create(_["ym"] = ym,
  //                     _["sigma"] = sigma
  //                     );
}



// [[Rcpp::export]]
List compute_param(const Eigen::Map<Eigen::MatrixXd>& y_t,
                     const Eigen::Map<Eigen::MatrixXd>& Ht,
                     const Eigen::Map<Eigen::MatrixXd>& y_t1,
                     const Eigen::Map<Eigen::MatrixXd>& RInv
                     ) {
   
  int k_outputs = y_t.cols(); // number of spatial locations
  int n = y_t.rows(); // number of model runs 60
  int q = Ht.cols(); // number of fixed basis functions  

  Eigen::MatrixXd A(q,q), RY(n,k_outputs), AHR(q, n), HY(q, k_outputs), 
                  K(n,n), Q(n,n);

  // compute intermediate terms
  Eigen::LDLT<Eigen::MatrixXd> ldlt;
  ldlt.compute(Ht.transpose()*RInv*Ht);
  A = ldlt.solve(Eigen::MatrixXd::Identity(q,q));
  RY = RInv*y_t;
  AHR = A*Ht.transpose()*RInv;
  HY = Ht.transpose()*RY; // q-by-k matrix
  K = RInv*Ht*AHR;
  Q = RInv - K;


  Eigen::VectorXd sigma(k_outputs);
  // create temp variables
  double omega; 
  Eigen::MatrixXd X(n, q+1);
  double ytQyt, ytKyt1, yKOKy, ytRyt1, yROKy, yRORy;
  // constant = (n-q-1) / (n-q-3);


  X.block(0,0,n,q) = Ht;


  Eigen::MatrixXd b = Eigen::MatrixXd::Identity(q+1, k_outputs);

  b.block(0,0,q,k_outputs) = A*HY;

  Eigen::MatrixXd eye(q+1,q+1), XRXInv(q+1,q+1);
  eye.setIdentity();

  for(int i=0; i<k_outputs; i++){
    // compute t(W_t1)%*%QH%*%W_t1
    omega = 1.0 / (y_t1.col(i).transpose() * Q * y_t1.col(i));

    // compute S^2
    ytQyt = y_t.col(i).transpose() * Q * y_t.col(i);
    ytKyt1 = y_t.col(i).transpose() * K * y_t1.col(i);
    yKOKy = ytKyt1 * ytKyt1 * omega;
    ytRyt1 = y_t.col(i).transpose()*RInv*y_t1.col(i);
    yROKy = ytRyt1 * ytKyt1 * omega;
    yRORy = ytRyt1 * ytRyt1 * omega;

    // sigma is k_output by 1 vector
    sigma(i) = (ytQyt - yKOKy + 2.0*yROKy - yRORy)/(n-q-1.0+2.0); 

    // compute XRXInv
    // X.block(0,q,n,1) = y_t1.col(i);
    X.col(q) = y_t1.col(i);

    // XRX = X.transpose() * RInv * X;
    ldlt.compute(X.transpose() * RInv * X);
    XRXInv = ldlt.solve(eye);  

    // compute beta 
    b.col(i) = XRXInv * (X.transpose() * (RY.col(i))); // it is the same as in R


  }


                
  
  return List::create(_["beta"] = b,
                      _["sigma2"] = sigma
                      );
}


