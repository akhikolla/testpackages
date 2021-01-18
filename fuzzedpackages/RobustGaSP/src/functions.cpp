
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath>
#include "ctools.h"
//#include "RobustGaSP.h"

using namespace Rcpp;
using namespace std;
using namespace Eigen;


//#define PI 3.14159265359

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//


///////kernel functions
// [[Rcpp::export]]
Eigen::MatrixXd matern_5_2_funct (const MapMat &d, double beta_i){
  //inline static Mat matern_5_2_funct (const Eigen::Map<Eigen::MatrixXd> & d, double beta_i){
  const double cnst = sqrt(5.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result +
	   result.array().pow(2.0).matrix()/3.0).cwiseProduct((-result).array().exp().matrix()));
  
}
// [[Rcpp::export]]
 Eigen::MatrixXd matern_3_2_funct (const Eigen::Map<Eigen::MatrixXd> & d, double beta_i){
  const double cnst = sqrt(3.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result ).cwiseProduct((-result).array().exp().matrix()));
  
}


// [[Rcpp::export]]
 Eigen::MatrixXd pow_exp_funct (const MapMat &d, double beta_i,double alpha_i){
  
  return (-(beta_i*d).array().pow(alpha_i)).exp().matrix();

}

// [[Rcpp::export]]
Eigen::MatrixXd periodic_gauss_funct(const MapMat &d, double beta_i){
  
  int Rnrow = d.rows();
  int Rncol = d.cols();
  
  Eigen::MatrixXd R=1.0/(2.0*sqrt(PI*beta_i))*Eigen::MatrixXd::Ones(Rnrow,Rncol);
    double two_beta_i=2*beta_i;
    int n_ti=std::min(std::max(11.0, two_beta_i),101.0);
    for(int ti=1; ti <n_ti; ti++){
      R=R+1.0/sqrt(PI*beta_i)*exp(-pow(ti,2.0)/(4.0*beta_i))* (ti*d).array().cos().matrix();
  }
  R=R/R(0,0);
  return R;
  
  


}


Eigen::MatrixXd periodic_gauss_funct_fixed_normalized_const(const MapMat &d, double beta_i,double perid_const_i){
  int Rnrow = d.rows();
  int Rncol = d.cols();
  
  Eigen::MatrixXd R=1.0/(2.0*sqrt(PI*beta_i))*Eigen::MatrixXd::Ones(Rnrow,Rncol);
  double two_beta_i=2*beta_i;
  int n_ti=std::min(std::max(11.0, two_beta_i),101.0);
  for(int ti=1; ti <n_ti; ti++){
    R=R+1.0/sqrt(PI*beta_i)*exp(-pow(ti,2.0)/(4.0*beta_i))* (ti*d).array().cos().matrix();
  }
  R=R/perid_const_i;
  return R;
  
}

// [[Rcpp::export]]
Eigen::MatrixXd periodic_exp_funct(const MapMat &d, double beta_i){
  
  int Rnrow = d.rows();
  int Rncol = d.cols();
  
  Eigen::MatrixXd R=1.0/(PI*beta_i)*Eigen::MatrixXd::Ones(Rnrow,Rncol);

  double five_beta_i=5*beta_i;
  int n_ti=std::min(std::max(21.0, five_beta_i),201.0);
  for(int ti=1; ti <n_ti; ti++){
    R=R+2.0*beta_i/((pow(beta_i,2.0)+pow(ti,2.0))*PI)*(ti*d).array().cos().matrix();
  }
  R=R/R(0,0);
  return R;
  

}




Eigen::MatrixXd periodic_exp_funct_fixed_normalized_const(const MapMat &d, double beta_i,double perid_const_i){
  int Rnrow = d.rows();
  int Rncol = d.cols();
  
  Eigen::MatrixXd R=1.0/(PI*beta_i)*Eigen::MatrixXd::Ones(Rnrow,Rncol);
  
  double five_beta_i=5*beta_i;
  int n_ti=std::min(std::max(21.0, five_beta_i),201.0);
  for(int ti=1; ti <n_ti; ti++){
    R=R+2.0*beta_i/((pow(beta_i,2.0)+pow(ti,2.0))*PI)*(ti*d).array().cos().matrix();
  }
  R=R/perid_const_i;
  return R;
  
}
       



//////derivative of kernel functions with regard to inverse range parameter
// [[Rcpp::export]]
 Eigen::MatrixXd  matern_5_2_deriv(const MapMat & R0_i,  const Mat R, double beta_i){
   
  const double sqrt_5 = sqrt(5.0);

  MatrixXd matOnes = Eigen::MatrixXd::Ones(R.rows(),R.cols());
  MatrixXd R0_i_2=R0_i.array().pow(2.0).matrix();

  MatrixXd part1= sqrt_5*R0_i+10.0/3*beta_i*R0_i_2;
  MatrixXd part2=matOnes+sqrt_5*beta_i*R0_i+5.0*pow(beta_i,2.0)*R0_i_2/3.0 ;
  return ((part1.cwiseQuotient(part2)  -sqrt_5*R0_i).cwiseProduct(R));
}

//inline static Eigen::MatrixXd

// [[Rcpp::export]]
 Eigen::MatrixXd  matern_3_2_deriv(const Eigen::Map<Eigen::MatrixXd> & R0_i,  const Eigen::MatrixXd R, double beta_i){
   
  const double sqrt_3 = sqrt(3.0);

  return(-sqrt(3)*R0_i.cwiseProduct(R)+sqrt_3*R0_i.cwiseProduct((-sqrt_3*beta_i*R0_i).array().exp().matrix()));
    
}


// [[Rcpp::export]]
Eigen::MatrixXd pow_exp_deriv(const MapMat & R0_i,  const Eigen::MatrixXd R, const double beta_i, const double alpha_i){
 return  -(R.array()*(R0_i.array().pow(alpha_i))).matrix()*alpha_i*pow(beta_i,alpha_i-1);
}

// [[Rcpp::export]]
Eigen::MatrixXd periodic_gauss_deriv(const MapMat &R0_i, const Eigen::MatrixXd & R, double beta_i){
  
  int Rnrow = R0_i.rows();
  int Rncol = R0_i.cols();
  
  Eigen::MatrixXd R_here=1.0/(2.0*sqrt(PI*beta_i))*Eigen::MatrixXd::Ones(Rnrow,Rncol);
  
  Eigen::MatrixXd R_partial=-pow(beta_i,-1.5)/(4*sqrt(PI))*R.Ones(Rnrow,Rncol);
  double two_beta_i=2*beta_i;
  int n_ti=std::min(std::max(11.0, two_beta_i),101.0);
  for(int ti=1; ti <n_ti; ti++){
    R_here=R_here+1.0/sqrt(PI*beta_i)*exp(-pow(ti,2.0)/(4.0*beta_i))* (ti*R0_i).array().cos().matrix();
    R_partial=R_partial+ exp(-pow(ti,2.0)/(4.0*beta_i) )*pow(beta_i,-1.5)/(2*sqrt(PI))*(pow(ti,2.0)/(beta_i*2.0)-1 )*(ti*R0_i).array().cos().matrix();
  }
  
  double c_norm=R_here(0,0);
  double c_norm_partial=R_partial(0,0);
  
  return (R.array()*( R_partial.array()*c_norm/R_here.array()- (c_norm_partial*R.Ones(Rnrow,Rncol)).array() )/c_norm).matrix();
  
  
}

// [[Rcpp::export]]
Eigen::MatrixXd periodic_exp_deriv(const MapMat &R0_i, const Eigen::MatrixXd & R, double beta_i){
  
  int Rnrow = R0_i.rows();
  int Rncol = R0_i.cols();
  
  Eigen::MatrixXd R_here=1.0/(PI*beta_i)*Eigen::MatrixXd::Ones(Rnrow,Rncol);
  
  Eigen::MatrixXd R_partial=-1.0/(PI*pow(beta_i,2.0) )*R.Ones(Rnrow,Rncol);
  double five_beta_i=5*beta_i;
  int n_ti=std::min(std::max(21.0, five_beta_i),201.0);
  for(int ti=1; ti <n_ti; ti++){
    R_here=R_here+2.0*beta_i/((pow(beta_i,2.0)+pow(ti,2.0))*PI)*(ti*R0_i).array().cos().matrix();
    R_partial=R_partial+ 2.0*(pow(ti,2.0)-pow(beta_i,2.0) )/(PI*(pow(beta_i,2.0)+pow(ti,2.0) ))*(ti*R0_i).array().cos().matrix();
  }
  
  double c_norm=R_here(0,0);
  double c_norm_partial=R_partial(0,0);
  
  return (R.array()*( R_partial.array()*c_norm/R_here.array()- (c_norm_partial*R.Ones(Rnrow,Rncol)).array() )/c_norm).matrix();
  
}



//this one is used in the old version and I still keep it here if other packages depends on it
// [[Rcpp::export]]
Eigen::MatrixXd separable_kernel (List R0, Eigen::VectorXd beta,String kernel_type, Eigen::VectorXd alpha ){
  Eigen::MatrixXd R0element = R0[0];
  int Rnrow = R0element.rows();
  int Rncol = R0element.cols();
  
  Eigen::MatrixXd R = R.Ones(Rnrow,Rncol);
  if(kernel_type=="matern_5_2"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (matern_5_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }else if(kernel_type=="matern_3_2"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (matern_3_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }
  else if(kernel_type=="pow_exp"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (pow_exp_funct(R0[i_ker],beta[i_ker],alpha[i_ker])).cwiseProduct(R);
    }
  }
  else if(kernel_type=="periodic_gauss"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (periodic_gauss_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }
  else if(kernel_type=="periodic_exp"){
    for (int i_ker = 0; i_ker < beta.size(); i_ker++){
      R = (periodic_exp_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
  }
  
  
  return R;
}

// [[Rcpp::export]]
Eigen::MatrixXd separable_multi_kernel (List R0, Eigen::VectorXd beta,Eigen::VectorXi kernel_type, Eigen::VectorXd alpha ){
  Eigen::MatrixXd R0element = R0[0];
  int Rnrow = R0element.rows();
  int Rncol = R0element.cols();

  Eigen::MatrixXd R = R.Ones(Rnrow,Rncol);
  //String kernel_type_i_ker;
  for (int i_ker = 0; i_ker < beta.size(); i_ker++){
   // kernel_type_i_ker=kernel_type[i_ker];
    if(kernel_type[i_ker]==3){
      R = (matern_5_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==2){
      R = (matern_3_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==1){
      R = (pow_exp_funct(R0[i_ker],beta[i_ker],alpha[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==4){
      R = (periodic_gauss_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==5){
      R = (periodic_exp_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }
      
  }
  return R;
}

////this is the for the normalizing constant in prediction
// [[Rcpp::export]]
Eigen::MatrixXd separable_multi_kernel_pred_periodic (const List R0, const Eigen::VectorXd beta,const Eigen::VectorXi kernel_type, const Eigen::VectorXd alpha, const Eigen::VectorXd perid_const){
  Eigen::MatrixXd R0element = R0[0];
  int Rnrow = R0element.rows();
  int Rncol = R0element.cols();
  
  Eigen::MatrixXd R = R.Ones(Rnrow,Rncol);
  //String kernel_type_i_ker;
  for (int i_ker = 0; i_ker < beta.size(); i_ker++){
    // kernel_type_i_ker=kernel_type[i_ker];
    if(kernel_type[i_ker]==3){
      R = (matern_5_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==2){
      R = (matern_3_2_funct(R0[i_ker],beta[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==1){
      R = (pow_exp_funct(R0[i_ker],beta[i_ker],alpha[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==4){
      R = (periodic_gauss_funct_fixed_normalized_const(R0[i_ker],beta[i_ker],perid_const[i_ker])).cwiseProduct(R);
    }else if(kernel_type[i_ker]==5){
      R = (periodic_exp_funct_fixed_normalized_const(R0[i_ker],beta[i_ker],perid_const[i_ker])).cwiseProduct(R);
    }
  }
  return R;
}



// [[Rcpp::export]]
Eigen::MatrixXd euclidean_distance(const MapMat & input1,const MapMat & input2){
  //input are n by p, where p is larger than n
  
  int num_obs1 = input1.rows();
  int num_obs2 = input2.rows();
  
  Eigen::MatrixXd R0=R0.Ones(num_obs1,num_obs2);
  
  for (int i = 0; i < num_obs1; i++){
    
    for (int j = 0; j < num_obs2; j++){
      R0(i,j)=sqrt((input1.row(i)-input2.row(j)).array().pow(2.0).sum());
    }
  }
  return R0;
}


// [[Rcpp::export]]
double log_marginal_lik(const Vec param,double nugget, const bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha ){
  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
    // nu=0;
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }

  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta, kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
	       
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition

  if(zero_mean=="Yes"){

   MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
  MatrixXd S_2= (yt_R_inv*output);
  double log_S_2=log(S_2(0,0));
    
   return (-(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
   
  }else{
    
  int q=X.cols();

  MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
  MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X

  LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
  MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
  MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solve
  MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
  MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
  double log_S_2=log(S_2(0,0));
  return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2);
  }
}


// [[Rcpp::export]]
double log_profile_lik(const Vec param,double nugget, const bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha ){
  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
    // nu=0;
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta, kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
  if(zero_mean=="Yes"){
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output);
    double log_S_2=log(S_2(0,0));
    
    return (-(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }else{
    
    //int q=X.cols();
    
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solve
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    double log_S_2=log(S_2(0,0));
    
    //change to profile likelihood seems only needs one step
    //return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2);
    return (-(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }
}


// [[Rcpp::export]]
double log_approx_ref_prior(const Vec param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b ){

  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  double t=CL.cwiseProduct(beta).sum()+nu;
  double part_I=-b*t;
  double part_II= a*log(t);
  return part_I+part_II;
}

// [[Rcpp::export]]
Eigen::VectorXd log_marginal_lik_deriv(const Eigen::VectorXd param,double nugget,  bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
    
  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
    
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
    
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  VectorXd ans=VectorXd::Ones(param_size);
  
  //String kernel_type_ti;
  
  if(zero_mean=="Yes"){
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output);
    //double log_S_2=log(S_2(0,0));

    
    MatrixXd dev_R_i;
    MatrixXd Vb_ti;
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0) ;  
    }
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
    }
    
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;

    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);

    MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
    MatrixXd dev_R_i;
    MatrixXd Wb_ti;
    //allow different choices of kernels
    
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      ans[ti]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      ans[p]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    

  }
      return ans;

}



// [[Rcpp::export]]
Eigen::VectorXd log_profile_lik_deriv(const Eigen::VectorXd param,double nugget,  bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  VectorXd ans=VectorXd::Ones(param_size);
  
  //String kernel_type_ti;
  MatrixXd Vb_ti;
  
  if(zero_mean=="Yes"){
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    
    
    MatrixXd dev_R_i;
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0) ;  
    }
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
    }
    
  }else{
   // int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    
    MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
    MatrixXd dev_R_i;
    MatrixXd Wb_ti;
    //allow different choices of kernels
    
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i)));
      Wb_ti=Vb_ti.transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i)));
      Wb_ti=Vb_ti.transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
      
      //Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      //ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
      
    }
    
    
  }
  return ans;
  
}


// [[Rcpp::export]]
Eigen::VectorXd log_approx_ref_prior_deriv(const Vec param,double nugget, bool nugget_est, const Eigen::VectorXd CL,const double a,const double b ){

  Eigen::VectorXd beta;
  Eigen::VectorXd return_vec;
  double nu=nugget;
  int param_size=param.size();
  if(!nugget_est){//not sure about the logical stuff. Previously (nugget_est==false)
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }

  //  double a=1/2.0;//let people specify
  // double b=(a+beta.size())/2.0;
  double t=CL.cwiseProduct(beta).sum()+nu;

  if(!nugget_est){
    return_vec=(a*CL/t- b*CL);
  }else{
    Eigen::VectorXd CL_1(param_size);
    CL_1.head(param_size-1)=CL;
    CL_1[param_size-1]=1;
    return_vec=(a*CL_1/t- b*CL_1);
  }
  return return_vec;

}

// [[Rcpp::export]]
double log_ref_marginal_post(const Eigen::VectorXd param,double nugget, bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
    
  Eigen::VectorXd beta;
  double nu=nugget;
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
    
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
    
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();

 // String kernel_type_ti;
  
  if(zero_mean=="Yes"){

    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output);
    double log_S_2=log(S_2(0,0));

    
    VectorXd ans=VectorXd::Ones(param_size);
    MatrixXd dev_R_i;
    List Vb(param_size);
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
    //  kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb[ti]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
    }

    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb[param_size-1]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
    }
  // int q=X.cols();
  MatrixXd IR(param_size+1,param_size+1);
  IR(0,0)=num_obs;
  
  for(int i=0;i<param_size;i++){
    MatrixXd Vb_i=Vb[i];
    IR(0,i+1)=IR(i+1,0)= Vb_i.trace();
    for(int j=0;j<param_size;j++){
      MatrixXd Vb_j=Vb[j];
      IR(i+1,j+1)=IR(j+1,i+1)=(Vb_i*Vb_j).trace();

    }
  }

  LLT<MatrixXd> lltOfIR(IR);
  MatrixXd LIR = lltOfIR.matrixL();
    
  return (-(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }else{
  int q=X.cols();
  MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
  MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;

  LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
  MatrixXd LX = lltOfXRinvX.matrixL();
  MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
  MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
  MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);

  MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
  MatrixXd dev_R_i;
  List Wb(param_size);
  
  
  for(int ti=0;ti<p;ti++){
  //  kernel_type_ti=kernel_type[ti];
    if(kernel_type[ti]==3){
      dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
    }else if(kernel_type[ti]==2){
      dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
    }else if(kernel_type[ti]==1){
      dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
    }else if(kernel_type[ti]==4){
      dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
    }else if(kernel_type[ti]==5){
      dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
    }
    Wb[ti]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
  }
  
  
  //the last one if the nugget exists
  if(nugget_est){
    dev_R_i=MatrixXd::Identity(num_obs,num_obs);
    Wb[param_size-1]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
  }
  MatrixXd IR(param_size+1,param_size+1);
  IR(0,0)=num_obs-q;
  for(int i=0;i<param_size;i++){
    MatrixXd Wb_i=Wb[i];
    IR(0,i+1)=IR(i+1,0)= Wb_i.trace();
    for(int j=0;j<param_size;j++){
      MatrixXd Wb_j=Wb[j];
      IR(i+1,j+1)=IR(j+1,i+1)=(Wb_i*Wb_j).trace();

    }
  }

  LLT<MatrixXd> lltOfIR(IR);
  MatrixXd LIR = lltOfIR.matrixL();

  double log_S_2=log(S_2(0,0));
    
  return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }
  //  return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+1/2.0*log(IR.determinant()) );
}



//this is a function to output a list, including theta_hat L (chlosky decomcoposition of R), LX(cholosky decomposition of Xt%*%R),  (the trend parameter), S2 (estimated sigma^2)
// [[Rcpp::export]]
List construct_rgasp(const Eigen::VectorXd beta,const double nu,  const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const  String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output,Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  List list_return(4);

  //similar to marginal likelihood
  //////// VectorXd beta= param.array().exp().matrix();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  // nu could be zero or nonzero

  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition

  list_return[0]=L; //first element to return
  if(zero_mean=="Yes"){
    list_return[1]=MatrixXd::Zero(1,1);
    list_return[2]= MatrixXd::Zero(1,1);
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd S_2= (yt_R_inv*output);
    list_return[3]=S_2(0,0)/(num_obs);

  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X

    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    list_return[1]=LX; //second element to return

    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd Xt_R_inv_y= X.transpose()*yt_R_inv.transpose();
    VectorXd theta_hat=LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(Xt_R_inv_y)); 
    list_return[2]=theta_hat;
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solver
    MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);

    list_return[3]=S_2(0,0)/(num_obs-q);
  }
  return list_return;

}


//this is a function to for prediction, including posterior mean, lower 95, upper 95 and standard deviation
// [[Rcpp::export]]
List pred_rgasp(const Eigen::VectorXd beta,const double nu, const  Eigen::Map<Eigen::MatrixXd> & input,  const Eigen::Map<Eigen::MatrixXd> & X,const  String zero_mean, const Eigen::Map<Eigen::MatrixXd> & output,const Eigen::Map<Eigen::MatrixXd> & testing_input, const Eigen::Map<Eigen::MatrixXd> & X_testing,
                const Eigen::Map<Eigen::MatrixXd> & L , Eigen::Map<Eigen::MatrixXd> & LX, Eigen::Map<Eigen::VectorXd> & theta_hat, double sigma2_hat,double q_025, double q_975, List r0,Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha,const String method, const bool interval_data){
  List pred(4);
    
  int num_testing_input=testing_input.rows();
  //int p=testing_input.cols();
  
  int p=beta.size();
  
  
  //int dim_inputs=input.cols();
  int num_obs=output.rows();
   
  //compute the vector for normalization for periodic 
  Eigen::VectorXd priodic_const=Eigen::VectorXd::Ones(p);
  for(int i_ker=0; i_ker<p;i_ker++){
    
    if(kernel_type[i_ker]==4){
      priodic_const[i_ker]=1.0/(2.0*sqrt(PI*beta[i_ker]));
      for(int ti=1; ti <11; ti++){
        priodic_const[i_ker]=priodic_const[i_ker]+1.0/sqrt(PI*beta[i_ker])*exp(-pow(ti,2.0)/(4.0*beta[i_ker]));
      }
    }else if(kernel_type[i_ker]==5){
      priodic_const[i_ker]=1.0/(PI*beta[i_ker]);
      for(int ti=1; ti<11; ti++){
        priodic_const[i_ker]=priodic_const[i_ker]+2.0*beta[i_ker]/((pow(beta[i_ker],2.0)+pow(ti,2.0))*PI);
      }
    }
  }
  
  //if(priodic_const.array().sum()==p ){
  //   MatrixXd r= separable_multi_kernel(r0,beta, kernel_type,alpha);
  //}else{
  MatrixXd r= separable_multi_kernel_pred_periodic(r0,beta, kernel_type,alpha,priodic_const);
  //}
    
  ////////I recalculate the following r
  /*
  MatrixXd  r(num_testing_input,num_obs);
  MatrixXd  diff1;
  List diff1_list(1); // 
  for(int qi=0; qi<num_testing_input;qi++){
    for(int i=0; i<dim_inputs;i++){
      diff1= (testing_input(qi,i)-input.col(i).array()).abs().matrix();
      diff1_list[0]=diff1;
      r0.col(i)=(matern_2_5_funct( diff1_list[0],beta[i])).col(0);                               
    }
    r.row(qi)=r0.col(0).transpose();
    if(dim_inputs>1){
      for(int i=1;i<dim_inputs;i++){
	r.row(qi)=(r.row(qi).array()*r0.col(i).transpose().array()).matrix();
      }
    }
  }
  */
  ///////////
  
  MatrixXd rt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(r.transpose()))).transpose();
  VectorXd c_star_star(num_testing_input);
  MatrixXd rtR_inv_r;

  

  if(zero_mean=="Yes"){
    if(interval_data){
      for(int i=0; i<num_testing_input;i++){
        rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
        c_star_star[i]=1+nu-rtR_inv_r(0,0);
      }
    }else{
      for(int i=0; i<num_testing_input;i++){
        rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
        c_star_star[i]=1-rtR_inv_r(0,0);
      }
    }
    VectorXd MU_testing=rt_R_inv*output;
    pred[0]=MU_testing;
    //VectorXd var=c_star_star*sigma2_hat;
    VectorXd var=c_star_star.array().abs().matrix()*sigma2_hat;  //when R is close to be singular, c_star_star can be very small negative
    pred[1]=MU_testing+var.array().sqrt().matrix()*q_025;
    pred[2]=MU_testing+var.array().sqrt().matrix()*q_975;
    if((method=="post_mode") || (method=="mmle")){
      pred[3]=var*(num_obs)/(num_obs-2);
    }else if(method=="mle"){
      pred[3]=var;
    }

  }else{

    //MatrixXd diff2;
    //MatrixXd X_testing_X_R_inv_r_i;
    if((method=="post_mode") || (method=="mmle")){
      int q=X.cols();
      MatrixXd  R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));  
      
      MatrixXd diff2;
      MatrixXd X_testing_X_R_inv_r_i;
      
      if(interval_data){
        for(int i=0; i<num_testing_input;i++){
          X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
    
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1+nu-rtR_inv_r(0,0)+diff2(0,0);
        }
      }else{
        for(int i=0; i<num_testing_input;i++){
          X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
          
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1-rtR_inv_r(0,0)+diff2(0,0);
        }
        
      }

      VectorXd MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);
      pred[0]=MU_testing;
      //VectorXd var=c_star_star*sigma2_hat;
      VectorXd var=c_star_star.array().abs().matrix()*sigma2_hat;  //when R is close to be singular, c_star_star can be very small negative
      pred[1]=MU_testing+var.array().sqrt().matrix()*q_025;
      pred[2]=MU_testing+var.array().sqrt().matrix()*q_975;
      pred[3]=var*(num_obs-q)/(num_obs-q-2);
    }else if(method=="mle"){
      if(interval_data){
        for(int i=0; i<num_testing_input;i++){
          //X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          //diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
          
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1+nu-rtR_inv_r(0,0);
        }
        
      }else{
        for(int i=0; i<num_testing_input;i++){
          //X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          //diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
          
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1-rtR_inv_r(0,0);
        }
      }
      
      VectorXd MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);
      pred[0]=MU_testing;
      //VectorXd var=c_star_star*sigma2_hat;
      VectorXd var=c_star_star.array().abs().matrix()*sigma2_hat;  //when R is close to be singular, c_star_star can be very small negative
      pred[1]=MU_testing+var.array().sqrt().matrix()*q_025;
      pred[2]=MU_testing+var.array().sqrt().matrix()*q_975;
      pred[3]=var;
      
    }
  }
  return pred;
   
}


//this is a function to for prediction, including posterior mean, lower 95, upper 95 and standard deviation
// [[Rcpp::export]]
List generate_predictive_mean_cov(const Eigen::VectorXd beta,const double nu, const  Eigen::Map<Eigen::MatrixXd> & input,  const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output,const Eigen::Map<Eigen::MatrixXd> & testing_input, const Eigen::Map<Eigen::MatrixXd> & X_testing, const Eigen::Map<Eigen::MatrixXd> & L , Eigen::Map<Eigen::MatrixXd> & LX, Eigen::Map<Eigen::VectorXd> & theta_hat,
                                  double sigma2_hat,List rr0, List r0,Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha,const String method,const bool sample_data){
  List mean_var(2);
    
  int num_testing_input=testing_input.rows();



  MatrixXd r= separable_multi_kernel(r0,beta, kernel_type,alpha); // looks this is num_testing_input x num_obs

  MatrixXd rr= separable_multi_kernel(rr0,beta, kernel_type,alpha);


  MatrixXd rt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(r.transpose()))).transpose();
  MatrixXd  rtR_inv_r= rt_R_inv*r.transpose();
  
  MatrixXd  C_star_star;
  if(zero_mean=="Yes"){
    if(sample_data){
         C_star_star= rr+nu*MatrixXd::Identity(num_testing_input,num_testing_input)-rtR_inv_r;
    }else{
        C_star_star= rr-rtR_inv_r;
    }
    VectorXd MU_testing=rt_R_inv*output;
    mean_var[0]=MU_testing;


    LLT<MatrixXd> lltOfC_star_star(C_star_star);
    MatrixXd LC_star_star = lltOfC_star_star.matrixL();

    mean_var[1]=sqrt(sigma2_hat)*LC_star_star;

  }else{
    MatrixXd  R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    VectorXd MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);
    
    if((method=="post_mode") || (method=="mmle")){
      MatrixXd  X_testing_X_R_inv_r=X_testing-r*R_inv_X;
      MatrixXd  diff2=X_testing_X_R_inv_r*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r.transpose())));
      if(sample_data){
         C_star_star= rr+nu*MatrixXd::Identity(num_testing_input,num_testing_input)-rtR_inv_r+diff2;
      }else{
        C_star_star= rr-rtR_inv_r+diff2;
      }
    }else if(method=="mle"){
      //MatrixXd  X_testing_X_R_inv_r=X_testing-r*R_inv_X;
      //MatrixXd  diff2=X_testing_X_R_inv_r*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r.transpose())));
      if(sample_data){
         C_star_star= rr+nu*MatrixXd::Identity(num_testing_input,num_testing_input)-rtR_inv_r;
      }else{
        C_star_star= rr-rtR_inv_r;
        
      }
    }
    
    LLT<MatrixXd> lltOfC_star_star(C_star_star);
    MatrixXd LC_star_star = lltOfC_star_star.matrixL();
    
    
    mean_var[0]=MU_testing;
    mean_var[1]=sqrt(sigma2_hat)*LC_star_star;
  
    //  mean_var[0]=r;
    //mean_var[1]=rr;
  }
  return mean_var;
   
}



//code for ppgasp



// [[Rcpp::export]]
double log_marginal_lik_ppgasp(const Vec param,double nugget, const bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha ){
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
    // nu=0;
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta, kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
  if(zero_mean=="Yes"){
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0));
    }
    
    //double log_S_2=log(S_2);
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }else{
    
    int q=X.cols();
    
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solve
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    //double log_S_2=((yt_R_inv.array()*output.array()).rowwise().sum().log()).sum();
    
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    }
    
    
    // double log_S_2=log(S_2);
    
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    return -k*(L.diagonal().array().log().matrix().sum())-k*(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2;
  }
  
  
}



// [[Rcpp::export]]
double log_profile_lik_ppgasp(const Vec param,double nugget, const bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha ){
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(!nugget_est){
    beta= param.array().exp().matrix();
    // nu=0;
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta, kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
  if(zero_mean=="Yes"){
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0));
    }
    
    //double log_S_2=log(S_2);
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }else{
    
    //int q=X.cols();
    
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solve
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    
    
    //double log_S_2=((yt_R_inv.array()*output.array()).rowwise().sum().log()).sum();
    
    
    double log_S_2=0;
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    }
     
    // double log_S_2=log(S_2);
    
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    
    //profile likelihood
    //return -k*(L.diagonal().array().log().matrix().sum())-k*(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2;
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2);
    
  }
  
  
}




// [[Rcpp::export]]
double log_ref_marginal_post_ppgasp(const Eigen::VectorXd param,double nugget, bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  
  // String kernel_type_ti;
  
  if(zero_mean=="Yes"){
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    //double log_S_2=log(S_2(0,0));
    
    double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0));
    }
    
   // double log_S_2=log(S_2);
    
    VectorXd ans=VectorXd::Ones(param_size);
    MatrixXd dev_R_i;
    List Vb(param_size);
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //  kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb[ti]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
    }
    
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb[param_size-1]=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
    }
    // int q=X.cols();
    MatrixXd IR(param_size+1,param_size+1);
    IR(0,0)=num_obs;
    
    for(int i=0;i<param_size;i++){
      MatrixXd Vb_i=Vb[i];
      IR(0,i+1)=IR(i+1,0)= Vb_i.trace();
      for(int j=0;j<param_size;j++){
        MatrixXd Vb_j=Vb[j];
        IR(i+1,j+1)=IR(j+1,i+1)=(Vb_i*Vb_j).trace();
        
      }
    }
    
    LLT<MatrixXd> lltOfIR(IR);
    MatrixXd LIR = lltOfIR.matrixL();
    
    return (-k*(L.diagonal().array().log().matrix().sum())-(num_obs)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    
    
    //double log_S_2=((yt_R_inv.array()*output.array()).rowwise().sum().log()).sum();
    
    
     
     double log_S_2=0;
    
    for(int loc_i=0;loc_i<k;loc_i++){
      log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    }
    
    
    //double log_S_2=log(S_2);
    
   // MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
    MatrixXd dev_R_i;
    List Wb(param_size);
    
    
    for(int ti=0;ti<p;ti++){
      //  kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Wb[ti]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
    }
    
    
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Wb[param_size-1]=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
    }
    MatrixXd IR(param_size+1,param_size+1);
    IR(0,0)=num_obs-q;
    for(int i=0;i<param_size;i++){
      MatrixXd Wb_i=Wb[i];
      IR(0,i+1)=IR(i+1,0)= Wb_i.trace();
      for(int j=0;j<param_size;j++){
        MatrixXd Wb_j=Wb[j];
        IR(i+1,j+1)=IR(j+1,i+1)=(Wb_i*Wb_j).trace();
        
      }
    }
    
    LLT<MatrixXd> lltOfIR(IR);
    MatrixXd LIR = lltOfIR.matrixL();
    
   // double log_S_2=log(S_2(0,0));
    
    return (-k*(L.diagonal().array().log().matrix().sum())-k*(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+ LIR.diagonal().array().log().matrix().sum());
  }
  //  return (-(L.diagonal().array().log().matrix().sum())-(LX.diagonal().array().log().matrix().sum())-(num_obs-q)/2.0*log_S_2+1/2.0*log(IR.determinant()) );
}


// [[Rcpp::export]]
Eigen::VectorXd log_marginal_lik_deriv_ppgasp(const Eigen::VectorXd param,double nugget,  bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  VectorXd ans=VectorXd::Ones(param_size);
  
  //String kernel_type_ti;
  
  if(zero_mean=="Yes"){
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    
    //double log_S_2=log(S_2(0,0));
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0);
      
    }
    MatrixXd dev_R_i;
    MatrixXd Vb_ti;
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i) )(0,0))/S_2_vec[loc_i];
      }
      ans[ti]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      
      //ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0) ;  
    }
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      double ratio=0;
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      //ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
      
    }
    
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    MatrixXd dev_R_i;
    MatrixXd Wb_ti;
    //allow different choices of kernels
    
    
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    }
    
    
    // double log_S_2=0;
    
    //for(int loc_i=0;loc_i<k;loc_i++){
    //  log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    //}
    
    
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      
      
      
      //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
      
      //MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
      
      ans[ti]=-0.5*k*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*ratio; 
      
      
      //ans[ti]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
    
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      //double S_2_dev=0;
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Wb_ti.diagonal().sum()  +(num_obs-q)/2.0*ratio; 
      
      //ans[p]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
  }
  return ans;
  
}


// [[Rcpp::export]]
Eigen::VectorXd log_profile_lik_deriv_ppgasp(const Eigen::VectorXd param,double nugget,  bool nugget_est, const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output, Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  
  Eigen::VectorXd beta;
  double nu=nugget;
  int k=output.cols();
  int param_size=param.size();
  if(nugget_est==false){//not sure about the logical stuff
    beta= param.array().exp().matrix();
  }else{
    beta=param.head(param_size-1).array().exp().matrix(); 
    nu=exp(param[param_size-1]); //nugget
  }
  int p=beta.size();
  int num_obs=output.rows();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  MatrixXd R_ori=  R;  // this is the one without the nugget
  
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  //not sure 
  
  LLT<MatrixXd> lltOfR(R);
  MatrixXd L = lltOfR.matrixL();
  VectorXd ans=VectorXd::Ones(param_size);
  
  //String kernel_type_ti;
  MatrixXd Vb_ti;
  
  if(zero_mean=="Yes"){
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    
    //double log_S_2=log(S_2(0,0));
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0);
      
    }
    MatrixXd dev_R_i;
    //allow different choices of kernels
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i) )(0,0))/S_2_vec[loc_i];
      }
      ans[ti]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      
      //ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0) ;  
    }
    //the last one if the nugget exists
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      Vb_ti=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i));
      
      double ratio=0;
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Vb_ti*(yt_R_inv.transpose()).col(loc_i))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Vb_ti.diagonal().sum()+num_obs/2.0*ratio;
      //ans[p]=-0.5*Vb_ti.diagonal().sum()+(num_obs/2.0*output.transpose()*Vb_ti*yt_R_inv.transpose()/ S_2(0,0))(0,0); 
      
    }
    
  }else{
    //int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X;
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X);
    MatrixXd LX = lltOfXRinvX.matrixL();
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose();
    MatrixXd dev_R_i;
    MatrixXd Wb_ti;
    //allow different choices of kernels
    
    
    VectorXd S_2_vec=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_vec[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    }
    
    
    // double log_S_2=0;
    
    //for(int loc_i=0;loc_i<k;loc_i++){
    //  log_S_2=log_S_2+log((yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0));
    //}
    
    
    for(int ti=0;ti<p;ti++){
      //kernel_type_ti=kernel_type[ti];
      if(kernel_type[ti]==3){
        dev_R_i=matern_5_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==2){
        dev_R_i=matern_3_2_deriv( R0[ti],R_ori,beta[ti]);  //now here I have R_ori instead of R
      }else if(kernel_type[ti]==1){
        dev_R_i=pow_exp_deriv( R0[ti],R_ori,beta[ti],alpha[ti]);   //now here I have R_ori instead of R
      }else if(kernel_type[ti]==4){
        dev_R_i = (periodic_gauss_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }else if(kernel_type[ti]==5){
        dev_R_i = (periodic_exp_deriv(R0[ti],R_ori,beta[ti])).cwiseProduct(R);
      }
      
      Vb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i)));
      Wb_ti=Vb_ti.transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      //ans[ti]=-0.5*Vb_ti.diagonal().sum()+(num_obs)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
      
      
     // Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      
      
      
      //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
      
      //MatrixXd Q_output= yt_R_inv.transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output;
      
      ans[ti]=-0.5*k*Vb_ti.diagonal().sum()+(num_obs)/2.0*ratio; 
      
      //ans[ti]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
    
    if(nugget_est){
      dev_R_i=MatrixXd::Identity(num_obs,num_obs);
      
      Vb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i)));
      Wb_ti=Vb_ti.transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      //Wb_ti=(L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(dev_R_i))).transpose()-dev_R_i*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv;
      
      //double S_2_dev=0;
      double ratio=0;
      
      for(int loc_i=0;loc_i<k;loc_i++){
        ratio=ratio+((output.col(loc_i).transpose()*Wb_ti.transpose()*(yt_R_inv.row(loc_i).transpose()-R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i)))(0,0))/S_2_vec[loc_i];
      }
      ans[p]=-0.5*k*Vb_ti.diagonal().sum()  +(num_obs)/2.0*ratio; 
      //ans[p]=-0.5*Wb_ti.diagonal().sum()+(num_obs-q)/2.0*(output.transpose()*Wb_ti.transpose()*Q_output/S_2(0,0))(0,0); 
    }
    
    
  }
  return ans;
  
}


//this is a function to output a list, including theta_hat L (chlosky decomcoposition of R), LX(cholosky decomposition of Xt%*%R),  (the trend parameter), S2 (estimated sigma^2)
// [[Rcpp::export]]
List construct_ppgasp(const Eigen::VectorXd beta,const double nu,  const List R0, const Eigen::Map<Eigen::MatrixXd> & X,const  String zero_mean,const Eigen::Map<Eigen::MatrixXd> & output,Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha){
  List list_return(4);
  
  //similar to marginal likelihood
  //////// VectorXd beta= param.array().exp().matrix();
  int num_obs=output.rows();
  int k=output.cols();
  MatrixXd R= separable_multi_kernel(R0,beta,kernel_type,alpha);
  R=R+nu*MatrixXd::Identity(num_obs,num_obs);  // nu could be zero or nonzero
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  
  list_return[0]=L; //first element to return
  if(zero_mean=="Yes"){
    list_return[1]=MatrixXd::Zero(1,1);
    list_return[2]= MatrixXd::Zero(1,1);
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    //MatrixXd S_2= (yt_R_inv*output);
    VectorXd S_2_all=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_all[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0);
    }
    
    list_return[3]= S_2_all/(num_obs);
    
  }else{
    int q=X.cols();
    MatrixXd R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X)); //one forward and one backward to compute R.inv%*%X
    MatrixXd Xt_R_inv_X=X.transpose()*R_inv_X; //Xt%*%R.inv%*%X
    
    LLT<MatrixXd> lltOfXRinvX(Xt_R_inv_X); // cholesky decomposition of Xt_R_inv_X called lltOfXRinvX
    MatrixXd LX = lltOfXRinvX.matrixL();  //  retrieve factor LX  in the decomposition 
    list_return[1]=LX; //second element to return
    
    MatrixXd yt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(output))).transpose(); 
    MatrixXd Xt_R_inv_y= X.transpose()*yt_R_inv.transpose();
    MatrixXd theta_hat=LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(Xt_R_inv_y)); 
    list_return[2]=theta_hat;
    MatrixXd R_inv_X_Xt_R_inv_X_inv_Xt_R_inv= R_inv_X*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(R_inv_X.transpose())));          //compute  R_inv_X_Xt_R_inv_X_inv_Xt_R_inv through one forward and one backward solver
    //MatrixXd S_2= (yt_R_inv*output-output.transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output);
    VectorXd S_2_all=VectorXd::Zero(k);
    
    for(int loc_i=0;loc_i<k;loc_i++){
      S_2_all[loc_i]=(yt_R_inv.row(loc_i)*output.col(loc_i))(0,0)-(output.col(loc_i).transpose()*R_inv_X_Xt_R_inv_X_inv_Xt_R_inv*output.col(loc_i))(0,0);
    }
    list_return[3]=S_2_all/(num_obs-q);
  }
  return list_return;
  
}


// [[Rcpp::export]]
List pred_ppgasp(const Eigen::VectorXd beta,const double nu, const  Eigen::Map<Eigen::MatrixXd> & input,  const Eigen::Map<Eigen::MatrixXd> & X,const  String zero_mean, const Eigen::Map<Eigen::MatrixXd> & output,const Eigen::Map<Eigen::MatrixXd> & testing_input, const Eigen::Map<Eigen::MatrixXd> & X_testing, 
                 const Eigen::Map<Eigen::MatrixXd> & L , Eigen::Map<Eigen::MatrixXd> & LX, Eigen::Map<Eigen::MatrixXd> & theta_hat,    const Eigen::Map<Eigen::VectorXd> &  sigma2_hat,double q_025, double q_975, List r0,Eigen::VectorXi kernel_type,const Eigen::VectorXd alpha, const  String method, const bool interval_data){
  List pred(4);
  
  int num_testing_input=testing_input.rows();
  //int p=testing_input.cols();
  int p=beta.size();
  
  //int dim_inputs=input.cols();
  int num_obs=output.rows();
  int k=output.cols();
  
  
  //compute the vector for normalization for periodic 
  Eigen::VectorXd priodic_const=Eigen::VectorXd::Ones(p);
  for(int i_ker=0; i_ker<p;i_ker++){
    
    if(kernel_type[i_ker]==4){
      priodic_const[i_ker]=1.0/(2.0*sqrt(PI*beta[i_ker]));
      for(int ti=1; ti <11; ti++){
        priodic_const[i_ker]=priodic_const[i_ker]+1.0/sqrt(PI*beta[i_ker])*exp(-pow(ti,2.0)/(4.0*beta[i_ker]));
      }
    }else if(kernel_type[i_ker]==5){
      priodic_const[i_ker]=1.0/(PI*beta[i_ker]);
      for(int ti=1; ti<11; ti++){
        priodic_const[i_ker]=priodic_const[i_ker]+2.0*beta[i_ker]/((pow(beta[i_ker],2.0)+pow(ti,2.0))*PI);
      }
    }
  }
  
  //if(priodic_const.array().sum()==p ){
    //   MatrixXd r= separable_multi_kernel(r0,beta, kernel_type,alpha);
    //}else{
  MatrixXd r= separable_multi_kernel_pred_periodic(r0,beta, kernel_type,alpha,priodic_const);
    //}
    


  MatrixXd rt_R_inv= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(r.transpose()))).transpose();
  VectorXd c_star_star(num_testing_input);
  MatrixXd rtR_inv_r;
  
  
  
  if(zero_mean=="Yes"){
    if(interval_data){
      for(int i_loc=0; i_loc<num_testing_input;i_loc++){
        rtR_inv_r=(rt_R_inv.row(i_loc)*r.row(i_loc).transpose());
        c_star_star[i_loc]=1+nu-rtR_inv_r(0,0);
      }
    }else{
      for(int i_loc=0; i_loc<num_testing_input;i_loc++){
        rtR_inv_r=(rt_R_inv.row(i_loc)*r.row(i_loc).transpose());
        c_star_star[i_loc]=1-rtR_inv_r(0,0);
      }
      
    }
    MatrixXd MU_testing=rt_R_inv*output;
    pred[0]=MU_testing;
    //VectorXd var=c_star_star*sigma2_hat;
    MatrixXd pred_var=MatrixXd::Zero(num_testing_input,k);

    for(int loc_i=0;loc_i<k;loc_i++){
      pred_var.col(loc_i)=  sigma2_hat[loc_i]*c_star_star.array().abs().matrix();
    }
    //VectorXd var=c_star_star.array().abs().matrix()*sigma2_hat;  //when R is close to be singular, c_star_star can be very small negative
    pred[1]=MU_testing+pred_var*q_025;
    pred[2]=MU_testing+pred_var*q_975;
    //pred[3]=pred_var*(num_obs)/(num_obs-2);
    
    if( (method=="post_mode") || (method=="mmle") ){
      pred[3]=pred_var*(num_obs)/(num_obs-2);
    }else if(method=="mle"){
      pred[3]=pred_var;
    }
    
  }else{
    
    if((method=="post_mode") || (method=="mmle")){
      int q=X.cols();
      MatrixXd diff2;
      MatrixXd  R_inv_X=L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(X));  
      MatrixXd X_testing_X_R_inv_r_i;
      if(interval_data){
        for(int i=0; i<num_testing_input;i++){
          X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
          
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1+nu-rtR_inv_r(0,0)+diff2(0,0);
        }
      }else{
        for(int i=0; i<num_testing_input;i++){
          X_testing_X_R_inv_r_i=X_testing.row(i)-r.row(i)*R_inv_X;
          diff2=X_testing_X_R_inv_r_i*(LX.transpose().triangularView<Upper>().solve(LX.triangularView<Lower>().solve(X_testing_X_R_inv_r_i.transpose())));
          
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1-rtR_inv_r(0,0)+diff2(0,0);
        }
      }
      MatrixXd MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);
      pred[0]=MU_testing;
      //VectorXd var=c_star_star*sigma2_hat;
      //VectorXd var=c_star_star.array().abs().matrix()*sigma2_hat;  //when R is close to be singular, c_star_star can be very small negative
      
      MatrixXd pred_var=MatrixXd::Zero(num_testing_input,k);
      
      for(int loc_i=0;loc_i<k;loc_i++){
        pred_var.col(loc_i)=  sigma2_hat[loc_i]*c_star_star.array().abs().matrix();
      }
      
      pred[1]=MU_testing+pred_var.array().sqrt().matrix()*q_025;
      pred[2]=MU_testing+pred_var.array().sqrt().matrix()*q_975;
      pred[3]=pred_var*(num_obs-q)/(num_obs-q-2);
    
    }else if(method=="mle"){
      if(interval_data){
        for(int i=0; i<num_testing_input;i++){
  
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1+nu-rtR_inv_r(0,0);
        }
      }else{
        for(int i=0; i<num_testing_input;i++){
          rtR_inv_r=(rt_R_inv.row(i)*r.row(i).transpose());
          c_star_star[i]=1-rtR_inv_r(0,0);
        }
      }
      MatrixXd MU_testing=X_testing*theta_hat+rt_R_inv*(output-X*theta_hat);
      pred[0]=MU_testing;
      
      MatrixXd pred_var=MatrixXd::Zero(num_testing_input,k);
      
      for(int loc_i=0;loc_i<k;loc_i++){
        pred_var.col(loc_i)=  sigma2_hat[loc_i]*c_star_star.array().abs().matrix();
      }
      pred[1]=MU_testing+pred_var.array().sqrt().matrix()*q_025;
      pred[2]=MU_testing+pred_var.array().sqrt().matrix()*q_975;
      pred[3]=pred_var;
      
    }

   //pred[3]=c_star_star; //test
  }
  return pred;
  
}


// [[Rcpp::export]]
bool test_const_column (const MapMat &d){
  
  int nrow = d.rows();
  int ncol = d.cols();
  bool res=false;//false means no constant column and true means yes
  double cur_value;
  for(int i=0; i<ncol; i++){
    cur_value=d(0,i);
    for(int j=1; j<nrow; j++){
      if(d(j,i)!=cur_value ){
        cur_value=d(j,i);
        break;
      }
    }
    if(cur_value==d(0,i)){
       res=true;
       break;
    }
    
  }
   return res;
}


