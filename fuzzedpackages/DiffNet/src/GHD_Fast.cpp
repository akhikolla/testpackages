#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
// If your compiler is to old, just disable / remove the following line
// [[Rcpp::plugins(cpp11)]]
#include<vector>
#include<cmath>
#include<map>
using namespace std;
//typedef Eigen::SparseMatrix<double> SpMat;
//typedef SpMat::InnerIterator InIterMat;

// [[Rcpp::export]]
double GHD_Fast(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
  
  int N = A.cols();
  double ghd_val = 0.0;
  double temp = 0.0;
  double global_A = 0.0 , global_B = 0.0;
  
  
  for (int j=0; j<N; ++j){
    for (Eigen::SparseMatrix<double>::InnerIterator i_(A,j); i_ ; ++i_){
      if (i_.index()!=j){
        global_A += i_.value();
      }
    }
    for (Eigen::SparseMatrix<double>::InnerIterator k_(B,j); k_; ++k_){
      if (k_.index()!=j){
        global_B += k_.value();
      }
    }
  }

  global_A = ((global_A)/N)/(N-1);
  global_B = ((global_B)/N)/(N-1);
  A -= B;

  double diff = (global_A-global_B);
  double part3 = (N*(N-1))*(diff*diff);
  
  for (int j=0; j<A.outerSize() ; j++)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator i_(A,j); i_; ++i_){
      if (i_.index()!=j)
      {
        ghd_val += (i_.value()*i_.value()) - (2*i_.value()*diff);
      }
    }
  }
  
  ghd_val = ((ghd_val + part3)/N)/(N-1);
  return(ghd_val);
}

// [[Rcpp::export]]
double MU_Fast(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
  
  int N = A.cols();
  double sa2 = 0.00, sb2 = 0.00, sa1 = 0.00, sb1 = 0.00, temp=0.00, temp1=0.00;
  double mu_permutation = 0.00;
  int nnz_A = A.nonZeros();
  int nnz_B = B.nonZeros();
  double thresh_val = (nnz_A+nnz_B)/(2*N);
  //For dense networks calculate directly
  if (thresh_val>(N/10)){
    sa1 = A.sum();
    sb1 = B.sum();
    sa2 = A.squaredNorm();
    sb2 = B.squaredNorm();
  }
  else{
    for (int j=0; j<N; ++j){
      for (Eigen::SparseMatrix<double>::InnerIterator i_(A,j); i_ ; ++i_){
        if (i_.index()!=j){
          sa1 += i_.value();
          sa2 += i_.value()*i_.value();
        }
      }
      for (Eigen::SparseMatrix<double>::InnerIterator k_(B,j); k_; ++k_){
        if (k_.index()!=j){
          sb1 += k_.value();
          sb2 += k_.value()*k_.value();
        }
      }
    }
  }
  mu_permutation = ((sa2+sb2)/N)/(N-1) - ((((2*sa1*sb1)/N)/N)/(N-1))/(N-1);
  return(mu_permutation);
}

// [[Rcpp::export]]
double STD_Fast(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B){
  
  int N = A.cols(), flag=0;
  double sa2 = 0.0, sb2 = 0.0, sa1 = 0.0, sb1 = 0.0, Ta = 0.0, Tb = 0.0, temp=0.0, temp1=0.0, temp2 = 0.0, temp3 = 0.0;
  double mu_permutation = 0.0;
  int nnz_A = A.nonZeros();
  int nnz_B = B.nonZeros();
  double thresh_val = (nnz_A+nnz_B)/(2*N);
  //For dense networks calculate directly
  if (thresh_val>(N/10)){
    sa1 = A.sum();
    sb1 = B.sum();
    sa2 = A.squaredNorm();
    sb2 = B.squaredNorm();
    flag=1;
  }
  for (int j=0; j<N; ++j){
    temp2 = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator i_(A,j); i_ ; ++i_){
      if (i_.index()!=j){
        temp = i_.value();
        temp2 += temp;
        if (flag==0)
        {
          sa1 += temp;
          sa2 += temp*temp;
        }
      }
    }
    Ta += (temp2)*(temp2);
    temp3 = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator k_(B,j); k_; ++k_){
      if (k_.index()!=j){
        temp1 = k_.value();
        temp3 += temp1;
        if (flag==0)
        {  
          sb1 += temp1;
          sb2 += temp1*temp1;
        }
      }
    }
    Tb += (temp3)*(temp3);
  }
  
  double Aa = (sa1*sa1);
  double Ab = (sb1*sb1);
  
  double Ba = (Ta - sa2);
  double Bb = (Tb - sb2);
  
  double Ca = (Aa + 2*sa2 - 4*Ta);
  double Cb = (Ab + 2*sb2 - 4*Tb);
    
  double D = (2*sa2*sb2);
  
  int p = N;
  
  float variance_permutation = 0.0; 
  variance_permutation = (4*(D + (4*Ba*Bb)/(p-2) + ((Ca*Cb)/(p-2))/(p-3) - (((Aa*Ab)/p)/(p-1))));
  return(variance_permutation);
}
