#include <Rcpp.h>
#include <iostream>
#include <RcppEigen.h>
//#include <random>
//#include <stdlib.h>
//[[Rcpp::depends(RcppEigen)]]
#include <math.h>
#include<algorithm>
//#include <stdio.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using namespace std;



// This is a simplRcppEigen.package.skeleton("test")
//e example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
double calculate(MatrixXd C,int select)
{
  double temp;
  double mean = 0;
  double distance = 0;
  MatrixXd Distance(C.rows(),C.rows());

  for(int Mrow = 0;Mrow<(C.rows()-1);Mrow++)
    for(int Mcol = Mrow+1;Mcol<C.rows();Mcol++)
    {
      temp = (C.row(Mrow) - C.row(Mcol)).squaredNorm();
      Distance(Mcol,Mrow) = sqrt(temp);
      Distance(Mrow,Mcol) = sqrt(temp);

    }
    for(int i =0;i<C.rows();i++)
      for(int j =0;j<C.rows ();j++)
      {
        if(i==j)
        {
          Distance(i,j)=0;
        }
      }
      for(int i = 0;i<Distance.rows();i++)
      {
        sort(Distance.col(i).data(),Distance.col(i).data()+Distance.rows());
      }

      if(Distance.rows()-1<= select)
      {
        distance = Distance.sum()/(Distance.rows()*(Distance.rows()-1));
        mean = distance;
      }
      else
      {
        for(int i = 1;i< select+1;i++)
          for(int j = 0;j<Distance.cols();j++)
          {
            distance = distance + Distance(i,j);
          }
          mean = distance/(select * Distance.cols());
      }
      return mean;
}

// [[Rcpp::export]]
List permute (Eigen::MatrixXd iTest,Eigen::MatrixXd jTest,SEXP times,SEXP selectvec)
{

  double cTemp = 0;
  double mean = 0;
  int ptime = Rcpp::as<int>(times);
  NumericVector pselect = Rcpp::as<NumericVector> (selectvec);
  int length = pselect.length();
  NumericVector pmean(length);
  NumericVector psd(length);
  int row1 = iTest.rows();
  int dim = pselect.length();
  int combine = iTest.cols()+jTest.cols();
  MatrixXd C(row1,combine);
  MatrixXd temp(ptime,dim);
  VectorXi indices = VectorXi::LinSpaced(row1,0,row1);

  for (int i = 0; i < ptime ;i++)
  {
    unsigned int rand = R::runif(1000000,10000000);
    for (int j = 0; j< row1; j++ )
    {
        indices[j] = (indices[j]+rand) % row1;
    }
    jTest = indices.asPermutation() * jTest;
    C<<iTest,jTest;
    for (int j = 0; j<dim;j++)
    {
      temp(i,j) = calculate(C,pselect[j]);
    }
  }
  for(int i = 0;i<temp.cols();i++)
  {
    for(int j = 0;j<temp.rows();j++)
    {
      mean = temp.col(i).mean();
      cTemp+=(temp(j,i)-mean)*(temp(j,i)-mean);
    }
    psd[i] = sqrt(cTemp/(temp.rows()-1));
    pmean[i] = temp.col(i).mean();
    cTemp = 0;

  }

  return Rcpp::List::create(Rcpp::Named("pmean") = pmean,Rcpp::Named("psd") = psd,Rcpp::Named("temp") = temp);

}

// [[Rcpp::export]]
 double  dis_value(Eigen::MatrixXd iTest,Eigen::MatrixXd jTest,SEXP select)
{

  int row = iTest.rows();
  int tcol1 = iTest.cols();
  int tcol2 = jTest.cols();
  int select1 = Rcpp::as<int>(select);
  Eigen::MatrixXd Distance(row,row);
  Eigen::MatrixXd C(row,tcol1+tcol2);
  double temp;
  double mean = 0;
  C<<iTest,jTest;

  for(int Mrow = 0;Mrow<(C.rows()-1);Mrow++)
    for(int Mcol = Mrow+1;Mcol<C.rows();Mcol++)
    {
      temp = (C.row(Mrow) - C.row(Mcol)).squaredNorm();
      Distance(Mcol,Mrow) = sqrt(temp);
      Distance(Mrow,Mcol) = sqrt(temp);
    }
    for(int i =0;i<C.rows();i++)
      for(int j =0;j<C.rows ();j++)
      {
        if(i==j)
        {
          Distance(i,j)=0;
        }

      }

      for(int i = 0;i<Distance.rows();i++)
      {
        sort(Distance.col(i).data(),Distance.col(i).data()+Distance.rows());

      }
      double distance=0;

  if(Distance.rows()-1<=select1)
  {
    mean = Distance.sum()/(Distance.rows()*(Distance.rows()-1));
  }
  else
  {
    for(int i = 1;i<select1+1;i++)
      for(int j = 0;j<Distance.cols();j++)
      {
        distance = distance + Distance(i,j);
      }
      mean = distance/(select1 * row);
  }
  return mean;

}

// calculate mean and sd



