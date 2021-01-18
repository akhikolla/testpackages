#include "curstat.h"

//////////////////////////////////////////////////////////////////////////////////
//
// Maximum Likellood Estimator (MLE)
//
// Export: Yes
// Date Created: 07.07.2017
// Created by: Kim Hendrickx
//////////////////////////////////////////////////////////////////////////////////
//'@title Maximum Likelihood Estimator
//'
//'@description The function ComputeMLE computes the Maximum Likelihood Estimator of the distribution function under current status data.
//'
//'@param data Dataframe with three variables:
//'\describe{
//'     \item{t}{Observation points t sorted in ascending order. All observations need to be positive. The total number of unique observation points equals \code{length(t)}.}
//'     \item{freq1}{Frequency of observation t satisfying \eqn{x \le t}.
//'                   The total number of observations with censoring indicator \eqn{\delta =1} equals \code{sum(freq1)}. }
//'     \item{freq2}{Frequency of observation t. The sample size equals \code{sum(freq2)}. If no tied observations are present in the data \code{length(t)} equals \code{sum(freq2)}. }
//'}
//'
//'@details In the current status model, the variable of interest \eqn{X} with distribution function \eqn{F} is not observed directly.
//'A censoring variable \eqn{T} is observed instead together with the indicator \eqn{\Delta = (X \le T)}.
//' ComputeMLE computes the MLE of \eqn{F} based on a sample of size \code{n <- sum(data$freq2)}.
//'
//'
//'
//'@return Dataframe with two variables :
//'\describe{
//'     \item{x}{jump locations of the MLE}
//'     \item{mle}{MLE evaluated at the jump locations}
//' }
//'
//'@references Groeneboom, P. and Hendrickx, K. (2017). The nonparametric bootstrap for the current status model. Electronic Journal of Statistics 11(2):3446-3848.
//'@seealso \code{\link{ComputeConfIntervals}}
//'
//'
//'@examples
//'library(Rcpp)
//'library(curstatCI)
//'
//'# sample size
//'n <- 1000
//'
//'# Uniform data  U(0,2)
//'set.seed(2)
//'y <- runif(n,0,2)
//'t <- runif(n,0,2)
//'delta <- as.numeric(y <= t)
//'
//'A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
//'mle <-ComputeMLE(A)
//'plot(mle$x, mle$mle,type ='s', ylim=c(0,1), main= "",ylab="",xlab="",las=1)
//'
//'@export
// [[Rcpp::export]]
DataFrame ComputeMLE(DataFrame data)
{
  double          *data0,*cumw,*cs,*F,*jumploc,*y;
  int             i,j,n,njumps,*freq1,*freq2;

  DataFrame DF = Rcpp::DataFrame(data);
  NumericVector xx = DF[0];
  IntegerVector freq01 = DF[1];
  IntegerVector freq02 = DF[2];

  n = (int)xx.size();

  data0= new double[n+1];
  freq1= new int[n+1];
  freq2= new int[n+1];

  data0[0]=0;

  for (i=1;i<=n;i++)
  {
    data0[i]=(double)xx[i-1];
    freq1[i]=(int)freq01[i-1];
    freq2[i]=(int)freq02[i-1];
  }

  F= new double[n+1];
  cumw= new double[n+1];
  cs= new double[n+1];
  y= new double[n+1];
  jumploc= new double[n+1];


  F[0]=0;
  cumw[0]=cs[0]=0;

  y[0]=0;

  for (i=1;i<=n;i++)
  {
    cs[i]=cs[i-1]+(double)freq1[i];
    cumw[i]=cumw[i-1]+(double)freq2[i];
  }

  convexmin(n,cumw,cs,y);

  j=0;

  jumploc[0]=0;

  for (i=1;i<=n;i++)
  {
    if (y[i]>y[i-1])
    {
      j++;
      F[j]=y[i];
      jumploc[j]=data0[i];
    }
  }

  njumps=j;

  NumericVector out1(njumps+1);
  NumericVector out2(njumps+1);

  for (i=0;i<=njumps;i++)
  {
    out1[i]=jumploc[i];
    out2[i] = F[i];
  }


  DataFrame out = DataFrame::create(Rcpp::Named("x")=out1, Rcpp::Named("mle")=out2);


  // free memory

  delete[] F, delete[] cumw, delete[] cs, delete[] y,
  delete[] jumploc,  delete[] data0;
  delete[] freq1; delete[] freq2;

  return out;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Smoothed Maximum Likellood Estimator (SMLE)
//
// Export: Yes
// Date Created: 07.07.2017
// Created by: Kim Hendrickx
//////////////////////////////////////////////////////////////////////////////////

//'@title Smoothed Maximum Likelihood Estimator
//'
//'@description The function ComputeSMLE computes the Smoothed Maximum Likelihood Estimator of the distribution function under current status data.
//'
//'@param data Dataframe with three variables:
//'\describe{
//'     \item{t}{Observation points t sorted in ascending order. All observations need to be positive. The total number of unique observation points equals \code{length(t)}.}
//'     \item{freq1}{Frequency of observation t satisfying \eqn{x \le t}.
//'                   The total number of observations with censoring indicator \eqn{\delta =1} equals \code{sum(freq1)}. }
//'     \item{freq2}{Frequency of observation t. The sample size equals \code{sum(freq2)}. If no tied observations are present in the data \code{length(t)} equals \code{sum(freq2)}. }
//'}
//'
//'@param bw numeric vector of size \code{length(x)}. This vector contains the pointwise bandwidth values for each point in the vector x.
//'
//'@param x numeric vector containing the points where the confidence intervals are computed.
//'
//'@details In the current status model, the variable of interest \eqn{X} with distribution function \eqn{F} is not observed directly.
//'A censoring variable \eqn{T} is observed instead together with the indicator \eqn{\Delta = (X \le T)}.
//' ComputeSMLE computes the SMLE of \eqn{F} based on a sample of size \code{n <- sum(data$freq2)}.
//' The bandwidth parameter vector that minimizes the pointwise Mean Squared Error using the subsampling principle in combination with undersmoothing is returned by the function \link{ComputeBW}.
//'
//'@return SMLE(x) Smoothed Maximum Likelihood Estimator. This is a vector of size \code{length(x)} containing the values of the SMLE for each point in the vector x.
//'
//'@references Groeneboom, P. and Hendrickx, K. (2017). The nonparametric bootstrap for the current status model. Electronic Journal of Statistics 11(2):3446-3848.
//'
//'@seealso \code{\link{ComputeConfIntervals}}
//'
//'@examples
//'library(Rcpp)
//'library(curstatCI)
//'
//'# sample size
//'n <- 1000
//'
//'# Uniform data  U(0,2)
//'set.seed(2)
//'y <- runif(n,0,2)
//'t <- runif(n,0,2)
//'delta <- as.numeric(y <= t)
//'
//'A<-cbind(t[order(t)], delta[order(t)], rep(1,n))
//'grid <-seq(0,2 ,by = 0.01)
//'
//'# bandwidth vector
//'h<-rep(2*n^-0.2,length(grid))
//'
//'smle <-ComputeSMLE(A,grid,h)
//'plot(grid, smle,type ='l', ylim=c(0,1), main= "",ylab="",xlab="",las=1)
//'
//'
//'@export
// [[Rcpp::export]]
NumericVector ComputeSMLE(DataFrame data, NumericVector x, NumericVector bw)
{
  double          *data0,*cumw,*cs,*F,*p,*jumploc,*y,*grid,A,B, *h;
  int             i,j,k,n,ngrid,njumps,*freq1,*freq2;

  DataFrame DF = Rcpp::DataFrame(data);
  NumericVector xx = DF[0];
  IntegerVector freq01 = DF[1];
  IntegerVector freq02 = DF[2];

  n = (int)xx.size();
  ngrid = (int)x.size();

  data0= new double[n+1];
  freq1= new int[n+1];
  freq2= new int[n+1];

  data0[0]=0;

  for (i=1;i<=n;i++)
  {
    data0[i]=(double)xx[i-1];
    freq1[i]=(int)freq01[i-1];
    freq2[i]=(int)freq02[i-1];
  }

  A = 0.0;
  B = data0[n];

  F= new double[n+1];
  p= new double[n+1];
  cumw= new double[n+1];
  cs= new double[n+1];
  y= new double[n+1];
  jumploc= new double[n+1];
  grid= new double[ngrid];
  h= new double[ngrid];

  for (i=0;i<ngrid;i++){
    grid[i]=(double)x[i];
    h[i] = (double)bw[i];
  }


  F[0]=0;
  cumw[0]=cs[0]=0;

  y[0]=0;

  for (i=1;i<=n;i++)
  {
    cs[i]=cs[i-1]+(double)freq1[i];
    cumw[i]=cumw[i-1]+(double)freq2[i];
  }

  convexmin(n,cumw,cs,y);

  j=0;

  jumploc[0]=0;

  for (i=1;i<=n;i++)
  {
    if (y[i]>y[i-1])
    {
      j++;
      p[j]=y[i]-y[i-1];
      F[j]=y[i];
      jumploc[j]=data0[i];
    }
  }

  njumps=j;

  double sum, t1,t2,t3;
  NumericVector out(ngrid);

  for(i=0; i<ngrid;i++)
  {
    sum=0;
    for (k=1;k<=njumps;k++)
    {
      t1=(grid[i]-jumploc[k])/h[i];
      t2=(grid[i]+jumploc[k]-2*A)/h[i];
      t3=(2*B-grid[i]-jumploc[k])/h[i];
      sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
    }
    out[i] =fmax(0,sum) ;
  }


  // free memory

  delete[] F, delete[] cumw, delete[] cs, delete[] y, delete[] p;
  delete[] jumploc,  delete[] data0, delete[] grid;
  delete[] freq1; delete[] freq2;

  return out;

}
