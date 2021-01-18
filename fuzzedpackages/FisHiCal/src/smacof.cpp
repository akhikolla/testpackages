#include <cstdlib>
#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;

/// functions that are not exported


// a: a matrix
// return a shallow copy of a
arma::mat shallow_copy(arma::mat a)
{
  int nrow = a.n_rows, ncol = a.n_cols;
  arma::mat b(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      b(i,j) = a(i,j);
    }
  }
  return b;
}

// x: a matrix providing a k-d configuration
// return the euclidean distance matrix induced by x
NumericMatrix calcDist(arma::mat x)
{
  //euclidian distances
  int nrow = x.n_rows, ncol = x.n_cols;
  double d = 0;
  NumericMatrix dist(nrow, nrow);
  
  for (int i = 0; i < nrow; i++)
  {
    for (int j = i+1; j < nrow; j++)
    {
       d = 0;
       for (int k = 0; k < ncol; k++)
       {
         d = d + ((x(i,k) - x(j,k))*(x(i,k) - x(j,k)));
       }
       dist(i,j) = sqrt(d);
       dist(j,i) = dist(i,j);
    }
    
  }
  return dist;
}

// dist: distances matrix
// diss: dissimilarities matrix
// w: weights matrix
// return B (see smacof documentation), where:
// B[i,j] = - diss[i,j]*w[i,j]/dist[i,j] if dist[i,j] != 0, i !=j
//       0 if dist[i,j] = 0, i != j
//       - sum(B[i,k])  k = 1..n, k!=i, i==j
arma::mat calcB(NumericMatrix dist, NumericMatrix diss, NumericMatrix w)
{
   int nrow = dist.nrow(), ncol = dist.ncol();
   double diagVal = 0.0;
   NumericMatrix B(nrow, ncol);
   NumericVector diag(nrow); 
   for (int i = 0; i < nrow; i++) 
   {
      diagVal = 0.0;
      for (int j = 0; j < ncol; j++) 
      {
        if(dist(i,j) != 0 && i != j )
        {
          B(i,j) = (-w(i,j)*diss(i,j))/dist(i,j);
          diagVal += B(i,j);
        }
        
      }
      diag(i) = diagVal;
   }
   for (int i = 0; i < nrow; i++)
   {
     B(i,i) = -1*diag(i);
   }
   
   // convert to arma::mat for matrix manipulation purposes
   arma::mat BMat(B.begin(), nrow, ncol, false);
   return BMat;
}

// w: a weight matrix
// return v (see smacof documentation), where
// v[i,j] = -w[i,j] if i!j
//          sum(w[i,k])  k = 1..n, k!=i, i==j
arma::mat calcV(NumericMatrix w)
{
  int nrow = w.nrow(), ncol = w.ncol();
  NumericMatrix v(nrow, ncol);
  NumericVector diag(nrow); 
  for (int i = 0; i < nrow; i++) 
   {
      double diagVal = 0.0;
      for (int j = 0; j < ncol; j++) 
      {
        if(i != j )
        {       
          diagVal += w(i,j);
          v(i,j) = -w(i,j);
        }
      }
      diag(i) = diagVal;
   }
   for (int i = 0; i < nrow; i++)
   {
     v(i,i) = diag(i);
   }
   
   // convert to arma::mat for matrix manipulation purposes
   arma::mat vMat(v.begin(), nrow, ncol, false);
   return vMat;
}

// calculates the stress function = sum(w[i,j](diss[i,j] - dist[i,j])^2) i < j
// dist: distances matrix
// diss: dissimilarities matrix
// w: weights matrix
// returns the stress value
double fstress(NumericMatrix dist, NumericMatrix diss, NumericMatrix w)
{
  double stress = 0.0;
  int nrow = dist.nrow(), ncol = dist.ncol();
  for(int i = 0; i < nrow; i++)
  {
    for (int j = i+1; j < ncol; j++)
    {
      stress = stress + w(i,j)*(diss(i,j)-dist(i,j))*(diss(i,j)-dist(i,j));
      
    }
  }
  return stress; 
}



// diss: n*n matrix of dissimilarities
// w: n*n matrix of weights - designed for local weights and 
//    assumed to be irreducible
// initConf: n*k initial configuration matrix
// infD: inifinte distance value - for unknown dissimilarities
// niter: number of iterations for smacof
// eps: difference value for convergence
// verbose: whetehr or not to print stress and iteration data
// implements the smacof (Scaling by Majorizing a Complicated Function)algorithm 
// returns the k-d configuration
// [[Rcpp::export]]
Rcpp::NumericMatrix smacof(NumericMatrix diss, NumericMatrix w,  NumericMatrix initConf, 
                            int niter, double eps, bool verbose)
{

     
        int nrow = diss.nrow();
        arma::mat init(initConf.begin(), nrow, initConf.ncol(), false);
        arma::mat z = shallow_copy(init);
        try
        {
          //prepare matrices
          arma::mat v = calcV(w);
          // Moore-Penrose inverse
          arma::mat mpV =  arma::pinv(v);
          // shallow copy is required is order not to change z when changing x
          arma::mat x = shallow_copy(init);
          NumericMatrix dist = calcDist(z);
          
          // initial stress
          double prevStress = fstress(dist, diss, w);
          double stressDiff = 2*eps;
          double stress = 0.0;
          if (verbose)
          {
            Rcout << "inital stress is: "<< prevStress <<std::endl;  
          }
          // start iteration
          int iter = 0;
          while(iter < niter && (stressDiff > eps))
          {
            iter = iter + 1;
            x = mpV*calcB(dist, diss, w)*z;
            dist = calcDist(x);
            stress = fstress(dist, diss, w);
            stressDiff = prevStress - stress;
            z = shallow_copy(x);
            prevStress = stress;
            if (verbose)
            {
              Rcout << "stress in iteration " << iter << " is: " << stress << std::endl;      
            }
          }
          //complete iterations or converged
          if (verbose && stressDiff > eps)
          {
            Rcout << "lsmacof didn't converge after " << iter << " iterations. Consider increasing number of iterations." << std::endl;    
          }
        } catch (const runtime_error& error)
        {
           Rcout << "ERROR IN SMACOF: failed to find the Moore-Penrose inverse" << std::endl;
        }
        
        return Rcpp::as<NumericMatrix>(wrap(z));
}
