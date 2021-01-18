#include <Rcpp.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <stdio.h>
#include "tf_cox.h"

extern "C" {
#include "tf.h"
#include "utils.h"
}


using namespace Rcpp;

// [[Rcpp::export]]
List tfCox_onelambda(int ord, double alpha, double lambda, IntegerVector discrete, int ndis,
            double stepSize, NumericMatrix X, NumericMatrix init_theta, 
            NumericVector Time, IntegerVector status, IntegerVector indx_time, 
            IntegerVector tie,  int n, int p, IntegerMatrix Perm, 
            IntegerMatrix Rank, IntegerVector thin, NumericVector vec_xtol, double tol, int niter,
            int backtracking) {
  
  NumericMatrix theta(n,p);
  NumericMatrix theta_new(n,p);
  NumericVector timeOrd(n);
  IntegerVector statusOrd(n), indices(n);
  NumericVector grad(n), beta_new(n), temp(n);
  int i, j, k, u, thinning, nbacktrack, flag, n_avg, start;
  int iter=0; int converge=0;
  double sum1=0; double sum2=0; double sumsq = 0; double dev=0;
  double diff1, diff2, new_loss, old_loss, avg; 
  double intercept = 0, lambda_new;
  double scalef, x_tol;
  double *x, *y, *w;

  x = (double *) malloc(n * sizeof(double));
  y = (double *) malloc(n * sizeof(double));
  w = (double *) malloc(n * sizeof(double));
  for (i=0; i<n; i++) w[i] = 1;

  for (int i=0; i<n; i++) {
    for (int j=0; j<p; j++) {
      theta(i,j) = init_theta(i,j);
    }
  }
  
  /* allocated space for return */
  double *beta, *weight;  int *m; 
  beta = (double *) malloc(n*sizeof(double));
  weight = (double *)malloc(n*sizeof(double));
  m = (int *) malloc(sizeof(int));
  
  while ((converge==0) & (iter<niter)) {

    iter++;
    grad = gradient(theta, n, p, status, indx_time, tie); 
    
    for (j=0; j<p; j++) {
      /* ------------------------------------------ */
      flag = 0;
      for (k=0; k<ndis; k++) {
        if (discrete[k] == j) {
          flag = 1;
          break;
        }
      }
      
      if (flag) {
        /*  gradient descent  */
        for (i=0; i<n; i++) temp[i] = theta(i,j) - stepSize*grad[i];
        for (i=0; i<n; i++) {
          y[i] = temp[Perm(i,j)-1];
          x[i] = X(Perm(i,j)-1,j);
        }
        
        avg = y[0];
        n_avg = 1;
        start = 0; 
        for (i=0; i<n-1; i++) {
          if (x[i] == x[i+1]) {
            avg += y[i+1];
            n_avg++;
          }
          else {
            for (k=start; k<start+n_avg; k++) beta_new[k] = avg/n_avg;
            avg=y[start+n_avg];
            start = start+n_avg;
            n_avg=1;
          }
        }
        for (k=start; k<start+n_avg; k++) beta_new[k] = avg/n_avg;
        
        for (i=0; i<n; i++) theta_new(i,j) = beta_new[Rank(i,j)-1]; 
        
      } else {
        /* proximal gradient descent  */
        thinning = thin[j];
        // Rcout << "j=" << j << "," << thinning <<"\n"; 
        for (i=0; i<n; i++) temp[i] = theta(i,j) - stepSize*grad[i];
        for (i=0; i<n; i++) {
          y[i] = temp[Perm(i,j)-1];
          x[i] = X(Perm(i,j)-1,j);
        }
        lambda_new = alpha*lambda*stepSize;

        /* trend filtering  */
        x_tol = vec_xtol[j]; 
        tf_main(x, y, w, n, ord, lambda_new, thinning, x_tol, beta, weight, m); 
        
        
        /* reorder the fitted values */
        u = 0; 
        for (i=0; i<m[0]; i++) {
          for (k=0; k < ((int)weight[i]); k++) {
            indices[u] = i;
            u++;
          }
        }
        for (i=0; i<n; i++) beta_new[i] = beta[indices[i]];
        for (i=0; i<n; i++) theta_new(i,j) = beta_new[Rank(i,j)-1]; 
      }
      /* ------------------------------------------ */
    }
      
    
    /* ------------------------------------------ */
    if (backtracking) {
      /* backtracking line search condition */
      new_loss = negloglik(theta_new, n, p, status, indx_time, tie);
      old_loss = negloglik(theta, n, p, status, indx_time, tie);
      diff1 = 0; diff2 = 0;
      for (j=0; j<p; j++) {
        for (i=0; i<n; i++) {
          diff1 += grad[i] * (theta_new(i,j)-theta(i,j));
          diff2 += pow(theta_new(i,j)-theta(i,j), 2.0);
        }
      }
      
      /* shorten stepSize and repeat previous steps */
      nbacktrack = 0;
      while (new_loss > old_loss + diff1 + 0.5*stepSize*diff2) {
        stepSize = stepSize * 0.8;
        nbacktrack++;
        if (nbacktrack > 100) break;
        
        /* descent step */
        for (j=0; j<p; j++) {
          /* ------------------------------------------ */
          flag = 0;
          for (k=0; k<ndis; k++) {
            if (discrete[k] == j) {
              flag = 1;
              break;
            }
          }
          
          if (flag) {
            /*  gradient descent  */
            for (i=0; i<n; i++) temp[i] = theta(i,j) - stepSize*grad[i];
            for (i=0; i<n; i++) {
              y[i] = temp[Perm(i,j)-1];
              x[i] = X(Perm(i,j)-1,j);
            }
            
            avg = y[0];
            n_avg = 1;
            start = 0; 
            for (i=0; i<n-1; i++) {
              if (x[i] == x[i+1]) {
                avg += y[i+1];
                n_avg++;
              }
              else {
                for (k=start; k<start+n_avg; k++) beta_new[k] = avg/n_avg;
                avg=y[start+n_avg];
                start = start+n_avg;
                n_avg=1;
              }
            }
            for (k=start; k<start+n_avg; k++) beta_new[k] = avg/n_avg;
            
            for (i=0; i<n; i++) theta_new(i,j) = beta_new[Rank(i,j)-1]; 
            
          } else {
            /* proximal gradient descent  */
            thinning = thin[j];
            for (i=0; i<n; i++) temp[i] = theta(i,j) - stepSize*grad[i];
            for (i=0; i<n; i++) {
              y[i] = temp[Perm(i,j)-1];
              x[i] = X(Perm(i,j)-1,j);
            }
            lambda_new = alpha*lambda*stepSize;
            
            /* trend filtering  */
            x_tol = vec_xtol[j]; 
            tf_main(x, y, w, n, ord, lambda_new, thinning, x_tol, beta, weight, m); 
            
            
            /* reorder the fitted values */
            u = 0; 
            for (i=0; i<m[0]; i++) {
              for (k=0; k < ((int)weight[i]); k++) {
                indices[u] = i;
                u++;
              }
            }
            for (i=0; i<n; i++) beta_new[i] = beta[indices[i]];
            for (i=0; i<n; i++) theta_new(i,j) = beta_new[Rank(i,j)-1]; 
          }
          /* ------------------------------------------ */
        }
        
        /* backtracking line search condition */
        new_loss = negloglik(theta_new, n, p, status, indx_time, tie);
        old_loss = negloglik(theta, n, p, status, indx_time, tie);
        diff1 = 0; diff2 = 0;
        for (j=0; j<p; j++) {
          for (i=0; i<n; i++) {
            diff1 += grad[i] * (theta_new(i,j)-theta(i,j));
            diff2 += pow(theta_new(i,j)-theta(i,j), 2.0);
          }
        }
      }

    }
    /* end of backtracking */
    /* ------------------------------------------ */
    
    
    for (j=0; j<p; j++) {
      /* centering and intercept  */
      sum1 = 0;
      for (i=0; i<n; i++) sum1 += theta_new(i,j);
      intercept += sum1/n;
      
      sumsq = 0;
      for (i=0; i<n; i++) {
        theta_new(i,j) = theta_new(i,j) - sum1/n;
        sumsq += pow(theta_new(i,j),2.0);
      }
      
      /* truncation for alpha < 1 */
      if (alpha < 1) {
        scalef = 1-(1-alpha)*lambda*stepSize/pow(sumsq,0.5);
        if (scalef <= 0) {
          for (i=0; i<n; i++) {
            theta_new(i,j) = 0;
          }
        }
        else {
          for (i=0; i<n; i++) {
            theta_new(i,j) = theta_new(i,j)*scalef;
          }
        }
      }
    }
    
    // convergence criteria
    sum1 = 0;
    sum2 = 0;
    for (i=0; i<n; i++) {
      for (j=0; j<p; j++) {
        sum1 += pow(theta_new(i,j) - theta(i,j), 2.0);
        sum2 += pow(theta(i,j), 2.0);
      }
    }
    
    for (i=0; i<n; i++) {
      for (j=0; j<p; j++) {
        theta(i,j) = theta_new(i,j);
      }
    }
    
    if (iter > 1) {
      if (sum2) {
        dev = sum1/sum2;
        if (dev < tol) converge = 1;
      }
      else 
        converge = 1;
    }

  }

  free(x);
  free(y);
  free(w);
  free(beta);
  free(weight);
  free(m);
  return List::create(
    Named("theta") = theta,
    Named("intercept") = intercept,
    Named("iter") = iter);
}


// ordered x and y but need thining 
// thinning = 1 means that it needs thinning

void tf_main(double *x, double *y, double *w, int n, int ord, double lambda, 
             int thinning, double x_tol, double *beta, double *weight, int *m)
{
  /* constant of trend filtering */
  int family = 0;
  int lam_flag = 1;
  int nlambda = 1;
  double lambda_min_ratio = 1e-5;
  int rho=1;
  double obj_tol=1e-5; 
  double obj_tol_newton=1e-5;
  int max_iter=200L;
  int max_iter_newton=50L;
  double alpha_ls=0.5;
  double gamma_ls=0.8;
  int max_iter_ls=30L; 
  int tridiag=0;
  int verbose=0;
  double *Lambda;
  double *obj;
  int *iter;
  int *opt_status;
  int *df;
  int i; 
  df = (int *)    malloc(sizeof(int));
  obj    = (double *) malloc(max_iter * sizeof(double));
  iter   = (int *)    malloc(sizeof(int));
  opt_status = (int *)    malloc(sizeof(int));
  Lambda = (double *)    malloc(sizeof(double)); 
  df[0] = 0;
  for (i = 0; i < max_iter; i++) obj[i] = 0;
  iter[0] = 0;
  opt_status[0] = 0;
  Lambda[0] = lambda;
  
  /* thining */
  double * xt;
  double * yt;
  double * wt;
  int nt;
  xt = yt = wt = NULL;
  if (thinning) {
    thin(x,y,w,n,ord,&xt,&yt,&wt,&nt,x_tol);
    if( xt != NULL )
    {
      x = xt;
      y = yt;
      w = wt;
      n = nt;
    }
  }
  
  /* Initalize output arrays with 0's */
  double *out;
  out = (double *) malloc(n * sizeof(double));
  for (i = 0; i < n; i++) out[i] = 0;
  double *beta0; beta0=NULL;

  /* trend filtering */
  tf_admm(x, y, w, n, ord, family, max_iter, beta0, lam_flag, Lambda,
          nlambda, lambda_min_ratio, tridiag, df, out, obj, iter, opt_status,
          rho, obj_tol, obj_tol_newton, alpha_ls, gamma_ls, max_iter_ls,
          max_iter_newton, verbose);
  
  /* prepare output of the allocated space */
  for (i=0; i<n; i++) {
    beta[i] = out[i];
    weight[i] = w[i];
  }
  m[0] = n;
  return; 
  /* free storage */
  free(df);
  free(obj);
  free(iter);
  free(opt_status);
  free(Lambda);
  free(out);
  if(xt != NULL) free(xt);
  if(yt != NULL) free(yt);
  if(wt != NULL) free(wt);
}


NumericVector gradient(NumericMatrix theta, int n, int p,
                       IntegerVector status, IntegerVector indx_time, 
                       IntegerVector tie) 
{
  NumericVector out(n);
  NumericVector grad(n);
  NumericVector exp_linear(n);
  NumericVector sum_riskset(n);
  IntegerVector statusOrd(n);
  double temp=0, prob_ui, sum_count=0; 
  int i, j, k, u, firstobs; 
  firstobs=n;
  
  for (i=0; i<n; i++) {
    statusOrd[i] = status[indx_time[i]-1];
    grad[i] = 0;
  }
  for (i=0; i<n; i++) {
    if (statusOrd[i]==1) {
      firstobs = i;
      break;
    }
  }
  for (i=0; i<n; i++) {
    temp = 0;
    for (j=0; j<p; j++) {
      temp += theta(indx_time[i]-1,j);
    }
    exp_linear[i] = exp(temp);
  }
  for (i=0; i<n; i++) {
    temp = 0;
    for (k=i; k<n; k++) {
      temp += exp_linear[k];
    }
    sum_riskset[i] = temp;
  }

  for (u=0; u<n; u++) {
    for (i=firstobs; i<=u; i++) {
      if (statusOrd[i]) {
        prob_ui = exp_linear[u]/sum_riskset[i]; 
        sum_count = 0; 
        for (k=i; k<=i+tie[i]; k++) sum_count += statusOrd[k];
        grad[u] += sum_count * prob_ui;
      }
      if (tie[i]) i += tie[i];
    }
    grad[u] = grad[u]/n;
    if (statusOrd[u]) grad[u] = grad[u] - 1.0/n;
  }
  
  for (i=0; i<n; i++) out[indx_time[i]-1] = grad[i]; 

  return out;
}

// negative log partial likelihood 
double negloglik(NumericMatrix theta, int n, int p, 
                 IntegerVector status, IntegerVector indx_time, 
                 IntegerVector tie)
{
  NumericVector out(n);
  IntegerVector statusOrd(n);
  NumericVector linear_term(n); 
  double temp=0, s=0;
  int i, j, k, count; 
  
  for (i=0; i<n; i++) {
    statusOrd[i] = status[indx_time[i]-1];
    temp = 0;
    for (j=0; j<p; j++) {
      temp += theta(indx_time[i]-1,j);
    }
    linear_term[i] = temp;
  }
  
  for (i=0; i<n; i++) {
    count = 0;
    temp = 0;
    if (statusOrd[i]) {
      for (k=i; k<=i+tie[i]; k++) {
        s += linear_term[k]; 
        count += statusOrd[k];
      }
      for (k=i; k<n; k++) {
        temp += exp(linear_term[k]);
      }
      s += count*log(temp);
    }
    if (tie[i]) i += tie[i]; 
  }
  s = s*(-1/n);
  return s;
}




