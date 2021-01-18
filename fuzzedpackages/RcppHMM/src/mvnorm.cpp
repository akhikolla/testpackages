
#include "mvnorm.h"

using namespace arma;

colvec dmvnormMultiple(mat x,  colvec mean,  mat sigma, bool logd) 
{ 
  int n = x.n_cols;
  int xdim = x.n_rows;
  vec out(n);
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    vec z = rooti * ( x.col(i) - mean) ;    
    out(i)      = constants - 0.5 * sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

mat rmvnormMultiple(int n, colvec mean, mat sigma)
{
    int ncols = sigma.n_cols;
    mat Y = randn(ncols, n);
    mat result = repmat(mean, 1, n) + chol(sigma) * Y;
    return result;
}

double dmvnormSingle(colvec x,  colvec mean,  mat sigma, bool logd) 
{   
  int xdim = x.n_rows;
  
  double out = 0.0;
  mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
  vec z = rooti * (x - mean) ;    
  out = constants - 0.5 * sum(z%z) + rootisum;       
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

colvec rmvnormSingle(colvec mean, mat sigma)
{
    int ncols = sigma.n_cols;
    vec Y = randn(ncols);
    vec result = mean + chol(sigma) * Y;
    return result;
}

bool isPositiveDefinite(mat matrix, double tol)
{        
    if (!matrix.is_square())
    // Rf_error("The matrix needs to be a square matrix");
      return false;
    if(!approx_equal(matrix, matrix.t(), "absdif", tol))
      return false;
    // Rf_error("The matrix needs to be a symmetric matrix");

    try{
      vec eigval = eig_sym(matrix);      
      int n = eigval.size();
      for(int i = 0; i < n; i++)
      if(eigval[i] < tol)
          return false;
      return true;
    }catch(const std::exception& e)
    {
      return false;
    }
    
}
