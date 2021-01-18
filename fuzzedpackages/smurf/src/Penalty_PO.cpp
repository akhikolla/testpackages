#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Soft thresholding operator S_y(x) where x is a vector
arma::vec soft_thresh(const arma::vec& x, const double y) {
  
  arma::vec z = zeros(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; i++) {
    //z(i) = std::max(0.0, x(i) - y) - std::max(0.0, -x(i) - y);
    // sign(x) * max(|x| - y, 0)
    z(i) = std::copysign(1, x(i)) * std::max(std::abs(x(i)) - y, 0.0);
  }
  return z;
}

// Soft thresholding operator S_y(x) where both x and y are vectors
arma::vec soft_thresh_vec(const arma::vec& x, const arma::vec& y) {
  
  arma::vec z = zeros(x.n_elem);
  for (unsigned int i = 0; i < x.n_elem; i++) {
    // sign(x) * max(|x| - y, 0)
    z(i) = std::copysign(1, x(i)) * std::max(std::abs(x(i)) - y(i), 0.0);
  }
  return z;
}

// Group soft thresholding operator S_y(x) where x is a vector
arma::vec group_soft_thresh(const arma::vec& x, const double y) {
  
  double norm_x = norm(x, 2);
  arma::vec z = zeros(x.n_elem);
  
  // Avoid problems if negative or zero norm
  if (norm_x > 0) {
    for (unsigned int i = 0; i < x.n_elem; i++) {
      z(i) = x(i) * std::max(1 - y / norm_x, 0.0);
    }
  }
  
  return z;
}


// Implementation of ADMM algorithm to compute proximal operator for Fused Lasso, Generalized Fused Lasso,
// 2D Fused Lasso and Graph-Guided Fused Lasso
//
// beta_tilde: Vector of coefficients after gradient update
// slambda: Penalty parameter (lambda) multiplied with step size (s)
// lambda1: Vector with lambda1 multiplied with penalty weights (for Lasso). The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
// lambda2: lambda2 multiplied with penalty weight (for Group Lasso). The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
// penmat: (Weighted) Penalty matrix
// Q: Matrix with eigenvectors of (weighted) penalty matrix in the columns
// eigval: Vector with eigenvalues of (weighted) penalty matrix
// fast: Boolean indicating if the fast version to compute inverses, in the ADMM update of x, is used
// maxiter: Maximum number of iterations
// rho: Initial value for augmented Lagrangian parameter
// beta_old: Vector with previous coefficient estimates
// [[Rcpp::export]]
arma::vec admm_po_cpp(const arma::vec& beta_tilde, const double slambda, const arma::vec& lambda1, const double lambda2, 
                      const arma::mat& penmat, const arma::mat& Q, const arma::vec& eigval, const bool fast,
                      const int maxiter, double rho, const arma::vec& beta_old) {
  
  arma::mat penmat_t = penmat.t();
  int m = penmat.n_rows;
  int d = penmat.n_cols;
  
  // Initialize values
  arma::vec x = zeros(d);
  arma::vec xhat = zeros(d);
  arma::vec z_old = zeros(m);
  // Use starting value for beta
  arma::vec z_new = penmat * beta_old;
  // Starting value for u is zero arma::vector
  arma::vec u = zeros(m);
  
  // Relative tolerance
  double eps_rel = pow(10, -10.0);
  // Absolute tolerance
  double eps_abs = pow(10, -12.0);
  // Tolerance for primal feasibility condition
  double eps_pri = sqrt((double) m) * eps_abs + eps_rel * std::max(norm(penmat * x, 2), norm(z_new, 2));
  // Tolerance for dual feasibility condition
  double eps_dual = sqrt((double) d) * eps_abs + eps_rel * rho * norm(penmat_t * u, 2);
  // Norm of primal residuals
  double r_norm = norm(penmat * x - z_new, 2);
  // Norm of dual residuals
  double s_norm = norm(- rho * penmat_t * (z_new - z_old), 2);
  
  
  // Relaxation parameter
  double xi = 1.5;
  // Iteration counter, note that we start from 1!
  int iter = 1;
  
  arma::mat ADMM_aux;
  arma::mat Qt;
  // Auxiliary matrix
  if (fast) {
    Qt = Q.t();
    // Fast version to compute inverse
    ADMM_aux = eye(d, d) - rho * Q * diagmat(1/(1/eigval + rho)) * Qt;
    
  } else {
    // Slower (standard) version to compute inverse
    ADMM_aux = inv(eye(d, d) + rho * penmat_t * penmat);
  }
  
  
  // Parameters to update rho, see Zhu (2017), JCGS
  double mu = 10;
  double eta = 2;
  double rho_old;
  
  while ((r_norm > eps_pri || s_norm > eps_dual || iter==1) && iter < maxiter) {
    
    // Check for interrupt every 1000 iterations
    if (iter % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
        
    
    z_old = z_new;
    
    // Update x
    x = ADMM_aux * (beta_tilde + rho * penmat_t * (z_old - u));
    // Relaxation
    xhat = xi * penmat * x + (1 - xi) * z_old;
    // Update z
    z_new = soft_thresh(xhat + u, slambda / rho);
    // Update u
    u = u + xhat - z_new;
    
    
    // Tolerance for primal feasibility condition
    eps_pri = sqrt((double) m) * eps_abs + eps_rel * std::max(norm(penmat * x, 2), norm(z_new, 2));
    // Tolerance for dual feasibility condition
    eps_dual = sqrt((double) d) * eps_abs + eps_rel * rho * norm(penmat_t * u, 2);
    // Norm of primal residuals
    r_norm = norm(penmat * x - z_new, 2);
    // Norm of dual residuals
    s_norm = norm(- rho * penmat_t * (z_new - z_old), 2);
    
    
    // Check if rho needs to be increased
    if (r_norm / eps_pri >= mu * s_norm / eps_dual) {
      
      rho_old = rho;
      // Update rho
      rho = eta * rho_old;
      // Update auxiliary matrix
      if (fast) {
        // Fast version to compute inverse
        ADMM_aux = eye(d, d) - rho * Q * diagmat(1/(1/eigval + rho)) * Qt;
        
      } else {
        // Slower (standard) version to compute inverse
        ADMM_aux = inv(eye(d, d) + rho * penmat_t * penmat);
      }
      // update u
      u = u * rho_old / rho;
      
      // Check if rho needs to be decreased
    } else if (s_norm / eps_dual >= mu * r_norm / eps_pri) {
      
      rho_old = rho;
      // Update rho
      rho = rho_old / eta;
      // Update auxiliary matrix
      if (fast) {
        // Fast version to compute inverse
        ADMM_aux = eye(d, d) - rho * Q * diagmat(1/(1/eigval + rho)) * Qt;
        
      } else {
        // Slower (standard) version to compute inverse
        ADMM_aux = inv(eye(d, d) + rho * penmat_t * penmat);
      }
      // update u
      u = u * rho_old / rho;
    }
    
    iter++;
  }

  if (lambda2 > 0) {
    // In case lambda2 is non-zero, an extra step is needed to obtain the proximal operator (see Liu et al. (2010))
    x = group_soft_thresh(x, slambda * lambda2);
  }
  
  // In case lambda1 is non-zero, an extra step is needed to obtain the proximal operator (see Liu et al. (2010))
  if (lambda1.n_elem > 1) {
    // Vector version
    x = soft_thresh_vec(x, slambda * lambda1);
    
  } else {
    // Only one element, check if non-zero first
    if (lambda1(0) > 0) {
      x = soft_thresh(x, slambda * lambda1(0));
    }
  }
  
  return x;
}
