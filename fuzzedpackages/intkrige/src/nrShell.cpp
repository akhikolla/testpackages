// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <cmath>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Remember that % denotes element-wise multiplication

// Basic calculation of the Newton-Rhapson step for simple kriging
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Vector of three coefficients for generalized L2 distance.
// - pen    - value of penalty term (used to impose constraint that
//            sum(abs(lam)) = 1)
// Outputs:
// - lamUp  - an updated version of the lambda vector.
// [[Rcpp::export]]
arma::colvec nrStep(const arma::colvec & lam, const arma::mat & pCovC,
                    const arma::colvec & lCovC, const arma::mat & pCovR,
                    const arma::colvec & lCovR, const arma::mat & pCovCR,
                    const arma::colvec & lCovCR, const arma::colvec & A,
                    const double & pen) {
  // Calculate absolute value and sign of lambas only once
  arma::colvec absLam = arma::abs(lam);
  arma::colvec sgnLam = arma::sign(lam);

  // Calculate CR/lambda summation (used in the hessian also)
  arma::colvec calcCR = (lam.t() * pCovCR).t();

  // Calculate the penalty
  int l = lam.size();
  double penalty = 0;

  for(int i = 0; i < l; i++){
    // Shorthand for sum = sum + x[i]
    penalty += absLam[i];
  }

  // Gradient Calculations
  // Calculate radius portion of gradient as we use it in the Hessian.
  arma::colvec rGrad = (absLam.t() * pCovR).t() - lCovR;

  // Calculate gradient
  // Line 1 - center, line 2 - radius, line 3,4 - center/radius
  arma::colvec grad = 2*(A(0)*((lam.t() * pCovC).t() - lCovC) +
    A(1)*sgnLam % rGrad + // Calculated previously
    A(2)*((absLam.t() * pCovCR).t() + sgnLam % calcCR -
    (sgnLam + 1) % lCovCR));
  // Penalty parameter added later.

  // Now calculate the hessian matrix
  arma::mat hessian(size(pCovC), arma::fill::zeros);

  // Penalty parameter tucked in with Radius calculation
  for(int i = 0; i < l; i++){
    for(int j = 0; j < l; j++){
      hessian(i, j) = 2*(A(0)*(pCovC(i, j)) +
        (sgnLam(i)*sgnLam(j)*(A(1)*pCovR(i, j) + pen)) +
        A(2)*((sgnLam(i)+sgnLam(j))*pCovCR(i, j)));
    }
  }

  // Update the gradient calculation with the penalty.
  // Update Hessian diagonals for quadratic adjustment
  for(int i = 0; i < l; i++){
    grad(i) += 2*sgnLam(i)*pen*(penalty - 1);
    hessian(i, i) += 2*( A(1)*rGrad(i) + A(2)*(calcCR(i) - lCovCR(i)) + pen*(penalty - 1) )/absLam(i);
  }

  // Update the values of lambda
  arma::colvec lamUp = arma::solve(hessian, hessian*lam - grad);

  return lamUp;
}

// Basic calculation of the Newton-Rhapson step for simple kriging
// Uses penalty to discourage lambda values near 0 initially.
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Vector of three coefficients for generalized L2 distance.
// - pen    - value of penalty term (used to impose constraint that
//            sum(abs(lam)) = 1)
// Outputs:
// - lamUp  - an updated version of the lambda vector.
// [[Rcpp::export]]
arma::colvec nrStep_long(const arma::colvec & lam, const arma::mat & pCovC,
                         const arma::colvec & lCovC, const arma::mat & pCovR,
                         const arma::colvec & lCovR, const arma::mat & pCovCR,
                         const arma::colvec & lCovCR, const arma::colvec & A,
                         const double & pen, const double & len) {
  // Calculate absolute value and sign of lambas only once
  arma::colvec absLam = arma::abs(lam);
  arma::colvec sgnLam = arma::sign(lam);

  // Calculate CR/lambda summation (used in the hessian also)
  arma::colvec calcCR = (lam.t() * pCovCR).t();

  // Calculate the penalty
  int l = lam.size();
  double penalty = 0;

  for(int i = 0; i < l; i++){
    // Shorthand for sum = sum + x[i]
    penalty += absLam[i];
  }

  // Gradient Calculations
  // Calculate radius portion of gradient as we use it in the Hessian.
  arma::colvec rGrad = (absLam.t() * pCovR).t() - lCovR;

  // Calculate gradient
  // Line 1 - center, line 2 - radius, line 3,4 - center/radius
  arma::colvec grad = 2*( A(0)*((lam.t() * pCovC).t() - lCovC) +
    A(1)*sgnLam % rGrad + // Calculated previously
    A(2)*((absLam.t() * pCovCR).t() + (sgnLam % calcCR) -
    (sgnLam + 1) % lCovCR) );
  // Penalty parameter added later.

  // Now calculate the hessian matrix
  arma::mat hessian(size(pCovC), arma::fill::zeros);

  for(int i = 0; i < l; i++){
    for(int j = 0; j < l; j++){
      hessian(i, j) = 2*( ( A(0)*pCovC(i, j) ) +
        ( sgnLam(i)*sgnLam(j)*A(1)*pCovR(i, j) ) +
        ( A(2)*(sgnLam(i)+sgnLam(j))*pCovCR(i, j) )  +
        (sgnLam(i)*sgnLam(j)/pen) );
    }
  }

  // Update the gradient calculation with the penalty.
  // Update Hessian diagonals for quadratic adjustment
  for(int i = 0; i < l; i++){
    grad(i) += ( 2*sgnLam(i)*(penalty - 1)/pen ) - ( 2*pen/(lam(i)*len*len) );
    hessian(i, i) += 2*pen/(lam(i)*lam(i)*len*len);
  }

  // Update the values of lambda
  arma::colvec lamUp = arma::solve(hessian, hessian*lam - grad);
  //arma::colvec lamUp = hessian*lam - grad;

  return lamUp;
}

// Basic calculation of the Newton-Rhapson step for ordinary kriging
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Vector of three coefficients for generalized L2 distance.
// - r      - value of penalty term (used to impose the n+1 constraints that
//            sum(lam) = 1 and each lam_i >= 0)
// - threshold - calculated in newRap_2. Determine the effect magnitude of 0.
//               Used to shift the lower bound of the barrier.
// - len       - length of the original lambda vector.
// Outputs:
// - lamUp  - an updated version of the lambda vector.
// [[Rcpp::export]]
arma::colvec nrStep_2(const arma::colvec & lam, const arma::mat & pCovC,
                      const arma::colvec & lCovC, const arma::mat & pCovR,
                      const arma::colvec & lCovR, const arma::mat & pCovCR,
                      const arma::colvec & lCovCR, const arma::colvec & A,
                      const double & r, const double & threshold,
                      const double & len) {

  // The absolute value of the lambdas SHOULD be equivalent to lambda itself
  // in this case. This condition should be ensured OUTSIDE of this function.

  // Calculate lambda penalty parameters
  double pen1 = arma::sum(lam) - 1;

  // Gradient Calculations
  // Calculate radius portion of gradient as we use it in the Hessian.
  arma::colvec grad = 2*(lam.t() * (A(0)*pCovC + A(1)*pCovR + 2*A(2)*pCovCR)).t() -
    2*(A(0)*lCovC + A(1)*lCovR + 2*A(2)*lCovCR) -
    ( r/( (lam + threshold) ) ) + ( 2.0*pen1/(r*len) ); // penalty parameter included

  // Now calculate the hessian matrix
  arma::mat hessian = 2*(A(0)*pCovC + A(1)*pCovR + 2*A(2)*pCovCR) + (2.0/(r*len));

  // Update the gradient calculation with the penalty.
  // Update Hessian diagonals for quadratic adjustment
  int l = lam.n_elem;
  for(int i = 0; i < l; i++){
    hessian(i, i) += r/( (lam(i) + threshold)*(lam(i) + threshold) );
  }

  // Update the values of lambda
  arma::colvec lamUp = arma::solve(hessian, hessian*lam - grad);

  return lamUp;
}


// Calculation of the Newton-Rhapson step for simple kriging
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Vector of three generalized L2 distance weights.
// - threshold - when |lambda_i| < threshold, value is set to 0.
// - tolq   - convergence tolerance for NR step
// - tolp   - convergence tolerance for penalty criteria
// - maxq   - max number of iterations for NR step
// - maxp   - max number of iterations to grow the penalty term.
// - eta    - growth parmameter for penalty criteria (eta > 1).
// Outputs:
// - lam    - a (hopefully) convergent solution for the lambda vector
// [[Rcpp::export]]
arma::colvec newRap(arma::colvec lam, const arma::mat & pCovC,
                    const arma::colvec & lCovC, const arma::mat & pCovR,
                    const arma::colvec & lCovR, const arma::mat & pCovCR,
                    const arma::colvec & lCovCR, const arma::colvec & A,
                    const double & threshold,
                    const double & tolq, const int & maxq,
                    const double & tolp, const int & maxp,
                    const double & eta) {

  // Determine initial length of lambda vector (size will likely change)
  // with each iteration.
  double len = lam.size();

  // Create a colvec with one entry that stores any warnings in the
  // iterative process.
  arma::colvec isWarn(1, arma::fill::zeros);

  // Initialize variables used in the following loops
  arma::colvec lamUp(len, arma::fill::zeros);
  arma::colvec tdiff(len, arma::fill::zeros);

  // Initialize vector of indices. + 1 ensures that all elements
  // are included initially.
  arma::uvec tindices = find( (arma::abs(lam) + 1) > 0 );

  // Update the values of lambda using no penalty initially.
  int i = 0;
  do{
    // Resize our updated lambda matrix
    lamUp.set_size( tindices.size() );

    lamUp = nrStep(lam.elem(tindices), pCovC.submat(tindices, tindices),
                   lCovC.elem(tindices), pCovR.submat(tindices, tindices),
                   lCovR.elem(tindices), pCovCR.submat(tindices, tindices),
                   lCovCR.elem(tindices), A, 0);

    // Compute the differences in the vector
    tdiff = arma::abs(lam.elem(tindices) - lamUp);

    // Overwrite lambda now that we have calculated the difference
    lam.elem(tindices) = lamUp;

    // Set elements of lambda close enough to 0, equal to 0
    //lam.elem( find(arma::abs(lam) <  threshold) ).zeros();
    lam.clean(threshold);

    tindices = find(lam);

    // If only one lambda remains, return this and suspend the search.
    if( tindices.n_elem < 2 ){
      if(tindices.n_elem < 1){
        Rcpp::Rcout << "No viable lambdas exist, returning prediction of 0..." << std::endl;
        isWarn(0) = 1;
      }
      lam = lam.zeros();
      lam.elem(tindices) += 1;

      arma::colvec lam2 = arma::join_vert(isWarn, lam);

      return lam2;
    }

    i++; // move along the counter
  } while(tdiff.max() > tolq && i < maxq); // End the do-loop

  // Now, determine the appropriate scale of a starting penalty parameter.
  // (Ignoring additive constants)
  // std:: for element by element abs() calculations.
  double penalty = arma::as_scalar(arma::sum(arma::abs(lam)) - 1);

  // Adjust the penalty function to be roughly the same scale as the
  // minimization equation. Do this by ensuring that the first
  // non-zero penalty is on the order of 1.
  double adj = 1.0/(penalty * penalty);

  int j = 0; // Reset counter parameters (must redeclare j in new scope)

  // Indicator variable that will end the penalty march for non-convergent solutions.
  bool lastTry = TRUE;

  while( (std::abs(penalty) > tolp || (i == maxq && lastTry)) && j < maxp ) {
    i=0; // reset the inner counter.

    // If we have found a non-convergent, yet feasible, solution. Then try one
    // more iteration of the penalty parameter.
    if(std::abs(penalty) <= tolp){
      lastTry = FALSE;
    }

    do{
      // Resize our updated lambda matrix
      lamUp.set_size( tindices.size() );

      lamUp = nrStep(lam.elem(tindices), pCovC.submat(tindices, tindices),
                     lCovC.elem(tindices), pCovR.submat(tindices, tindices),
                     lCovR.elem(tindices), pCovCR.submat(tindices, tindices),
                     lCovCR.elem(tindices), A, adj); // adj is the penalty term

      // Compute the differences in the vector
      tdiff = arma::abs(lam.elem(tindices) - lamUp);

      // Overwrite lambda now that we have calculated the difference
      lam.elem(tindices) = lamUp;

      // Set elements of lambda close enough to 0, equal to 0
      //lam.elem( find(arma::abs(lam) <  threshold) ).zeros();
      lam.clean(threshold);

      tindices = find(lam);

      // If only one lambda remains, return this and suspend the search.
      if( tindices.n_elem < 2 ){
        if(tindices.n_elem < 1){
          Rcpp::Rcout << "No viable lambdas exist, returning prediction of 0..." << std::endl;
          isWarn(0) = 1;
        }
        lam = lam.zeros();
        lam.elem(tindices) += 1;

        arma::colvec lam2 = arma::join_vert(isWarn, lam);

        return lam2;
      }


      i++; // move along the counter
    } while(tdiff.max() > tolq && i < maxq); // end the inner do-loop

    // Calculate the penalty (will determine if we stay in the loop)
    penalty = arma::as_scalar(arma::sum(arma::abs(lam)) - 1);

    // Grow the adjustment parameter according to the growth parameter eta.
    adj /= eta;

    j++; // move along the j counter

    // Warn the user if the returned solution was the result on a non-convergent
    // optimization.
    if(tdiff.max() > tolq && !lastTry){
      Rcpp::Rcout << "Feasible solution obtained from a non-convergent optimization step " <<
        "with max(diff) = " << tdiff.max() << "." << std::endl;
      isWarn(0) = 1;
    }

  }

  if(std::abs(penalty) > tolp){
    Rcpp::Rcout <<
      "Convergent, feasible solution not obtained with sum(abs(lam)) = " <<
        penalty + 1 << "," << " and max(diff) = " << tdiff.max() << "." << std::endl;
    isWarn(0) = 1;
  }

  arma::colvec lam2 = arma::join_vert(isWarn, lam);

  return lam2;
}

// Calculation of the Newton-Rhapson step for simple kriging
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Vector of three generalized L2 distance weights.
// - threshold - when |lambda_i| < threshold, value is set to 0.
// - tolq   - convergence tolerance for NR step
// - tolp   - convergence tolerance for penalty criteria
// - maxq   - max number of iterations for NR step
// - maxp   - max number of iterations to grow the penalty term.
// - eta    - growth parmameter for penalty criteria (eta > 1).
// Outputs:
// - lam    - a (hopefully) convergent solution for the lambda vector
// [[Rcpp::export]]
arma::colvec newRap_long(arma::colvec lam, const arma::mat & pCovC,
                         const arma::colvec & lCovC, const arma::mat & pCovR,
                         const arma::colvec & lCovR, const arma::mat & pCovCR,
                         const arma::colvec & lCovCR, const arma::colvec & A,
                         const double & threshold,
                         const double & tolq, const int & maxq,
                         const double & tolp, const int & maxp,
                         const double & eta, const double & r) {

  // Determine initial length of lambda vector (size will likely change)
  // with each iteration.
  double len = lam.size();

  // Create a colvec with one entry that stores any warnings in the
  // iterative process.
  arma::colvec isWarn(1, arma::fill::zeros);

  // Initialize variables used in the following loops
  arma::colvec lamUp(len, arma::fill::zeros);
  arma::colvec tdiff(len, arma::fill::zeros);

  // Initialize vector of indices. + 1 ensures that all elements
  // are included initially.
  arma::uvec tindices = find( (arma::abs(lam) + 1) > 0 );

  int i = 0;
  int j = 0; // Reset counter parameters (must redeclare j in new scope)
  double penalty = tolp + 1;
  double adj = r;

  bool lastTry = TRUE;

  while( (std::abs(penalty) > tolp || (i == maxq && lastTry)) && j < maxp ) {
    i=0; // reset the inner counter.

    // If we have found a non-convergent, yet feasible, solution. Then try one
    // more iteration of the penalty parameter.
    if(std::abs(penalty) <= tolp){
      lastTry = FALSE;
    }

    do{
      // Resize our updated lambda matrix
      lamUp.set_size( tindices.size() );

      lamUp = nrStep_long(lam.elem(tindices), pCovC.submat(tindices, tindices),
                          lCovC.elem(tindices), pCovR.submat(tindices, tindices),
                          lCovR.elem(tindices), pCovCR.submat(tindices, tindices),
                          lCovCR.elem(tindices), A, adj, len); // adj is the penalty term

      // Compute the differences in the vector
      tdiff = arma::abs(lam.elem(tindices) - lamUp);

      // Overwrite lambda now that we have calculated the difference
      lam.elem(tindices) = lamUp;

      // Set elements of lambda close enough to 0, equal to 0
      //lam.elem( find(arma::abs(lam) <  threshold) ).zeros();
      lam.clean(threshold);

      tindices = find(lam);

      // If only one lambda remains, return this and suspend the search.
      if( tindices.n_elem < 2 ){
        if(tindices.n_elem < 1){
          Rcpp::Rcout << "No viable lambdas exist, returning prediction of 0..." << std::endl;
          isWarn(0) = 1;
        }

        lam = lam.zeros();
        lam.elem(tindices) += 1;

        arma::colvec lam2 = arma::join_vert(isWarn, lam);

        return lam2;
      }


      i++; // move along the counter
    } while(tdiff.max() > tolq && i < maxq); // end the inner do-loop

    // Calculate the penalty (will determine if we stay in the loop)
    penalty = arma::as_scalar(arma::sum(arma::abs(lam)) - 1);

    // Shrink the adjustment parameter according to the growth parameter eta.
    adj *= eta;

    j++; // move along the j counter

    // Warn the user if the returned solution was the result on a non-convergent
    // optimization.
    if(tdiff.max() > tolq && !lastTry){
      Rcpp::Rcout << "Feasible solution obtained from a non-convergent optimization step " <<
        "with max(diff) = " << tdiff.max() << "." << std::endl;
      isWarn(0) = 1;
    }

  } // End the outer-loop.

  if(std::abs(penalty) > tolp){
    Rcpp::Rcout <<
      "Convergent, feasible solution not obtained with sum(abs(lam)) = " <<
        penalty + 1 << "," << " and max(diff) = " << tdiff.max() << "." << std::endl;
    isWarn(0) = 1;
  }

  arma::colvec lam2 = arma::join_vert(isWarn, lam);

  return lam2;
}


// Calculation of the Newton-Rhapson step for ordinary kriging
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - A      - Generalized L2 distance weights.
// - threshold - when |lambda_i| < threshold, value is set to 0.
// - tolq   - convergence tolerance for NR step
// - tolp   - convergence tolerance for penalty criteria
// - maxq   - max number of iterations for NR step
// - maxp   - max number of iterations for penalty criteria
// - eta    - shrink parmameter for penalty criteria (eta < 1).
// - r      - initial (and somewhat large) penalty parmater that ensures intitial
//            solutions resides in feasible region.
// Outputs:
// - lam    - a (hopefully) convergent solution for the lambda vector
// [[Rcpp::export]]
arma::colvec newRap_2(arma::colvec lam, const arma::mat & pCovC,
                      const arma::colvec & lCovC, const arma::mat & pCovR,
                      const arma::colvec & lCovR, const arma::mat & pCovCR,
                      const arma::colvec & lCovCR, const arma::colvec & A,
                      const double & threshold,
                      const double & tolq, const int & maxq,
                      const double & tolp, const int & maxp,
                      const double & eta, double r) {

  // Determine initial length of lambda vector (size will likely change)
  // with each iteration.
  double len = lam.size();

  // Create a colvec with one entry that stores any warnings in the
  // iterative process.
  arma::colvec isWarn(1, arma::fill::zeros);

  // Initialize variables used in the following loops
  arma::colvec lamUp(len, arma::fill::zeros);
  arma::colvec tdiff(len, arma::fill::zeros);

  arma::uvec tindices = find(lam);

  // Initialize penalty parameter
  double penalty = tolp + 1;

  // Initialize counter parameters
  int i = 0;
  int j = 0;

  // Ensures algorithm does not leave feasible region without warning.
  bool feasible = TRUE;
  bool lastTry = TRUE;

  while( (std::abs(penalty) > tolp || (i == maxq && lastTry)) && j < maxp && feasible){
    i=0; // reset the inner counter.

    // If we have found a non-convergent, yet feasible, solution. Then try one
    // more iteration of the penalty parameter.
    if(std::abs(penalty) <= tolp){
      lastTry = FALSE;
    }

    do{
      // Resize our updated lambda matrix
      lamUp.set_size( tindices.size() );

      lamUp = nrStep_2(lam.elem(tindices), pCovC.submat(tindices, tindices),
                       lCovC.elem(tindices), pCovR.submat(tindices, tindices),
                       lCovR.elem(tindices), pCovCR.submat(tindices, tindices),
                       lCovCR.elem(tindices), A, r, threshold, len); // adj is the penalty term

      // solution MUST be non-negative. It this is violated, recommend a slower
      // shrink parameter.
      if(lamUp.min() < -threshold){
        Rcpp::Rcout <<
          "Left feasible region. Consider using a slower shrink parameter or larger intital penalty." <<
            std::endl;
        feasible = 0;
        isWarn(0) = 1;
        break;
      }

      // Compute the differences in the vector
      tdiff = arma::abs(lam.elem(tindices) - lamUp);

      // Overwrite lambda now that we have calculated the difference
      lam.elem(tindices) = lamUp;

      // Set elements of lambda close enough to 0, equal to 0
      //lam.elem( find(arma::abs(lam) < threshold) ).zeros();
      lam.clean(threshold);

      tindices = find(lam);

      // If only one lambda remains, return this and suspend the search.
      if( tindices.n_elem < 2 ){
        if(tindices.n_elem < 1){
          Rcpp::Rcout << "No viable lambdas exist, returning prediction of 0..." << std::endl;
          isWarn(0) = 1;
        }
        lam = lam.zeros();
        lam.elem(tindices) += 1;

        arma::colvec lam2 = arma::join_vert(isWarn, lam);

        return lam2;
      }

      i++; // move along the counter
    } while(tdiff.max() > tolq && i < maxq); // end the inner do-loop

    // Calculate the penalty (will determine if we stay in the loop)
    penalty = sum(arma::abs(lam)) - 1;

    // Shirnk the penalty parameter according to the growth parameter eta.
    r *= eta;

    j++; // move along the j counter

    // Warn the user if the returned solution was the result on a non-convergent
    // optimization.
    if(tdiff.max() > tolq && !lastTry){
      Rcpp::Rcout << "Feasible solution obtained from a non-convergent optimization step " <<
        "with max(diff) = " << tdiff.max() << "." << std::endl;
      isWarn(0) = 1;
    }
  } // end the outer loop

  if(std::abs(penalty) > tolp){
    Rcpp::Rcout <<
      "Convergent, feasible solution not obtained with sum(abs(lam)) = " <<
        penalty + 1 << "," << " and max(diff) = " << tdiff.max() << "." << std::endl;
    isWarn(0) = 1;
  }

  arma::colvec lam2 = arma::join_vert(isWarn, lam);

  return lam2;
}


// Wrapper to the Newton - Rhapson algorithm that returns a two column matrix
// of interval values predictions at the locations of interest.
// Inputs:
// - lam - vector of lambdas from previous iteration
// - pCovC  - pairwise Covariance matrix for centers
// - pCovR  - pairwise Covariance matrix for Radii
// - pCovCR - pairwise Covariance matrix for Center Radius interaction
// - lCovC  - vector of covariances (centers) at location of interest
// - lCovR  - vector of covariances (radius) at location of interest
// - lCovCR - vector of covariances (center/radius) at location of interest
// - values - A two column matrix holding the centers and radii at the
//            measurement locations.
// - A      - Vector of three generalized L2 distance weights.
// - thresh - when |lambda_i| < abs(1/n*thresh), value is set to 0.
// - tolq   - convergence tolerance for NR step
// - tolp   - convergence tolerance for penalty criteria
// - maxq   - max number of iterations for NR step
// - maxp   - max number of iterations for penalty criteria
// - eta    - growth/shrink parmameter for penalty criteria (eta < 1).
// - r      - initial (and somewhat large) penalty parmater that ensures
//            intitial solutions resides in feasible region. Used for
//            ordinary kriging only.
// - trend  - When defined, used to detrend the centers in the simple
//            kriging calculation. When trend = -9999.99, ordinary kriging
//            is used instead.
// Outputs:
// - predicts - A two column matrix representing the centers and radii at each
//            - resective location of interest.
// [[Rcpp::export]]
arma::mat nrShell(const arma::mat & pCovC, const arma::mat & pCovR,
                  const arma::mat & pCovCR,
                  const arma::mat & lCovC, const arma::mat & lCovR,
                  const arma::mat & lCovCR,
                  const arma::mat & values, // Observations at locations.
                  const arma::colvec & A,
                  const double & thresh,
                  const double & tolq, const int & maxq,
                  const double & tolp, const int & maxp,
                  const double & eta, double r,
                  double trend, const bool & fast){

  int len1 = pCovC.n_rows;
  int len2 = lCovC.n_cols;

  // Calculate the threshold for which values of lambda will be set to 0.
  double len = values.n_rows;
  double threshold = std::abs( 1/(len*thresh) );
  double predVar = 0;

  // Pre-allocate variables we will either fill or re-use over and over.
  arma::mat predicts(len2, 4, arma::fill::zeros);

  for(int i = 0; i < len2; i++){

    // Make initial guess using naive simple kriging estimate of the centers.
    arma::colvec lam = arma::solve(pCovC, lCovC.col(i));

    // Preserve sign information of initial guess.
    arma::colvec lamSgn = arma::sign(lam);
    lamSgn.elem( find(lamSgn == 0) ) += 1; // Elimate any zeros in sign vector.

    // Don't let initial guesses be too close to 0.
    // This should help to avoid singularities.
    lam = arma::abs(lam);
    //lam.elem( find(lam < threshold) ).zeros();
    lam.clean(threshold); // replaces previous line
    lam.elem( find(lam == 0) ) += threshold;

    // If trend is not specified, use ordinary kriging.
    if(trend == -9999.99){

      lam = newRap_2(lam, pCovC, lCovC.col(i), pCovR, lCovR.col(i),
                     pCovCR, lCovCR.col(i), A,
                     threshold, tolq, maxp, tolp, maxq, eta, r);

      // Extract the warning and subset lambda accordingly
      predicts(i, 3) = lam(0);
      lam = lam.subvec(1, lam.n_elem-1);

      // Make predictions for the center.
      predicts(i, 0) = arma::as_scalar(lam.t() * values.col(0));
    }else{

      // Re-introduce negative lambdas in the initial guess for simple kriging.
      lam = lam % lamSgn;

      if(fast){
        lam = newRap(lam, pCovC, lCovC.col(i), pCovR, lCovR.col(i), pCovCR, lCovCR.col(i), A,
                     threshold, tolq, maxp, tolp, maxq, eta);
      }else{
        lam = newRap_long(lam, pCovC, lCovC.col(i), pCovR, lCovR.col(i), pCovCR, lCovCR.col(i), A,
                          threshold, tolq, maxp, tolp, maxq, eta, r);
      }

      // Extract the warning and subset lambda accordingly
      predicts(i, 3) = lam(0);
      lam = lam.subvec(1, lam.n_elem-1);

      // Make predictions for the center, removing the trend.
      predicts(i, 0) = arma::as_scalar(lam.t() * (values.col(0) - trend)) + trend;
    }

    predicts(i, 1) = arma::as_scalar(arma::abs(lam).t() * values.col(1));

    predVar = A(0)*pCovC(0, 0) + A(1)*pCovR(0, 0) + A(2)*pCovCR(0, 0);

    for(int k = 0; k < len1; k++){
      for(int q = 0; q < len1; q++){
        predVar += A(0)*lam(k)*lam(q)*pCovC(k, q) + A(1)*std::abs(lam(k)*lam(q))*pCovR(k, q) +
          A(2)*std::abs(lam(k))*lam(q)*pCovCR(k, q);
      }
      predVar += -2*( A(0)*lam(k)*lCovC(k, i) + A(1)*std::abs(lam(k))*lCovR(k, i) ) -
        A(2)*( lam(k) + std::abs(lam(k)) )*lCovCR(k, i);
    }

    predicts(i, 2) = predVar;

    // Check for suer interrupts
    Rcpp::checkUserInterrupt();
  }

  return predicts;
}







