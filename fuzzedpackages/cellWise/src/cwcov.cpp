// 
// #include "cwcov.h"
// 
// 
// // [[Rcpp::export]]
// Rcpp::List cwcov_cpp(arma::mat & Z,
//                      arma::vec & m0,
//                      arma::mat & S0,
//                      const double & crit = 0.01,
//                      const int & maxits = 5,
//                      const double & quant = 0.99,
//                      const double & maxImp = 0.25) {
//   
//   
//   
//   
//   try
//   {
//     const int n = Z.n_rows;
//     const int d = Z.n_cols;
//     
//     int nbits    = 0;
//     double convcrit = 1;
//     
//     // containters for iterative step
//     arma::cube betamat(n, d + 1, d,arma::fill::zeros);
//     arma::cube Bmat(n, d + 1, std::pow(d, 2),arma::fill::zeros); // matrix containing bias terms, we unfold the bias matrix
//     arma::mat  orderings(n, d, arma::fill::zeros);
//     arma::mat  distances(n, d + 1, arma::fill::zeros);
//     arma::mat  pvals(n, d, arma::fill::ones);
//     
//     arma::vec mu = m0;
//     invMat invOut = mpinv(S0);
//     arma::mat Sigmai = invOut.Inv;
//     arma::mat Sigmaisqrt =  invOut.InvSqrt;
//     
//     
//     while ((nbits < maxits) && (convcrit > crit)) { //start of main loop
//       
//       
//       for (unsigned int i = 0; i < n; i++) {
//         arma::vec z = Z.row(i);
//         arma::vec response = Sigmaisqrt * (z - mu);
//         arma::vec weights = huberweights(z - mu, 1.5);
//         arma::vec signs    <-   1  -2 * ((z - mu) < 0);
//         
//         larOut   <- flagCells(predictors, response,
//                               weights = weights,
//                               signs = signs,
//                               Sigmai = Sigmai)
//         diffRSS  <- abs(diff(larOut$RSS))
//         
//         pvals[i, larOut$ordering] <- 1 - pchisq(diffRSS, 1)
//         pvals[i, larOut$ordering] <- rev(cummin(rev(pvals[i,larOut$ordering])))
//         
//         betamat[i, , ] <- larOut$beta
//         Bmat[i, , , ]  <- larOut$biasmat
//         distances[i, ] <- larOut$RSS
//         orderings[i, ] <- larOut$ordering
//       }
//       
//     }//end of main loop
//     
//     
//   } catch( std::exception& __ex__ )
//   {
//     forward_exception_to_r( __ex__ );
//   } catch(...)
//   {
//     ::Rf_error( "c++ exception " "(unknown reason)" );
//   }
//   return Rcpp::wrap(NA_REAL);
// }
